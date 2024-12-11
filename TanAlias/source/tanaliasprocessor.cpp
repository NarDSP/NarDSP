//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "tanaliasprocessor.h"
#include "tanaliascids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

#include <fftw3.h>

#include <array>

using namespace Steinberg;

const double pi = 3.141592653589793;
const double half_pi = 1.5707963267948966;

namespace nardsp {

constexpr double fft_resolution = 1.0 / 512.0;

static double estimateHighestFrequencyFFTW(
	const std::array<double, 512>& inputSignal,
	const double sampleRate,
	const int fftSize,
	double* fftwInput,
	fftw_complex* fftwOutput,
	const fftw_plan fftwPlan
) {

	if (fftwInput == nullptr || fftwOutput == nullptr)
		return 0.0;

	for (int i = 0; i < 512; i++) {
		fftwInput[i] = inputSignal[i];
	}

	fftw_execute(fftwPlan);

	// Find the highest frequency
	double maxSquaredMagnitude = 0.0;
	int highestIndex = 0;
	int halfSize = 257;

	for (int i = 0; i < halfSize; ++i) {
		double real = fftwOutput[i][0];
		double imaginary = fftwOutput[i][1];
		double squaredMagnitude = real * real + imaginary * imaginary;

		if (squaredMagnitude > maxSquaredMagnitude) {
			maxSquaredMagnitude = squaredMagnitude;
			highestIndex = i;
		}
	}

	return highestIndex * sampleRate * fft_resolution;
}

static double attenFactor(double freq, double sampleRate) {
	if (freq < 0.425 * sampleRate) return 1.0;
	return 0.0;
}

//------------------------------------------------------------------------
// TanAliasProcessor
//------------------------------------------------------------------------
TanAliasProcessor::TanAliasProcessor()
{

	A = 0.0;
	B = 0.5;
	C = 1.0;

	for (int x = 0; x < fftSize; x++) {
		in1Fixed[x] = 0.0;
		in2Fixed[x] = 0.0;
	}

	for (int x = 0; x < 30; x++) {
		prevAtten[x] = 0.0;
		atten[x] = 0.0;
	}

	fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

	fftwInput = (double*)fftw_malloc(sizeof(double) * fftSize);
	fftwOutput = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 257);
	fftwPlan = fftw_plan_dft_r2c_1d(fftSize, fftwInput, fftwOutput, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	//--- set the wanted controller for our processor
	setControllerClass(kTanAliasControllerUID);
}

//------------------------------------------------------------------------
TanAliasProcessor::~TanAliasProcessor() = default;

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::initialize(FUnknown* context)
{
	//---always initialize the parent-------
	tresult result = AudioEffect::initialize(context);
	// if everything Ok, continue
	if (result != kResultOk)
	{
		return result;
	}

	//--- create Audio IO ------
	addAudioInput(STR16("Stereo In"), Steinberg::Vst::SpeakerArr::kStereo);
	addAudioOutput(STR16("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

	/* If you don't need an event bus, you can remove the next line */
	addEventInput(STR16("Event In"), 1);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::terminate()
{
	fftw_destroy_plan(fftwPlan);
	fftw_free(fftwInput);
	fftw_free(fftwOutput);

	//---do not forget to call parent ------
	return AudioEffect::terminate();
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::setActive(TBool state)
{
	//--- called when the Plug-in is enabled/disabled (On/Off) -----
	return AudioEffect::setActive(state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::process(Vst::ProcessData& data)
{

	if (data.inputParameterChanges)
	{
		int32 numParamsChanged = data.inputParameterChanges->getParameterCount();
		for (int32 index = 0; index < numParamsChanged; index++)
		{
			Vst::IParamValueQueue* paramQueue = data.inputParameterChanges->getParameterData(index);
			if (paramQueue)
			{
				Vst::ParamValue value;
				int32 sampleOffset;
				int32 numPoints = paramQueue->getPointCount();

				switch (paramQueue->getParameterId()) {

				case TanAliasParams::kParamAId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
						A = value;
					drive = 0.88 * (0.05 + 0.95 * A); // from 0.05 to pi/6
					break;

				case TanAliasParams::kParamBId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
						B = value;
					break;

				case TanAliasParams::kParamCId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
						C = value;
					break;

				}
			}
		}
	}

	int32 samples_to_process = data.numSamples;
	double inv_samples_to_process = 1.0 / samples_to_process;
	int32 chunkSize = std::min(samples_to_process, fftSize);

	double sampleRate = processSetup.sampleRate;

	double process_bits = processSetup.symbolicSampleSize;

	if (process_bits == Vst::kSample64) {

		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

		double* in1 = (double*)in[0];
		double* in2 = (double*)in[1];
		double* out1 = (double*)out[0];
		double* out2 = (double*)out[1];

		int remaining = fftSize;
		// Filling fixed-size buffers
		for (int i = 0; i < chunkSize; ++i)
		{
			in1Fixed[i] = in1[i];
			in2Fixed[i] = in2[i];
			remaining--;
		}
		if (chunkSize < fftSize) {

			// Zero padding
			for (int j = fftSize - remaining; j < fftSize; j++) {
				in1Fixed[j] = 0.0;
				in2Fixed[j] = 0.0;
			}
		}

		// Left channel coefficients
		double atten0tan3 = tan3;
		double atten2tan5 = tan5;
		double atten4tan7 = tan7;
		double atten6tan9 = tan9;
		double atten8tan11 = tan11;
		double atten10tan13 = tan13;
		double atten12tan15 = tan15;
		double atten14tan17 = tan17;
		double atten16tan19 = tan19;
		double atten18tan21 = tan21;
		double atten20tan23 = tan23;
		double atten22tan25 = tan25;
		double atten24tan27 = tan27;
		double atten26tan29 = tan29;
		double atten28tan31 = tan31;

		// Right channel coefficients
		double atten1tan3 = tan3;
		double atten3tan5 = tan5;
		double atten5tan7 = tan7;
		double atten7tan9 = tan9;
		double atten9tan11 = tan11;
		double atten11tan13 = tan13;
		double atten13tan15 = tan15;
		double atten15tan17 = tan17;
		double atten17tan19 = tan19;
		double atten19tan21 = tan21;
		double atten21tan23 = tan23;
		double atten23tan25 = tan25;
		double atten25tan27 = tan27;
		double atten27tan29 = tan29;
		double atten29tan31 = tan31;

		double t = 0.0;

		double distL = 0.0;
		double distR = 0.0;

		double inL = 0.0;
		double inR = 0.0;

		double dryL = 0.0;
		double dryR = 0.0;

		double outL = 0.0;
		double outR = 0.0;

		highestFrequencyL = estimateHighestFrequencyFFTW(in1Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);
		highestFrequencyR = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);

		for (int i = 0; i < 30; i++) {
			prevAtten[i] = atten[i];
		}

		atten[0] = attenFactor(highestFrequencyL * 3, sampleRate);
		atten[2] = attenFactor(highestFrequencyL * 5, sampleRate);
		atten[4] = attenFactor(highestFrequencyL * 7, sampleRate);
		atten[6] = attenFactor(highestFrequencyL * 9, sampleRate);
		atten[8] = attenFactor(highestFrequencyL * 11, sampleRate);
		atten[10] = attenFactor(highestFrequencyL * 13, sampleRate);
		atten[12] = attenFactor(highestFrequencyL * 15, sampleRate);
		atten[14] = attenFactor(highestFrequencyL * 17, sampleRate);
		atten[16] = attenFactor(highestFrequencyL * 19, sampleRate);
		atten[18] = attenFactor(highestFrequencyL * 21, sampleRate);
		atten[20] = attenFactor(highestFrequencyL * 23, sampleRate);
		atten[22] = attenFactor(highestFrequencyL * 25, sampleRate);
		atten[24] = attenFactor(highestFrequencyL * 27, sampleRate);
		atten[26] = attenFactor(highestFrequencyL * 29, sampleRate);
		atten[28] = attenFactor(highestFrequencyL * 31, sampleRate);

		atten[1] = attenFactor(highestFrequencyR * 3, sampleRate);
		atten[3] = attenFactor(highestFrequencyR * 5, sampleRate);
		atten[5] = attenFactor(highestFrequencyR * 7, sampleRate);
		atten[7] = attenFactor(highestFrequencyR * 9, sampleRate);
		atten[9] = attenFactor(highestFrequencyR * 11, sampleRate);
		atten[11] = attenFactor(highestFrequencyR * 13, sampleRate);
		atten[13] = attenFactor(highestFrequencyR * 15, sampleRate);
		atten[15] = attenFactor(highestFrequencyR * 17, sampleRate);
		atten[17] = attenFactor(highestFrequencyR * 19, sampleRate);
		atten[19] = attenFactor(highestFrequencyR * 21, sampleRate);
		atten[21] = attenFactor(highestFrequencyR * 23, sampleRate);
		atten[23] = attenFactor(highestFrequencyR * 25, sampleRate);
		atten[25] = attenFactor(highestFrequencyR * 27, sampleRate);
		atten[27] = attenFactor(highestFrequencyR * 29, sampleRate);
		atten[29] = attenFactor(highestFrequencyR * 31, sampleRate);

		for (size_t i = 0; i < samples_to_process; i++) {

			double position = (double)i * inv_samples_to_process;
			double common_factor = sin(half_pi * position) * sin(half_pi * position);

			double atten0 = prevAtten[0] + (atten[0] - prevAtten[0]) * common_factor;
			double atten1 = prevAtten[1] + (atten[1] - prevAtten[1]) * common_factor;

			double atten2 = prevAtten[2] + (atten[2] - prevAtten[2]) * common_factor;
			double atten3 = prevAtten[3] + (atten[3] - prevAtten[3]) * common_factor;

			double atten4 = prevAtten[4] + (atten[4] - prevAtten[4]) * common_factor;
			double atten5 = prevAtten[5] + (atten[5] - prevAtten[5]) * common_factor;

			double atten6 = prevAtten[6] + (atten[6] - prevAtten[6]) * common_factor;
			double atten7 = prevAtten[7] + (atten[7] - prevAtten[7]) * common_factor;

			double atten8 = prevAtten[8] + (atten[8] - prevAtten[8]) * common_factor;
			double atten9 = prevAtten[9] + (atten[9] - prevAtten[9]) * common_factor;

			double atten10 = prevAtten[10] + (atten[10] - prevAtten[10]) * common_factor;
			double atten11 = prevAtten[11] + (atten[11] - prevAtten[11]) * common_factor;

			double atten12 = prevAtten[12] + (atten[12] - prevAtten[12]) * common_factor;
			double atten13 = prevAtten[13] + (atten[13] - prevAtten[13]) * common_factor;

			double atten14 = prevAtten[14] + (atten[14] - prevAtten[14]) * common_factor;
			double atten15 = prevAtten[15] + (atten[15] - prevAtten[15]) * common_factor;

			double atten16 = prevAtten[16] + (atten[16] - prevAtten[16]) * common_factor;
			double atten17 = prevAtten[17] + (atten[17] - prevAtten[17]) * common_factor;

			double atten18 = prevAtten[18] + (atten[18] - prevAtten[18]) * common_factor;
			double atten19 = prevAtten[19] + (atten[19] - prevAtten[19]) * common_factor;

			double atten20 = prevAtten[20] + (atten[20] - prevAtten[20]) * common_factor;
			double atten21 = prevAtten[21] + (atten[21] - prevAtten[21]) * common_factor;

			double atten22 = prevAtten[22] + (atten[22] - prevAtten[22]) * common_factor;
			double atten23 = prevAtten[23] + (atten[23] - prevAtten[23]) * common_factor;

			double atten24 = prevAtten[24] + (atten[24] - prevAtten[24]) * common_factor;
			double atten25 = prevAtten[25] + (atten[25] - prevAtten[25]) * common_factor;

			double atten26 = prevAtten[26] + (atten[26] - prevAtten[26]) * common_factor;
			double atten27 = prevAtten[27] + (atten[27] - prevAtten[27]) * common_factor;

			double atten28 = prevAtten[28] + (atten[28] - prevAtten[28]) * common_factor;
			double atten29 = prevAtten[29] + (atten[29] - prevAtten[29]) * common_factor;

			// Left channel coefficients
			atten0tan3 = atten0 * tan3;
			atten2tan5 = atten2 * tan5;
			atten4tan7 = atten4 * tan7;
			atten6tan9 = atten6 * tan9;
			atten8tan11 = atten8 * tan11;
			atten10tan13 = atten10 * tan13;
			atten12tan15 = atten12 * tan15;
			atten14tan17 = atten14 * tan17;
			atten16tan19 = atten16 * tan19;
			atten18tan21 = atten18 * tan21;
			atten20tan23 = atten20 * tan23;
			atten22tan25 = atten22 * tan25;
			atten24tan27 = atten24 * tan27;
			atten26tan29 = atten26 * tan29;
			atten28tan31 = atten28 * tan31;

			// Right channel coefficients
			atten1tan3 = atten1 * tan3;
			atten3tan5 = atten3 * tan5;
			atten5tan7 = atten5 * tan7;
			atten7tan9 = atten7 * tan9;
			atten9tan11 = atten9 * tan11;
			atten11tan13 = atten11 * tan13;
			atten13tan15 = atten13 * tan15;
			atten15tan17 = atten15 * tan17;
			atten17tan19 = atten17 * tan19;
			atten19tan21 = atten19 * tan21;
			atten21tan23 = atten21 * tan23;
			atten23tan25 = atten23 * tan25;
			atten25tan27 = atten25 * tan27;
			atten27tan29 = atten27 * tan29;
			atten29tan31 = atten29 * tan31;

			inL = in1[i];
			inR = in2[i];

			// denormal
			if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
			if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

			dryL = inL;
			dryR = inR;

			// inv_sinA
			t = drive;
			double T[32]{};
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			double tanA = tan1 * t;
			tanA += atten0tan3 * T[3];
			tanA += atten2tan5 * T[5];
			tanA += atten4tan7 * T[7];
			tanA += atten6tan9 * T[9];
			tanA += atten8tan11 * T[11];
			tanA += atten10tan13 * T[13];
			tanA += atten12tan15 * T[15];
			tanA += atten14tan17 * T[17];
			tanA += atten16tan19 * T[19];
			tanA += atten18tan21 * T[21];
			tanA += atten20tan23 * T[23];
			tanA += atten22tan25 * T[25];
			tanA += atten24tan27 * T[27];
			tanA += atten26tan29 * T[29];
			tanA += atten28tan31 * T[31];
			double inv_tanA = 1.0 / tanA;

			// Left
			t = drive * inL;
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			distL = tan1 * t;
			distL += atten0tan3 * T[3];
			distL += atten2tan5 * T[5];
			distL += atten4tan7 * T[7];
			distL += atten6tan9 * T[9];
			distL += atten8tan11 * T[11];
			distL += atten10tan13 * T[13];
			distL += atten12tan15 * T[15];
			distL += atten14tan17 * T[17];
			distL += atten16tan19 * T[19];
			distL += atten18tan21 * T[21];
			distL += atten20tan23 * T[23];
			distL += atten22tan25 * T[25];
			distL += atten24tan27 * T[27];
			distL += atten26tan29 * T[29];
			distL += atten28tan31 * T[31];

			t = drive * inR;
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			distR = tan1 * t;
			distR += atten1tan3 * T[3];
			distR += atten3tan5 * T[5];
			distR += atten5tan7 * T[7];
			distR += atten7tan9 * T[9];
			distR += atten9tan11 * T[11];
			distR += atten11tan13 * T[13];
			distR += atten13tan15 * T[15];
			distR += atten15tan17 * T[17];
			distR += atten17tan19 * T[19];
			distR += atten19tan21 * T[21];
			distR += atten21tan23 * T[23];
			distR += atten23tan25 * T[25];
			distR += atten25tan27 * T[27];
			distR += atten27tan29 * T[29];
			distR += atten29tan31 * T[31];

			outL = distL * inv_tanA;
			outR = distR * inv_tanA;

			outL = ((B + 0.5) * outL - dryL) * C + dryL;
			outR = ((B + 0.5) * outR - dryR) * C + dryR;

			fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
			fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;

			out1[i] = std::clamp(outL, -1.4, 1.4);
			out2[i] = std::clamp(outR, -1.4, 1.4);

		}

	}
	else {

		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

		float* in1 = (float*)in[0];
		float* in2 = (float*)in[1];
		float* out1 = (float*)out[0];
		float* out2 = (float*)out[1];

		int remaining = fftSize;
		// Filling fixed-size buffers
		for (int i = 0; i < chunkSize; ++i)
		{
			in1Fixed[i] = in1[i];
			in2Fixed[i] = in2[i];
			remaining--;
		}
		if (chunkSize < fftSize) {

			// Zero padding
			for (int j = fftSize - remaining; j < fftSize; j++) {
				in1Fixed[j] = 0.0;
				in2Fixed[j] = 0.0;
			}
		}

		// Left channel coefficients
		double atten0tan3 = tan3;
		double atten2tan5 = tan5;
		double atten4tan7 = tan7;
		double atten6tan9 = tan9;
		double atten8tan11 = tan11;
		double atten10tan13 = tan13;
		double atten12tan15 = tan15;
		double atten14tan17 = tan17;
		double atten16tan19 = tan19;
		double atten18tan21 = tan21;
		double atten20tan23 = tan23;
		double atten22tan25 = tan25;
		double atten24tan27 = tan27;
		double atten26tan29 = tan29;
		double atten28tan31 = tan31;

		// Right channel coefficients
		double atten1tan3 = tan3;
		double atten3tan5 = tan5;
		double atten5tan7 = tan7;
		double atten7tan9 = tan9;
		double atten9tan11 = tan11;
		double atten11tan13 = tan13;
		double atten13tan15 = tan15;
		double atten15tan17 = tan17;
		double atten17tan19 = tan19;
		double atten19tan21 = tan21;
		double atten21tan23 = tan23;
		double atten23tan25 = tan25;
		double atten25tan27 = tan27;
		double atten27tan29 = tan29;
		double atten29tan31 = tan31;

		double t = 0.0;

		double distL = 0.0;
		double distR = 0.0;

		double inL = 0.0;
		double inR = 0.0;

		double dryL = 0.0;
		double dryR = 0.0;

		double outL = 0.0;
		double outR = 0.0;

		highestFrequencyL = estimateHighestFrequencyFFTW(in1Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);
		highestFrequencyR = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);

		for (int i = 0; i < 30; i++) {
			prevAtten[i] = atten[i];
		}

		atten[0] = attenFactor(highestFrequencyL * 3, sampleRate);
		atten[2] = attenFactor(highestFrequencyL * 5, sampleRate);
		atten[4] = attenFactor(highestFrequencyL * 7, sampleRate);
		atten[6] = attenFactor(highestFrequencyL * 9, sampleRate);
		atten[8] = attenFactor(highestFrequencyL * 11, sampleRate);
		atten[10] = attenFactor(highestFrequencyL * 13, sampleRate);
		atten[12] = attenFactor(highestFrequencyL * 15, sampleRate);
		atten[14] = attenFactor(highestFrequencyL * 17, sampleRate);
		atten[16] = attenFactor(highestFrequencyL * 19, sampleRate);
		atten[18] = attenFactor(highestFrequencyL * 21, sampleRate);
		atten[20] = attenFactor(highestFrequencyL * 23, sampleRate);
		atten[22] = attenFactor(highestFrequencyL * 25, sampleRate);
		atten[24] = attenFactor(highestFrequencyL * 27, sampleRate);
		atten[26] = attenFactor(highestFrequencyL * 29, sampleRate);
		atten[28] = attenFactor(highestFrequencyL * 31, sampleRate);

		atten[1] = attenFactor(highestFrequencyR * 3, sampleRate);
		atten[3] = attenFactor(highestFrequencyR * 5, sampleRate);
		atten[5] = attenFactor(highestFrequencyR * 7, sampleRate);
		atten[7] = attenFactor(highestFrequencyR * 9, sampleRate);
		atten[9] = attenFactor(highestFrequencyR * 11, sampleRate);
		atten[11] = attenFactor(highestFrequencyR * 13, sampleRate);
		atten[13] = attenFactor(highestFrequencyR * 15, sampleRate);
		atten[15] = attenFactor(highestFrequencyR * 17, sampleRate);
		atten[17] = attenFactor(highestFrequencyR * 19, sampleRate);
		atten[19] = attenFactor(highestFrequencyR * 21, sampleRate);
		atten[21] = attenFactor(highestFrequencyR * 23, sampleRate);
		atten[23] = attenFactor(highestFrequencyR * 25, sampleRate);
		atten[25] = attenFactor(highestFrequencyR * 27, sampleRate);
		atten[27] = attenFactor(highestFrequencyR * 29, sampleRate);
		atten[29] = attenFactor(highestFrequencyR * 31, sampleRate);

		for (size_t i = 0; i < samples_to_process; i++) {

			double position = (double)i * inv_samples_to_process;
			double common_factor = sin(half_pi * position) * sin(half_pi * position);

			double atten0 = prevAtten[0] + (atten[0] - prevAtten[0]) * common_factor;
			double atten1 = prevAtten[1] + (atten[1] - prevAtten[1]) * common_factor;

			double atten2 = prevAtten[2] + (atten[2] - prevAtten[2]) * common_factor;
			double atten3 = prevAtten[3] + (atten[3] - prevAtten[3]) * common_factor;

			double atten4 = prevAtten[4] + (atten[4] - prevAtten[4]) * common_factor;
			double atten5 = prevAtten[5] + (atten[5] - prevAtten[5]) * common_factor;

			double atten6 = prevAtten[6] + (atten[6] - prevAtten[6]) * common_factor;
			double atten7 = prevAtten[7] + (atten[7] - prevAtten[7]) * common_factor;

			double atten8 = prevAtten[8] + (atten[8] - prevAtten[8]) * common_factor;
			double atten9 = prevAtten[9] + (atten[9] - prevAtten[9]) * common_factor;

			double atten10 = prevAtten[10] + (atten[10] - prevAtten[10]) * common_factor;
			double atten11 = prevAtten[11] + (atten[11] - prevAtten[11]) * common_factor;

			double atten12 = prevAtten[12] + (atten[12] - prevAtten[12]) * common_factor;
			double atten13 = prevAtten[13] + (atten[13] - prevAtten[13]) * common_factor;

			double atten14 = prevAtten[14] + (atten[14] - prevAtten[14]) * common_factor;
			double atten15 = prevAtten[15] + (atten[15] - prevAtten[15]) * common_factor;

			double atten16 = prevAtten[16] + (atten[16] - prevAtten[16]) * common_factor;
			double atten17 = prevAtten[17] + (atten[17] - prevAtten[17]) * common_factor;

			double atten18 = prevAtten[18] + (atten[18] - prevAtten[18]) * common_factor;
			double atten19 = prevAtten[19] + (atten[19] - prevAtten[19]) * common_factor;

			double atten20 = prevAtten[20] + (atten[20] - prevAtten[20]) * common_factor;
			double atten21 = prevAtten[21] + (atten[21] - prevAtten[21]) * common_factor;

			double atten22 = prevAtten[22] + (atten[22] - prevAtten[22]) * common_factor;
			double atten23 = prevAtten[23] + (atten[23] - prevAtten[23]) * common_factor;

			double atten24 = prevAtten[24] + (atten[24] - prevAtten[24]) * common_factor;
			double atten25 = prevAtten[25] + (atten[25] - prevAtten[25]) * common_factor;

			double atten26 = prevAtten[26] + (atten[26] - prevAtten[26]) * common_factor;
			double atten27 = prevAtten[27] + (atten[27] - prevAtten[27]) * common_factor;

			double atten28 = prevAtten[28] + (atten[28] - prevAtten[28]) * common_factor;
			double atten29 = prevAtten[29] + (atten[29] - prevAtten[29]) * common_factor;

			// Left channel coefficients
			atten0tan3 = atten0 * tan3;
			atten2tan5 = atten2 * tan5;
			atten4tan7 = atten4 * tan7;
			atten6tan9 = atten6 * tan9;
			atten8tan11 = atten8 * tan11;
			atten10tan13 = atten10 * tan13;
			atten12tan15 = atten12 * tan15;
			atten14tan17 = atten14 * tan17;
			atten16tan19 = atten16 * tan19;
			atten18tan21 = atten18 * tan21;
			atten20tan23 = atten20 * tan23;
			atten22tan25 = atten22 * tan25;
			atten24tan27 = atten24 * tan27;
			atten26tan29 = atten26 * tan29;
			atten28tan31 = atten28 * tan31;

			// Right channel coefficients
			atten1tan3 = atten1 * tan3;
			atten3tan5 = atten3 * tan5;
			atten5tan7 = atten5 * tan7;
			atten7tan9 = atten7 * tan9;
			atten9tan11 = atten9 * tan11;
			atten11tan13 = atten11 * tan13;
			atten13tan15 = atten13 * tan15;
			atten15tan17 = atten15 * tan17;
			atten17tan19 = atten17 * tan19;
			atten19tan21 = atten19 * tan21;
			atten21tan23 = atten21 * tan23;
			atten23tan25 = atten23 * tan25;
			atten25tan27 = atten25 * tan27;
			atten27tan29 = atten27 * tan29;
			atten29tan31 = atten29 * tan31;

			inL = in1[i];
			inR = in2[i];

			// denormal
			if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
			if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

			dryL = inL;
			dryR = inR;

			// inv_sinA
			t = drive;
			double T[32]{};
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			double tanA = tan1 * t;
			tanA += atten0tan3 * T[3];
			tanA += atten2tan5 * T[5];
			tanA += atten4tan7 * T[7];
			tanA += atten6tan9 * T[9];
			tanA += atten8tan11 * T[11];
			tanA += atten10tan13 * T[13];
			tanA += atten12tan15 * T[15];
			tanA += atten14tan17 * T[17];
			tanA += atten16tan19 * T[19];
			tanA += atten18tan21 * T[21];
			tanA += atten20tan23 * T[23];
			tanA += atten22tan25 * T[25];
			tanA += atten24tan27 * T[27];
			tanA += atten26tan29 * T[29];
			tanA += atten28tan31 * T[31];
			double inv_tanA = 1.0 / tanA;

			// Left
			t = drive * inL;
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			distL = tan1 * t;
			distL += atten0tan3 * T[3];
			distL += atten2tan5 * T[5];
			distL += atten4tan7 * T[7];
			distL += atten6tan9 * T[9];
			distL += atten8tan11 * T[11];
			distL += atten10tan13 * T[13];
			distL += atten12tan15 * T[15];
			distL += atten14tan17 * T[17];
			distL += atten16tan19 * T[19];
			distL += atten18tan21 * T[21];
			distL += atten20tan23 * T[23];
			distL += atten22tan25 * T[25];
			distL += atten24tan27 * T[27];
			distL += atten26tan29 * T[29];
			distL += atten28tan31 * T[31];

			t = drive * inR;
			T[0] = 1.0;
			T[1] = t;
			for (int n = 2; n < 32; n++) {
				T[n] = 2.0 * t * T[n - 1] - T[n - 2];
			}

			distR = tan1 * t;
			distR += atten1tan3 * T[3];
			distR += atten3tan5 * T[5];
			distR += atten5tan7 * T[7];
			distR += atten7tan9 * T[9];
			distR += atten9tan11 * T[11];
			distR += atten11tan13 * T[13];
			distR += atten13tan15 * T[15];
			distR += atten15tan17 * T[17];
			distR += atten17tan19 * T[19];
			distR += atten19tan21 * T[21];
			distR += atten21tan23 * T[23];
			distR += atten23tan25 * T[25];
			distR += atten25tan27 * T[27];
			distR += atten27tan29 * T[29];
			distR += atten29tan31 * T[31];

			outL = distL * inv_tanA;
			outR = distR * inv_tanA;

			outL = ((B + 0.5) * outL - dryL) * C + dryL;
			outR = ((B + 0.5) * outR - dryR) * C + dryR;

			//begin 32 bit stereo floating point dither
			int expon;
			uint32_t bits;
			bits = *reinterpret_cast<uint64_t*>(&outL);
			expon = ((bits >> 23) & 0xFF) - 127;
			fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
			uint64_t factor = 1ULL << (expon + 62);
			outL += ((double(fpdL) - uint32_t(0x7fffffff)) * 6e-37l * factor);
			bits = *reinterpret_cast<uint64_t*>(&outR);
			expon = ((bits >> 23) & 0xFF) - 127;
			fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;
			factor = 1ULL << (expon + 62);
			outR += ((double(fpdR) - uint32_t(0x7fffffff)) * 6e-37l * factor);
			//end 32 bit stereo floating point dither

			out1[i] = std::clamp(outL, -1.4, 1.4);
			out2[i] = std::clamp(outR, -1.4, 1.4);

		}

	}

	return kResultOk;

}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::setupProcessing(Vst::ProcessSetup& newSetup)
{

	//--- called before any processing ----
	return AudioEffect::setupProcessing(newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::canProcessSampleSize(int32 symbolicSampleSize)
{
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue;

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::setState(IBStream* state)
{
	if (!state)
		return kResultFalse;

	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer(state, kLittleEndian);

	float savedParam1 = 0.0f;
	float savedParam2 = 0.0f;
	float savedParam3 = 0.0f;

	if (streamer.readFloat(savedParam1) == false)
		return kResultFalse;
	if (streamer.readFloat(savedParam2) == false)
		return kResultFalse;
	if (streamer.readFloat(savedParam3) == false)
		return kResultFalse;

	A = savedParam1;
	B = savedParam2;
	C = savedParam3;

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasProcessor::getState(IBStream* state)
{
	auto toSaveParam1 = static_cast<float>(A);
	auto toSaveParam2 = static_cast<float>(B);
	auto toSaveParam3 = static_cast<float>(C);

	// here we need to save the model
	IBStreamer streamer(state, kLittleEndian);

	streamer.writeFloat(toSaveParam1);
	streamer.writeFloat(toSaveParam2);
	streamer.writeFloat(toSaveParam3);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace nardsp
