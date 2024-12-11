//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "sinaliasprocessor.h"
#include "sinaliascids.h"

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
// SinAliasProcessor
//------------------------------------------------------------------------
SinAliasProcessor::SinAliasProcessor()
{

	A = 0.0;
	B = 0.5;
	C = 1.0;

	for (int x = 0; x < fftSize; x++) {
		in1Fixed[x] = 0.0;
		in2Fixed[x] = 0.0;
	}

	for (int x = 0; x < 8; x++) {
		prevAtten[x] = 0.0;
		atten[x] = 0.0;
	}

	fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

	fftwInput = (double*)fftw_malloc(sizeof(double) * fftSize);
	fftwOutput = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 257);
	fftwPlan = fftw_plan_dft_r2c_1d(fftSize, fftwInput, fftwOutput, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	//--- set the wanted controller for our processor
	setControllerClass(kSinAliasControllerUID);
}

//------------------------------------------------------------------------
SinAliasProcessor::~SinAliasProcessor() = default;

//------------------------------------------------------------------------
tresult PLUGIN_API SinAliasProcessor::initialize(FUnknown* context)
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
tresult PLUGIN_API SinAliasProcessor::terminate()
{
	fftw_destroy_plan(fftwPlan);
	fftw_free(fftwInput);
	fftw_free(fftwOutput);

	//---do not forget to call parent ------
	return AudioEffect::terminate();
}

//------------------------------------------------------------------------
tresult PLUGIN_API SinAliasProcessor::setActive(TBool state)
{
	//--- called when the Plug-in is enabled/disabled (On/Off) -----
	return AudioEffect::setActive(state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API SinAliasProcessor::process(Vst::ProcessData& data)
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

				case SinAliasParams::kParamAId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
						A = value;
						drive = 0.05 + 0.95 * A; // from 0.05 to pi/6
					break;

				case SinAliasParams::kParamBId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
						B = value;
					break;

				case SinAliasParams::kParamCId:
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
		double atten0sin3 = sin3;
		double atten2sin5 = sin5;
		double atten4sin7 = sin7;
		double atten6sin9 = sin9;

		// Right channel coefficients
		double atten1sin3 = sin3;
		double atten3sin5 = sin5;
		double atten5sin7 = sin7;
		double atten7sin9 = sin9;

		double t = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		double t5 = 0.0;
		double t7 = 0.0;
		double t9 = 0.0;

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

		prevAtten[0] = atten[0];
		prevAtten[1] = atten[1];
		prevAtten[2] = atten[2];
		prevAtten[3] = atten[3];
		prevAtten[4] = atten[4];
		prevAtten[5] = atten[5];
		prevAtten[6] = atten[6];
		prevAtten[7] = atten[7];

		atten[0] = attenFactor(highestFrequencyL * 3, sampleRate);
		atten[2] = attenFactor(highestFrequencyL * 5, sampleRate);
		atten[4] = attenFactor(highestFrequencyL * 7, sampleRate);
		atten[6] = attenFactor(highestFrequencyL * 9, sampleRate);

		atten[1] = attenFactor(highestFrequencyR * 3, sampleRate);
		atten[3] = attenFactor(highestFrequencyR * 5, sampleRate);
		atten[5] = attenFactor(highestFrequencyR * 7, sampleRate);
		atten[7] = attenFactor(highestFrequencyR * 9, sampleRate);

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

			// Left channel coefficients
			atten0sin3 = atten0 * sin3;
			atten2sin5 = atten2 * sin5;
			atten4sin7 = atten4 * sin7;
			atten6sin9 = atten6 * sin9;

			// Right channel coefficients
			atten1sin3 = atten1 * sin3;
			atten3sin5 = atten3 * sin5;
			atten5sin7 = atten5 * sin7;
			atten7sin9 = atten7 * sin9;

			inL = in1[i];
			inR = in2[i];

			// denormal
			if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
			if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

			dryL = inL;
			dryR = inR;

			// inv_sinA
			t = drive;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;
			double sinA = sin1 * t;
			sinA += atten0sin3 * (4.0 * t3 - 3.0 * t);
			sinA += atten2sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			sinA += atten4sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			sinA += atten6sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);
			double inv_sinA = 1.0 / sinA;

			// Left
			t = drive * inL;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;

			distL = sin1 * t;
			distL += atten0sin3 * (4.0 * t3 - 3.0 * t);
			distL += atten2sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			distL += atten4sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			distL += atten6sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);

			// Right
			t = drive * inR;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;

			distR = sin1 * t;
			distR += atten1sin3 * (4.0 * t3 - 3.0 * t);
			distR += atten3sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			distR += atten5sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			distR += atten7sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);

			outL = distL * inv_sinA;
			outR = distR * inv_sinA;

			outL = ((B + 0.5) * outL - dryL) * C + dryL;
			outR = ((B + 0.5) * outR - dryR) * C + dryR;

			// 64 bit dither is overkill
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
		double atten0sin3 = sin3;
		double atten2sin5 = sin5;
		double atten4sin7 = sin7;
		double atten6sin9 = sin9;

		// Right channel coefficients
		double atten1sin3 = sin3;
		double atten3sin5 = sin5;
		double atten5sin7 = sin7;
		double atten7sin9 = sin9;

		double t = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		double t5 = 0.0;
		double t7 = 0.0;
		double t9 = 0.0;

		double distL = 0.0;
		double distR = 0.0;

		double inL = 0.0;
		double inR = 0.0;

		double dryL = 0.0;
		double dryR = 0.0;

		double outL = 0.0;
		double outR = 0.0;

		highestFrequencyL = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);
		highestFrequencyR = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);

		prevAtten[0] = atten[0];
		prevAtten[1] = atten[1];
		prevAtten[2] = atten[2];
		prevAtten[3] = atten[3];
		prevAtten[4] = atten[4];
		prevAtten[5] = atten[5];
		prevAtten[6] = atten[6];
		prevAtten[7] = atten[7];

		atten[0] = attenFactor(highestFrequencyL * 3, sampleRate);
		atten[2] = attenFactor(highestFrequencyL * 5, sampleRate);
		atten[4] = attenFactor(highestFrequencyL * 7, sampleRate);
		atten[6] = attenFactor(highestFrequencyL * 9, sampleRate);

		atten[1] = attenFactor(highestFrequencyR * 3, sampleRate);
		atten[3] = attenFactor(highestFrequencyR * 5, sampleRate);
		atten[5] = attenFactor(highestFrequencyR * 7, sampleRate);
		atten[7] = attenFactor(highestFrequencyR * 9, sampleRate);

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

			// Left channel coefficients
			atten0sin3 = atten0 * sin3;
			atten2sin5 = atten2 * sin5;
			atten4sin7 = atten4 * sin7;
			atten6sin9 = atten6 * sin9;

			// Right channel coefficients
			atten1sin3 = atten1 * sin3;
			atten3sin5 = atten3 * sin5;
			atten5sin7 = atten5 * sin7;
			atten7sin9 = atten7 * sin9;

			inL = in1[i];
			inR = in2[i];

			// denormal
			if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
			if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

			dryL = inL;
			dryR = inR;

			// inv_sinA
			t = drive;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;
			double sinA = sin1 * t;
			sinA += atten0sin3 * (4.0 * t3 - 3.0 * t);
			sinA += atten2sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			sinA += atten4sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			sinA += atten6sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);
			double inv_sinA = 1.0 / sinA;

			// Left
			t = drive * inL;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;

			distL = sin1 * t;
			distL += atten0sin3 * (4.0 * t3 - 3.0 * t);
			distL += atten2sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			distL += atten4sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			distL += atten6sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);

			// Right
			t = drive * inR;
			t2 = t * t;
			t3 = t * t2;
			t5 = t3 * t2;
			t7 = t5 * t2;
			t9 = t7 * t2;

			distR = sin1 * t;
			distR += atten1sin3 * (4.0 * t3 - 3.0 * t);
			distR += atten3sin5 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);
			distR += atten5sin7 * (64.0 * t7 - 112.0 * t5 + 56.0 * t3 - 7.0 * t);
			distR += atten7sin9 * (256.0 * t9 - 576.0 * t7 + 432.0 * t5 - 120.0 * t3 + 9.0 * t);

			outL = distL * inv_sinA;
			outR = distR * inv_sinA;

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
tresult PLUGIN_API SinAliasProcessor::setupProcessing(Vst::ProcessSetup& newSetup)
{

	//--- called before any processing ----
	return AudioEffect::setupProcessing(newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API SinAliasProcessor::canProcessSampleSize(int32 symbolicSampleSize)
{
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue;

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API SinAliasProcessor::setState(IBStream* state)
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
tresult PLUGIN_API SinAliasProcessor::getState(IBStream* state)
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
