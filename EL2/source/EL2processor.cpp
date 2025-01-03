//------------------------------------------------------------------------
// Copyright(c) 2025 NarDSP.
//------------------------------------------------------------------------

#include "EL2processor.h"
#include "EL2cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include <pluginterfaces/vst/ivsteditcontroller.h>

#include <fftw3.h>
#include <array>

using namespace Steinberg;

namespace nardsp {

const double half_pi = 1.5707963267948966;

constexpr double fft_resolution = 1.0 / 512.0;

static double estimateHighestFrequencyFFTW(
	const std::array<double, 512>& inputSignal,
	const double sampleRate,
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
// EL2Processor
//------------------------------------------------------------------------
EL2Processor::EL2Processor ()
{
	A = 0.0;
	B = 0.5;
	C = 0.5;

	for (size_t x = 0; x < 16; x++) {
		atten[x] = 1.0;
		prevAtten[x] = 1.0;
	}

	for (size_t x = 0; x < fftSize; x++) {
		in1Fixed[x] = 0.0;
		in2Fixed[x] = 0.0;
	}

	fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

	fftwInput = (double*)fftw_malloc(sizeof(double) * fftSize);
	fftwOutput = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 257);
	fftwPlan = fftw_plan_dft_r2c_1d(fftSize, fftwInput, fftwOutput, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	//--- set the wanted controller for our processor
	setControllerClass (kEL2ControllerUID);
}

//------------------------------------------------------------------------
EL2Processor::~EL2Processor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated
	
	//---always initialize the parent-------
	tresult result = AudioEffect::initialize (context);
	// if everything Ok, continue
	if (result != kResultOk)
	{
		return result;
	}

	//--- create Audio IO ------
	addAudioInput (STR16 ("Stereo In"), Steinberg::Vst::SpeakerArr::kStereo);
	addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

	/* If you don't need an event bus, you can remove the next line */
	addEventInput (STR16 ("Event In"), 1);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::terminate ()
{
	fftw_destroy_plan(fftwPlan);
	fftw_free(fftwInput);
	fftw_free(fftwOutput);

	//---do not forget to call parent ------
	return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::setActive (TBool state)
{
	//--- called when the Plug-in is enable/disable (On/Off) -----
	return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::process (Vst::ProcessData& data)
{
	if (data.inputParameterChanges)
	{
		int32 numCrashParamsChanged = data.inputParameterChanges->getParameterCount();
		for (int32 index = 0; index < numCrashParamsChanged; index++)
		{
			Vst::IParamValueQueue* paramQueue = data.inputParameterChanges->getParameterData(index);
			if (paramQueue)
			{
				Vst::ParamValue value;
				int32 sampleOffset;
				int32 numPoints = paramQueue->getPointCount();

				switch (paramQueue->getParameterId())
				{
				case EL2Params::kParamAId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
						A = value;
					}
					break;

				case EL2Params::kParamBId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
						B = value;
						drive = 0.01 + 0.95 * B;
					}
					break;

				case EL2Params::kParamCId:
					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
						C = value;
					}
					break;

				}

			}
		}
	}

	int32 samples_to_process = data.numSamples;
	double inv_samples_to_process = 1.0 / samples_to_process;
	int32 chunkSize = std::min(samples_to_process, fftSize);

	double sampleRate = processSetup.sampleRate;

	void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
	void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

	float* in1 = (float*)in[0];
	float* in2 = (float*)in[1];
	float* out1 = (float*)out[0];
	float* out2 = (float*)out[1];

	// Filling fixed-size buffers
	for (size_t i = 0; i < chunkSize; ++i)
	{
		in1Fixed[i] = (double)in1[i];
		in2Fixed[i] = (double)in2[i];
	}

	double t = 0.0;

	double inL = 0.0;
	double inR = 0.0;

	double distL = 0.0;
	double distR = 0.0;

	double outL = 0.0;
	double outR = 0.0;

	int mode = static_cast<int>(2.99 * A + 0.01);
	if (mode == 0) {
		c0 = p1_0;
		c1 = p1_1;
		c2 = p1_2;
		c3 = p1_3;
		c4 = p1_4;
		c5 = p1_5;
		c6 = p1_6;
		c7 = p1_7;
		c8 = p1_8;
		c9 = p1_9;
	}
	else if(mode == 1) {
		c0 = p2_0;
		c1 = p2_1;
		c2 = p2_2;
		c3 = p2_3;
		c4 = p2_4;
		c5 = p2_5;
		c6 = p2_6;
		c7 = p2_7;
		c8 = p2_8;
		c9 = p2_9;
	}
	else if (mode == 2) {
		c0 = t1_0;
		c1 = t1_1;
		c2 = t1_2;
		c3 = t1_3;
		c4 = t1_4;
		c5 = t1_5;
		c6 = t1_6;
		c7 = t1_7;
		c8 = t1_8;
		c9 = t1_9;
	}
	else if (mode == 3) {
		c0 = t2_0;
		c1 = t2_1;
		c2 = t2_2;
		c3 = t2_3;
		c4 = t2_4;
		c5 = t2_5;
		c6 = t2_6;
		c7 = t2_7;
		c8 = t2_8;
		c9 = t2_9;
	}

	highestFrequencyL = estimateHighestFrequencyFFTW(in1Fixed, sampleRate, fftwInput, fftwOutput, fftwPlan);
	highestFrequencyR = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftwInput, fftwOutput, fftwPlan);

	for (size_t i = 0; i < 16; i++) {
		prevAtten[i] = atten[i];
	}

	atten[0] = attenFactor(highestFrequencyL * 2, sampleRate);
	atten[2] = attenFactor(highestFrequencyL * 3, sampleRate);
	atten[4] = attenFactor(highestFrequencyL * 4, sampleRate);
	atten[6] = attenFactor(highestFrequencyL * 5, sampleRate);
	atten[8] = attenFactor(highestFrequencyL * 6, sampleRate);
	atten[10] = attenFactor(highestFrequencyL * 7, sampleRate);
	atten[12] = attenFactor(highestFrequencyL * 8, sampleRate);
	atten[14] = attenFactor(highestFrequencyL * 9, sampleRate);

	atten[1] = attenFactor(highestFrequencyR * 2, sampleRate);
	atten[3] = attenFactor(highestFrequencyR * 3, sampleRate);
	atten[5] = attenFactor(highestFrequencyR * 4, sampleRate);
	atten[7] = attenFactor(highestFrequencyR * 5, sampleRate);
	atten[9] = attenFactor(highestFrequencyR * 6, sampleRate);
	atten[11] = attenFactor(highestFrequencyR * 7, sampleRate);
	atten[13] = attenFactor(highestFrequencyR * 8, sampleRate);
	atten[15] = attenFactor(highestFrequencyR * 9, sampleRate);

	// inv_cA
	t = drive;
	double T[10]{};
	T[0] = 1.0;
	T[1] = t;
	for (int n = 2; n < 10; n++) {
		T[n] = 2.0 * t * T[n - 1] - T[n - 2];
	}

	double cA = c0 + c1 * t;
	cA += c2 * T[2];
	cA += c3 * T[3];
	cA += c4 * T[4];
	cA += c5 * T[5];
	cA += c6 * T[6];
	cA += c7 * T[7];
	cA += c8 * T[8];
	cA += c9 * T[9];
	double inv_cA = 1.0 / cA;

	for (size_t i = 0; i < samples_to_process; i++) {

		double position = (double)i * inv_samples_to_process;
		double common_factor_odd = sin(half_pi * position) * sin(half_pi * position);

		double common_factor_even = sin(half_pi * position);
		common_factor_even = sin(half_pi * common_factor_even) * sin(half_pi * common_factor_even);
		common_factor_even = sin(half_pi * common_factor_even) * sin(half_pi * common_factor_even);

		double atten0 = prevAtten[0] + (atten[0] - prevAtten[0]) * common_factor_even;
		double atten1 = prevAtten[1] + (atten[1] - prevAtten[1]) * common_factor_even;

		double atten2 = prevAtten[2] + (atten[2] - prevAtten[2]) * common_factor_odd;
		double atten3 = prevAtten[3] + (atten[3] - prevAtten[3]) * common_factor_odd;

		double atten4 = prevAtten[4] + (atten[4] - prevAtten[4]) * common_factor_even;
		double atten5 = prevAtten[5] + (atten[5] - prevAtten[5]) * common_factor_even;

		double atten6 = prevAtten[6] + (atten[6] - prevAtten[6]) * common_factor_odd;
		double atten7 = prevAtten[7] + (atten[7] - prevAtten[7]) * common_factor_odd;

		double atten8 = prevAtten[8] + (atten[8] - prevAtten[8]) * common_factor_even;
		double atten9 = prevAtten[9] + (atten[9] - prevAtten[9]) * common_factor_even;

		double atten10 = prevAtten[10] + (atten[10] - prevAtten[10]) * common_factor_odd;
		double atten11 = prevAtten[11] + (atten[11] - prevAtten[11]) * common_factor_odd;

		double atten12 = prevAtten[12] + (atten[12] - prevAtten[12]) * common_factor_even;
		double atten13 = prevAtten[13] + (atten[13] - prevAtten[13]) * common_factor_even;

		double atten14 = prevAtten[14] + (atten[14] - prevAtten[14]) * common_factor_odd;
		double atten15 = prevAtten[15] + (atten[15] - prevAtten[15]) * common_factor_odd;

		// Left channel coefficients
		double atten0c2 = atten0 * c2;
		double atten2c3 = atten2 * c3;
		double atten4c4 = atten4 * c4;
		double atten6c5 = atten6 * c5;
		double atten8c6 = atten8 * c6;
		double atten10c7 = atten10 * c7;
		double atten12c8 = atten12 * c8;
		double atten14c9 = atten14 * c9;

		// Right channel coefficients
		double atten1c2 = atten1 * c2;
		double atten3c3 = atten3 * c3;
		double atten5c4 = atten5 * c4;
		double atten7c5 = atten7 * c5;
		double atten9c6 = atten9 * c6;
		double atten11c7 = atten11 * c7;
		double atten13c8 = atten13 * c8;
		double atten15c9 = atten15 * c9;

		inL = in1[i];
		inR = in2[i];

		// denormal
		if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
		if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

		// Noise
		t = 0.0;
		T[0] = 1.0;
		T[1] = t;
		for (int n = 2; n < 10; n++) {
			T[n] = 2.0 * t * T[n - 1] - T[n - 2];
		}
		double noiseL = c0;
		noiseL += atten0c2 * T[2];
		noiseL += atten2c3 * T[3];
		noiseL += atten4c4 * T[4];
		noiseL += atten6c5 * T[5];
		noiseL += atten8c6 * T[6];
		noiseL += atten10c7 * T[7];
		noiseL += atten12c8 * T[8];
		noiseL += atten14c9 * T[9];

		double noiseR = c0;
		noiseR += atten1c2 * T[2];
		noiseR += atten3c3 * T[3];
		noiseR += atten5c4 * T[4];
		noiseR += atten7c5 * T[5];
		noiseR += atten9c6 * T[6];
		noiseR += atten11c7 * T[7];
		noiseR += atten13c8 * T[8];
		noiseR += atten15c9 * T[9];

		// Left
		t = drive * inL;
		T[0] = 1.0;
		T[1] = t;
		for (int n = 2; n < 10; n++) {
			T[n] = 2.0 * t * T[n - 1] - T[n - 2];
		}

		distL = c0 + c1 * t;
		distL += atten0c2 * T[2];
		distL += atten2c3 * T[3];
		distL += atten4c4 * T[4];
		distL += atten6c5 * T[5];
		distL += atten8c6 * T[6];
		distL += atten10c7 * T[7];
		distL += atten12c8 * T[8];
		distL += atten14c9 * T[9];

		t = drive * inR;
		T[0] = 1.0;
		T[1] = t;
		for (int n = 2; n < 10; n++) {
			T[n] = 2.0 * t * T[n - 1] - T[n - 2];
		}

		distR = c0 + c1 * t;
		distR += atten1c2 * T[2];
		distR += atten3c3 * T[3];
		distR += atten5c4 * T[4];
		distR += atten7c5 * T[5];
		distR += atten9c6 * T[6];
		distR += atten11c7 * T[7];
		distR += atten13c8 * T[8];
		distR += atten15c9 * T[9];

		outL = inv_cA * (distL - noiseL);
		outR = inv_cA * (distR - noiseR);

		outL *= (C + 0.5);
		outR *= (C + 0.5);

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

		out1[i] = std::clamp(outL, -1.2, 1.2);
		out2[i] = std::clamp(outR, -1.2, 1.2);

	}

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----
	return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
	// by default kSample32 is supported
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	// disable the following comment if your processing support kSample64
	/* if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue; */

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API EL2Processor::setState (IBStream* state)
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
tresult PLUGIN_API EL2Processor::getState (IBStream* state)
{
	float toSaveParam1 = A;
	float toSaveParam2 = B;
	float toSaveParam3 = C;

	// here we need to save the model
	IBStreamer streamer(state, kLittleEndian);

	streamer.writeFloat(toSaveParam1);
	streamer.writeFloat(toSaveParam2);
	streamer.writeFloat(toSaveParam3);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace nardsp
