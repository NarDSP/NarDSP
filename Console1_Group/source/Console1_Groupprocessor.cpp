//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "Console1_Groupprocessor.h"
#include "Console1_Groupcids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

#include <fftw3.h>

#include <array>
#include <vector>

using namespace Steinberg;

namespace nardsp {

constexpr double pi = 3.141592653589793;
constexpr double half_pi = 0.5 * pi;

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
// Console1_GroupProcessor
//------------------------------------------------------------------------
Console1_GroupProcessor::Console1_GroupProcessor()
{

	for (int x = 0; x < fftSize; x++) {
		in1Fixed[x] = 0.0;
		in2Fixed[x] = 0.0;
	}

	for (int x = 0; x < 8; x++) {
		atten[x] = 0.0;
	}

	fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

	fftwInput = (double*)fftw_malloc(sizeof(double) * fftSize);
	fftwOutput = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 257);
	fftwPlan = fftw_plan_dft_r2c_1d(fftSize, fftwInput, fftwOutput, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	//--- set the wanted controller for our processor
	setControllerClass(kConsole1_GroupControllerUID);
}

//------------------------------------------------------------------------
Console1_GroupProcessor::~Console1_GroupProcessor()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::initialize(FUnknown* context)
{
	// Here the Plug-in will be instantiated

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
tresult PLUGIN_API Console1_GroupProcessor::terminate()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return AudioEffect::terminate();
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::setActive(TBool state)
{
	//--- called when the Plug-in is enable/disable (On/Off) -----
	return AudioEffect::setActive(state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::process(Vst::ProcessData& data)
{

	int32 numSamples = data.numSamples;
	double inv_samples_to_process = 1.0 / numSamples;

	int32 chunkSize = std::min(numSamples, fftSize);

	double sampleRate = processSetup.sampleRate;

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

	double t = 0.0;
	double t2 = 0.0;
	double t3 = 0.0;
	double t5 = 0.0;

	double distL = 0.0;
	double distR = 0.0;

	double inL = 0.0;
	double inR = 0.0;

	double outL = 0.0;
	double outR = 0.0;

	highestFrequencyL = estimateHighestFrequencyFFTW(in1Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);
	highestFrequencyR = estimateHighestFrequencyFFTW(in2Fixed, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);

	atten[4] = atten[0];
	atten[5] = atten[1];
	atten[6] = atten[2];
	atten[7] = atten[3];

	atten[0] = attenFactor(highestFrequencyL * 3.0, sampleRate);
	atten[1] = attenFactor(highestFrequencyR * 3.0, sampleRate);
	atten[2] = attenFactor(highestFrequencyL * 5.0, sampleRate);
	atten[3] = attenFactor(highestFrequencyR * 5.0, sampleRate);

	double attenL3 = 1.0;
	double attenR3 = 1.0;
	double attenL5 = 1.0;
	double attenR5 = 1.0;

	for (size_t i = 0; i < numSamples; i++) {

		double position = (double)i * inv_samples_to_process;

		double sqrt_alpha = sin(half_pi * position);
		double alpha = sqrt_alpha * sqrt_alpha;

		inL = in1[i];
		inR = in2[i];

		if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
		if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;

		attenL3 = atten[4] + alpha * (atten[0] - atten[4]);
		attenR3 = atten[5] + alpha * (atten[1] - atten[5]);
		attenL5 = atten[6] + alpha * (atten[2] - atten[6]);
		attenR5 = atten[7] + alpha * (atten[3] - atten[7]);

		t = 0.7142857142857143 * inL;
		t2 = t * t;
		t3 = t * t2;
		t5 = t3 * t2;

		distL = 1.6910658468027662 * t;
		distL += attenL3 * 0.13547110619736719 * (4.0 * t3 - 3.0 * t);
		distL += attenL5 * -0.06962696354433665 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);

		t = 0.7142857142857143 * inR;
		t2 = t * t;
		t3 = t * t2;
		t5 = t3 * t2;

		distR = 1.6910658468027662 * t;
		distR += attenR3 * 0.13547110619736719 * (4.0 * t3 - 3.0 * t);
		distR += attenR5 * -0.06962696354433665 * (16.0 * t5 - 20.0 * t3 + 5.0 * t);

		outL = std::min(std::max(0.85 * distL, -1.1), 1.1);
		outR = std::min(std::max(0.85 * distR, -1.1), 1.1);

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

		out1[i] = outL;
		out2[i] = outR;

	}

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::setupProcessing(Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----
	return AudioEffect::setupProcessing(newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::canProcessSampleSize(int32 symbolicSampleSize)
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
tresult PLUGIN_API Console1_GroupProcessor::setState(IBStream* state)
{
	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer(state, kLittleEndian);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupProcessor::getState(IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer(state, kLittleEndian);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace nardsp
