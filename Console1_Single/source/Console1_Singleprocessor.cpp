//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "Console1_Singleprocessor.h"
#include "Console1_Singlecids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

#include <fftw3.h>

#include <array>
#include <vector>

#include <immintrin.h>

using namespace Steinberg;

namespace nardsp {

const double pi = 3.141592653589793;

constexpr double fft_resolution = 1.0 / 1024.0;

static double estimateHighestFrequencyFFTW(
	const std::array<double, 1024>& inputSignal,
	const double sampleRate,
	const int fftSize,
	double* fftwInput,
	fftw_complex* fftwOutput,
	const fftw_plan fftwPlan
) {

	if (fftwInput == nullptr || fftwOutput == nullptr)
		return 0.0;

	std::memcpy(fftwInput, inputSignal.data(), fftSize * sizeof(double));

	fftw_execute(fftwPlan);

	// Find the highest frequency
	double maxSquaredMagnitude = 0.0;
	int highestIndex = 0;
	int halfSize = 513; // 1024/2 + 1

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

static double attenFactor(const double freq, const double sampleRate) {
	
	if (freq < 0.5 * 0.9 * sampleRate) return 1.0;
	return 0.0;

	/*
	if (freq <= 0.9 * 0.5 * sampleRate) {
		double y = cos(pi * freq / (0.9 * sampleRate));
		return y * y;
	}
	return 0.0;
	*/
}

//------------------------------------------------------------------------
// Console1_SingleProcessor
//------------------------------------------------------------------------
Console1_SingleProcessor::Console1_SingleProcessor ()
{

	for (size_t n = 0; n < fftSize; n++) {
		in1[n] = 0.0;
		in2[n] = 0.0;
	}

	for (size_t n = 0; n < 8; n++) {
		atten[n] = 0.0;
	}

	x[0] = 0.0;
	x[1] = 0.0;
	y[0] = 0.0;
	y[1] = 0.0;

	b[0] = 0.0;
	b[1] = 0.0;
	b[2] = 0.0;

	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;

	fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
	fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

	fftwInput = (double*)fftw_malloc(sizeof(double) * fftSize);
	fftwOutput = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 513); // 1024/2 + 1
	fftwPlan = fftw_plan_dft_r2c_1d(fftSize, fftwInput, fftwOutput, FFTW_MEASURE);
	std::fill(fftwInput, fftwInput + fftSize, 0.0);

	//--- set the wanted controller for our processor
	setControllerClass (kConsole1_SingleControllerUID);
}

//------------------------------------------------------------------------
Console1_SingleProcessor::~Console1_SingleProcessor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::initialize (FUnknown* context)
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
tresult PLUGIN_API Console1_SingleProcessor::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!
	
	//---do not forget to call parent ------
	return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::setActive (TBool state)
{
	//--- called when the Plug-in is enable/disable (On/Off) -----
	return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::process (Vst::ProcessData& data)
{

	int32 numSamples = data.numSamples;

	int32 chunkSize = std::min(numSamples, fftSize);

	double sampleRate = processSetup.sampleRate;

	void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
	void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

	float* in1Variable = (float*)in[0];
	float* in2Variable = (float*)in[1];

	// Filling fixed-size buffers
	for (size_t i = 0; i < chunkSize; i++) {
		in1[i] = in1Variable[i];
		in2[i] = in2Variable[i];
	}

	float* out1 = (float*)out[0];
	float* out2 = (float*)out[1];

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

	highestFrequencyL = estimateHighestFrequencyFFTW(in1, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);
	highestFrequencyR = estimateHighestFrequencyFFTW(in2, sampleRate, fftSize, fftwInput, fftwOutput, fftwPlan);

	const double L3 = attenFactor(highestFrequencyL * 3, sampleRate);
	const double R3 = attenFactor(highestFrequencyR * 3, sampleRate);
	const double L5 = attenFactor(highestFrequencyL * 5, sampleRate);
	const double R5 = attenFactor(highestFrequencyR * 5, sampleRate);

	int32 samples_to_process = data.numSamples;

	const double half_pi = 1.5707963267948966;

	// Setup EQ

	double f0 = 8000.0; // 7000.0
	double Q = 0.7071067811865476;

	double wc = 2.0 * pi * f0 / 44100.0;
	double alpha_filter = Q * sin(wc);

	double cwc = 1.0 - cos(wc);
	b[0] = 0.5 * cwc;
	b[1] = 1.0 - cwc;
	b[2] = b[0];
	a[0] = 1.0 + alpha_filter;
	a[1] = -2.0 * cwc;
	a[2] = 1.0 - alpha_filter;

	double inv_a0 = 1.0 / a[0];
	b[0] *= inv_a0;
	b[1] *= inv_a0;
	b[2] *= inv_a0;
	a[1] *= inv_a0;
	a[2] *= inv_a0;

	double gain = 0.105;
	constexpr double output_pad = 1.0 - 0.5 * 0.105;

	// End setup EQ

	for (size_t i = 0; i < numSamples; ++i) {

		// Set up attenuation factors

		double position = (double)i / samples_to_process;
		double alpha = sin(half_pi * position) * sin(half_pi * position);

		atten[0] = atten[0] + alpha * (L3 - atten[0]);
		atten[1] = atten[1] + alpha * (R3 - atten[1]);
		atten[2] = atten[2] + alpha * (L5 - atten[2]);
		atten[3] = atten[3] + alpha * (R5 - atten[3]);

		// Begin processing

		inL = in1[i];
		inR = in2[i];

		if (fabs(inL) < 1.18e-23) inL = fpdL * 1.18e-17;
		if (fabs(inR) < 1.18e-23) inR = fpdR * 1.18e-17;
		fpdL ^= fpdL << 13; fpdL ^= fpdL >> 17; fpdL ^= fpdL << 5;
		fpdR ^= fpdR << 13; fpdR ^= fpdR >> 17; fpdR ^= fpdR << 5;

		t = 0.7142857142857143 * inL;
		t2 = t * t;
		t3 = t * t2;
		t5 = t3 * t2;

		distL = 1.0758135261590016 * t;
		distL += atten[0] * (-0.46446567770981484 * t3 + 0.34834925828236113 * t);
		distL += atten[2] * (0.3432396617563104 * t5 - 0.42904957719538794 * t3 + 0.10726239429884699 * t);

		t = 0.7142857142857143 * inR;
		t2 = t * t;
		t3 = t * t2;
		t5 = t3 * t2;

		distR = 1.0758135261590016 * t;
		distR += atten[1] * (-0.46446567770981484 * t3 + 0.34834925828236113 * t);
		distR += atten[3] * (0.3432396617563104 * t5 - 0.42904957719538794 * t3 + 0.10726239429884699 * t);

		// Left channel EQ

		double ap_output_L = distL - alpha_filter * (ap_x1_L + ap_x2_L);
		ap_x2_L = ap_x1_L;
		ap_x1_L = ap_output_L;
		ap_y2_L = ap_y1_L;
		ap_y1_L = ap_output_L;

		double mix_L = 0.5 * (ap_output_L + distL);

		double bp_output_L = b[0] * mix_L + b[1] * x1_L + b[2] * x2_L
			- a[1] * y1_L - a[2] * y2_L;

		x2_L = x1_L;
		x1_L = mix_L;
		y2_L = y1_L;
		y1_L = bp_output_L;

		outL = distL + (gain * bp_output_L);

		// Right channel EQ

		double ap_output_R = distR - alpha_filter * (ap_x1_R + ap_x2_R);
		ap_x2_R = ap_x1_R;
		ap_x1_R = ap_output_R;
		ap_y2_R = ap_y1_R;
		ap_y1_R = ap_output_R;

		double mix_R = 0.5 * (ap_output_R + distR);

		double bp_output_R = b[0] * mix_R + b[1] * x1_R + b[2] * x2_R
			- a[1] * y1_R - a[2] * y2_R;

		x2_R = x1_R;
		x1_R = mix_R;
		y2_R = y1_R;
		y1_R = bp_output_R;

		outR = distR + (gain * bp_output_R);

		outL = std::min(std::max(1.16 * output_pad * outL, -1.05), 1.05);
		outR = std::min(std::max(1.16 * output_pad * outR, -1.05), 1.05);

		out1[i] = static_cast<float>(outL);
		out2[i] = static_cast<float>(outR);

	}

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::setupProcessing (Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----
	return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::canProcessSampleSize (int32 symbolicSampleSize)
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
tresult PLUGIN_API Console1_SingleProcessor::setState (IBStream* state)
{
	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer (state, kLittleEndian);
	
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_SingleProcessor::getState (IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer (state, kLittleEndian);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace nardsp
