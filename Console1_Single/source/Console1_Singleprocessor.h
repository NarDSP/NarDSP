//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <fftw3.h>

#include <array>

namespace nardsp {

//------------------------------------------------------------------------
//  Console1_SingleProcessor
//------------------------------------------------------------------------
class Console1_SingleProcessor : public Steinberg::Vst::AudioEffect
{
public:
	Console1_SingleProcessor ();
	~Console1_SingleProcessor () SMTG_OVERRIDE;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new Console1_SingleProcessor; 
	}

	//--- ---------------------------------------------------------------------
	// AudioEffect overrides:
	//--- ---------------------------------------------------------------------
	/** Called at first after constructor */
	Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
	
	/** Called at the end before destructor */
	Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;
	
	/** Switch the Plug-in on/off */
	Steinberg::tresult PLUGIN_API setActive (Steinberg::TBool state) SMTG_OVERRIDE;

	/** Will be called before any process call */
	Steinberg::tresult PLUGIN_API setupProcessing (Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;
	
	/** Asks if a given sample size is supported see SymbolicSampleSizes. */
	Steinberg::tresult PLUGIN_API canProcessSampleSize (Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

	/** Here we go...the process call */
	Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
		
	/** For persistence */
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;

//------------------------------------------------------------------------
protected:

	double highestFrequencyL = 20.0;
	double highestFrequencyR = 20.0;

	std::array<double, 1024> in1;
	std::array<double, 1024> in2;

	std::array<double, 8> atten;

	double x[2];
	double y[2];

	double b[3];
	double a[3];

	double x1_L = 0.0;
	double x2_L = 0.0;
	double y1_L = 0.0;
	double y2_L = 0.0;

	double ap_x1_L = 0.0;
	double ap_x2_L = 0.0;
	double ap_y1_L = 0.0;
	double ap_y2_L = 0.0;

	double x1_R = 0.0;
	double x2_R = 0.0;
	double y1_R = 0.0;
	double y2_R = 0.0;

	double ap_x1_R = 0.0;
	double ap_x2_R = 0.0;
	double ap_y1_R = 0.0;
	double ap_y2_R = 0.0;

	double* fftwInput;
	fftw_complex* fftwOutput;
	fftw_plan fftwPlan;
	const int fftSize = 1024;

	uint32_t fpdL;
	uint32_t fpdR;

};

//------------------------------------------------------------------------
} // namespace nardsp
