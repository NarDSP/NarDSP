//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <fftw3.h>

#include <array>

namespace nardsp {

//------------------------------------------------------------------------
//  SinAliasProcessor
//------------------------------------------------------------------------
class SinAliasProcessor : public Steinberg::Vst::AudioEffect
{
public:
	SinAliasProcessor ();
	~SinAliasProcessor () SMTG_OVERRIDE;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new SinAliasProcessor; 
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

	Steinberg::Vst::ParamValue A;
	Steinberg::Vst::ParamValue B;
	Steinberg::Vst::ParamValue C;

	double highestFrequencyL = 20.0;
	double highestFrequencyR = 20.0;

	double prevInL = 0.0;
	double prevInR = 0.0;
	double prevOutL = 0.0;
	double prevOutR = 0.0;

	double one_over_five_A = 1.0;

	const int fftSize = 512;
	double* fftwInput;
	fftw_complex* fftwOutput;
	fftw_plan fftwPlan;

	// Fixed-size buffers
	std::array<double, 512> in1Fixed;
	std::array<double, 512> in2Fixed;

	std::array<double, 8> atten;
	std::array<double, 8> prevAtten;

	double sin1 = 1.133470896065218;
	double sin3 = -0.13788414624347689;
	double sin5 = 0.0044798166408466459;
	double sin7 = -6.7466721314357603e-05;
	double sin9 = 5.898035473961283e-07;

	double output_pad = 1.0;
	double drive = 0.05;

	uint32_t fpdL;
	uint32_t fpdR;

};

//------------------------------------------------------------------------
} // namespace nardsp
