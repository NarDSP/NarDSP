//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "public.sdk/source/vst/vstaudioeffect.h"

namespace nardsp {

//------------------------------------------------------------------------
//  LowSinc2Processor
//------------------------------------------------------------------------
class LowSinc2Processor : public Steinberg::Vst::AudioEffect
{
public:
	LowSinc2Processor ();
	~LowSinc2Processor () SMTG_OVERRIDE;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new LowSinc2Processor; 
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

	double bL[22];
	double bR[22];
	double f[22];

	double previousSampleL = 0.0;
	double previousSampleR = 0.0;

	double alpha = 0.50264;
	double s = 0.05; // smoothing factor
	double smoothA = 0.0;
	double smoothB = 0.5;
	double smoothC = 1.0;

	uint32_t fpdL;
	uint32_t fpdR;

};

//------------------------------------------------------------------------
} // namespace nardsp
