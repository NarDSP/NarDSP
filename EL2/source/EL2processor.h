//------------------------------------------------------------------------
// Copyright(c) 2025 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <fftw3.h> 
#include <array>

namespace nardsp {

//------------------------------------------------------------------------
//  EL2Processor
//------------------------------------------------------------------------
class EL2Processor : public Steinberg::Vst::AudioEffect
{
public:
	EL2Processor ();
	~EL2Processor () SMTG_OVERRIDE;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new EL2Processor; 
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
	double drive = 0.5;

	double highestFrequencyL = 20.0;
	double highestFrequencyR = 20.0;

	double atten[16];
	double prevAtten[16];

	const int fftSize = 512;
	double* fftwInput;
	fftw_complex* fftwOutput;
	fftw_plan fftwPlan;

	// Fixed-size buffers
	std::array<double, 512> in1Fixed;
	std::array<double, 512> in2Fixed;

	double p1_0 = -0.08677868;
	double p1_1 = -0.89591807;
	double p1_2 = -0.0631735;
	double p1_3 = 0.06175486;
	double p1_4 = 0.01628709;
	double p1_5 = -0.01002569;
	double p1_6 = -0.00836354;
	double p1_7 = -0.00281562;
	double p1_8 = -0.00104544;
	double p1_9 = 0.0;

	double p2_0 = 1.43203656e-02;
	double p2_1 = -1.07054432e+00;
	double p2_2 = 2.00901443e-02;
	double p2_3 = 1.33182705e-01;
	double p2_4 = 1.19903839e-02;
	double p2_5 = -1.09978313e-02;
	double p2_6 = -9.84602629e-04;
	double p2_7 = -7.00181071e-03;
	double p2_8 = -7.20520787e-03;
	double p2_9 = -1.57938985e-03;

	double t1_0 = -0.10068046;
	double t1_1 = -0.75860694;
	double t1_2 = -0.11586637;
	double t1_3 = -0.00524329;
	double t1_4 = -0.01477255;
	double t1_5 = 0.00638286;
	double t1_6 = 0.00511599;
	double t1_7 = -0.00163751;
	double t1_8 = 0.00470263;
	double t1_9 = -0.00361433;

	double t2_0 = -0.11690122;
	double t2_1 = -0.72815378;
	double t2_2 = -0.13228972;
	double t2_3 = -0.00526361;
	double t2_4 = -0.01062647;
	double t2_5 = 0.0075646;
	double t2_6 = 0.00608398;
	double t2_7 = -0.00271511;
	double t2_8 = 0.00132195;
	double t2_9 = -0.00239476;

	double c0 = p1_0;
	double c1 = p1_1;
	double c2 = p1_2;
	double c3 = p1_3;
	double c4 = p1_4;
	double c5 = p1_5;
	double c6 = p1_6;
	double c7 = p1_7;
	double c8 = p1_8;
	double c9 = p1_9;

	uint32_t fpdL;
	uint32_t fpdR;

};

//------------------------------------------------------------------------
} // namespace nardsp
