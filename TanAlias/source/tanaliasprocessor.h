//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <fftw3.h>

#include <array>

namespace nardsp {

	//------------------------------------------------------------------------
	//  TanAliasProcessor
	//------------------------------------------------------------------------
	class TanAliasProcessor : public Steinberg::Vst::AudioEffect
	{
	public:
		TanAliasProcessor();
		~TanAliasProcessor() SMTG_OVERRIDE;

		// Create function
		static Steinberg::FUnknown* createInstance(void* /*context*/)
		{
			return (Steinberg::Vst::IAudioProcessor*)new TanAliasProcessor;
		}

		//--- ---------------------------------------------------------------------
		// AudioEffect overrides:
		//--- ---------------------------------------------------------------------
		/** Called at first after constructor */
		Steinberg::tresult PLUGIN_API initialize(Steinberg::FUnknown* context) SMTG_OVERRIDE;

		/** Called at the end before destructor */
		Steinberg::tresult PLUGIN_API terminate() SMTG_OVERRIDE;

		/** Switch the Plug-in on/off */
		Steinberg::tresult PLUGIN_API setActive(Steinberg::TBool state) SMTG_OVERRIDE;

		/** Will be called before any process call */
		Steinberg::tresult PLUGIN_API setupProcessing(Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;

		/** Asks if a given sample size is supported see SymbolicSampleSizes. */
		Steinberg::tresult PLUGIN_API canProcessSampleSize(Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

		/** Here we go...the process call */
		Steinberg::tresult PLUGIN_API process(Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;

		/** For persistence */
		Steinberg::tresult PLUGIN_API setState(Steinberg::IBStream* state) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API getState(Steinberg::IBStream* state) SMTG_OVERRIDE;

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

		std::array<double, 30> atten;
		std::array<double, 30> prevAtten;

		double tan1 = 1.2073855733208536;
		double tan3 = -0.28401846109755496;
		double tan5 = 0.097513958890027341;
		double tan7 = -0.035250937833235164;
		double tan9 = 0.012873546866224071;
		double tan11 = -0.0047118578854305343;
		double tan13 = 0.0017254603555171694;
		double tan15 = -0.00063192864473000238;
		double tan17 = 0.00023144232846261659;
		double tan19 = -8.4765848930710985e-05;
		double tan21 = 3.1045942009186311e-05;
		double tan23 = -1.1371755605005252e-05;
		double tan25 = 4.1680976686983562e-06;
		double tan27 = -1.5352705573512537e-06;
		double tan29 = 5.8605602260702883e-07;
		double tan31 = -2.7952383654722066e-07;

		double output_pad = 1.0;
		double drive = 0.05;

		uint32_t fpdL;
		uint32_t fpdR;

	};

	//------------------------------------------------------------------------
} // namespace nardsp
