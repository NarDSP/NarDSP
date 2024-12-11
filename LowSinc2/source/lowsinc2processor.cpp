//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "lowsinc2processor.h"
#include "lowsinc2cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

constexpr double pi = 3.141592653589793;

using namespace Steinberg;

namespace nardsp {

	static double sinc(double x) {
		if (x == 0) return 1.0;
		return std::sin(pi * x) / (pi * x);
	}
	//------------------------------------------------------------------------
	// LowSinc2Processor
	//------------------------------------------------------------------------
	LowSinc2Processor::LowSinc2Processor()
	{

		A = 0.0;
		B = 0.5;
		C = 1.0;

		for (int i = 0; i < 22; i++) {
			f[i] = 0.0;
			bL[i] = 0.0;
			bR[i] = 0.0;
		}

		fpdL = 1.0; while (fpdL < 16386) fpdL = rand() * UINT32_MAX;
		fpdR = 1.0; while (fpdR < 16386) fpdR = rand() * UINT32_MAX;

		//--- set the wanted controller for our processor
		setControllerClass(kLowSinc2ControllerUID);
	}

	//------------------------------------------------------------------------
	LowSinc2Processor::~LowSinc2Processor()
	{}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::initialize(FUnknown* context)
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
	tresult PLUGIN_API LowSinc2Processor::terminate()
	{
		// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

		//---do not forget to call parent ------
		return AudioEffect::terminate();
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::setActive(TBool state)
	{
		//--- called when the Plug-in is enable/disable (On/Off) -----

		previousSampleL = 0.0;
		previousSampleR = 0.0;

		for (int i = 0; i < 22; i++) {
			bL[i] = 0.0;
			bR[i] = 0.0;
		}

		return AudioEffect::setActive(state);
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::process(Vst::ProcessData& data)
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
					case LowSinc2Params::kParamAId:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
							A = value;
						break;

					case LowSinc2Params::kParamBId:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
							B = value;
						break;

					case LowSinc2Params::kParamCId:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
							C = value;
						break;

					}
				}
			}
		}

		// Flush case: only update the parameters, no processing is possible
		if (data.numInputs == 0 || data.numOutputs == 0)
			return kResultOk;

		int32 numSamples = data.numSamples;

		// Return if there are no samples to process
		if (numSamples <= 0)
			return kResultOk;

		int TAPS = int(smoothB * 19.99 + 0.01) + 1;

		double sum = 0.0;
		for (int j = 0; j < TAPS; j++) {
			f[j] = sinc(smoothA * j);
			sum += f[j];
		}
		for (int n = 0; n < TAPS; ++n) {
			f[n] /= sum;
		}

		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

		if (processSetup.symbolicSampleSize == Vst::kSample32) {

			Vst::Sample32* in1 = (Vst::Sample32*)in[0];
			Vst::Sample32* in2 = (Vst::Sample32*)in[1];
			Vst::Sample32* out1 = (Vst::Sample32*)out[0];
			Vst::Sample32* out2 = (Vst::Sample32*)out[1];

			double wet = smoothC * 2.0 - 1.0;

			for (int i = 0; i < numSamples; i++) {

				smoothA = s * (0.99 * A + 0.01) + (1.0 - s) * smoothA;
				smoothB = s * B + (1.0 - s) * smoothB;
				smoothC = s * C + (1.0 - s) * smoothC;

				double inL = in1[i];
				double inR = in2[i];

				double dryL = inL;
				double dryR = inR;

				for (int x = TAPS - 1; x >= 0; x--) {
					bL[x + 1] = bL[x];
					bR[x + 1] = bR[x];
				}
				bL[0] = inL;
				bR[0] = inR;

				double FXL = 0.0;
				double FXR = 0.0;
				for (int n = 0; n < TAPS; n++) {
					FXL += f[n] * bL[n];
					FXR += f[n] * bR[n];
				}

				inL = FXL;
				inR = FXR;

				out1[i] = inL * wet + dryL * (1.0 - wet);
				out2[i] = inR * wet + dryR * (1.0 - wet);

				out1[i] = out1[i] * alpha + previousSampleL * (1.0 - alpha);
				out2[i] = out2[i] * alpha + previousSampleR * (1.0 - alpha);

				previousSampleL = out1[i];
				previousSampleR = out2[i];
			}

			return kResultOk;

		}

		else if (processSetup.symbolicSampleSize == Vst::kSample64) {

			Vst::Sample64* in1 = (Vst::Sample64*)in[0];
			Vst::Sample64* in2 = (Vst::Sample64*)in[1];
			Vst::Sample64* out1 = (Vst::Sample64*)out[0];
			Vst::Sample64* out2 = (Vst::Sample64*)out[1];

			double wet = smoothC * 2.0 - 1.0;

			for (int i = 0; i < numSamples; i++) {

				smoothA = s * (0.99 * A + 0.01) + (1.0 - s) * smoothA;
				smoothB = s * B + (1.0 - s) * smoothB;
				smoothC = s * C + (1.0 - s) * smoothC;

				double inL = in1[i];
				double inR = in2[i];

				double dryL = inL;
				double dryR = inR;

				for (int x = TAPS - 1; x >= 0; x--) {
					bL[x + 1] = bL[x];
					bR[x + 1] = bR[x];
				}
				bL[0] = inL;
				bR[0] = inR;

				double FXL = 0.0;
				double FXR = 0.0;
				for (int n = 0; n < TAPS; n++) {
					FXL += f[n] * bL[n];
					FXR += f[n] * bR[n];
				}

				inL = FXL;
				inR = FXR;

				out1[i] = inL * wet + dryL * (1.0 - wet);
				out2[i] = inR * wet + dryR * (1.0 - wet);

				out1[i] = out1[i] * alpha + previousSampleL * (1.0 - alpha);
				out2[i] = out2[i] * alpha + previousSampleR * (1.0 - alpha);

				previousSampleL = out1[i];
				previousSampleR = out2[i];
			}

			return kResultOk;

		}

		return kResultOk;

	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::setupProcessing(Vst::ProcessSetup& newSetup)
	{
		//--- called before any processing ----
		return AudioEffect::setupProcessing(newSetup);
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::canProcessSampleSize(int32 symbolicSampleSize)
	{
		// by default kSample32 is supported
		if (symbolicSampleSize == Vst::kSample32)
			return kResultTrue;

		// disable the following comment if your processing support kSample64
		if (symbolicSampleSize == Vst::kSample64)
			return kResultTrue;

		return kResultFalse;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API LowSinc2Processor::setState(IBStream* state)
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
	tresult PLUGIN_API LowSinc2Processor::getState(IBStream* state)
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
