//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "lowsinc2controller.h"
#include "lowsinc2cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include <base/source/fstreamer.h>

using namespace Steinberg;

namespace nardsp {

//------------------------------------------------------------------------
// LowSinc2Controller Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API LowSinc2Controller::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated

	//---do not forget to call parent ------
	tresult result = EditControllerEx1::initialize (context);
	if (result != kResultOk)
	{
		return result;
	}

	// Here you could register some parameters
	parameters.addParameter(STR16("Cutoff"), STR16(""), 0, 0.0, Vst::ParameterInfo::kCanAutomate, LowSinc2Params::kParamAId);
	parameters.addParameter(STR16("Ripple"), STR16(""), 0, 0.5, Vst::ParameterInfo::kCanAutomate, LowSinc2Params::kParamBId);
	parameters.addParameter(STR16("Inv/Wet"), STR16(""), 0, 1.0, Vst::ParameterInfo::kCanAutomate, LowSinc2Params::kParamCId);

	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API LowSinc2Controller::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API LowSinc2Controller::setComponentState (IBStream* state)
{
	// Here you get the state of the component (Processor part)
	if (!state)
		return kResultFalse;

	float savedParam1 = 0.0f;
	float savedParam2 = 0.0f;
	float savedParam3 = 0.0f;

	IBStreamer streamer(state, kLittleEndian);

	if (streamer.readFloat(savedParam1) == false)
		return kResultFalse;
	if (streamer.readFloat(savedParam2) == false)
		return kResultFalse;
	if (streamer.readFloat(savedParam3) == false)
		return kResultFalse;

	if (auto param = parameters.getParameter(LowSinc2Params::kParamAId))
		param->setNormalized(savedParam1);
	if (auto param = parameters.getParameter(LowSinc2Params::kParamBId))
		param->setNormalized(savedParam2);
	if (auto param = parameters.getParameter(LowSinc2Params::kParamCId))
		param->setNormalized(savedParam3);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API LowSinc2Controller::setState (IBStream* state)
{
	// Here you get the state of the controller

	return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API LowSinc2Controller::getState (IBStream* state)
{
	// Here you are asked to deliver the state of the controller (if needed)
	// Note: the real state of your plug-in is saved in the processor

	return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API LowSinc2Controller::createView (FIDString name)
{
	// Here the Host wants to open your editor (if you have one)
	if (FIDStringsEqual (name, Vst::ViewType::kEditor))
	{
		// create your editor here and return a IPlugView ptr of it
		auto* view = new VSTGUI::VST3Editor (this, "view", "lowsinc2editor.uidesc");
		return view;
	}
	return nullptr;
}

//------------------------------------------------------------------------
} // namespace nardsp
