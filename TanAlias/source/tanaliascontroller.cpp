//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "tanaliascontroller.h"
#include "tanaliascids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include <base/source/fstreamer.h>

using namespace Steinberg;

namespace nardsp {

//------------------------------------------------------------------------
// TanAliasController Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasController::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated

	//---do not forget to call parent ------
	tresult result = EditControllerEx1::initialize (context);
	if (result != kResultOk)
	{
		return result;
	}

	// Here you could register some parameters
	parameters.addParameter(STR16("A"), STR16(""), 0, 0.0, Vst::ParameterInfo::kCanAutomate, TanAliasParams::kParamAId);
	parameters.addParameter(STR16("B"), STR16(""), 0, 0.5, Vst::ParameterInfo::kCanAutomate, TanAliasParams::kParamBId);
	parameters.addParameter(STR16("C"), STR16(""), 0, 1.0, Vst::ParameterInfo::kCanAutomate, TanAliasParams::kParamCId);

	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasController::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasController::setComponentState (IBStream* state)
{
	// Here, you get the state of the component (Processor part)
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

	if (auto param = parameters.getParameter(TanAliasParams::kParamAId))
		param->setNormalized(savedParam1);
	if (auto param = parameters.getParameter(TanAliasParams::kParamBId))
		param->setNormalized(savedParam2);
	if (auto param = parameters.getParameter(TanAliasParams::kParamCId))
		param->setNormalized(savedParam3);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasController::setState (IBStream* state)
{
	// Here you get the state of the controller

	return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API TanAliasController::getState (IBStream* state)
{
	// Here you are asked to deliver the state of the controller (if needed)
	// Note: the real state of your plug-in is saved in the processor

	return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API TanAliasController::createView (FIDString name)
{
	// Here the Host wants to open your editor (if you have one)
	if (FIDStringsEqual (name, Vst::ViewType::kEditor))
	{
		// create your editor here and return a IPlugView ptr of it
		auto* view = new VSTGUI::VST3Editor (this, "view", "tanaliaseditor.uidesc");
		return view;
	}
	return nullptr;
}

//------------------------------------------------------------------------
} // namespace nardsp
