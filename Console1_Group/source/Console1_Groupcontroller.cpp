//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#include "Console1_Groupcontroller.h"
#include "Console1_Groupcids.h"
#include "vstgui/plugin-bindings/vst3editor.h"

using namespace Steinberg;

namespace nardsp {

//------------------------------------------------------------------------
// Console1_GroupController Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupController::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated

	//---do not forget to call parent ------
	tresult result = EditControllerEx1::initialize (context);
	if (result != kResultOk)
	{
		return result;
	}

	// Here you could register some parameters

	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupController::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupController::setComponentState (IBStream* state)
{
	// Here you get the state of the component (Processor part)
	if (!state)
		return kResultFalse;

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupController::setState (IBStream* state)
{
	// Here you get the state of the controller

	return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API Console1_GroupController::getState (IBStream* state)
{
	// Here you are asked to deliver the state of the controller (if needed)
	// Note: the real state of your plug-in is saved in the processor

	return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API Console1_GroupController::createView (FIDString name)
{
	// No GUI
	return nullptr;

	/*
	// Here the Host wants to open your editor (if you have one)
	if (FIDStringsEqual (name, Vst::ViewType::kEditor))
	{
		// create your editor here and return a IPlugView ptr of it
		auto* view = new VSTGUI::VST3Editor (this, "view", "Console1_Groupeditor.uidesc");
		return view;
	}
	return nullptr;
	*/
}

//------------------------------------------------------------------------
} // namespace nardsp
