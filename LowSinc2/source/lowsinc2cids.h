//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace nardsp {
//------------------------------------------------------------------------
static const Steinberg::FUID kLowSinc2ProcessorUID (0x4FD7ED46, 0x3E175268, 0x9A32D402, 0x06984F7A);
static const Steinberg::FUID kLowSinc2ControllerUID (0x65D63FE6, 0x483058E5, 0x9D648500, 0x98A75490);

#define LowSinc2VST3Category "NarDSP"

enum LowSinc2Params {

	kParamAId = 1,
	kParamBId = 2,
	kParamCId = 3

};

//------------------------------------------------------------------------
} // namespace nardsp
