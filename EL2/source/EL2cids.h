//------------------------------------------------------------------------
// Copyright(c) 2025 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace nardsp {
//------------------------------------------------------------------------
static const Steinberg::FUID kEL2ProcessorUID (0x836B53F3, 0x7E2E5EA6, 0x9432F0D7, 0xB49575BE);
static const Steinberg::FUID kEL2ControllerUID (0x3A160508, 0x1767519C, 0x9E5F6435, 0xF063B1EA);

#define EL2VST3Category "Fx"

enum EL2Params {

	kParamAId = 1,
	kParamBId = 2,
	kParamCId = 3,

};

//------------------------------------------------------------------------
} // namespace nardsp
