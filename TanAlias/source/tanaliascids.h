//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace nardsp {
//------------------------------------------------------------------------
static const Steinberg::FUID kTanAliasProcessorUID (0xD4AB1DD7, 0x5DC3586F, 0x9242767E, 0x720D22E9);
static const Steinberg::FUID kTanAliasControllerUID (0x79F00898, 0xE4BC507A, 0x93C5B156, 0x73A60189);

#define TanAliasVST3Category "NarDSP"

enum TanAliasParams {

	kParamAId = 1,
	kParamBId = 2,
	kParamCId = 3,

};

//------------------------------------------------------------------------
} // namespace nardsp
