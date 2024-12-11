//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace nardsp {
//------------------------------------------------------------------------
static const Steinberg::FUID kSinAliasProcessorUID (0x829E13F3, 0x7D845CCB, 0xBE6B7C17, 0xB44283BC);
static const Steinberg::FUID kSinAliasControllerUID (0x5DD4CFAD, 0x52575CB7, 0x8A1EDDF9, 0x0043DF3A);

#define SinAliasVST3Category "NarDSP"

enum SinAliasParams {

	kParamAId = 1,
	kParamBId = 2,
	kParamCId = 3

};

//------------------------------------------------------------------------
} // namespace nardsp
