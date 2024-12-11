//------------------------------------------------------------------------
// Copyright(c) 2024 NarDSP.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace nardsp {
//------------------------------------------------------------------------
static const Steinberg::FUID kConsole1_SingleProcessorUID (0x70249425, 0x21985A04, 0x8FEC1356, 0x83D52568);
static const Steinberg::FUID kConsole1_SingleControllerUID (0x6AE3460F, 0x39D4543C, 0x85B3BD9E, 0x6312736A);

#define Console1_SingleVST3Category "NarDSP"

//------------------------------------------------------------------------
} // namespace nardsp
