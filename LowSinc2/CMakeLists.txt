cmake_minimum_required(VERSION 3.14.0)

option(SMTG_ENABLE_VST3_PLUGIN_EXAMPLES "Enable VST 3 Plug-in Examples" OFF)
option(SMTG_ENABLE_VST3_HOSTING_EXAMPLES "Enable VST 3 Hosting Examples" OFF)

set(SMTG_RUN_VST_VALIDATOR OFF)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

set(vst3sdk_SOURCE_DIR "C:/DSP/vst-sdk_3.7.12_build-20_2024-07-23/VST_SDK/vst3sdk")
if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(LowSinc2
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "LowSinc2 VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(LowSinc2
    source/version.h
    source/lowsinc2cids.h
    source/lowsinc2processor.h
    source/lowsinc2processor.cpp
    source/lowsinc2controller.h
    source/lowsinc2controller.cpp
    source/lowsinc2entry.cpp
)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(LowSinc2
        PRIVATE
            resource/lowsinc2editor.uidesc
    )
    target_link_libraries(LowSinc2
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(LowSinc2
        RESOURCES
            "resource/lowsinc2editor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (LowSinc2
    RESOURCES
        resource/4FD7ED463E1752689A32D40206984F7A_snapshot.png
        resource/4FD7ED463E1752689A32D40206984F7A_snapshot_2.0x.png
)

target_link_libraries(LowSinc2
    PRIVATE
        sdk
)

smtg_target_configure_version_file(LowSinc2)

if(SMTG_MAC)
    smtg_target_set_bundle(LowSinc2
        BUNDLE_IDENTIFIER com.nardsp.lowsinc2
        COMPANY_NAME "NarDSP"
    )
    smtg_target_set_debug_executable(LowSinc2
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(LowSinc2 PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT LowSinc2)

        smtg_target_set_debug_executable(LowSinc2
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)
