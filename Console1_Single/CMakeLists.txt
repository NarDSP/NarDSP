cmake_minimum_required(VERSION 3.14.0)

option(SMTG_ENABLE_VST3_PLUGIN_EXAMPLES "Enable VST 3 Plug-in Examples" OFF)
option(SMTG_ENABLE_VST3_HOSTING_EXAMPLES "Enable VST 3 Hosting Examples" OFF)

set(SMTG_RUN_VST_VALIDATOR OFF)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

set(vst3sdk_SOURCE_DIR "C:/DSP/vst-sdk_3.7.12_build-20_2024-07-23/VST_SDK/vst3sdk")
if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Console1_Single
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "Console1_Single VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Console1_Single
    source/version.h
    source/Console1_Singlecids.h
    source/Console1_Singleprocessor.h
    source/Console1_Singleprocessor.cpp
    source/Console1_Singlecontroller.h
    source/Console1_Singlecontroller.cpp
    source/Console1_Singleentry.cpp
)

### FFTW
target_link_libraries(Console1_Single PRIVATE C:/DSP/NarDSP/Console1_Single/libs/fftw/fftw3.lib)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(Console1_Single
        PRIVATE
            resource/Console1_Singleeditor.uidesc
    )
    target_link_libraries(Console1_Single
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Console1_Single
        RESOURCES
            "resource/Console1_Singleeditor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (Console1_Single
    RESOURCES
        resource/7024942521985A048FEC135683D52568_snapshot.png
        resource/7024942521985A048FEC135683D52568_snapshot_2.0x.png
)

target_link_libraries(Console1_Single
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Console1_Single)

if(SMTG_MAC)
    smtg_target_set_bundle(Console1_Single
        BUNDLE_IDENTIFIER com.nardsp.console1_single
        COMPANY_NAME "NarDSP"
    )
    smtg_target_set_debug_executable(Console1_Single
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Console1_Single PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Console1_Single)

        smtg_target_set_debug_executable(Console1_Single
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)
