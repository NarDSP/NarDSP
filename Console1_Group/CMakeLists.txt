cmake_minimum_required(VERSION 3.14.0)

option(SMTG_ENABLE_VST3_PLUGIN_EXAMPLES "Enable VST 3 Plug-in Examples" OFF)
option(SMTG_ENABLE_VST3_HOSTING_EXAMPLES "Enable VST 3 Hosting Examples" OFF)

set(SMTG_RUN_VST_VALIDATOR OFF)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

set(vst3sdk_SOURCE_DIR "C:/DSP/vst-sdk_3.7.12_build-20_2024-07-23/VST_SDK/vst3sdk")
if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Console1_Group
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "Console1_Group VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Console1_Group
    source/version.h
    source/Console1_Groupcids.h
    source/Console1_Groupprocessor.h
    source/Console1_Groupprocessor.cpp
    source/Console1_Groupcontroller.h
    source/Console1_Groupcontroller.cpp
    source/Console1_Groupentry.cpp
)

### FFTW
target_link_libraries(Console1_Group PRIVATE C:/DSP/NarDSP/Console1_Group/libs/fftw/fftw3.lib)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(Console1_Group
        PRIVATE
            resource/Console1_Groupeditor.uidesc
    )
    target_link_libraries(Console1_Group
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Console1_Group
        RESOURCES
            "resource/Console1_Groupeditor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (Console1_Group
    RESOURCES
        resource/E17D08BAC7E55E9FBBF779A061CD1C7B_snapshot.png
        resource/E17D08BAC7E55E9FBBF779A061CD1C7B_snapshot_2.0x.png
)

target_link_libraries(Console1_Group
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Console1_Group)

if(SMTG_MAC)
    smtg_target_set_bundle(Console1_Group
        BUNDLE_IDENTIFIER com.nardsp.console1_group
        COMPANY_NAME "NarDSP"
    )
    smtg_target_set_debug_executable(Console1_Group
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Console1_Group PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Console1_Group)

        smtg_target_set_debug_executable(Console1_Group
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)
