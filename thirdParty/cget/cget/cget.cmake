set(CGET_PREFIX "/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_sync/SAIGE/thirdParty/cget")
set(CMAKE_PREFIX_PATH "/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_sync/SAIGE/thirdParty/cget")
if (${CMAKE_VERSION} VERSION_LESS "3.6.0")
    include_directories(SYSTEM ${CGET_PREFIX}/include)
    else ()
        set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "${CGET_PREFIX}/include")
        set(CMAKE_C_STANDARD_INCLUDE_DIRECTORIES "${CGET_PREFIX}/include")
endif()
if (CMAKE_CROSSCOMPILING)
    list(APPEND CMAKE_FIND_ROOT_PATH "/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_sync/SAIGE/thirdParty/cget")
endif()
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_sync/SAIGE/thirdParty/cget")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_ENABLE_PARALLEL_BUILD_FLAG "/MP")
endif()
if (BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "ON" CACHE BOOL "")
endif()
set(CMAKE_FIND_FRAMEWORK "LAST" CACHE STRING "")
set(CMAKE_INSTALL_RPATH "${CGET_PREFIX}/lib" CACHE STRING "")
