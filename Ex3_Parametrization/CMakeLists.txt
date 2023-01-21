cmake_minimum_required(VERSION 3.10)

set(targetName Ex3_Parametrization)

project(${targetName})

# OpenMesh
set(OpenMesh_DIR ../OpenMesh/src)
include_directories(${OpenMesh_DIR})

# Gmm++
set(Gmm_DIR ../gmm/include)
include_directories(${Gmm_DIR})

## Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines 

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wno-error  -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    message("Setting gcc-specific options")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG(-march=native  COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(MARCH_NATIVE "-march=native")
  else()
    set(MARCH_NATIVE "")
  endif()
  SET(CMAKE_CXX_FLAGS_RELEASE        "${MARCH_NATIVE} -O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (-DNOMINMAX)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


set(OpenMesh_SRC_DIR
        ${OpenMesh_DIR}/OpenMesh
        ${OpenMesh_DIR}/OpenMesh/Core
        ${OpenMesh_DIR}/OpenMesh/Core/Geometry
        ${OpenMesh_DIR}/OpenMesh/Core/IO
        ${OpenMesh_DIR}/OpenMesh/Core/IO/exporter
        ${OpenMesh_DIR}/OpenMesh/Core/IO/importer
        ${OpenMesh_DIR}/OpenMesh/Core/IO/reader
        ${OpenMesh_DIR}/OpenMesh/Core/IO/writer
        ${OpenMesh_DIR}/OpenMesh/Core/Mesh
        ${OpenMesh_DIR}/OpenMesh/Core/Mesh/gen
        ${OpenMesh_DIR}/OpenMesh/Core/System
        ${OpenMesh_DIR}/OpenMesh/Core/Templates
        ${OpenMesh_DIR}/OpenMesh/Core/Utils
        )

file(GLOB PROJECT_SRC *.cpp
        ${OpenMesh_DIR}/OpenMesh/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/Geometry/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/IO/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/IO/exporter/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/IO/importer/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/IO/reader/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/IO/writer/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/Mesh/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/Mesh/gen/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/System/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/Templates/*.cc
        ${OpenMesh_DIR}/OpenMesh/Core/Utils/*.cc
       )

# Create an executable
add_executable(${targetName} ${PROJECT_SRC})
