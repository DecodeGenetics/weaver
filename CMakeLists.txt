cmake_minimum_required (VERSION 3.2)
project (weaver LANGUAGES C CXX)

include(ExternalProject)

# Build weaver in release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build." FORCE)
endif()

# Betatest build type is release but also checks for assertions
string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_CXX_FLAGS_BETATEST "${CMAKE_CXX_FLAGS_RELEASE}" )
string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_C_FLAGS_BETATEST "${CMAKE_C_FLAGS_RELEASE}" )

# Profile build type is release with gprof flags
set(CMAKE_CXX_FLAGS_PROFILE "-pg ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "Profiling flags" FORCE)

# Build paw (dependency) in release by default
if(NOT PAW_CMAKE_BUILD_TYPE)
  set(PAW_CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build for paw." FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

##############
## Git info ##
##############
# Get the current working git branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest abbreviated SHA commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_SHORT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_LONG_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git diff COMMAND wc -l
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_NUM_DIRTY_LINES
  OUTPUT_STRIP_TRAILING_WHITESPACE)


############
## Weaver ##
############
set (weaver_VERSION_MAJOR 0)
set (weaver_VERSION_MINOR 2)
set (weaver_VERSION_PATCH 0)

add_subdirectory(src) # Exposes "weaver_sources", which contains all source files of weaver
add_library(weaver_objects OBJECT ${weaver_sources})

target_compile_features(weaver_objects PUBLIC cxx_std_17)
target_compile_options(weaver_objects PUBLIC -Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros)
add_dependencies(weaver_objects project_paw) # weaver depends on building paw

target_include_directories(weaver_objects PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}/include/weaver)

# Add weaver executable
add_executable(weaver src/main.cpp $<TARGET_OBJECTS:weaver_objects>)

# Inherit stuff from weaver_objects
target_include_directories(weaver PRIVATE $<TARGET_PROPERTY:weaver_objects,INTERFACE_INCLUDE_DIRECTORIES>)
target_compile_features(weaver PRIVATE $<TARGET_PROPERTY:weaver_objects,COMPILE_FEATURES>)
target_compile_options(weaver PRIVATE $<TARGET_PROPERTY:weaver_objects,COMPILE_OPTIONS>)

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  ${PROJECT_SOURCE_DIR}/src/in.constants.hpp
  ${PROJECT_BINARY_DIR}/include/weaver/constants.hpp)

target_include_directories(weaver_objects PUBLIC ${PROJECT_BINARY_DIR}/include)

########################
## Build dependencies ##
########################

set (STATIC_DIR "" CACHE STRING "Build in 'static' mode and include libraries in this directory.")

### gfatools ###
ExternalProject_Add(
  project_gfatools
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/gfatools
  BUILD_IN_SOURCE 0
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gfatools
  CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/submodules/gfatools ${CMAKE_CURRENT_BINARY_DIR}/
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/gfatools libgfa1.a
  INSTALL_COMMAND "")

add_library(gfatools STATIC IMPORTED)
target_include_directories(weaver_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/gfatools)
set(gfatools_location ${CMAKE_CURRENT_BINARY_DIR}/gfatools/libgfa1.a)
message(STATUS "gfatools target location is ${gfatools_location}")
set_property(TARGET gfatools PROPERTY IMPORTED_LOCATION ${gfatools_location})
add_dependencies(gfatools project_gfatools)
add_dependencies(weaver_objects gfatools)
target_link_libraries(weaver PUBLIC gfatools)

### paw ###
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  message (STATUS "paw in debug mode")

  ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw
  BUILD_IN_SOURCE 0
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paw
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -H${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw -B${CMAKE_CURRENT_BINARY_DIR}/paw -DCMAKE_BUILD_TYPE=${PAW_CMAKE_BUILD_TYPE} -DFORCE_AVX_FAST=1
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/paw static
  INSTALL_COMMAND "")
else()
  message (STATUS "paw in non-debug mode")

  ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw
  BUILD_IN_SOURCE 0
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paw
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -H${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw -B${CMAKE_CURRENT_BINARY_DIR}/paw -DCMAKE_BUILD_TYPE=${PAW_CMAKE_BUILD_TYPE} -DFORCE_AVX512=1
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/paw static
  INSTALL_COMMAND "")
endif()

add_library(paw STATIC IMPORTED)
target_include_directories(weaver_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw/include)
target_include_directories(weaver_objects SYSTEM PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/paw/include)
set(paw_location ${CMAKE_CURRENT_BINARY_DIR}/paw/lib/libpaw.a)
message(STATUS "paw target location is ${paw_location}")
set_property(TARGET paw PROPERTY IMPORTED_LOCATION ${paw_location})
add_dependencies(paw project_paw)
add_dependencies(weaver_objects paw)
target_link_libraries(weaver PUBLIC paw)

### libdeflate ###
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/submodules/libdeflate/Makefile)
    message(FATAL_ERROR "libdeflate not found. Run 'git submodule update --init' to fetch all submodules.")
endif ()

ExternalProject_Add(
    project_libdeflate
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
    CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/submodules/libdeflate ${CMAKE_CURRENT_BINARY_DIR}/
    BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
        PREFIX="${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "CC=${CMAKE_C_COMPILER}"
        "CFLAGS=-fPIC -O3 ${CMAKE_C_FLAGS}" libdeflate.a
    INSTALL_COMMAND ""
)

add_library(libdeflate STATIC IMPORTED)
set(libdeflate_location ${CMAKE_CURRENT_BINARY_DIR}/libdeflate/libdeflate.a)
message(STATUS "libdeflate target location is ${libdeflate_location}")

set_property(TARGET libdeflate PROPERTY IMPORTED_LOCATION ${libdeflate_location})
add_dependencies(libdeflate project_libdeflate)

### htslib ###
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/submodules/htslib/Makefile)
    message(FATAL_ERROR "htslib not found. Run 'git submodule update --init' to fetch all submodules.")
endif()

ExternalProject_Add(
    project_htslib
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib
    CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/submodules/htslib ${CMAKE_CURRENT_BINARY_DIR}/ COMMAND autoreconf -i COMMAND ${CMAKE_CURRENT_BINARY_DIR}/htslib/configure --disable-libcurl --disable-gcs --disable-lzma --disable-bz2 --with-libdeflate
        "CFLAGS=${MYCFLAGS} -g -Wall -O3 ${CMAKE_C_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "LDFLAGS=${MYLDFLAGS} -L${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "CC=${CMAKE_C_COMPILER}"
    BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/htslib libhts.a
    INSTALL_COMMAND ""
)

add_library(htslib STATIC IMPORTED)


#target_include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/htslib/htslib)

set(htslib_location ${CMAKE_CURRENT_BINARY_DIR}/htslib/libhts.a)
message(STATUS "htslib target location is ${htslib_location}")
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
add_dependencies(htslib project_htslib)
add_dependencies(project_htslib libdeflate)
target_include_directories(weaver_objects PUBLIC SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/submodules/htslib)
target_link_libraries(weaver PUBLIC htslib)
target_link_libraries(weaver PUBLIC libdeflate)

### parallel_hashmap ###
target_include_directories(weaver_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/parallel-hashmap)

### threads ###
find_package(Threads)
target_link_libraries(weaver PUBLIC ${CMAKE_THREAD_LIBS_INIT})

### rt and filesystem ###
if (STATIC_DIR STREQUAL "")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Using GCC")
    target_link_libraries(weaver PUBLIC "rt") # GCC requires linking to rt
    target_link_libraries(weaver PUBLIC "stdc++fs")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(STATUS "Using Clang")
  else()
    message(WARNING "Unsupported compiler")
  endif ()
else()
  target_link_libraries(weaver PUBLIC "${STATIC_DIR}/libpthread.a")

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(weaver PUBLIC "${STATIC_DIR}/librt.a")
    target_link_libraries(weaver PUBLIC "${STATIC_DIR}/libstdc++fs.a")
  endif ()
endif()

### cereal
target_include_directories(weaver_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/cereal/include)

### zstd
if (STATIC_DIR STREQUAL "")
    find_package(ZSTD REQUIRED)
    target_link_libraries(weaver PUBLIC ${ZSTD_LIBRARIES})
else()
    message (STATUS "Using prebuilt zstd")
    target_link_libraries(weaver PUBLIC "${STATIC_DIR}/libzstd.a")
endif()

### zlib ###
if (STATIC_DIR STREQUAL "")
  message (STATUS "Checking for zlib")
  find_package(ZLIB REQUIRED)
  target_include_directories(weaver_objects SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
  target_link_libraries(weaver PUBLIC ${ZLIB_LIBRARIES})
else()
  message (STATUS "Using prebuilt zlib")
  target_link_libraries(weaver PUBLIC "${STATIC_DIR}/libz.a")
endif()

# LOCAL binaries have static GCC, PREBUILT are all static
if (STATIC_DIR STREQUAL "")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

#############
## Testing ##
#############
set(WEAVER_TESTS
  edit
  gfa
  gfa_location
  icu
  index_io
  make_mmi
  region
  sketch
  segment
  sequence_utils
  small_variant
  sr_align
  sr_cigar
  sr_seed
  sr_seed_pair
  stable_contigs)

enable_testing(true)
add_subdirectory(test)

## Add unit tests test
foreach(TEST IN LISTS WEAVER_TESTS)
  add_test(NAME test/test_${TEST} COMMAND test_${TEST})
endforeach()

## Add black box test
add_test(NAME WeaverBlackBox COMMAND sh -c "./weaver index ${CMAKE_CURRENT_SOURCE_DIR}/test/data/test_human_10k.v2.gfa.gz --threads 4 -k 17 -w 5 --vcf=${CMAKE_CURRENT_SOURCE_DIR}/test/data/truth.vcf.gz --log=${CMAKE_CURRENT_BINARY_DIR}/black_box_index.log --vverbose && ./weaver map ${CMAKE_CURRENT_SOURCE_DIR}/test/data/test_human_10k.v2.gfa.gz ${CMAKE_CURRENT_SOURCE_DIR}/test/data/small.read1.fq.gz --fq2=${CMAKE_CURRENT_SOURCE_DIR}/test/data/small.read2.fq.gz -H ${CMAKE_CURRENT_SOURCE_DIR}/test/data/extra_header_lines.tsv -k 17 -w 5 --log=${CMAKE_CURRENT_BINARY_DIR}/black_box_map.log --vverbose --rta3-quals --fastq-data-buffer-size=650 --max-open-files=12 -r chr20 > ${CMAKE_CURRENT_BINARY_DIR}/test/output.sam && cat ${CMAKE_CURRENT_BINARY_DIR}/test/output.sam | md5sum > ${CMAKE_CURRENT_BINARY_DIR}/test/output.sam.md5")

# add_test(NAME WeaverBlackBoxDifferentSoff COMMAND sh -c "./weaver map ${CMAKE_CURRENT_SOURCE_DIR}/test/data/test_b38_chr20.gfa ${CMAKE_CURRENT_SOURCE_DIR}/test/data/small.read1.fq.gz --fq2=${CMAKE_CURRENT_SOURCE_DIR}/test/data/small.read2.fq.gz -k 19 -w 5 --log=${CMAKE_CURRENT_BINARY_DIR}/black_box_different_soff.log -R \"@RG\\tID:foo\\tPL:ILLUMINA\\tSM:bar\" --vverbose --fastq-data-buffer-size=650 --max-open-files=12 | samtools markdup -u - - | samtools view --remove-tag ms -b -o ${CMAKE_CURRENT_BINARY_DIR}/test/output_different_soff.bam && samtools index ${CMAKE_CURRENT_BINARY_DIR}/test/output_different_soff.bam")

################
# clang-format #
################
# Add clang-format check target
find_program(CLANG_FORMAT "clang-format-12")

if (CLANG_FORMAT STREQUAL "" OR CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "" AND NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/test" "-name" "'*.[ch]pp'" "-exec" "${CLANG_FORMAT}" "--style=file" "--dry-run" "-Werror" "{}" "+" COMMENT "Checking all files format with clang-format.")
endif()

if (NOT CLANG_FORMAT STREQUAL "" AND NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  add_custom_target (do_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/test" "-name" "'*.[ch]pp'" "-exec" "${CLANG_FORMAT}" "--style=file" "--i" "{}" "+" COMMENT "Formatting all files with clang-format.")
endif()

#################
# Documentation #
#################
## Adds a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  #install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()
