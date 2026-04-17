cmake_minimum_required(VERSION 3.2)

set(CATCH_BUILD_STATIC_LIBRARY 1 CACHE BOOL "Build Catch2 as static" FORCE)
add_subdirectory(../submodules/catch catch)


############
# catch-wv #
############
# catch-wv is an object library that contains the entire catch library. This should be the only time
# that catch is compiled.

## Compile test.cpp as the main Catch file
add_library(catch-wv OBJECT test.cpp)

## Specify includes
target_include_directories(catch-wv PUBLIC
  ${PROJECT_SOURCE_DIR}/submodules/catch/single_include/
  ${PROJECT_SOURCE_DIR}/submodules/gfatools/
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include)

## Inherit compile options and features from weaver_objects
target_compile_features(catch-wv PUBLIC $<TARGET_PROPERTY:weaver_objects,COMPILE_FEATURES>)
target_compile_options(catch-wv PUBLIC $<TARGET_PROPERTY:weaver_objects,COMPILE_OPTIONS>)

###############
# test_weaver #
###############
# test_weaver contains all unit tests for weaver.
## Unit test files
foreach(TEST IN LISTS WEAVER_TESTS)
  add_executable(test_${TEST}
    test_${TEST}.cpp
    $<TARGET_OBJECTS:catch-wv>
    $<TARGET_OBJECTS:weaver_objects>)

  ## Inherit compile options and includes from weaver_objects
  target_include_directories(test_${TEST} PUBLIC
    $<TARGET_PROPERTY:weaver_objects,INTERFACE_INCLUDE_DIRECTORIES>
    ${PROJECT_SOURCE_DIR}/submodules/catch/single_include/
    ${PROJECT_SOURCE_DIR}/include)

  target_compile_features(test_${TEST} PUBLIC $<TARGET_PROPERTY:weaver_objects,COMPILE_FEATURES>)
  target_compile_options(test_${TEST} PUBLIC $<TARGET_PROPERTY:weaver_objects,COMPILE_OPTIONS>)

  ## Inherit linking libraries from weaver
  target_link_libraries(test_${TEST} PUBLIC $<TARGET_PROPERTY:weaver,INTERFACE_LINK_LIBRARIES>)
endforeach()
