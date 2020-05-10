file(GLOB_RECURSE sources "*.cpp" "*.h")

add_library(bmdp ${sources})

# find SYSTEM libraries
find_library(ARMA_LIB armadillo)
find_library(GINAC_LIB ginac)

# setup scaffolding for third-party non-cmake libraries
set(IMDP_LIB_PATH "${CMAKE_SOURCE_DIR}/lib/IMDP-synthesis")

add_custom_target(IMDP_LIB_TARGET)
add_custom_command(
  TARGET IMDP_LIB_TARGET
  COMMAND make lib
  WORKING_DIRECTORY "${IMDP_LIB_PATH}"
  COMMENT "Building IMDP library"
  )

# 1.0: create an library target out of the library compilation result
add_library(IMDP_LIB STATIC IMPORTED GLOBAL)
add_dependencies(IMDP_LIB IMDP_LIB_TARGET)

# 2.0: specify where the library is and where to find the headers
set_target_properties(IMDP_LIB
    PROPERTIES
    IMPORTED_LOCATION "${IMDP_LIB_PATH}/lib/libsynth.a"
    INTERFACE_INCLUDE_DIRECTORIES "${IMDP_LIB_PATH}")

# 3.0: finally, link
add_dependencies(bmdp IMDP_LIB)

# link with third-party SYSTEM libraries
target_link_libraries(bmdp 
  PUBLIC common
  PUBLIC ${ARMA_LIB} 
  PUBLIC ${GINAC_LIB}
  PUBLIC nlopt
  PUBLIC IMDP_LIB)

add_custom_command(
 TARGET bmdp POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/phi2.txt" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
 COMMENT "copying CS1.mat => ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/phi2.txt" DESTINATION bin)

add_custom_command(
 TARGET bmdp POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/phi1.txt" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
 COMMENT "copying CS1.mat => ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/phi1.txt" DESTINATION bin)


generic_includes(bmdp)
