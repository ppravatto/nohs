cmake_minimum_required(VERSION 3.12.0)

project(nohs VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g>"
)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
endif()

add_library(nohs SHARED src/hermite.cpp src/solver.cpp src/optimizer.cpp)
set_target_properties(nohs PROPERTIES
  PUBLIC_HEADER include/nohs.h include/nohs_exceptions.h include/nohs_optimizer.h
)
set_target_properties(nohs PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(nohs PRIVATE include)
target_include_directories(nohs PRIVATE src)
target_link_libraries(nohs ${LIBS})


install(
    TARGETS nohs
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    )

install(
  FILES include/nohs.h include/nohs_exceptions.h include/nohs_optimizer.h
  DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
  )

option(COMPILE_EXAMPLES "Set whether the examples are compiled" OFF)

if(COMPILE_EXAMPLES)

    set(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Set example installation path")

    add_executable(hermite_test examples/hermite_test.cpp)
    target_include_directories(hermite_test PRIVATE include)
    target_link_libraries(hermite_test PRIVATE nohs)
    
    add_executable(quartic examples/quartic.cpp)
    target_include_directories(quartic PRIVATE include)
    target_link_libraries(quartic PRIVATE nohs)

    install(
        TARGETS hermite_test quartic
        RUNTIME DESTINATION "${EXAMPLES_INSTALL_PATH}"
        )
    
endif()