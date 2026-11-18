cmake_minimum_required(VERSION 3.10)
project(fourier-signal-lab C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

include_directories(src)

add_executable(fourier-lab
  src/main.c
  src/complex.h
  src/fft.c
  src/fft.h
  src/window.c
  src/window.h
  src/signal.c
  src/signal.h
  src/spectrum.c
  src/spectrum.h
  src/io.c
  src/io.h
)
target_link_libraries(fourier-lab m)

enable_testing()

add_executable(run_tests
  tests/test_main.c
  tests/test_fft.c
  tests/test_window.c
  tests/test_signal.c
  tests/test_spectrum.c
  tests/test_io.c
  src/fft.c
  src/window.c
  src/signal.c
  src/spectrum.c
  src/io.c
)
target_include_directories(run_tests PRIVATE tests)
target_link_libraries(run_tests m)

add_test(NAME unit_tests COMMAND run_tests)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  add_test(NAME numpy_validation
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_numpy.py
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
