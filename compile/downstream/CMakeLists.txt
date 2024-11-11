cmake_minimum_required(VERSION 3.12)

project(
  JetReconstruction_downstream
  LANGUAGES C
  VERSION 0.1.0)

find_package(JetReconstruction 0.4.3 REQUIRED)

add_executable(jetreconstruction_test jetreconstruction_test.c)
target_link_libraries(jetreconstruction_test
                      PUBLIC JetReconstruction::JetReconstruction)

enable_testing()
add_test(NAME JetReconstructionTest COMMAND jetreconstruction_test)
