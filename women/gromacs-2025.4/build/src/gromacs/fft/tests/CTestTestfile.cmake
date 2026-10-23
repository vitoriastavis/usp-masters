# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/fft/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/fft/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[FFTUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/fft-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/FFTUnitTests.xml")
set_tests_properties([=[FFTUnitTests]=] PROPERTIES  LABELS "GTest;SlowTest;QuickGpuTest" PROCESSORS "1" TIMEOUT "480" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/fft/tests/CMakeLists.txt;41;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/fft/tests/CMakeLists.txt;0;")
