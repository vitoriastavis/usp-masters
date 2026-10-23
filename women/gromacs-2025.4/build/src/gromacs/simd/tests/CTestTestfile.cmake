# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/simd/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/simd/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[SimdUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/simd-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/SimdUnitTests.xml")
set_tests_properties([=[SimdUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;373;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/simd/tests/CMakeLists.txt;34;gmx_add_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/simd/tests/CMakeLists.txt;0;")
