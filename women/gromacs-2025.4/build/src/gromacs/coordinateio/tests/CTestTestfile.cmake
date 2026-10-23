# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/coordinateio/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/coordinateio/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[CoordinateIOTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/coordinateio-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/CoordinateIOTests.xml")
set_tests_properties([=[CoordinateIOTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;373;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/coordinateio/tests/CMakeLists.txt;34;gmx_add_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/coordinateio/tests/CMakeLists.txt;0;")
