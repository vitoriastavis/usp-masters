# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/util/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/nblib/util/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[NbLibUtilTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/nblib-util-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/NbLibUtilTests.xml")
set_tests_properties([=[NbLibUtilTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/util/tests/CMakeLists.txt;50;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/util/tests/CMakeLists.txt;0;")
