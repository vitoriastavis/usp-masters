# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/listed_forces/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/nblib/listed_forces/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[NbLibListedForcesTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/nblib-listed-forces-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/NbLibListedForcesTests.xml")
set_tests_properties([=[NbLibListedForcesTests]=] PROPERTIES  LABELS "GTest;IntegrationTest" PROCESSORS "1" TIMEOUT "120" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/listed_forces/tests/CMakeLists.txt;64;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/nblib/listed_forces/tests/CMakeLists.txt;0;")
