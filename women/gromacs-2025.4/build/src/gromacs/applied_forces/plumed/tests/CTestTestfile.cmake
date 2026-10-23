# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/applied_forces/plumed/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/applied_forces/plumed/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[PlumedAppliedForcesUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/plumed_applied_forces-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/PlumedAppliedForcesUnitTests.xml")
set_tests_properties([=[PlumedAppliedForcesUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;373;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/applied_forces/plumed/tests/CMakeLists.txt;34;gmx_add_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/applied_forces/plumed/tests/CMakeLists.txt;0;")
add_test([=[PlumedMDTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/plumed_md-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/PlumedMDTests.xml")
set_tests_properties([=[PlumedMDTests]=] PROPERTIES  LABELS "GTest;IntegrationTest;QuickGpuTest" PROCESSORS "1" TIMEOUT "120" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/applied_forces/plumed/tests/CMakeLists.txt;57;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/applied_forces/plumed/tests/CMakeLists.txt;0;")
