# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/utility/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/utility/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[UtilityUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/utility-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/UtilityUnitTests.xml")
set_tests_properties([=[UtilityUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;373;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/utility/tests/CMakeLists.txt;34;gmx_add_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/utility/tests/CMakeLists.txt;0;")
add_test([=[UtilityMpiUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/utility-mpi-test" "-ntmpi" "2" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/UtilityMpiUnitTests.xml")
set_tests_properties([=[UtilityMpiUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest;MpiTest" PROCESSORS "2" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;385;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/utility/tests/CMakeLists.txt;66;gmx_add_mpi_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/utility/tests/CMakeLists.txt;0;")
