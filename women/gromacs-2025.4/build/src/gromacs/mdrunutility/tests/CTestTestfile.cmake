# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/mdrunutility/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/mdrunutility/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[MdrunUtilityUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/mdrunutility-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/MdrunUtilityUnitTests.xml")
set_tests_properties([=[MdrunUtilityUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest" PROCESSORS "1" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;373;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/mdrunutility/tests/CMakeLists.txt;38;gmx_add_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/mdrunutility/tests/CMakeLists.txt;0;")
add_test([=[MdrunUtilityMpiUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/mdrunutility-mpi-test" "-ntmpi" "4" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/MdrunUtilityMpiUnitTests.xml")
set_tests_properties([=[MdrunUtilityMpiUnitTests]=] PROPERTIES  LABELS "GTest;UnitTest;MpiTest" PROCESSORS "4" TIMEOUT "30" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;385;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/mdrunutility/tests/CMakeLists.txt;45;gmx_add_mpi_unit_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/mdrunutility/tests/CMakeLists.txt;0;")
