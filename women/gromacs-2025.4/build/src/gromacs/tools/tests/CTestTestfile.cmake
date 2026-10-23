# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/tools/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/src/gromacs/tools/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[ToolUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/tool-test" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/ToolUnitTests.xml")
set_tests_properties([=[ToolUnitTests]=] PROPERTIES  LABELS "GTest;SlowTest" PROCESSORS "1" TIMEOUT "480" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/tools/tests/CMakeLists.txt;43;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/tools/tests/CMakeLists.txt;0;")
add_test([=[ToolWithLeaksUnitTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/tool-test-with-leaks" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/ToolWithLeaksUnitTests.xml")
set_tests_properties([=[ToolWithLeaksUnitTests]=] PROPERTIES  LABELS "GTest;SlowTest" PROCESSORS "1" TIMEOUT "480" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/tools/tests/CMakeLists.txt;50;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/gromacs/tools/tests/CMakeLists.txt;0;")
