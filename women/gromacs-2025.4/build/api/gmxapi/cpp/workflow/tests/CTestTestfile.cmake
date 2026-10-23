# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/workflow/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/gmxapi/cpp/workflow/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[GmxapiInternalInterfaceTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/workflow-details-test" "-ntomp" "2" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/GmxapiInternalInterfaceTests.xml")
set_tests_properties([=[GmxapiInternalInterfaceTests]=] PROPERTIES  LABELS "GTest;IntegrationTest;QuickGpuTest" PROCESSORS "2" TIMEOUT "120" WORKING_DIRECTORY "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/gmxapi/cpp/workflow/tests" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/workflow/tests/CMakeLists.txt;51;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/workflow/tests/CMakeLists.txt;0;")
