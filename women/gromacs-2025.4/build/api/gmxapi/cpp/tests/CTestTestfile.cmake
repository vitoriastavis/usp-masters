# CMake generated Testfile for 
# Source directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/tests
# Build directory: /home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/gmxapi/cpp/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[GmxapiExternalInterfaceTests]=] "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/bin/gmxapi-test" "-ntomp" "2" "--gtest_output=xml:/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/Testing/Temporary/GmxapiExternalInterfaceTests.xml")
set_tests_properties([=[GmxapiExternalInterfaceTests]=] PROPERTIES  LABELS "GTest;IntegrationTest;QuickGpuTest" PROCESSORS "2" TIMEOUT "120" WORKING_DIRECTORY "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/api/gmxapi/cpp/tests" _BACKTRACE_TRIPLES "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/testutils/TestMacros.cmake;355;add_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/tests/CMakeLists.txt;69;gmx_register_gtest_test;/home/gbiuser/Documents/vitoria/gromacs-2025.4/api/gmxapi/cpp/tests/CMakeLists.txt;0;")
