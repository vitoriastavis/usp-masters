# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/external/muparser")
  file(MAKE_DIRECTORY "/home/gbiuser/Documents/vitoria/gromacs-2025.4/src/external/muparser")
endif()
file(MAKE_DIRECTORY
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-build"
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix"
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/tmp"
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/src/muparser-populate-stamp"
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/src"
  "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/src/muparser-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/src/muparser-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/gbiuser/Documents/vitoria/gromacs-2025.4/build/_deps/muparser-subbuild/muparser-populate-prefix/src/muparser-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
