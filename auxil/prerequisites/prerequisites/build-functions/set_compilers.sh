# Set CC, CXX, and FC
# shellcheck disable=SC2154
set_compilers()
{
  # If FC is already set, use the set value.
  # Else if FC is empty, use the value specified via the -f or --with-fortran command-line arguments.
  # Else use the default specifed in ../build.sh-usage.
  FC=${FC:-${arg_f}}
  info "Fortran compiler: ${FC}"

  if [[ $(uname) == "Darwin" && "${package_to_build}" == "cmake"  ]]; then
    if [[ -x "/usr/bin/gcc" ]]; then
      [ ! -z "${CC:-}" ] && info "Overriding CC: cmake build requires Apple LLVM gcc, which XCode command-line tools puts in /usr/bin"
      CC=/usr/bin/gcc
    else
      info      "OS X detected.  Please install XCode command-line tools and "
      emergency "ensure that /usr/bin/gcc exists and is executable. Aborting."
    fi
  else
    # If CC is already set, use the set value.
    # Else if CC is empty, use the value specified via the -c or --with-c command-line arguments.
    # Else use the default specifed in ../build.sh-usage.
    CC=${CC:-${arg_c}}
  fi
  info "C compiler: ${CC}"

  if [[ $(uname) == "Darwin" && ${package_to_build} == "cmake"  ]]; then
    if [[ -x "/usr/bin/g++" ]]; then
      [ ! -z "${CXX:-}" ] && info "Overriding CXX: cmake build requires Apple LLVM g++, which XCode command-line tools puts in /usr/bin"
      CXX=/usr/bin/g++
    else
      info      "OS X detected.  Please install XCode command-line tools and "
      emergency "ensure that /usr/bin/g++ exists and is executable. Aborting."
    fi
  else
    # If CXX is already set, use the set value.
    # Else if CXX is empty, use the value specified via the -C or --with-cxx command-line arguments.
    # Else use the default specifed in ../build.sh-usage.
    CXX=${CXX:-${arg_C}}
  fi
  info "C++ compiler: ${CXX}"
}
