# shellcheck shell=bash disable=SC2154,SC2034,SC2148
find_or_install()
{
  package="$1"
  # This is a bash 3 hack standing in for a bash 4 hash (bash 3 is the lowest common
  # denominator because, for licensing reasons, OS X only has bash 3 by default.)
  # See http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash
  package_executable_array=(
    "gcc:gfortran"
    "cmake:cmake"
    "mpich:mpifort"
    "flex:flex"
    "bison:yacc"
    "m4:m4"
    "_unknown:0"
  )
  for element in "${package_executable_array[@]}" ; do
     KEY="${element%%:*}"
     VALUE="${element##*:}"
     if [[ "$KEY" == "_unknown" ]]; then
       # No recognizeable argument passed so print usage information and exit:
       printf "%s: Package name (%s) not recognized in find_or_install function [exit 40].\n" "$this_script" "$package"
       exit 40
     elif [[ $package == "$KEY" ]]; then
       executable=$VALUE
       break
     fi
  done

  if [[ "$package" == "$executable" ]]; then
    printf "%s: Checking whether %s is in the PATH..." "$this_script" "$executable"
  else
    printf "%s: Checking whether %s executable $executable is in the PATH..." "$this_script" "$package"
  fi
  if type "$executable" >& /dev/null; then
    printf "yes.\n"
    package_in_path=true
    package_version_in_path=$( ("${executable}" --version 2>/dev/null || "${executable}" -V) | head -1)
  else
    printf "no.\n"
    package_in_path=false
  fi

  export default_package_install_path=$(./build.sh -P "$package")

  printf "Checking whether %s is in %s..." "$executable" "$default_package_install_path"
  if [[ -x "$default_package_install_path/bin/$executable" ]]; then
    printf "yes.\n"
    script_installed_package=true
    stack_push script_installed "$package" "$executable"
  else
    script_installed_package=false
    printf "no.\n"
  fi

  minimum_version=$(./build.sh -V "$package")

  if [[ "$package" == "cmake" ]]; then

    # We arrive here only by the explicit, direct call 'find_or_install cmake' inside
    # the build_opencoarrays function. Because there is no possibility of arriving here
    # by recursion (no packages depend on cmake except OpenCoarrays, which gets built
    # after all dependencies have been found or installed), cmake must add itself to
    # the dependency stack if no acceptable cmake is found.

    # Every branch that discovers an acceptable pre-existing installation must set the
    # CMAKE environment variable. Every branch must also manage the dependency stack.

    # If the user specified a CMake binary, use it
    if [[ ! -z "${arg_m:-}" ]]; then

      echo -e "$this_script: Using the $package specified by -m or --with-cmake: ${arg_m}\n"
      export CMAKE="${arg_m}"
      # Halt the recursion
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$script_installed_package" == true ]]; then
      echo -e "$this_script: Using the $package installed by $this_script\n"
      export CMAKE=$default_package_install_path/bin/$executable
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then
      echo -e "$this_script: Checking whether $package in PATH is version < $minimum_version... "

      if ! ./check_version.sh "$package" "$(./build.sh -V "$package")"; then
        printf "yes.\n"
        # Here we place $package on the dependency stack to trigger the build of the above file:
        stack_push dependency_pkg "$package" "none"
        stack_push dependency_exe "$package" "none"
        stack_push dependency_path "$(./build.sh -P cmake)" "none"

      else
        printf "no.\n"
        echo -e "$this_script: Using the $executable found in the PATH.\n"
        export CMAKE=$executable
        # Prevent recursion
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      fi

    else # Build package ($package has no prerequisites)
      stack_push dependency_pkg "$package" "none"
      stack_push dependency_exe "$package" "none"
      stack_push dependency_path "$(./build.sh -P "$package")" "none"
    fi

  elif [[ $package == "mpich" ]]; then

    # We arrive here only by the explicit, direct call 'find_or_install mpich' inside
    # the build_opencoarrays function. Because there is no possibility of arriving here
    # by recursion (no packages depend on mpich except OpenCoarrays, which gets built
    # after all dependencies have been found or installed), mpich must add itself to
    # the dependency stack if no acceptable mpich is found.

    # Every branch that discovers an acceptable pre-existing installation must set the
    # MPIFC, MPICC, and MPICXX environment variables. Every branch must also manage the
    # dependency stack.

    # If the user specified a Fortran compiler, verify that mpifort wraps the specified compiler
    if [[ ! -z "${arg_M:-}" ]]; then

      echo -e "$this_script: Using the $package specified by -M or --with-mpi: ${arg_M}\n"
      export MPIFC="${arg_M}"/bin/mpifort
      export MPICC="${arg_M}"/bin/mpicc
      export MPICXX="${arg_M}"/bin/mpicxx
      # Halt the recursion
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$script_installed_package" == true ]]; then

      echo -e "$this_script: Using the $package installed by $this_script\n"
      export MPIFC=$default_package_install_path/bin/mpifort
      export MPICC=$default_package_install_path/bin/mpicc
      export MPICXX=$default_package_install_path/bin/mpicxx
      # Halt the recursion
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then

      echo -e "$this_script: Checking whether $executable in PATH wraps gfortran... "
      mpifort_version_header=$(mpifort --version | head -1)
      first_three_characters=$(echo "$mpifort_version_header" | cut -c1-3)
      if [[ "$first_three_characters" != "GNU" ]]; then
        printf "no.\n"
        # Trigger 'find_or_install gcc' and subsequent build of $package
        stack_push dependency_pkg "none" "$package" "gcc"
        stack_push dependency_exe "none" "$executable" "gfortran"
        stack_push dependency_path "none" "$(./build.sh -P "$package")" "$(./build.sh -P gcc)"
      else
        printf "yes.\n"

        if [[ ! -z "${arg_f:-}" ]]; then

          info "-f (or --with-fortran) argument detected with value ${arg_f}"
          printf "yes.\n %s: Using the specified %s.\n" "$this_script" "$executable"
          export MPIFC=mpifort
          export MPICC=mpicc
          export MPICXX=mpicxx

          # Halt the recursion
          stack_push dependency_pkg "none"
          stack_push dependency_exe "none"
          stack_push dependency_path "none"

        else

          export minimum_compiler_version=$(./build.sh -V gcc)
          info "$this_script: Checking whether $executable in PATH wraps gfortran version >= $minimum_compiler_version... "
          $executable acceptable_compiler.f90 -o acceptable_compiler || true;
          $executable print_true.f90 -o print_true || true;
          if [[ -f ./acceptable_compiler && -f ./print_true ]]; then
            acceptable=$(./acceptable_compiler $minimum_compiler_version)
            is_true=$(./print_true)
            rm acceptable_compiler print_true
          else
            acceptable=false
          fi
          if [[ "$acceptable" == "${is_true:-}" ]]; then
            printf "yes.\n %s: Using the $executable found in the PATH.\n" "$this_script"
            export MPIFC=mpifort
            export MPICC=mpicc
            export MPICXX=mpicxx

            # Halt the recursion
            stack_push dependency_pkg "none"
            stack_push dependency_exe "none"
            stack_push dependency_path "none"
          else
            printf "no\n"
            # Trigger 'find_or_install gcc' and subsequent build of $package
            stack_push dependency_pkg "none" "$package" "gcc"
            stack_push dependency_exe "none" "$executable" "gfortran"
            stack_push dependency_path "none" "$(./build.sh -P "$package")" "$(./build.sh -P gcc)"
          fi
        fi
      fi

    else # $package not in PATH and not yet installed by this script
      # Trigger 'find_or_install gcc' and subsequent build of $package
      stack_push dependency_pkg  "none" "$package" "gcc"
      stack_push dependency_exe  "none" "$executable" "gfortran"
      stack_push dependency_path "none" "$(./build.sh -P "$package")" "$(./build.sh -P gcc)"
    fi

    # Check consistency of MPIFC, if set, and user-specified Fortran compiler
    if [[ ! -z ${MPIFC:-} && ! -z "${arg_f:-}" ]]; then
      MPIFC_wraps=$("${MPIFC}" --version)
      compiler=$(${arg_f} --version)
      if [[ "${MPIFC_wraps}" != "${compiler}"   ]]; then
        emergency "Specified MPI ${MPIFC_wraps} wraps a compiler other than the specified Fortran compiler ${compiler}"
      fi
    fi

  elif [[ $package == "gcc" ]]; then

    # We arrive here when the 'elif [[ $package == "mpich" ]]' block pushes "gcc" onto the
    # the dependency_pkg stack, resulting in the recursive call 'find_or_install gcc'

    # Every branch that discovers an acceptable pre-existing installation must set the
    # FC, CC, and CXX environment variables. Every branch must also manage the dependency stack.

    if [[ ! -z "${arg_f:-}" ]]; then

      info "-f (or --with-fortran) argument detected with value ${arg_f}"
      [ -z "${arg_c:-}" ] && emergency "-f (--with-fortran) specifies Fortran compiler; Please also specify C/C++ compilers"
      [ -z "${arg_C:-}" ] && emergency "-f (--with-fortran) specifies Fortran compiler; Please also specify C/C++ compilers"

      export FC="${arg_f}"
      export CC="${arg_c}"
      export CXX="${arg_C}"

      # Remove $package from the dependency stack
      stack_pop dependency_pkg package_done
      stack_pop dependency_exe executable_done
      stack_pop dependency_path package_done_path
      # Halt the recursion and signal that none of $package's prerequisites need to be built
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$script_installed_package" == true ]]; then
      echo -e "$this_script: Using the $package executable $executable installed by $this_script\n"
      export FC=$default_package_install_path/bin/gfortran
      export CC=$default_package_install_path/bin/gcc
      export CXX=$default_package_install_path/bin/g++
      gfortran_lib_paths="$default_package_install_path/lib64/:$default_package_install_path/lib"
      if [[ -z "${LD_LIBRARY_PATH:-}" ]]; then
        echo "$this_script: export LD_LIBRARY_PATH=\"$gfortran_lib_paths\""
                            export LD_LIBRARY_PATH="$gfortran_lib_paths"
      else
        echo "$this_script: export LD_LIBRARY_PATH=\"$gfortran_lib_paths:$LD_LIBRARY_PATH\""
                            export LD_LIBRARY_PATH="$gfortran_lib_paths:$LD_LIBRARY_PATH"
      fi
      # Remove $package from the dependency stack
      stack_pop dependency_pkg package_done
      stack_pop dependency_exe executable_done
      stack_pop dependency_path package_done_path
      # Put $package onto the script_installed log
      stack_push script_installed package_done
      stack_push script_installed executable_done
      stack_push script_installed package_done_path
      # Halt the recursion and signal that none of $package's prerequisites need to be built
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then
      info "$this_script: Checking whether $executable in PATH is version $minimum_version or later..."
      $executable -o acceptable_compiler acceptable_compiler.f90 || true;
      $executable -o print_true print_true.f90 || true;
      if [[ -f ./acceptable_compiler && -f ./print_true ]]; then
        is_true=$(./print_true)
        info "Executing './acceptable_compiler $minimum_version'"
        acceptable=$(./acceptable_compiler $minimum_version)
        rm acceptable_compiler print_true
      else
        acceptable=false
      fi

      if [[ "$acceptable" == "${is_true:-}" ]]; then
        printf "yes.\n"
        echo -e "$this_script: Using the $executable found in the PATH.\n"
        export FC=gfortran
        export CC=gcc
        export CXX=g++
        # Remove $package from the dependency stack
        stack_pop dependency_pkg package_done
        stack_pop dependency_exe executable_done
        stack_pop dependency_path package_done_path
        # Put $package onto the script_installed log
        stack_push script_installed package_done
        stack_push script_installed executable_done
        stack_push script_installed package_done_path
        # Halt the recursion and signal that none of $package's prerequisites need to be built
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      else
        printf "no.\n"
        # Trigger 'find_or_install flex' and subsequent build of $package
        stack_push dependency_pkg "flex"
        stack_push dependency_exe "flex"
        stack_push dependency_path "$(./build.sh -P flex)"
      fi

    else # $package is not in PATH and not yet installed by this script
      # Trigger 'find_or_install flex' and subsequent build of $package
      stack_push dependency_pkg "flex"
      stack_push dependency_exe "flex"
      stack_push dependency_path "$(./build.sh -P flex)"
    fi

  elif [[ $package == "flex" ]]; then

    # We arrive here only if the 'elif [[ $package == "gcc" ]]' block has pushed "flex"
    # onto the dependency_pkg stack, resulting in the recursive call 'find_or_install flex'.
    # flex therefore does not need to add itself to the stack.

    # Every branch that discovers an acceptable pre-existing installation must set the
    # FLEX environment variable. Every branch must also manage the dependency stack.

    if [[ "$script_installed_package" == true ]]; then
      echo -e "$this_script: Using the $executable installed by $this_script\n"
      export FLEX=$default_package_install_path/bin/$executable
      # Remove flex from the dependency stack
      stack_pop dependency_pkg package_done
      stack_pop dependency_exe executable_done
      stack_pop dependency_path package_done_path
      # Put $package onto the script_installed log
      stack_push script_installed package_done
      stack_push script_installed executable_done
      stack_push script_installed package_done_path
      # Halt the recursion and signal that no prerequisites need to be built
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then

      echo -e "$this_script: Checking whether $package in PATH is version < $minimum_version... "
      if ! ./check_version.sh "$package" "$(./build.sh -V "$package")"; then
        printf "yes\n"

        export FLEX="$default_package_install_path/bin/$executable"
        # Trigger 'find_or_install bison' and subsequent build of $package
        stack_push dependency_pkg "bison"
        stack_push dependency_exe "yacc"
        stack_push dependency_path "$(./build.sh -P bison)"
      else
        printf "no.\n"
        echo -e "$this_script: Using the $executable found in the PATH.\n"
        export FLEX=$executable
        # Remove $package from the dependency stack
        stack_pop dependency_pkg package_done
        stack_pop dependency_exe executable_done
        stack_pop dependency_path package_done_path
        # Put $package onto the script_installed log
        stack_push script_installed package_done
        stack_push script_installed executable_done
        stack_push script_installed package_done_path
        # Halt the recursion and signal that none of $package's prerequisites need to be built
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      fi

    else  # $package is not in the PATH and not yet installed by $this_script
      # Trigger 'find_or_install bison' and subsequent build of $package
      stack_push dependency_pkg "bison"
      stack_push dependency_exe "yacc"
      stack_push dependency_path "$(./build.sh -P bison)"
    fi

  elif [[ $package == "bison" ]]; then

    # We arrive when the 'elif [[ $package == "flex" ]]' block pushes "bison" onto the
    # the dependency_pkg stack, resulting in the recursive call 'find_or_install bison'

    # Every branch that discovers an acceptable pre-existing installation must set the
    # YACC environment variable. Every branch must also manage the dependency stack.

    if [[ "$script_installed_package" == true ]]; then
      echo -e "$this_script: Using the $package executable $executable installed by $this_script\n"
      export YACC=$default_package_install_path/bin/yacc
      # Remove bison from the dependency stack
      stack_pop dependency_pkg package_done
      stack_pop dependency_exe executable_done
      stack_pop dependency_path package_done_path
      # Put $package onto the script_installed log
      stack_push script_installed package_done
      stack_push script_installed executable_done
      stack_push script_installed package_done_path
      # Halt the recursion and signal that there are no prerequisites to build
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then
      echo -e "$this_script: Checking whether $package executable $executable in PATH is version < $minimum_version... "
      if ! ./check_version.sh "$package" "$(./build.sh -V "$package")"; then
        printf "yes.\n"
        export YACC="$default_package_install_path/bin/$executable"
        # Trigger 'find_or_install m4' and subsequent build of $package
        stack_push dependency_pkg "m4"
        stack_push dependency_exe "m4"
        stack_push dependency_path "$(./build.sh -P m4)"
      else
        printf "no.\n"
        echo -e "$this_script: Using the $package executable $executable found in the PATH.\n"
        YACC=yacc
        # Remove bison from the dependency stack
        stack_pop dependency_pkg package_done
        stack_pop dependency_exe executable_done
        stack_pop dependency_path package_done_path
        # Put $package onto the script_installed log
        stack_push script_installed package_done
        stack_push script_installed executable_done
        stack_push script_installed package_done_path
        # Halt the recursion and signal that there are no prerequisites to build
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      fi

    else # $package not in PATH and not yet installed by this script
      # Trigger 'find_or_install m4' and subsequent build of $package
      stack_push dependency_pkg "m4"
      stack_push dependency_exe "m4"
      stack_push dependency_path "$(./build.sh -P m4)"
    fi

  elif [[ $package == "m4" ]]; then

    # We arrive when the 'elif [[ $package == "bison" ]]' block pushes "m4" onto the
    # the dependency_pkg stack, resulting in the recursive call 'find_or_install m4'

    # Every branch that discovers an acceptable pre-existing installation must set the
    # M4 environment variable. Every branch must also manage the dependency stack.

    if [[ "$script_installed_package" == true ]]; then
      echo -e "$this_script: Using the $package executable $executable installed by $this_script\n"
      export M4=$default_package_install_path/bin/m4
      # Remove m4 from the dependency stack
      stack_pop dependency_pkg package_done
      stack_pop dependency_exe executable_done
      stack_pop dependency_path package_done_path
      # Put $package onto the script_installed log
      stack_push script_installed package_done
      stack_push script_installed executable_done
      stack_push script_installed package_done_path
      # Halt the recursion and signal that there are no prerequisites to build
      stack_push dependency_pkg "none"
      stack_push dependency_exe "none"
      stack_push dependency_path "none"

    elif [[ "$package_in_path" == "true" ]]; then
      echo -e "$this_script: Checking whether $package executable $executable in PATH is version < $minimum_version... "
      if ! ./check_version.sh "$package" "$(./build.sh -V "$package")"; then
        printf "yes.\n"
        export M4="$default_package_install_path/bin/m4"
        # Halt the recursion and signal that there are no prerequisites to build
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      else
        printf "no.\n"
        echo -e "$this_script: Using the $package executable $executable found in the PATH.\n"
        M4=m4
        # Remove m4 from the dependency stack
        stack_pop dependency_pkg package_done
        stack_pop dependency_exe executable_done
        stack_pop dependency_path package_done_path
        # Put $package onto the script_installed log
        stack_push script_installed package_done
        stack_push script_installed executable_done
        stack_push script_installed package_done_path
        # Halt the recursion and signal that there are no prerequisites to build
        stack_push dependency_pkg "none"
        stack_push dependency_exe "none"
        stack_push dependency_path "none"
      fi

    else # $package not in PATH and not yet installed by this script
      # Halt the recursion and signal that there are no prerequisites to build
      export M4="$default_package_install_path/bin/m4"
      stack_push dependency_pkg  "none"
      stack_push dependency_exe  "none"
      stack_push dependency_path "none"
    fi

  else
    if [[ -z "${package:-}" ]]; then
      echo -e "$this_script: empty package name passed to find_or_install function. [exit 50]\n"
      exit 50
    else
      echo -e "$this_script: unknown package name ($package) passed to find_or_install function. [exit 55]\n"
      exit 55
    fi
  fi

  echo "$this_script: Updated dependency stack (top to bottom = left to right):"
  stack_print dependency_pkg

  stack_size dependency_pkg num_stacked
  (( num_dependencies=num_stacked-1 )) || true

  if [[ $num_dependencies -lt 0 ]]; then
    emergency "The procedure named in the external call to find_or_install is not on the dependency stack. [exit 60]\n"

  elif [[ $num_dependencies -gt 0 ]]; then
    stack_pop  dependency_pkg  prerequisite_pkg
    stack_pop  dependency_exe  prerequisite_exe
    stack_pop  dependency_path prerequisite_path

    if [[ $prerequisite_pkg != "none" ]]; then
      stack_push dependency_pkg  "$prerequisite_pkg"
      stack_push dependency_exe  "$prerequisite_exe"
      stack_push dependency_path "$prerequisite_path"
      echo -e "$this_script: Building $package from source requires $prerequisite_pkg.\n"
      find_or_install "$prerequisite_pkg"
    fi
  fi

  echo "$this_script: Remaining $package dependency stack (top to bottom = left to right):"
  stack_print dependency_pkg

  stack_pop dependency_pkg package
  stack_pop dependency_exe executable
  stack_pop dependency_path default_package_install_path

  if [[ $package != "none" ]]; then

    export default_package_version=$(./build.sh -V "${package}")
    export package_version_to_install=${arg_I:-${default_package_version}}

    prefix_root="${prefix_root:-${default_package_install_path%${package}/${package_version_to_install}}}"

    if [[ ${package} == "gcc" ]]; then
      export package_directory_name="gnu"
    else
      export package_directory_name=${package}
    fi

    package_install_prefix="${prefix_root}/${package_directory_name}/${package_version_to_install}"

    if [[ "$package" == "$executable" ]]; then
      echo "$this_script: Ready to install $executable in $package_install_prefix"
    else
      echo "$this_script: Ready to install $package executable $executable in $package_install_prefix"
    fi


    if [[ "${arg_y}" == "${__flag_present}" ]]; then
      info "-y or --yes-to-all flag present. Proceeding with non-interactive build."
    else
      echo -e "$this_script: Ok to download (if necessary), build, and install $package from source? (Y/n) "
      read -r proceed_with_build

      if [[ "$proceed_with_build" == "n" || "$proceed_with_build" == "no" ]]; then
        printf "n\n"
        echo -e "$this_script: OpenCoarrays installation requires $package. Aborting. [exit 70]\n"
        exit 70
      else # permission granted to build
        printf "Y\n"
      fi
    fi

    # On OS X, CMake must be built with Apple LLVM gcc, which XCode command-line tools puts in /usr/bin
    if [[ $(uname) == "Darwin" && $package == "cmake"  ]]; then
      if [[ -x "/usr/bin/gcc" ]]; then
        CC=/usr/bin/gcc
      else
        echo -e "$this_script: OS X detected.  Please install XCode command-line tools and \n"
        echo -e "$this_script: ensure that /usr/bin/gcc exists and is executable. Aborting. [exit 75]\n"
        exit 75
      fi
    # Otherwise, if no CC has been defined yet, use the gcc in the user's PATH
    elif [[ -z "${CC:-}" ]]; then
      CC=gcc
    fi

    # On OS X, CMake must be built with Apple LLVM g++, which XCode command-line tools puts in /usr/bin
    if [[ $(uname) == "Darwin" && $package == "cmake"  ]]; then
      if [[ -x "/usr/bin/g++" ]]; then
        CXX=/usr/bin/g++
      else
        echo -e "$this_script: OS X detected.  Please install XCode command-line tools \n"
        echo -e "$this_script: and ensure that /usr/bin/g++ exists and is executable. Aborting. [exit 76]\n"
        exit 76
      fi
    # Otherwise, if no CXX has been defined yet, use the g++ in the user's PATH
    elif [[ -z "${CXX:-}" ]]; then
      CXX=g++
    fi

    # If no FC has been defined yet, use the gfortran in the user's PATH
    if [[ -z "${FC:-}" ]]; then
      FC=gfortran
    fi

    if [[ "${arg_y}" == "${__flag_present}" ]]; then
      yes_to_all="-y"
    fi

    if [[ "${arg_Z}" == "${__flag_present}" ]]; then
      bootstrap="-Z"
    fi

    echo -e "$this_script: Downloading, building, and installing $package \n"
    echo "$this_script: Build command: FC=$FC CC=$CC CXX=$CXX ./build.sh -p $package -i $package_install_prefix -j $num_threads ${yes_to_all:-} ${bootstrap:-}"

    FC="$FC" CC="$CC" CXX="$CXX" ./build.sh -p "$package" -i "$package_install_prefix" -j "$num_threads" "${yes_to_all:-}" "${bootstrap:-}"

    if [[ ! -x "$package_install_prefix/bin/$executable" ]]; then

      echo -e "$this_script: Installation unsuccessful. "
      emergency "$executable is not in the path ($package_install_prefix/bin) or the user lacks executable permission."

    else
      echo -e "$this_script: Installation successful.\n"
      if [[ "$package" == "$executable" ]]; then
        echo -e "$this_script: $executable is in $package_install_prefix/bin \n"
      else
        echo -e "$this_script: $package executable $executable is in $package_install_prefix/bin \n"
      fi

      if [[ $package == "cmake" ]]; then
        echo "$this_script: export CMAKE=$package_install_prefix/bin/$executable"
                            export CMAKE="$package_install_prefix/bin/$executable"
      elif [[ $package == "bison" ]]; then
        echo "$this_script: export YACC=$package_install_prefix/bin/$executable"
                            export YACC="$package_install_prefix/bin/$executable"
      elif [[ $package == "flex" ]]; then
        echo "$this_script: export FLEX=$package_install_prefix/bin/$executable"
                            export FLEX="$package_install_prefix/bin/$executable"
      elif [[ $package == "m4" ]]; then
        echo "$this_script: export M4=$package_install_prefix/bin/$executable"
                            export M4="$package_install_prefix/bin/$executable"
      elif [[ $package == "gcc" ]]; then
        echo "$this_script: export FC=$package_install_prefix/bin/gfortran"
                            export FC="$package_install_prefix/bin/gfortran"
        echo "$this_script: export CC=$package_install_prefix/bin/gcc"
                            export CC="$package_install_prefix/bin/gcc"
        echo "$this_script: export CXX=$package_install_prefix/bin/g++"
                            export CXX="$package_install_prefix/bin/g++"
        gfortran_lib_paths="$package_install_prefix/lib64/:$package_install_prefix/lib"
        if [[ -z "${LD_LIBRARY_PATH:-}" ]]; then
          export LD_LIBRARY_PATH="$gfortran_lib_paths"
          else
            export LD_LIBRARY_PATH="$gfortran_lib_paths:$LD_LIBRARY_PATH"
          fi
        elif [[ $package == "mpich" ]]; then
          echo "$this_script: export MPIFC=$package_install_prefix/bin/mpifort"
                              export MPIFC="$package_install_prefix/bin/mpifort"
          echo "$this_script: export MPICC= $package_install_prefix/bin/mpicc"
                              export MPICC="$package_install_prefix/bin/mpicc"
          echo "$this_script: export MPICXX=$package_install_prefix/bin/mpicxx"
                              export MPICXX="$package_install_prefix/bin/mpicxx"
        else
          echo -e "$this_script: WARNING: $package executable $executable installed correctly but the \n"
          echo -e "$this_script:          corresponding environment variable(s) have not been set. This \n"
          echo -e "$this_script:          could prevent a complete build of OpenCoarrays. Please report this\n"
          echo -e "$this_script:          issue at https://github.com/sourceryinstitute/opencoarrays/issues\n"
        fi
        if [[ -z "${PATH:-}" ]]; then
          export PATH="$package_install_prefix/bin"
        else
          export PATH="$package_install_prefix/bin:$PATH"
        fi
      fi
    fi # End 'if [[ ! -x "$package_install_prefix/bin/$executable" ]]; then'
}
