# shellcheck shell=bash disable=SC2154,SC2129,SC2148
report_results()
{
  fully_qualified_FC="$(type -P "${FC}")"
  if [[ ${fully_qualified_FC} != *gfortran* ]]; then
    emergency "report_results.sh: non-gfortran compiler: \${fully_qualified_FC}=${fully_qualified_FC}"
  fi
  # Set path_to_FC fully-qualified gfortran location
  compiler_install_root="${fully_qualified_FC%bin/gfortran*}"

  fully_qualified_MPIFC="$(type -P "${MPIFC}")"
  mpi_install_root="${fully_qualified_MPIFC%bin/mpifort*}"

  fully_qualified_CMAKE="$(type -P "${CMAKE}")"
  cmake_install_path="${fully_qualified_CMAKE%/cmake*}"

  # Report installation success or failure and record locations for software stack:
  if [[ -x "${install_path%/}/bin/caf" && -x "${install_path%/}/bin/cafrun" ]]; then

    # Installation succeeded
    info "$this_script: Done."
    info ""
    info "*** The OpenCoarrays compiler wrapper (caf) and program  ***"
    info "*** launcher (cafrun) are in the following directory:    ***"
    info ""
    info "${install_path%/}/bin."
    info ""
    if [[ -f setup.sh ]]; then
      ${SUDO:-} rm setup.sh
    fi
    if [[ -f setup.csh ]]; then
      ${SUDO:-} rm setup.csh
    fi
    # Prepend the OpenCoarrays license to the setup.sh script:
    while IFS='' read -r line || [[ -n "${line}" ]]; do
        echo "# ${line}" >> setup.sh
    done < "${OPENCOARRAYS_SRC_DIR}/LICENSE"
    while IFS='' read -r line || [[ -n "${line}" ]]; do
        echo "# ${line}" >> setup.csh
    done < "${OPENCOARRAYS_SRC_DIR}/LICENSE"
    echo "#                                                                               " | tee -a setup.csh setup.sh
    echo "# Execute this script via the following command:                                " | tee -a setup.csh setup.sh
    echo "# source ${install_path%/}/setup.sh                                             " | tee -a setup.csh setup.sh
    echo "                                                                                " | tee -a setup.csh setup.sh
    if [[ -x "${cmake_install_path}/cmake" ]]; then
      echo "# Prepend the CMake path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                       " >> setup.sh
      echo "  export PATH=\"${cmake_install_path%/}\"                            " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  if ! [[ \"\${PATH}\" =~ \"${cmake_install_path}\" ]] ; then        " >> setup.sh
      echo "    export PATH=\"${cmake_install_path%/}\":\${PATH}                 " >> setup.sh
      echo "  fi                                                                 " >> setup.sh
      echo "fi                                                                   " >> setup.sh
      echo "set path = (\"${cmake_install_path%/}\" \$path)                      " >> setup.csh
    fi
    if [[ -x "${fully_qualified_FC}" ]]; then
      echo "# Prepend the compiler path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                                " >> setup.sh
      echo "  export PATH=\"${compiler_install_root%/}/bin\"                              " >> setup.sh
      echo "else                                                                          " >> setup.sh
      echo "  export PATH=\"${compiler_install_root%/}/bin:\${PATH}\"                     " >> setup.sh
      echo "fi                                                                            " >> setup.sh
      echo "set path = (\"${compiler_install_root%/}\"/bin \$path)                        " >> setup.csh
    fi
    LD_LIB_P_VAR=LD_LIBRARY_PATH
    if [[ "${OSTYPE:-}" =~ [Dd]arwin ]]; then
      LD_LIB_P_VAR=DYLD_LIBRARY_PATH
    fi
    if [[ -d "${compiler_install_root%/}/lib" || -d "${compiler_install_root%/}/lib64" ]]; then
      echo "# Prepend the compiler library paths to the ${LD_LIB_P_VAR} environment variable:" | tee -a setup.sh setup.csh
      compiler_lib_paths="${compiler_install_root%/}/lib64/:${compiler_install_root%/}/lib"
      echo "if [[ -z \"\${${LD_LIB_P_VAR}}\" ]]; then                                    " >> setup.sh
      echo "  export ${LD_LIB_P_VAR}=\"${compiler_lib_paths%/}\"                          " >> setup.sh
      echo "else                                                                          " >> setup.sh
      echo "  export ${LD_LIB_P_VAR}=\"${compiler_lib_paths%/}:\${${LD_LIB_P_VAR}}\"     " >> setup.sh
      echo "fi                                                                            " >> setup.sh
      echo "setenv LD_LIBRARY_PATH \"${compiler_lib_paths%/}:\${LD_LIBRARY_PATH}\"        " >> setup.csh
    fi
    echo "                                                                                " >> setup.sh
    if [[ -x "${mpi_install_root}/bin/mpifort" ]]; then
      echo "# Prepend the MPI path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                       " >> setup.sh
      echo "  export PATH=\"${mpi_install_root%/}/bin\"                          " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  export PATH=\"${mpi_install_root%/}/bin\":\${PATH}                 " >> setup.sh
      echo "fi                                                                   " >> setup.sh
      echo "set path = (\"${mpi_install_root%/}/bin\" \$path)                    " >> setup.csh
    fi
    # In all likelihood, the following paths are only needed if OpenCoarrays built them,
    # In by far the most common such use case, they would have been built in a recursive
    # build of all the OpenCoarrays dependency tree (rather than built indvidually via
    # ./install --package) so we only need check the default location in which OpenCoarrays
    # would install them.  If they are not there, then it is very likely the case that the
    # the system versions of these packages are present and in the user's path or that the
    # user doesn't need them at all (e.g. there was no need to build gfortran from source).
    flex_install_path=$("${build_script}" -P flex)
    if [[ -x "${flex_install_path}/bin/flex" ]]; then
      echo "# Prepend the flex path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                         " >> setup.sh
      echo "  export PATH=\"${flex_install_path}/bin\"                             " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  export PATH=\"${flex_install_path}/bin\":\${PATH}                      " >> setup.sh
      echo "set path = (\"$flex_install_path\"/bin \$path)                      " >> setup.csh
      echo "fi                                                                   " >> setup.sh
    fi
    bison_install_path=$("${build_script}" -P bison)
    if [[ -x "${bison_install_path}/bin/yacc" ]]; then
      echo "# Prepend the bison path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                         " >> setup.sh
      echo "  export PATH=\"${bison_install_path}/bin\"                            " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  export PATH=\"${bison_install_path}/bin\":\${PATH}                     " >> setup.sh
      echo "fi                                                                   " >> setup.sh
      echo "set path = (\"$bison_install_path\"/bin \$path)                      " >> setup.csh
    fi
    m4_install_path=$("${build_script}" -P m4)
    if [[ -x "${m4_install_path}/bin/m4" ]]; then
      echo "# Prepend the m4 path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                         " >> setup.sh
      echo "  export PATH=\"${m4_install_path}/bin\"                               " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  export PATH=\"${m4_install_path}/bin\":\${PATH}                        " >> setup.sh
      echo "fi                                                                   " >> setup.sh
      echo "set path = (\"$m4_install_path\"/bin \$path)                      " >> setup.csh
    fi
    opencoarrays_install_path="${install_path}"
    if [[ -x "${opencoarrays_install_path}/bin/caf" ]]; then
      echo "# Prepend the OpenCoarrays path to the PATH environment variable:" | tee -a setup.sh setup.csh
      echo "if [[ -z \"\${PATH}\" ]]; then                                         " >> setup.sh
      echo "  export PATH=\"${opencoarrays_install_path%/}/bin\"                     " >> setup.sh
      echo "else                                                                 " >> setup.sh
      echo "  export PATH=\"${opencoarrays_install_path%/}/bin\":\${PATH}              " >> setup.sh
      echo "fi                                                                   " >> setup.sh
      echo "set path = (\"${opencoarrays_install_path%/}\"/bin \$path)                      " >> setup.csh
    fi
    if ${SUDO:-} mv setup.sh "${opencoarrays_install_path}"; then
       setup_sh_location=${opencoarrays_install_path}
    else
       setup_sh_location=${PWD}
    fi
    if ${SUDO:-} mv setup.csh "${opencoarrays_install_path}"; then
       setup_csh_location=${opencoarrays_install_path}
    else
       setup_csh_location=${PWD}
    fi
    echo "*** To set up your environment for using caf and cafrun, please   ***"
    echo "*** source the installed setup.sh file in a bash shell or source  ***"
    echo "*** setup.csh in a C-shell or add one of the following statements ***"
    echo "*** to your login file:                                           ***"
    echo ""
    echo " source ${setup_sh_location%/}/setup.sh"
    echo " source ${setup_csh_location%/}/setup.csh"
    echo ""
    echo "*** Installation complete.                                        ***"

  else # Installation failed

    echo "Something went wrong. Either the user lacks executable permissions for the"
    echo "OpenCoarrays compiler wrapper (caf), program launcher (cafrun), or prerequisite"
    echo "package installer (build), or these programs are not in the following, expected"
    echo "location:"
    echo "${install_path}/bin."
    echo "Please review the following file for more information:"
    echo "${install_path}/${installation_record}"
    echo "and submit an bug report at https://github.com/sourceryinstitute/opencoarrays/issues/new"
    echo "[exit 100]"
    exit 100

  fi # Ending check for caf, cafrun, build not in expected path
}
