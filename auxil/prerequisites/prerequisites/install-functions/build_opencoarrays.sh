#shellcheck shell=bash
# shellcheck disable=SC2154
build_opencoarrays()
{
  print_header
  info "Invoking find_or_install mpich"
  find_or_install mpich
  info "Invoking find_or_install cmake"
  find_or_install cmake
  build_path="${build_path}"/opencoarrays/$("${OPENCOARRAYS_SRC_DIR}"/install.sh -V opencoarrays)
  if [[ -d "${build_path}" ]]; then
    rm -rf "${build_path}"
  fi
  mkdir -p "$build_path"
  pushd "$build_path"
  if [[ -z ${MPIFC:-} || -z ${MPICC:-} ]]; then
    emergency "build_opencoarrays.sh: empty \${MPIFC}=${MPIFC:-} or \${MPICC}=${MPICC:-}"
  fi
  MPIFC_show=($("$MPIFC" -show))
  MPICC_show=($("$MPICC" -show))
  if [[ "${MPIFC_show[0]}" != *gfortran* ]]; then
    emergency "build_opencoarrays.sh: MPI doesn't wrap gfortran/gcc: \${MPIFC_show}=${MPIFC_show[*]}"
  fi
  if [[ -z "${OPENCOARRAYS_DEVELOPER:-}" ]]; then
      # We should examine the value too, but CMake has many ways of saying "true"
      WDEVFLAG=-Wno-dev
  else
      WDEVFLAG=-Wdev
  fi
    # Set FC to the MPI implementation's gfortran command with any preceding path but without any subsequent arguments:
  FC="${MPIFC_show[0]}"
  # Set CC to the MPI implementation's gcc command...
  CC="${MPICC_show[0]}"
  # try to find mpiexec
  MPI_BIN_DIR="$(type -P "${MPICC}")"
  MPIEXEC_CANDIDATES=($(find "${MPI_BIN_DIR%/*}" -name 'mpiexec' -o -name 'mpirun' -o -name 'lamexec' -o -name 'srun'))
  if ! ((${#MPIEXEC_CANDIDATES[@]} >= 1)); then
    emergency "Could not find a suitable \`mpiexec\` in directory containing mpi wrappers (${MPICC%/*})"
  else
    MPIEXEC="${MPIEXEC_CANDIDATES[0]}"
  fi

  source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_installation_path.sh"
  set_or_print_installation_path

  info "Configuring OpenCoarrays in ${PWD} with the command:"
  info "CC=\"${CC}\" FC=\"${FC}\" $CMAKE \"${OPENCOARRAYS_SRC_DIR}\" \"${WDEVFLAG}\" -DCMAKE_INSTALL_PREFIX=\"${install_path}\" -DMPIEXEC=\"${MPIEXEC}\" -DMPI_C_COMPILER=\"${MPICC}\" -DMPI_Fortran_COMPILER=\"${MPIFC}\""
  CC="${CC}" FC="${FC}" $CMAKE "${OPENCOARRAYS_SRC_DIR}" "${WDEVFLAG}" -DCMAKE_INSTALL_PREFIX="${install_path}" -DMPIEXEC="${MPIEXEC}" -DMPI_C_COMPILER="${MPICC}" -DMPI_Fortran_COMPILER="${MPIFC}"
  info "Building OpenCoarrays in ${PWD} with the command make -j${num_threads}"
  make "-j${num_threads}"
  if [[ ! -z ${SUDO:-} ]]; then
    printf "\nThe chosen installation path requires sudo privileges. Please enter password if prompted.\n"
  fi
  info "Installing OpenCoarrays in ${install_path} with the command ${SUDO:-} make install"
  ${SUDO:-} make install
}
