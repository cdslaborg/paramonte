# shellcheck shell=bash
# Make the build directory, configure, and build
# shellcheck disable=SC2154

# shellcheck source=prerequisites/build-functions/edit_GCC_download_prereqs_file_if_necessary.sh
source "${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/edit_GCC_download_prereqs_file_if_necessary.sh"

function version { echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }'; }

build_and_install()
{
  num_threads=${arg_j}
  build_path="${OPENCOARRAYS_SRC_DIR}/prerequisites/builds/${package_to_build}-${version_to_build}"

  info "Building ${package_to_build} ${version_to_build}"
  info "Build path: ${build_path}"
  info "Installation path: ${install_path}"

  set_SUDO_if_needed_to_write_to_directory "${build_path}"
  set_SUDO_if_needed_to_write_to_directory "${install_path}"
  if [[ -d "${build_path}" ]]; then
    rm -rf "${build_path}"
  fi
  mkdir -p "${build_path}"
  info "pushd ${build_path}"
  pushd "${build_path}" || emergency "build_and_install.sh: pushd failed"

  if [[ "${package_to_build}" != "gcc" ]]; then

    if [[ "${package_to_build}" == "mpich" ]]; then
      if [[ "${version_to_build}" == "3.2" ]]; then
        info "Patching MPICH 3.2 on Mac OS due to segfault bug (see http://lists.mpich.org/pipermail/discuss/2016-May/004764.html)."
        sed 's/} MPID_Request ATTRIBUTE((__aligned__(32)));/} ATTRIBUTE((__aligned__(32))) MPID_Request;/g' \
          "${download_path}/${package_source_directory}/src/include/mpiimpl.h" > "${download_path}/${package_source_directory}/src/include/mpiimpl.h.patched"
        cp "${download_path}/${package_source_directory}/src/include/mpiimpl.h.patched" "${download_path}/${package_source_directory}/src/include/mpiimpl.h"
      fi

      FC_version=$($FC --version)
      text_before_dot="${FC_version%%.*}" # grab text before first dot
      major_version="${text_before_dot##* }" # grab text after final space
      if (( ${major_version} >= 10 )); then
       export FFLAGS="-w -fallow-argument-mismatch"
      fi
    fi

    if [[ "${package_to_build}" == "cmake"  && $(uname) == "Linux" ]]; then

      export cmake_binary_installer="${download_path}/cmake-${version_to_build}-Linux-x86_64.sh"
      ${SUDO:-} mkdir -p "$install_path"
      chmod u+x "${cmake_binary_installer}"
      if [[ -n "${SUDO:-}" ]]; then
        info "You do not have write permissions to the installation path ${install_path}"
        info "If you have administrative privileges, enter your password to install ${package_to_build}"
      fi
      info "Installing Cmake with the following command: "
      info "${SUDO:-} \"${cmake_binary_installer}\" --prefix=\"$install_path\" --exclude-subdir"
      ${SUDO:-} "${cmake_binary_installer}" --prefix="$install_path" --exclude-subdir

    else # build from source

      info "Configuring ${package_to_build} ${version_to_build} with the following command:"
      info "FFLAGS=${FFLAGS:-} FC=\"${FC:-'gfortran'}\" CC=\"${CC:-'gcc'}\" CXX=\"${CXX:-'g++'}\" \"${download_path}/${package_source_directory}\"/configure --prefix=\"${install_path}\""
      FFLAGS=${FFLAGS:-} FC="${FC:-'gfortran'}" CC="${CC:-'gcc'}" CXX="${CXX:-'g++'}" "${download_path}/${package_source_directory}"/configure --prefix="${install_path}"
      info "Building with the following command:"
      info "FC=\"${FC:-'gfortran'}\" CC=\"${CC:-'gcc'}\" CXX=\"${CXX:-'g++'}\" make -j\"${num_threads}\""
      FC="${FC:-'gfortran'}" CC="${CC:-'gcc'}" CXX="${CXX:-'g++'}" make "-j${num_threads}"
      info "Installing ${package_to_build} in ${install_path}"
      if [[ -n "${SUDO:-}" ]]; then
        info "You do not have write permissions to the installation path ${install_path}"
        info "If you have administrative privileges, enter your password to install ${package_to_build}"
      fi
      info "Installing with the following command: ${SUDO:-} make install"
      ${SUDO:-} make install

    fi

  elif [[ ${package_to_build} == "gcc" ]]; then

    # Warn about header prerequisite on macOS Mojave or subsequent versions
    if [[ $(uname) == "Darwin" ]]; then
      export kernel=$(uname -r)
      export Mojave="18.7.0"
      if [ $(version $kernel) -ge $(version $Mojave) ]; then
        info ""
        info "______________________________________________________________________________"
        info "Detected Darwin $kernel >= $Mojave (Mojave).  If $package_to_build build fails"
        info "due to a missing header (*.h) file, please try something like the following bash command:"
        info "open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg"
        info "Follow the prompts to install the missing headers. Then restart this $this_script."
        info "See https://bit.ly/build-gcc-on-mojave for more details."
        if [[ "${arg_y}" == "${__flag_present}" ]]; then
          info "-y or --yes-to-all flag present. Proceeding with non-interactive build."
        else
          info "Would you like to proceed anyway? (Y/n)"
          read -r proceed
          if [[ "${proceed}" == "n" || "${proceed}" == "N" || "${proceed}" == "no"  ]]; then
            info "n"
            emergency "Aborting. [exit 80]"
          else
            info "y"
          fi
        fi
      fi
    fi

    info "pushd ${download_path}/${package_source_directory} "
    pushd "${download_path}/${package_source_directory}" || emergency "build_and_install.sh: pushd failed"

    # Patch gfortran if necessary
    export patches_dir="${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/patches/${package_to_build}/${version_to_build}"
    if [[ -d "${patches_dir:-}" ]]; then
      for patch in "${patches_dir%/}"/*.diff ; do
	info "Applying patch ${patch##*/} to $package_to_build ${version_to_build}."
	patch -p1 < "$patch"
      done
    fi

    # Use GCC's contrib/download_prerequisites script after modifying it, if necessary, to use the
    # the preferred download mechanism set in prerequisites/build-functions/set_or_print_downloader.sh

    # Switch download mechanism if wget is not available
    edit_GCC_download_prereqs_file_if_necessary

    # Download GCC prerequisities
    "${PWD}"/contrib/download_prerequisites

    info "popd"
    popd || emergency "build_and_install.sh: popd failed"
    info "Configuring gcc/g++/gfortran builds with the following command:"
    info "${download_path}/${package_source_directory}/configure --prefix=${install_path} --enable-languages=c,c++,fortran,lto --disable-multilib --disable-werror ${bootstrap_configure}"
    "${download_path}/${package_source_directory}/configure" --prefix="${install_path}" --enable-languages=c,c++,fortran,lto --disable-multilib --disable-werror "${bootstrap_configure}"
    info "Building with the following command: make -j ${num_threads} ${bootstrap_build}"
    make -j ${num_threads} ${bootstrap_build:-}
    if [[ -n "${SUDO:-}" ]]; then
      info "You do not have write permissions to the installation path ${install_path}"
      info "If you have administrative privileges, enter your password to install ${package_to_build}"
    fi
    info "Installing with the following command: ${SUDO:-} make install"
    ${SUDO:-} make install

  else
     emergency "This branch should never be reached."
  fi # end if [[ "${package_to_build}" != "gcc"  && "${package_to_build}" != "cmake"  ]]; then


  info "popd"
  popd || emergency "build_and_install.sh: popd failed"
}
