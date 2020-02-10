# shellcheck shell=bash disable=SC2154,SC2148
print_header() {
  pushd "${__dir}" >/dev/null
  local gccver mpichver cmakever flexver bisonver m4ver
  gccver="$("${__file}" -V gcc)"
  mpichver="$("${__file}" -V mpich)"
  cmakever="$("${__file}" -V cmake)"
  flexver="$("${__file}" -V flex)"
  bisonver="$("${__file}" -V bison)"
  m4ver="$("${__file}" -V m4)"
  popd >/dev/null
  clear
  echo "__________ Please read these instructions carefully before proceeding. ___________"
  echo "__________________________________________________________________________________"
  echo "*** Most end users will find it easier to install OpenCoarrays via the native  ***"
  echo "*** package manager for the target operating system.  Advanced developers will ***"
  echo "*** find it best to install via CMake.  Please see the OpenCoarrays README.md  ***"
  echo "*** for a list of supported package installers or for instructions on using    ***"
  echo "*** CMake to install.  The current script offers a fallback when other methods ***"
  echo "*** prove infeasible. This script builds any missing prerequisites from source ***"
  echo "*** and then builds OpenCoarrays from source.  The process might take a few    ***"
  echo "*** minutes if the script finds all prerequisites or several hours if the      ***"
  echo "*** script finds no prerequisites. (Run './install.sh --help' without quotes   ***"
  echo "*** to see the installation options.)  This script traverses the following     ***"
  echo "*** dependency tree and asks permission to install packages that are not in    ***"
  echo "*** your PATH, not in the prerequisites/installations subdirectory of the      ***"
  echo "*** OpenCoarrays source tree, or have insufficiently recent version numbers:   ***"
  echo ""
  # Generate the following text using the `tree` command w/ dummy directory structure
  cat <<-EOF
    opencoarrays
    |__ cmake-${cmakever}
    |__ mpich-${mpichver}
        |__ gcc-${gccver}
            |__ flex-${flexver}
            |   |__ bison-${bisonver}
            |       |__ m4-${m4ver}
            |__ gmp
            |__ mpc
            |__ mpfr

EOF
  echo "*** For a non-interactive installation, relaunch this script with the command  ***"
  echo "*** './install.sh --yes-to-all' without quotes to instruct the script to       ***"
  echo "*** assume affirmative answers to all user queries.                            ***"
  echo ""
  printf "%s will be installed in %s\n" "${arg_p}" "${prefix_root:-${install_path}}"
  echo ""
  if [[ "${arg_y}" == "${__flag_present}" ]]; then
    info "-y or --yes-to-all flag present. Proceeding with non-interactive build."
  else
    printf "Ready to rock and roll? (Y/n)"
    read -r install_now
    printf '%s\n' "${install_now}"
    if [[ "${install_now}" == "n" || "${install_now}" == "no" ]]; then
      emergency "${this_script}: Aborting. [exit 85]"
    fi
  fi
}
