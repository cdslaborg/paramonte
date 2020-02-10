# If -P specifies a package, print the installation path and exit with normal status
# Otherwise, set install_path
# shellcheck disable=SC2154
set_or_print_installation_path()
{
  # Verify requirements
  [ ! -z "${arg_P}" ] && [ ! -z "${arg_p:-${arg_D:-${arg_U:-${arg_V}}}}" ] &&
    emergency "Please pass only one of {-D, -p, -P, -U, -V} or a longer equivalent (multiple detected)."

  if [[ $arg_p == "opencoarrays" ]]; then
    export version_to_build=$(${OPENCOARRAYS_SRC_DIR}/install.sh -V opencoarrays)
    export package_name="opencoarrays"
  fi
  arg_i_default="${OPENCOARRAYS_SRC_DIR}/prerequisites/installations/"
  if [[ "${arg_i%/}" == "${arg_i_default}" ]]; then
    if [[ -z "${arg_r}" ]]; then
      export install_path="${arg_i%/}/${arg_p:-${arg_D:-${arg_P:-${arg_U:-${arg_V}}}}}/${version_to_build}"
    else
      export install_path="${arg_r%/}/${arg_p:-${arg_D:-${arg_P:-${arg_U:-${arg_V}}}}}/${version_to_build}"
    fi
  else
    if [[ ! "${arg_i%/}" != "${arg_r}" ]]; then
      emergency "Please pass only one of {-i,-r} or a longer equivalent (multiple detected). [exit 105]"
    fi
    install_path="${arg_i%/}"
  fi

  # If -P is present, print ${install_path} and exit with normal status
  if [[ ! -z "${arg_P:-}" ]]; then
    printf "%s\n" "${install_path}"
    exit 0
  fi

  info "${package_name} ${version_to_build} installation path: ${install_path}"
}
