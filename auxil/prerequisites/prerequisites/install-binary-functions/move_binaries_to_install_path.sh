install_or_skip()
{
  SUDO=${1}
  install_path="${2}"
  package_to_install="${3}"
  download_path="${4}"
  if [[ -d "${install_path}/${package_to_install}" ]]; then
    info "The following installation path exists:"
    info "${install_path}/${package_to_install}"
    info "If you want to replace it, please remove it and restart this script."
    info "Skipping ${package_to_install} installation."
  else
    info "Installing ${package_to_install} with the following command:"
    info "${SUDO} mv ${download_path}/${package_to_install} ${install_path}"
    ${SUDO:-} mv "${download_path}/${package_to_install}" "${install_path}"
    if [[ -x "${install_path}/${package_to_install}"  ]]; then
      info "Installation complete for ${package_to_install} in the following location:"
      info "${install_path}/${package_to_install}"
    else
      info "Something went wrong.  Either ${package_to_install} is not in the following expected location"
      info "or the user lacks executable permissions for the directory:"
      emergency "${install_path}/${package_to_install}"
    fi
  fi
}

function move_binaries_to_install_path()
{
  info "Installation package names: ${install_names}"
  package_to_install="${install_names%%,*}" # remove longest back-end match for ,*
  remaining_packages="${install_names#*,}"  # remove shortest front-end match for *,
  if [[ ! -d "${install_path}" ]]; then
    ${SUDO:-} mkdir -p "${install_path}"
  fi
  while [[ "${package_to_install}" != "${remaining_packages}" ]]; do
    install_or_skip  "${SUDO:-}" "${install_path}" "${package_to_install}" "${download_path}"
    info "Remaining installation package names:  ${remaining_packages}"
    install_names="${remaining_packages}"
    package_to_install="${install_names%%,*}" # remove longest back-end match for ,*
    remaining_packages="${install_names#*,}"  # remove shortest front-end for *,
    info "Installing ${package_to_install} binary with the following command:"
  done
  install_or_skip  "${SUDO:-}" "${install_path}" "${package_to_install}" "${download_path}"
}
