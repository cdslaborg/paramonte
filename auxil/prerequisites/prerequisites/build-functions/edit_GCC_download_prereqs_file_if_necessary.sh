# shellcheck disable=SC2154
replace_wget()
{
  # Define a file extension for the download_prerequisites backup
  backup_extension=".original"
  backup_file="${download_prereqs_file}${backup_extension}"
  if [[ -f ${backup_file}  ]]; then
    # Prevent overwriting an existing backup:
    backup_extension=""
  fi

  if [[ "${gcc_prereqs_fetch}" == "ftp_url" ]]; then
    # Insert a new line after line 2 to include ftp_url.sh as a download option
    sed -i${backup_extension} -e '2 a\'$'\n'". ${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/ftp_url.sh"$'\n' "${download_prereqs_file}"
  fi

  arg_string="${gcc_prereqs_fetch_args[*]:-}"

  info "Using the following command to replace wget in the GCC download_prerequisites file:"
  info "sed -i${backup_extension} s/\"${wget_command}\"/\"${gcc_prereqs_fetch} ${arg_string} \"/ \"${download_prereqs_file}\""
  if [[ "${__operating_system}" == "Linux" ]]; then
    sed -i"${backup_extension}" s/wget/"${gcc_prereqs_fetch} ${arg_string} "/ "${download_prereqs_file}"
  else
    sed -i "${backup_extension}" s/wget/"${gcc_prereqs_fetch} ${arg_string} "/ "${download_prereqs_file}"
  fi
}

edit_GCC_download_prereqs_file_if_necessary()
{
  __operating_system=$(uname)
  download_prereqs_file="${PWD}/contrib/download_prerequisites"

  # Grab the line with the first occurence of 'wget'
  wget_line=$(grep wget "${download_prereqs_file}" | head -1) || true
  wget_command="${wget_line%%ftp*}" # grab everything before ftp

  # Check for wget format used before GCC 7
  if [[ ! -z "${wget_command}" ]]; then
    # Check whether wget is available on this system
    if ! type wget &> /dev/null; then
      replace_wget
    fi
  fi
}
