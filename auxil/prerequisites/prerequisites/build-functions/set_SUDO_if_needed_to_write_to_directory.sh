# Define the sudo command to be used if the installation path requires administrative permissions
set_SUDO_if_needed_to_write_to_directory()
{
  if [[ $# != "1" ]]; then
    emergency "set_SUDO_if_needed_to_write_to_directory takes exactly one argument"
  else
    directory_to_create=$1
  fi
  if [[ -z "${LD_LIBRARY_PATH:-}" ]]; then
     info "\${LD_LIBRARY_PATH} is empty. Try setting it if the compiler encounters linking problems."
  fi
  SUDO_COMMAND="sudo env LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-} PATH=${PATH:-}"
  info "Checking whether the directory ${directory_to_create} exists... "
  if [[ -d "${directory_to_create}" ]]; then
    info "yes"
    info "Checking whether I have write permissions to ${directory_to_create} ... "
    if [[ -w "${directory_to_create}" ]]; then
      info "yes"
    else
      info "no"
      SUDO="${SUDO_COMMAND}"
    fi
  else
    info "no"
    info "Checking whether I can create ${directory_to_create} ... "
    if mkdir -p "${directory_to_create}" >& /dev/null; then
      info "yes."
    else
      info "no."
      # shellcheck disable=SC2034
      SUDO="${SUDO_COMMAND}"
    fi
  fi
}
