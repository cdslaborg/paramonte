# shellcheck shell=bash disable=SC2148
# shellcheck source=./ftp_url.sh
source "${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/ftp_url.sh"
# shellcheck source=./set_SUDO_if_needed_to_write_to_directory.sh
source "${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/set_SUDO_if_needed_to_write_to_directory.sh"

# Download package if the tar ball is not already in the present working directory
# shellcheck disable=SC2154
download_if_necessary()
{
  download_path="${OPENCOARRAYS_SRC_DIR}/prerequisites/downloads"
  set_SUDO_if_needed_to_write_to_directory "${download_path}"

  # We set args regardless of whether this function performs a download because
  # GCC builds will need this to modify GCC's contrib/download_prerequisites script

  case "${fetch}" in
    "wget" )
      args=("--no-check-certificate")
      ;;
    "ftp_url" )
      args=("-n")
      ;;
    "git" )
      args=("clone")
      ;;
    "svn" )
       case "${arg_B:-}" in
         "gcc")
           args=("ls")
           ;;
         *)
           args=("checkout")
         ;;
       esac
       ;;
    "curl" )
      first_three_characters=$(echo "${package_url}" | cut -c1-3)
      case "${first_three_characters}" in
        "ftp"  )
           args=("-LO" "-u" "anonymous: ")
        ;;
        "htt"  )
           args=("-LO")
        ;;
        *)
          emergency "download_if_necessary.sh: Unrecognized URL."
        ;;
      esac
      ;;
    *)
      emergency "download_if_necessary.sh: Unrecognized \${fetch}=${fetch}."
      ;;
  esac

  case "${gcc_prereqs_fetch}" in
    "wget" )
      gcc_prereqs_fetch_args=("--no-check-certificate")
      ;;
    "ftp_url" )
      gcc_prereqs_fetch_args=("-n")
      ;;
    "curl" )
      gcc_prereqs_fetch_args=("-LO" "-u" "anonymous: ")
      ;;
    *)
      emergency "download_if_necessary.sh: Unrecognized \${gcc_prereqs_fetch}=${gcc_prereqs_fetch}."
      ;;
   esac


  if  [[ -f "${download_path}/${url_tail}" || -d "${download_path}/${url_tail##*branches/}" && ! -z ${url_tail##*branches/} ]]; then
    info "Found '${url_tail##*branches/}' in ${download_path}."
    info "If it resulted from an incomplete download, building ${package_name} could fail."
    if [[ "${arg_y}" == "${__flag_present}" ]]; then
      info "-y or --yes-to-all flag present. Proceeding with non-interactive build."
    else
      info "Would you like to proceed anyway? (Y/n)"
      read -r proceed
      if [[ "${proceed}" == "n" || "${proceed}" == "N" || "${proceed}" == "no"  ]]; then
        info "n"
        info "Please remove ${url_tail##*branches/} and restart the installation to to ensure a fresh download." 1>&2
        emergency "Aborting. [exit 80]"
      else
        info "y"
      fi
    fi
  elif ! type "${fetch}" &> /dev/null; then
    # The download mechanism is missing
    info "The default download mechanism for ${package_name} is ${fetch}."
    info "Please either ensure that ${fetch} is installed and in your PATH"
    info "or download the ${package_name} source from "
    info "${package_url}"
    info "Place the downloaded file in ${download_path} and restart this script."
    emergency "Aborting [exit 90]"
  else

    if [[ "${fetch}" == "svn" || "${fetch}" == "git" ]]; then
      package_source_directory="${url_tail}"
    else
      package_source_directory="${package_name}-${version_to_build}"
    fi
    info "Downloading ${package_name} ${version_to_build-} to the following location:"
    info "${download_path}/${package_source_directory}"
    info "Download command: \"${fetch}\" ${args[@]:-} ${package_url}"
    info "Depending on the file size and network bandwidth, this could take several minutes or longer."
    pushd "${download_path}"
    "${fetch}" ${args[@]:-} "${package_url}"
    popd
    if [[ ! -z "${arg_B:-}" ]]; then
      return
    else
      if [[ "${fetch}" == "svn" ]]; then
        search_path="${download_path}/${version_to_build}"
      else
        search_path="${download_path}/${url_tail}"
      fi
      if [[ -f "${search_path}" || -d "${search_path}" ]]; then
        info "Download succeeded. The ${package_name} source is in the following location:"
        info "${search_path}"
      else
        info "Download failed. The ${package_name} source is not in the following, expected location:"
        info "${search_path}"
        emergency "Aborting. [exit 110]"
      fi
    fi
  fi
}
