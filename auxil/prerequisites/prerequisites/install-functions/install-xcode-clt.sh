# shellcheck shell=bash disable=SC2154,SC2148
need_xcodeCLT() {
  if [[ "${OSTYPE}" != [Dd][Aa][Rr][Ww][Ii][Nn]* ]]; then
    echo "false"
  else
    local xcode_dir
    xcode_dir="$(/usr/bin/xcode-select -print-path 2>/dev/null)" || true
    if [[ ! -d "${xcode_dir}" || ! -x "${xcode_dir}/usr/bin/make" ]]; then
      echo "true"
    else
      echo "false"
    fi
  fi
}

xcode_clt_install () {
  app_store_label="$(softwareupdate -l | grep -B 1 -E "Command Line (Developer|Tools)" | awk -F"*" '/^ +\\*/ {print $2}' | sed 's/^ *//' | tail -n1)" || return 0
  if [[ "${arg_y}" == "${__flag_present}" ]]; then
    info "-y or --yes-to-all flag present. Proceeding with non-interactive download and install."
  else
    info "If you proceed with the Xcode-CLT installation you may be prompted for your password."
    printf "Ready to install %s? (Y/n)" "${app_store_label}"
    read -r install_xcodeclt
    printf '%s\n' "${install_xcodeclt}"
    if [[ "${install_xcodeclt}" == [nN]* ]]; then
      emergency "${this_script}: Aborting: cannot proceed without XCode CLT."
    fi
  fi
  info "install.sh will attempt to download and install XCode CLT for you now."
  info "Note: sudo privileges required. Please enter your password if prompted."
  place_holder="/tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress"
  sudo /usr/bin/touch "${place_holder}" || true
  sudo /usr/sbin/softwareupdate -i app_store_label || return 0
  sudo /bin/rm -f "${place_holder}" || true
  sudo /usr/bin/xcode-select --switch "/Library/Developer/CommandLineTools" || return 1
}

headless_xcode_clt_install () {
  # Fallback here if headless install/upgrade failed. This is the usual path through this code.
  if [[ "${arg_y}" == "${__flag_present}" ]]; then
    info "-y or --yes-to-all flag present. Proceeding with non-interactive download and install."
  else
    info "If you proceed with the Xcode-CLT installation you may be prompted for your password."
    printf "Ready to install Xcode-CLT using xcode-select? (Y/n)"
    read -r install_xcodeclt
    printf '%s\n' "${install_xcodeclt}"
    if [[ "${install_xcodeclt}" == [nN]* ]]; then
      emergency "${this_script}: Aborting: cannot proceed without XCode CLT."
    fi
  fi
  info "Installing Xcode Command Line Tools (CLT), using \`xcode-select --install\`."
  info "Note: sudo privileges required. Please enter your password if prompted."
  sudo /usr/bin/xcode-select --install || emergency "${this_script}: unable to run \`sudo xcode-select --install\`"
  printf "Please press <enter> once installation of Xcode command line tools has completed."
  read -r
  sudo /usr/bin/xcode-select --switch "/Library/Developer/CommandLineTools" || \
    emergency "${this_script}: Xcode-CLT installation and activation failed, unable to continue"
}

maybe_install_xcodeCLT () {
  if [[ "$(need_xcodeCLT)" == "true" ]]; then
    info "It appears that you are on Mac OS and do not have the Xcode command line tools (CLT) installed,"
    info "or they need to be upgraded.install.sh will now attempt to install/upgrade the Xcode-CLT using \`softwareupdate\`"
    info "This may take some time, please be patient."
    xcode_clt_install || true # This usually fails since `softwareupdate -i` doesn't return CLTs needing update
  fi
  if [[ "$(need_xcodeCLT)" == "true" ]]; then
    info "\`softwareupdate\` installation/upgrade failed. This is normal. Now trying with \`xcode-select --install\`"
    headless_xcode_clt_install || emergency "${this_script}: Could not install Xcode command line tools, unable to proceed. Please install them via the app store before retrying this installation."
  fi

  clang_output="$(/usr/bin/xzrun clang 2>&1)" || true
  if [[ "${clang_output}" =~ license ]]; then
    emergency "${this_script}: It appears you have not agreed to the Xcode license. Please do so before attempting to run this script again. This may be achieved by opening Xcode.app or running \`sudo xcodebuild -license\`"
  fi
}
