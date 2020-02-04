function download_all_prerequisites() 
{
  pushd ${OPENCOARRAYS_SRC_DIR}

  download_list=( "m4" "bison" "flex" "mpich" "cmake" )

  for package_to_download in "${download_list[@]}" ;
  do
    ./install.sh --package ${package_to_download} --only-download
  done
  
  if [[ ! -z ${arg_b:-}  ]]; then
    ./install.sh --package gcc  --only-download --install-branch "${arg_b}"
    cd prerequisites/downloads/${arg_b:-}
  else  # Download default version
    ./install.sh --package gcc  --only-download
    gcc_version=$(./install.sh -V gcc)
    cd prerequisites/downloads/
    tar xf gcc-${gcc_version}.tar.[bg]z* || emergency "tar didn't work"
    listing=$(ls -lt)
    cd gcc-${gcc_version}
  fi

  source ${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/set_or_print_downloader.sh
  set_or_print_downloader
  
  source ${OPENCOARRAYS_SRC_DIR}/prerequisites/build-functions/edit_GCC_download_prereqs_file_if_necessary.sh
  edit_GCC_download_prereqs_file_if_necessary
  
  ./contrib/download_prerequisites
  
  popd
}
