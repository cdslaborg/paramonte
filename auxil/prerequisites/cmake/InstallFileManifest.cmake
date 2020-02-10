
# This script will create and install an install manifest, including SHA256 hashes of each installed file
# Variables passed from CMake must be set with `install(CODE "set(...)")`

message(STATUS "Generating SHA256 checksumed receipt for all installed assets")

# Set receipt install destination
set(RECEIPT_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")

# Mimic cmake_install.cmake's handling of components
if(CMAKE_INSTALL_COMPONENT)
  set(HASHED_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.SHA256.txt")
else()
  set(HASHED_INSTALL_MANIFEST "install_manifest.SHA256.txt")
endif()

# Clean out any old install manifest on re-installation, if it exists
file(REMOVE "${CMAKE_BINARY_DIR}/${HASHED_INSTALL_MANIFEST}")

if(DEFINED ENV{DESTDIR})
  get_filename_component(ABS_DESTDIR "$ENV{DESTDIR}" ABSOLUTE)
else()
  set(ABS_DESTDIR "")
endif()
# Loop over files computing hashes
foreach(file IN LISTS CMAKE_INSTALL_MANIFEST_FILES)
  file(SHA256 "$ENV{DESTDIR}${file}" _file_sha256)
  file(RELATIVE_PATH FILE_REL_PATH "${ABS_DESTDIR}${RECEIPT_INSTALL_DIR}" "${ABS_DESTDIR}${file}")
  file(APPEND "${CMAKE_BINARY_DIR}/${HASHED_INSTALL_MANIFEST}" "${_file_sha256}  ${FILE_REL_PATH}\n")
endforeach()
file(APPEND "${CMAKE_BINARY_DIR}/${HASHED_INSTALL_MANIFEST}"
  "# Paths relative to \${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}\n")

file(INSTALL DESTINATION "${RECEIPT_INSTALL_DIR}" TYPE FILE PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ FILES
  "${CMAKE_BINARY_DIR}/${HASHED_INSTALL_MANIFEST}")

file(SHA256 "${CMAKE_BINARY_DIR}/${HASHED_INSTALL_MANIFEST}" MANIFEST_SHA256)
message(STATUS
  "Global checksum for OpenCoarrays installation:
       ${MANIFEST_SHA256}  ${HASHED_INSTALL_MANIFEST}")
message(STATUS "${PROJECT_NAME} was configured with SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH}")
message(STATUS "The current environment has SOURCE_DATE_EPOCH set to: $ENV{SOURCE_DATE_EPOCH}")
