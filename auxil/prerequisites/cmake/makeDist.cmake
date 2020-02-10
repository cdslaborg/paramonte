# CMake file to be called in script mode (${CMAKE_COMMAND} -P <file>) to
# Generate a source archive release asset from add_custom_command
#
# See SourceDistTarget.cmake

if(NOT CMAKE_ARGV3)
  message(FATAL_ERROR "Must pass the top level src dir to ${CMAKE_ARGV2} as the first argument")
endif()

if(NOT CMAKE_ARGV4)
  message(FATAL_ERROR "Must pass the top level src dir to ${CMAKE_ARGV2} as the second argument")
endif()

find_package(Git)
if(NOT GIT_FOUND)
  message( FATAL_ERROR "You can't create a source archive release asset without git!")
endif()

execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always
  RESULT_VARIABLE git_status
  OUTPUT_VARIABLE git_version
  WORKING_DIRECTORY "${CMAKE_ARGV3}"
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT (git_status STREQUAL "0"))
  message( FATAL_ERROR "git describe --always failed with exit status: ${git_status} and message:
${git_version}")
endif()

set(archive "OpenCoarrays-${git_version}")
set(l_archive "opencoarrays-${git_version}")
set(release_asset "${CMAKE_ARGV4}/${archive}.tar.gz")
execute_process(
  COMMAND "${GIT_EXECUTABLE}" archive "--prefix=${archive}/" -o "${release_asset}" "${git_version}"
  RESULT_VARIABLE git_status
  OUTPUT_VARIABLE git_output
  WORKING_DIRECTORY "${CMAKE_ARGV3}"
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT (git_status STREQUAL "0"))
  message( FATAL_ERROR "git archive ... failed with exit status: ${git_status} and message:
${git_output}")
else()
  message( STATUS "Source code release asset created from `git archive`: ${release_asset}")
endif()

file(SHA256 "${release_asset}" tarball_sha256)
set(sha256_checksum "${tarball_sha256}  ${archive}.tar.gz")
configure_file("${CMAKE_ARGV3}/cmake/opencoarrays-VER-SHA256.txt.in"
  "${CMAKE_ARGV4}/${l_archive}-SHA256.txt"
  @ONLY)
message( STATUS
  "SHA 256 checksum of release tarball written out as: ${CMAKE_ARGV4}/${l_archive}-SHA256.txt" )

find_program(GPG_EXECUTABLE
  gpg
  DOC "Location of GnuPG (gpg) executable")

if(GPG_EXECUTABLE)
  execute_process(
    COMMAND "${GPG_EXECUTABLE}" --armor --detach-sign --comment "@gpg_comment@" "${CMAKE_ARGV4}/${l_archive}-SHA256.txt"
    RESULT_VARIABLE gpg_status
    OUTPUT_VARIABLE gpg_output
    WORKING_DIRECTORY "${CMAKE_ARGV4}")
  if(NOT (gpg_status STREQUAL "0"))
    message( WARNING "GPG signing of ${CMAKE_ARGV4}/${l_archive}-SHA256.txt appears to have failed
with status: ${gpg_status} and output: ${gpg_output}")
  else()
    configure_file("${CMAKE_ARGV3}/cmake/opencoarrays-VER-SHA256.txt.asc.in"
      "${CMAKE_ARGV4}/${l_archive}-GPG.comment"
      @ONLY)
    file(READ "${CMAKE_ARGV4}/${l_archive}-GPG.comment" gpg_comment)
    configure_file("${CMAKE_ARGV4}/${l_archive}-SHA256.txt.asc"
      "${CMAKE_ARGV4}/${l_archive}-SHA256.txt.asc.out"
      @ONLY)
    file(RENAME "${CMAKE_ARGV4}/${l_archive}-SHA256.txt.asc.out"
      "${CMAKE_ARGV4}/${l_archive}-SHA256.txt.asc")
    message(STATUS "GPG signed SHA256 checksum created: ${CMAKE_ARGV4}/${l_archive}-SHA256.txt.asc")
  endif()
  execute_process(
    COMMAND "${GPG_EXECUTABLE}" --armor --detach-sign --comment "Detached signature of ${archive}.tar.gz." "${release_asset}"
    RESULT_VARIABLE gpg_status
    RESULT_VARIABLE gpg_output
    WORKING_DIRECTORY "${CMAKE_ARGV4}"
    )
    if(NOT (gpg_status STREQUAL "0"))
      message( WARNING "GPG signing of ${release_asset} appears to have failed
with status: ${gpg_status} and output: ${gpg_output}")
    endif()
endif()
