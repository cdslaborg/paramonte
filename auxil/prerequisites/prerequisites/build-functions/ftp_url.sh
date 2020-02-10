#! /bin/sh
# Download a file from an anonymous ftp site
#
# Usage:
#    ftp_url  <ftp-mode>  ftp://<fully-qualified-domain>:/<path-to-file>/<file-name>
#
# Example:
#    ftp_url -n ftp://ftp.gnu.org:/gnu/gcc/gcc-6.1.0/gcc-6.1.0.tar.bz2
ftp_url()
{
  ftp_mode="${1}"
  url="${2}"

  if [ "${ftp_mode}" != "-n" ]; then
    echo "Unexpected ftp mode received by ftp_url.sh: ${ftp_mode}"
  fi

  protocol="${url%%:*}" # grab text_before_first_colon
  if [ "${protocol}" != "ftp" ]; then
    echo "URL with unexpected protocol received by ftp_url.sh: ${text_before_first_colon}"
  fi

  text_after_double_slash="${url##*//}"
  FTP_SERVER="${text_after_double_slash%%/*}" # grab remaining text before first slash
  FILE_NAME="${url##*/}" # grab text after final slash

  text_after_final_colon="${url##*:}"
  FILE_PATH="${text_after_final_colon#*//}"
  FILE_PATH="${FILE_PATH#*/}" 
  FILE_PATH="${FILE_PATH%/*}" 

  USERNAME=anonymous
  PASSWORD="noone@nowhere.com"
  echo "starting anonymous download: ${protocol} ${ftp_mode} ${FTP_SERVER};... cd ${FILE_PATH}; ...; get ${FILE_NAME}"

ftp "${ftp_mode}" "${FTP_SERVER}" <<Done-ftp
user "${USERNAME}" "${PASSWORD}"
cd "${FILE_PATH}"
passive
binary
get "${FILE_NAME}"
bye
Done-ftp

echo "finished anonymous ftp"
}
