./install.sh --mem heap --lib shared --par "none mpi" -s gnu >install.sh.gnu.out 2>&1
./install.sh --mem heap --lib shared --par "none" -s intel >install.sh.intel.out 2>&1
rm -rf ./bin/*.tar.gz
./auxil/btar.sh --dir ./bin/