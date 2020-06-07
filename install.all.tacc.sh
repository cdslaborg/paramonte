compilerSuiteList="intel gcc"
for cs in $compilerSuiteList; do
    ml ${cs}
    ml matlab
    ./install.sh --mem heap --lib dynamic --par "none mpi" >install.sh.${cs}.out 2>&1
done