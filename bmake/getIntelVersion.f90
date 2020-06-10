program getIntelVersion 
use iso_fortran_env, only: compiler_version 
integer :: fileunit, startindex, endindex 
character(:), allocatable :: string 
string = trim(adjustl(compiler_version())) 
startindex = index(string,"Version") + 8 
endindex = index(string,"Build") - 2 
string = trim(adjustl(string(startindex:endindex))) 
open(newunit=fileunit,file="getIntelVersion.tmp") 
write(fileunit,"(A)") trim(adjustl(string)) 
end program getIntelVersion 
