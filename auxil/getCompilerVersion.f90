!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

program getCompilerVersion
use iso_fortran_env, only: compiler_version
implicit none
integer :: fileunit, startindex, endindex
logical :: isIntel = .false.
logical :: isGNU = .false.
character(:), allocatable :: string, isParaMonteCompatibleCompiler

isParaMonteCompatibleCompiler = "false"
string = trim(adjustl(compiler_version()))

if ( index(string,"GCC") > 0 ) then ! it's Gfortran
    startindex = index(string,"version") + 8
    endindex = len(string)
    isGNU = .true.
elseif ( index(string,"Intel") > 0 ) then ! it's Intel
    startindex = index(string,"Version") + 8
    endindex = index(string,"Build") - 2
    isIntel = .true.
end if

string = trim(adjustl(string(startindex:endindex)))
open(newunit=fileunit,file="getCompilerVersion.tmp")
write(fileunit,"(A)") string
close(fileunit)

! check for ParaMonteCompatibleCompiler

if (isIntel) then
    if ( lge(string(1:6),"18.0.0") ) isParaMonteCompatibleCompiler = "true"
elseif (isGNU) then
    if ( lge(string(1:3),"7.3") ) isParaMonteCompatibleCompiler = "true"
end if

open(newunit=fileunit,file="isParaMonteCompatibleCompiler.tmp")
write(fileunit,"(A)") isParaMonteCompatibleCompiler
close(fileunit)

end program getCompilerVersion