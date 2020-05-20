!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

program main

!use Test_Decoration_mod
!use Test_Err_mod
!use Test_DateTime_mod
!use Test_String_mod
!use Test_System_mod
!use Test_FileList_mod
!use Test_Path_mod
!use Test_Timer_mod
!use Test_TimerCPU_mod
!use Test_RandomSeed_mod
!use Test_FileContents_mod
!use Test_File_mod
!use Test_CrossCorr_mod
!use Test_Matrix_mod
!use Test_CorrCoef_mod
!use Test_TranGaus_mod
!use Test_Math_mod
!use Test_Batse_mod
!use Test_BandSpectrum_mod
!use Test_Statistics_mod
use Test_ParaMonte_mod

use iso_fortran_env, only: compiler_options

!use mpi
!implicit none
!integer :: ierrMPI

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!block
!use iso_fortran_env, only: IK=>int32, RK=>real64
!use Timer_mod, only: Timer_type
!integer(IK) :: i
!real(RK) :: unifrnd(10**7)
!type(Timer_type) :: Timer
!call Timer%tic()
!do i = 1, 10**7
!    call random_number(Unifrnd(i))
!end do
!call Timer%toc()
!write(*,*) sum(Unifrnd), Timer%Time%delta
!end block
!stop

! result:  5000444.27932245 0.147000074386597 in debug mode

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!call test_Decoration()
!call test_Err()
!call test_DateTime()
!call test_String()
!call test_System()
!call test_FileList()
!call test_Path()
!call test_Timer()
!call test_TimerCPU()
!call test_RandomSeed()
!call test_File()
!call test_CrossCorr()
!call test_Matrix()
!call test_CorrCoef()
!call test_TranGaus()
!call test_Math()
!call test_Batse()
!call test_BandSpectrum()
!call test_Statistics()

write(*,"(A)") compiler_options()

call test_ParaDRAM()

end program main