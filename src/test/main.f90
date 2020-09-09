!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!use Test_Misc_mod
!use Test_Batse_mod
!use Test_Statistics_mod
use Test_Optimization_mod
!!use Test_EconomicsToolbox_mod
!use Test_BandSpectrum_mod
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
!call test_Misc()
!!call test_EconomicsToolbox()
!call test_Batse()
!call test_BandSpectrum()
!call test_Statistics()
call test_Optimization()
call test_ParaMonte()

!write(*,"(A)") compiler_options()

#if defined MPI_ENABLED
block
    use Test_mod, only: finalizeMPI
    call finalizeMPI()
end block
#endif

end program main