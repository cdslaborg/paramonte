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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This is main entry to the tests of the ParaMonte kernel library.
!>  \author Amir Shahmoradi

program main

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!block
!    use iso_fortran_env, only: IK=>int32, RK=>real64
!    use Timer_mod, only: Timer_type
!    integer(IK) :: i
!    real(RK) :: unifrnd(10**7)
!    type(Timer_type) :: Timer
!    call Timer%tic()
!    do i = 1, 10**7
!        call random_number(Unifrnd(i))
!    end do
!    call Timer%toc()
!    write(*,*) sum(Unifrnd), Timer%Time%delta
!end block
!stop

! result:  5000444.27932245 0.147000074386597 in debug mode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use Test_mod, only: setup, finalize

call setup()

#if defined BASIC_TEST_ENABLED
!block; use Test_BandSpectrum_mod; call test_BandSpectrum(); end block
!block; use Test_Batse_mod; call test_Batse(); end block
!block; use Test_Constants_mod; call test_Constants(); end block
!block; use Test_CorrCoef_mod; call test_CorrCoef(); end block
!block; use Test_Cosmology_mod; call test_Cosmology(); end block
!block; use Test_CrossCorr_mod; call test_CrossCorr(); end block
!block; use Test_DateTime_mod; call test_DateTime(); end block
!block; use Test_Decoration_mod; call test_Decoration(); end block
!block; use Test_Err_mod; call test_Err(); end block
!block; use Test_File_mod; call test_File(); end block
!block; use Test_FileContents_mod; call test_FileContents(); end block
!block; use Test_FileList_mod; call test_FileList(); end block
!block; use Test_Integration_mod; call test_Integration(); end block
block; use Test_Kmeans_mod; call test_Kmeans(); end block
!block; use Test_KmeansOOP_mod; call test_KmeansOOP(); end block
!block; use Test_Math_mod; call test_Math(); end block
!block; use Test_Matrix_mod; call test_Matrix(); end block
!block; use Test_Misc_mod; call test_Misc(); end block
block; use Test_Knn_mod; call test_Knn(); end block
block; use Test_MultiSkewNorm_mod; call test_MultiSkewNorm(); end block
block; use Test_PartitionMaxDen_mod; call Test_Partition(); end block
block; use Test_PartitionMinVol_mod; call Test_Partition(); end block
!block; use Test_PartitionOptDen_mod; call Test_Partition(); end block
!block; use Test_PartitionBenchm_mod; call Test_Partition(); end block
!block; use Test_Optimization_mod; call test_Optimization(); end block
!block; use Test_Parallelism_mod; call test_Parallelism(); end block
!block; use Test_Path_mod; call test_Path(); end block
!block; use Test_RandomSeed_mod; call test_RandomSeed(); end block
!block; use Test_Sort_mod; call test_Sort(); end block
!block; use Test_StarFormation_mod; call test_StarFormation(); end block
!block; use Test_Statistics_mod; call test_Statistics(); end block
!block; use Test_GeoCyclicFit_mod; call test_GeoCyclicFit(); end block
!block; use Test_String_mod; call test_String(); end block
!block; use Test_System_mod; call test_System(); end block
!block; use Test_Timer_mod; call test_Timer(); end block
!block; use Test_TimerCPU_mod; call test_TimerCPU(); end block
!block; use Test_TranGaus_mod; call test_TranGaus(); end block
!block; use Test_Unique_mod; call test_Unique(); end block
#endif

#if defined SAMPLER_TEST_ENABLED
block; use Test_ParaDRAM_mod; call test_ParaDRAM(); end block
block; use Test_ParaDISE_mod; call test_ParaDISE(); end block
block; use Test_ParaDRAM_RefinedChain_mod; call test_RefinedChain(); end block
block; use Test_ParaDISE_RefinedChain_mod; call test_RefinedChain(); end block
block; use Test_ParaDRAM_ChainFileContents_mod; call test_ChainFileContents(); end block
block; use Test_ParaDISE_ChainFileContents_mod; call test_ChainFileContents(); end block
!block; use Test_ParaNest_ChainFileContents_mod; call test_ChainFileContents(); end block
#endif

call finalize()

end program main
