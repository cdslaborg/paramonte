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

!>  \brief This module contains tests of the module [ParaMonteChainFileContents_mod](@ref ParaMonteChainFileContents_mod).
!>  \author Amir Shahmoradi

module Test_ParaMonteChainFileContents_mod

    use Test_mod, only: Test_type, getLogFuncMVN
    use ParaMonteChainFileContents_mod
    implicit none

    private
    public :: test_ParaMonteChainFileContents

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_ParaMonteChainFileContents()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_nullifyChainFileContents_1, "test_nullifyChainFileContents_1")
        call Test%finalize()
    end subroutine test_ParaMonteChainFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_nullifyChainFileContents_1() result(assertion)
        use Constants_mod, only: IK, RK
        use ParaDRAM_mod, only: ParaDRAM_type
        implicit none
        logical                     :: assertion
        type(ParaDRAM_type)         :: PD
        type(ChainFileContents_type):: CFC
        real(RK)    , parameter     :: tolerance = 1.e-10_RK
        character(*), parameter     :: DELIM = "delim"
       !character(9), parameter     :: variableNameList(*) = ["var1", "var2"]
        integer(IK) , parameter     :: NDIM = 2_IK
        assertion = .true.
#if defined CODECOV_ENABLED
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_15" &
                            , mpiFinalizeRequested = .false. &
                            , outputRealPrecision = 15_IK &
                            , outputDelimiter = DELIM &
                            , sampleSize = 10_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        ! read the chain file, first call without the optional input file, in which case no error should happen.

        CFC = ChainFileContents_type( ndim = NDIM )
        assertion = assertion .and. .not. CFC%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        CFC = ChainFileContents_type( ndim = NDIM &
                                    !, variableNameList = variableNameList &
                                    , chainFilePath = PD%ChainFile%Path%original &
                                    , targetChainSize = 2 * PD%Chain%Count%compact &
                                    , chainSize = PD%Chain%Count%compact &
                                    )
        assertion = assertion .and. .not. CFC%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE
#endif
    end function test_nullifyChainFileContents_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_ParaMonteChainFileContents_mod