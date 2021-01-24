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

!>  \brief This include file contains the body of the submodules
!>  [ParaDRAM_mod@Test_runSampler_smod](@ref paradram_mod@test_runsampler_smod) and
!>  [ParaDISE_mod@@Test_runSampler_smod](@ref paradise_mod@test_runsampler_smod).
!>  \author Amir Shahmoradi

#if defined PARADRAM

#define ParaXXXX_mod ParaDRAM_mod
#define ParaXXXX_type ParaDRAM_type
#define test_ParaXXXX test_ParaDRAM
#define ParaXXXX_NML "&ParaDRAM"
#define ParaXXXX ParaDRAM
#define ParaXXXX_RefinedChain_mod ParaDRAM_RefinedChain_mod

#elif defined PARADISE

#define ParaXXXX_mod ParaDISE_mod
#define ParaXXXX_type ParaDISE_type
#define test_ParaXXXX test_ParaDISE
#define ParaXXXX_NML "&ParaDISE"
#define ParaXXXX ParaDISE
#define ParaXXXX_RefinedChain_mod ParaDISE_RefinedChain_mod

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with no input arguments or input file.
    module function test_runSampler_1() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_1" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with an internal input file.
    module function test_runSampler_2() result(assertion)
        use String_mod, only: num2str
        implicit none
        logical                 :: assertion
        type(ParaXXXX_type)     :: PD
        integer(IK) , parameter :: chainSize_ref = 100_IK
        integer(IK) , parameter :: userSeed_ref = 1111_IK

        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_2" &
                            , inputFile = ParaXXXX_NML//" randomSeed = "//num2str(userSeed_ref)//" chainSize = "//num2str(chainSize_ref)//" /" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred .and. PD%SpecBase%RandomSeed%userSeed==userSeed_ref .and. PD%SpecMCMC%ChainSize%val==chainSize_ref
        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Image%id, PD%Err%occurred                   = ", PD%Image%id, PD%Err%occurred
            write(Test%outputUnit,"(*(g0))") "Image%id, PD%Err%msg                        = ", PD%Image%id, PD%Err%msg
            write(Test%outputUnit,"(*(g0))") "Image%id, userSeed_ref                      = ", PD%Image%id, userSeed_ref
            write(Test%outputUnit,"(*(g0))") "Image%id, PD%SpecBase%RandomSeed%userSeed   = ", PD%Image%id, PD%SpecBase%RandomSeed%userSeed
            write(Test%outputUnit,"(*(g0))") "Image%id, chainSize_ref                     = ", PD%Image%id, chainSize_ref
            write(Test%outputUnit,"(*(g0))") "Image%id, PD%SpecMCMC%ChainSize%val         = ", PD%Image%id, PD%SpecMCMC%ChainSize%val
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
#endif
    end function test_runSampler_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with an external input file.
    module function test_runSampler_3() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaXXXX_mod@test_runSampler_3.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_3" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a path to a non-existing external input file.
    module function test_runSampler_4() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_4" &
                            , inputFile = " " &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with wrong SpecBase input values.
    module function test_runSampler_5() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_5" &
                            , inputFile = Test%inDir//"/Test_ParaXXXX_mod@test_runSampler_5.in" &
                            )
        assertion = assertion .and. PD%Err%occurred
#endif
    end function test_runSampler_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file.
    module function test_runSampler_6() result(assertion)
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD1, PD2
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_6" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_6" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP
#endif
    end function test_runSampler_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    module function test_runSampler_7() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&paramonte randomSeed = 1111 /" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecBase/test_runSampler_7" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with a wrong internal input namelist group:
    !> +    Infinity values for `domainLowerLimitVec` and `domainUpperLimitVec`.
    module function test_runSampler_8() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        type(ParaXXXX_type) :: PD
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaXXXX_mod@test_runSampler_8.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"@SpecBase/test_runSampler_8" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
#endif
    end function test_runSampler_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file.
    module function test_runSampler_9() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9.ref" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(2)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_9" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "process, PD2%Err%occurred", Test%Image%id, PD2%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START

            assertion = assertion .and. all(shape(PD2%RefinedChain%LogFuncState) == shape(PD1%RefinedChain%LogFuncState))

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "shape(PD1%RefinedChain%LogFuncState)", shape(PD1%RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))") "shape(PD2%RefinedChain%LogFuncState)", shape(PD2%RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if

            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! by default, the output precision is only 8 digits

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if

        end if
        ! LCOV_EXCL_STOP

        end block
#endif
    end function test_runSampler_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, with an ASCII output file.
    module function test_runSampler_10() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block

        type(ParaXXXX_type) :: PD1, PD2

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10.ref" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 50_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 50_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_10" &
                            , restartFileFormat = "ascii" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 50_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (PD2%Image%isLeader) assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! LCOV_EXCL_LINE ! by default, the output precision is only 8 digits

        end block
#endif
    end function test_runSampler_10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an
    !> output sample file, with an ASCII output file but with binary chain file.
    module function test_runSampler_11() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11.ref" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_11" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "verbose" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (PD2%Image%isLeader) assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! LCOV_EXCL_LINE ! by default, the output precision is only 8 digits

        end block
#endif
    end function test_runSampler_11

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an
    !> output sample file, with an ASCII output file but with binary chain file.
    module function test_runSampler_12() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2

        ! Run the fresh simulation as a reference run

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12.ref" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_12" &
                            , restartFileFormat = "ascii" &
                            , chainFileFormat = "binary" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD2%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (PD2%Image%isLeader) assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-6_RK ) ! LCOV_EXCL_LINE ! by default, the output precision is only 8 digits

        end block
#endif
    end function test_runSampler_12

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with an empty external input file.
    !> This should first search for the specific namelist, then for the generic namelist,
    !> then issues a warning and runs the simulation without setting the simulation specifications.
    module function test_runSampler_13() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED
        block

        type(ParaXXXX_type) :: PD

        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = Test%inDir//"/Test_ParaXXXX_mod@test_runSampler_13.in" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_13" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred

        end block
#endif
    end function test_runSampler_13

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler with an empty external input file.
    !> This should first search for the specific namelist, then for the generic namelist,
    !> then issues a warning and runs the simulation without setting the simulation specifications.
    module function test_runSampler_14() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD

        call PD%runSampler  ( ndim = 1_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , inputFile = "&nonsense /" &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_14" &
                            )
        assertion = assertion .and. .not. PD%Err%occurred

        end block

#endif
    end function test_runSampler_14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the read method of the `ParaMCMCRefinedChain_type` class can successfully read an external sample file.
    module function test_runSampler_15() result(assertion)
        use ParaXXXX_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
        use Statistics_mod, only: flatten
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        real(RK)    , parameter :: tolerance = 1.e-10_RK
        character(*), parameter :: DELIM = "delim"
        integer(IK) , parameter :: NDIM = 2_IK
        type(ParaXXXX_type)     :: PD
        type(RefinedChain_type) :: RefinedChain
        real(RK), allocatable   :: Difference(:,:)
        real(RK), allocatable   :: FlattenedLogFuncState(:,:)

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

        if (PD%Image%isLeader) then
        ! LCOV_EXCL_START

            RefinedChain = readRefinedChain( sampleFilePath = PD%SampleFile%Path%original, delimiter = PD%SpecBase%OutputDelimiter%val, ndim = PD%nd%val )

            assertion = assertion .and. RefinedChain%numRefinement == 0_IK

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "RefinedChain%numRefinement", RefinedChain%numRefinement
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if

            ! NOTE: Keep in mind that `PD%RefinedChain` is a weighted chain internally, but unweighted when read from the external file.
            FlattenedLogFuncState = flatten(nd = PD%nd%val + 1, np = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact, Point = PD%RefinedChain%LogFuncState, Weight = PD%RefinedChain%Weight)
            assertion = assertion .and. all( shape(FlattenedLogFuncState) == shape(RefinedChain%LogFuncState) )

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "shape(PD%RefinedChain%LogFuncState) ", shape(PD%RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))") "shape(RefinedChain%LogFuncState)    ", shape(RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "PD%RefinedChain%LogFuncState        ", PD%RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "RefinedChain%LogFuncState           ", RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if

            Difference = abs(FlattenedLogFuncState - RefinedChain%LogFuncState)
            assertion = assertion .and. all(Difference < tolerance)

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "PD%RefinedChain%LogFuncState", PD%RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "RefinedChain%LogFuncState   ", RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "Difference                  ", Difference
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if

        end if
        ! LCOV_EXCL_STOP

        end block

#endif
    end function test_runSampler_15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, in multichain parallelism.
    module function test_runSampler_16() result(assertion)
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_16" &
                            , parallelizationModel = "multi chain" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_16" &
                            , parallelizationModel = "multi chain" &
                            , outputRealPrecision = 16_IK &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP

        end block
#endif
    end function test_runSampler_16

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, in singlechain parallelism.
    !> Also, ensure the delayed rejection sampling is activated by setting a low value for `adaptiveUpdatePeriod` and
    !> a large value for `scaleFactor`.
    module function test_runSampler_17() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2
        real(RK), parameter :: targetAcceptanceRate(*) = [0.2_RK, 0.23_RK]

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_17" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "single chain" &
                            , greedyAdaptationCount = 30_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , delayedRejectionCount = 4_IK &
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , scaleFactor = "5 * gelman" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_17" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "single chain" &
                            , greedyAdaptationCount = 30_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , delayedRejectionCount = 4_IK &
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , scaleFactor = "5 * gelman" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP

        end block

#endif
    end function test_runSampler_17

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, in singlechain parallelism.
    !> Also, ensure the delayed rejection sampling is activated by setting a low value for `adaptiveUpdatePeriod` and
    !> a large value for `scaleFactor`.
    !> Also, use uniform proposal, as opposed to Normal in `test_runSampler_17()`.
    module function test_runSampler_18() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2
        real(RK), parameter :: targetAcceptanceRate(*) = [0.2_RK, 0.23_RK]

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_18" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "single chain" &
                            , greedyAdaptationCount = 30_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , delayedRejectionCount = 4_IK &
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_18" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "singlechain" &
                            , greedyAdaptationCount = 30_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , delayedRejectionCount = 4_IK &
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = 250_IK &
                            , chainSize = 500_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP

        end block

#endif
    end function test_runSampler_18

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, in multichain parallelism.
    !> Also, ensure the delayed rejection sampling is activated by setting a low value for `adaptiveUpdatePeriod` and
    !> a large value for `scaleFactor`.
    !> Also, avoid delayed rejection, as opposed to what is done in `test_runSampler_17()`.
    !> Also, request the `maxCumSum` sample refinement and count method.
    module function test_runSampler_19() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2
        real(RK), parameter :: targetAcceptanceRate(*) = [0.2_RK, 0.23_RK]

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_19" &
                            , sampleRefinementMethod = "MaxCumSumAutoCorr-median" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "multi chain" &
                            , greedyAdaptationCount = 2_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , chainFileFormat = "verbose" &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = -50_IK &
                            , chainSize = 700_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_19" &
                            , sampleRefinementMethod = "MaxCumSumAutoCorr-median" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "multi chain" &
                            , greedyAdaptationCount = 2_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , chainFileFormat = "verbose" &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = -50_IK &
                            , chainSize = 700_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP

        end block

#endif
    end function test_runSampler_19

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ParaXXXX sampler to restart a complete simulation without an output sample file, in multichain parallelism.
    !> Also, ensure the delayed rejection sampling is activated by setting a low value for `adaptiveUpdatePeriod` and
    !> a large value for `scaleFactor`.
    !> Also, avoid delayed rejection, as opposed to what is done in `test_runSampler_17()`.
    !> Also, request the `cutoffAutoCorr` sample refinement and count method.
    module function test_runSampler_20() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical             :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        type(ParaXXXX_type) :: PD1, PD2
        real(RK), parameter :: targetAcceptanceRate(*) = [0.2_RK, 0.23_RK]

        ! Run the fresh simulation

        call PD1%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_20" &
                            , sampleRefinementMethod = "cutoffAutoCorr-maximum" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "multi chain" &
                            , greedyAdaptationCount = 2_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , chainFileFormat = "binary" &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = -50_IK &
                            , chainSize = 700_IK &
                            )
        assertion = assertion .and. .not. PD1%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD1%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        call Test%Image%sync()

        block
            use System_mod, only: removeFile
            call removeFile(PD1%SampleFile%Path%original, PD1%Err) ! delete the sample file
        end block

        call Test%Image%sync()

        ! restart the simulation with the same configuration

        call PD2%runSampler ( ndim = 2_IK &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_20" &
                            , sampleRefinementMethod = "cutoffAutoCorr-maximum" &
                            , targetAcceptanceRate = targetAcceptanceRate &
                            , parallelizationModel = "multi chain" &
                            , greedyAdaptationCount = 2_IK & ! This must remain larger than adaptiveUpdatePeriod
                            , adaptiveUpdatePeriod = 1_IK &
                            , outputRealPrecision = 16_IK &
                            , chainFileFormat = "binary" &
                            , scaleFactor = "5 * gelman" &
                            , proposalModel = "uniform" &
                            , randomSeed = 12345_IK &
                            , sampleSize = -50_IK &
                            , chainSize = 700_IK &
                            )

        assertion = assertion .and. .not. PD2%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isVerboseMode) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))")   "process, PD2%Err%occurred(1)", Test%Image%id, PD1%Err%occurred
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (PD2%Image%isLeader) then
        ! LCOV_EXCL_START
            assertion = assertion .and. all( abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState) < 1.e-12_RK ) ! by default, the output precision is only 8 digits
            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(Test%outputUnit,"(*(g0,:,' '))")
                    write(Test%outputUnit,"(*(g0,:,' '))")   "process, Difference:", Test%Image%id, abs(PD2%RefinedChain%LogFuncState - PD1%RefinedChain%LogFuncState)
                    write(Test%outputUnit,"(*(g0,:,' '))")
                end if
                return
            end if
        end if
        ! LCOV_EXCL_STOP

        end block
#endif
    end function test_runSampler_20

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the read method of the `ParaMCMCRefinedChain_type` class can successfully read an external sample file.
    !> \remark
    !> This is similar to test #15, except that it also verifies the functionality of the optional argument `tenabled` to `readRefinedChain()`.
    module function test_runSampler_21() result(assertion)
        use ParaXXXX_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
        use Statistics_mod, only: flatten
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        real(RK)    , parameter :: tolerance = 1.e-10_RK
        character(*), parameter :: DELIM = "delim"
        integer(IK) , parameter :: NDIM = 2_IK
        type(ParaXXXX_type)     :: PD
        type(RefinedChain_type) :: RefinedChain
        real(RK), allocatable   :: Difference(:,:)
        real(RK), allocatable   :: FlattenedLogFuncState(:,:)

        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_21" &
                            , sampleRefineMentMethod = "batchMeans-med" &
                            , parallelizationModel = "multichain" &
                            , outputRealPrecision = 15_IK &
                            , outputDelimiter = DELIM &
                            , sampleSize = 10_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (PD%Image%isLeader) then

            RefinedChain = readRefinedChain ( sampleFilePath = PD%SampleFile%Path%original & ! LCOV_EXCL_LINE
                                            , delimiter = PD%SpecBase%OutputDelimiter%val & ! LCOV_EXCL_LINE
                                            , tenabled = .true. & ! LCOV_EXCL_LINE
                                            , ndim = PD%nd%val & ! LCOV_EXCL_LINE
                                            )

            RefinedChain = readRefinedChain ( sampleFilePath = PD%SampleFile%Path%original & ! LCOV_EXCL_LINE
                                            , delimiter = PD%SpecBase%OutputDelimiter%val & ! LCOV_EXCL_LINE
                                            , tenabled = .false. & ! LCOV_EXCL_LINE
                                            , ndim = PD%nd%val & ! LCOV_EXCL_LINE
                                            )

            assertion = assertion .and. RefinedChain%numRefinement == 0_IK

            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "RefinedChain%numRefinement", RefinedChain%numRefinement
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if
            ! LCOV_EXCL_STOP

            ! NOTE: Keep in mind that `PD%RefinedChain` is a weighted chain internally, but unweighted when read from the external file.
            FlattenedLogFuncState = flatten(nd = PD%nd%val + 1, np = PD%RefinedChain%Count(PD%RefinedChain%numRefinement)%compact, Point = PD%RefinedChain%LogFuncState, Weight = PD%RefinedChain%Weight)
            assertion = assertion .and. all( shape(FlattenedLogFuncState) == shape(RefinedChain%LogFuncState) )

            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "shape(PD%RefinedChain%LogFuncState) ", shape(PD%RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))") "shape(RefinedChain%LogFuncState)    ", shape(RefinedChain%LogFuncState)
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "PD%RefinedChain%LogFuncState        ", PD%RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "RefinedChain%LogFuncState           ", RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if
            ! LCOV_EXCL_STOP

            Difference = abs(FlattenedLogFuncState - RefinedChain%LogFuncState)
            assertion = assertion .and. all(Difference < tolerance)

            if (.not. assertion) then
            ! LCOV_EXCL_START
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "PD%RefinedChain%LogFuncState", PD%RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "RefinedChain%LogFuncState   ", RefinedChain%LogFuncState
                    write(*,"(10(g0,:,', '))") "Difference                  ", Difference
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if
            ! LCOV_EXCL_STOP

        end if

        end block

#endif
    end function test_runSampler_21

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether an object of [ChainFileContents_type](@ref paramcmcrefinedchain_mod::chainfilecontents_type)
    !> class can handle zero-count input MCMC chain for refinement.
    !> \remark
    module function test_runSampler_22() result(assertion)
        use ParaXXXX_RefinedChain_mod, only: readRefinedChain, RefinedChain_type
        implicit none
        logical                 :: assertion
        assertion = .true.
#if defined CODECOV_ENABLED || defined SAMPLER_TEST_ENABLED

        block

        real(RK)    , parameter :: tolerance = 1.e-10_RK
        character(*), parameter :: DELIM = "delim"
        integer(IK) , parameter :: NDIM = 2_IK
        type(RefinedChain_type) :: RefinedChain
        type(ParaXXXX_type)     :: PD

        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFuncMVN &
                            , mpiFinalizeRequested = .false. &
                            , outputFileName = Test%outDir//"/"//MODULE_NAME//"/test_runSampler_22" &
                            , sampleRefineMentMethod = "maxCumSumAutoCorr-minimum" &
                            , parallelizationModel = "multi chain" &
                            , chainSize = NDIM + 1_IK &
                            )
        assertion = assertion .and. .not. PD%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (PD%Image%isLeader) then
        ! LCOV_EXCL_START

            ! Now, set PD%Chain%ndim and PD%Chain%Count%Compact and to zero. 
            ! This should lead to the use of shape(PD%Chain%State(:,:)) to infer 
            ! ndim and Compact chain length.
            ! To test the behavior of the procedure in the presence of 1 point in the chain,
            ! we will also resize the input PD%Chain%State(1:ndim,1:PD%Chain%Count%Compact) to PD%Chain%State(1:ndim,1:1)

            PD%Chain%ndim = 0
            PD%Chain%Count%compact = 0
            PD%Chain%State = PD%Chain%State(1:NDIM, 1:1)

            call PD%RefinedChain%get(CFC = PD%Chain, Err = PD%Err)

            assertion = assertion .and. .not. PD%Err%occurred .and. all(PD%RefinedChain%IAC==0_IK) .and. size(PD%RefinedChain%Count)==1_IK .and. PD%RefinedChain%numRefinement==0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            if (.not. assertion) then
                if (Test%isVerboseMode) then
                    write(*,"(10(g0,:,', '))")
                    write(*,"(10(g0,:,', '))") "Test%Image%id, RefinedChain%numRefinement", RefinedChain%numRefinement
                    write(*,"(10(g0,:,', '))") "Test%Image%id, RefinedChain%IAC", RefinedChain%IAC
                    write(*,"(10(g0,:,', '))")
                end if
                return
            end if

        end if
        ! LCOV_EXCL_STOP

        end block

#endif
    end function test_runSampler_22

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaXXXX_RefinedChain_mod
#undef ParaXXXX_mod
#undef ParaXXXX_type
#undef test_ParaXXXX
#undef ParaXXXX_NML
#undef ParaXXXX
