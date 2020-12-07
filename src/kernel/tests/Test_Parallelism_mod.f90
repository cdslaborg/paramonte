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

!>  \brief This module contains tests of the module [Parallelism_mod](@ref parallelism_mod).
!>  \author Amir Shahmoradi

module Test_Parallelism_mod

    use Test_mod, only: Test_type
    use Parallelism_mod
    implicit none

    private
    public :: test_Parallelism

    type(Test_type) :: Test

    real(RK)   , parameter :: comSecTime = 5.E-007_RK
    real(RK)   , parameter :: parSecTime = .7E-03_RK
    real(RK)   , parameter :: successProb = 0.234_RK
    integer(IK), parameter :: lenProcessID = 100_IK
    integer(IK), parameter :: processCount = 8_IK
    integer(IK), parameter :: ProcessID(lenProcessID) = [ 1_IK, 3_IK, 7_IK, 6_IK, 6_IK, 1_IK, 2_IK, 8_IK, 8_IK, 2_IK &
                                                        , 8_IK, 3_IK, 7_IK, 6_IK, 3_IK, 5_IK, 8_IK, 4_IK, 8_IK, 4_IK &
                                                        , 1_IK, 5_IK, 3_IK, 3_IK, 7_IK, 6_IK, 2_IK, 2_IK, 4_IK, 3_IK &
                                                        , 2_IK, 2_IK, 4_IK, 4_IK, 7_IK, 3_IK, 8_IK, 4_IK, 6_IK, 4_IK &
                                                        , 3_IK, 2_IK, 8_IK, 3_IK, 8_IK, 1_IK, 1_IK, 1_IK, 4_IK, 2_IK &
                                                        , 8_IK, 4_IK, 6_IK, 5_IK, 3_IK, 1_IK, 8_IK, 4_IK, 1_IK, 3_IK &
                                                        , 4_IK, 3_IK, 4_IK, 8_IK, 3_IK, 6_IK, 2_IK, 7_IK, 8_IK, 4_IK &
                                                        , 7_IK, 5_IK, 4_IK, 5_IK, 1_IK, 2_IK, 1_IK, 3_IK, 3_IK, 6_IK &
                                                        , 5_IK, 6_IK, 4_IK, 6_IK, 1_IK, 5_IK, 3_IK, 8_IK, 8_IK, 7_IK &
                                                        , 7_IK, 6_IK, 8_IK, 1_IK, 2_IK, 3_IK, 1_IK, 3_IK, 8_IK, 4_IK ]
    real(RK)   , parameter :: ForkJoinSpeedupScaling(*)=[ 1.0000000000000000_RK, 1.9709401945470562_RK, 2.9130986783067012_RK, 3.8267945495787576_RK &
                                                        , 4.7123840519965787_RK, 5.5702568576898672_RK, 6.4008325027127082_RK, 7.2045569830497680_RK &
                                                        , 7.9818995172384897_RK, 8.7333494795853746_RK, 9.4594135061080085_RK, 10.160612773695840_RK &
                                                        , 10.837480451544774_RK, 11.490559322673228_RK, 12.120399572259624_RK, 12.727556738641406_RK &
                                                        , 13.312589822070557_RK, 13.876059545717652_RK, 14.418526762942477_RK, 14.940551004491214_RK &
                                                        , 15.442689159026127_RK, 15.925494280231174_RK, 16.389514513655225_RK, 16.835292136442131_RK &
                                                        , 17.263362703145123_RK, 17.674254290921780_RK, 18.068486837547017_RK, 18.446571565858278_RK &
                                                        , 18.809010488451573_RK, 19.156295986674149_RK, 19.488910458202938_RK, 19.807326027754005_RK &
                                                        , 20.112004315731770_RK, 20.403396259894631_RK, 20.681941985383180_RK, 20.948070718725017_RK &
                                                        , 21.202200741694575_RK, 21.444739381165583_RK, 21.676083031346092_RK, 21.896617205030289_RK &
                                                        , 22.106716610736861_RK, 22.306745252829629_RK, 22.497056551931891_RK, 22.677993483151678_RK &
                                                        , 22.849888729829892_RK, 23.013064850708005_RK, 23.167834458585791_RK, 23.314500408703378_RK &
                                                        , 23.453355995235302_RK, 23.584685154428104_RK, 23.708762673047069_RK, 23.825854400923014_RK &
                                                        , 23.936217466505919_RK, 24.040100494440729_RK, 24.137743824279756_RK, 24.229379729539652_RK &
                                                        , 24.315232636394896_RK, 24.395519341379647_RK, 24.470449227540744_RK, 24.540224478552073_RK &
                                                        , 24.605040290360275_RK, 24.665085079988089_RK, 24.720540691171845_RK, 24.771582596555930_RK &
                                                        , 24.818380096209165_RK, 24.861096512266098_RK, 24.899889379530705_RK, 24.934910631911251_RK &
                                                        , 24.966306784583264_RK, 24.994219111802845_RK, 25.018783820315196_RK, 25.040132218323745_RK &
                                                        , 25.058390880003248_RK, 25.073681805556479_RK, 25.086122576828277_RK, 25.095826508503546_RK &
                                                        , 25.102902794926649_RK, 25.107456652589516_RK, 25.109589458344086_RK, 25.109398883401962_RK &
                                                        , 25.106979023190462_RK, 25.102420523139443_RK, 25.095810700477681_RK, 25.087233662121250_RK &
                                                        , 25.076770418739297_RK, 25.064498995084946_RK, 25.050494536680663_RK, 25.034829412948952_RK &
                                                        , 25.017573316879893_RK, 24.998793361327269_RK, 24.978554172025603_RK, 24.956917977419540_RK &
                                                        , 24.933944695397066_RK, 24.909692017016823_RK, 24.884215487319214_RK, 24.857568583309476_RK &
                                                        , 24.829802789199849_RK, 24.800967668996421_RK, 24.771110936514656_RK, 24.740278522906117_RK &
                                                        , 24.708514641776901_RK, 24.675861851976773_RK, 24.642361118136034_RK, 24.608051869025356_RK &
                                                        , 24.572972053811721_RK, 24.537158196282036_RK, 24.500645447103722_RK, 24.463467634189893_RK &
                                                        , 24.425657311234712_RK, 24.387245804482649_RK, 24.348263257793459_RK, 24.308738676062795_RK &
                                                        , 24.268699967056719_RK, 24.228173981716250_RK, 24.187186552986518_RK, 24.145762533223461_RK &
                                                        , 24.103925830228839_RK, 24.061699441963224_RK, 24.019105489984572_RK, 23.976165251658554_RK &
                                                        , 23.932899191185239_RK, 23.889326989485234_RK, 23.845467572986806_RK, 23.801339141354195_RK &
                                                        , 23.756959194195876_RK, 23.712344556790100_RK, 23.667511404863891_RK, 23.622475288460180_RK ]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Parallelism()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_constructForkJoin_1, "test_constructForkJoin_1")
        call Test%run(test_constructForkJoin_2, "test_constructForkJoin_2")
        call Test%run(test_constructForkJoin_3, "test_constructForkJoin_3")
        call Test%run(test_constructForkJoin_4, "test_constructForkJoin_4")
        call Test%run(test_constructForkJoin_5, "test_constructForkJoin_5")
        call Test%run(test_constructForkJoin_6, "test_constructForkJoin_6")
        call Test%run(test_constructForkJoin_7, "test_constructForkJoin_7")
        call Test%finalize()
    end subroutine test_Parallelism

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with valid input.
    function test_constructForkJoin_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        real(RK), parameter :: tolerance = 1.e-10_RK
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = successProb &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = .not. ForkJoin%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        assertion = assertion .and. ForkJoin%UniqueProcess%count == 8_IK
        assertion = assertion .and. all(ForkJoin%UniqueProcess%Identity == [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK, 7_IK, 8_IK])
        assertion = assertion .and. all(ForkJoin%UniqueProcess%Frequency == [13_IK, 11_IK, 18_IK, 16_IK, 7_IK, 11_IK, 8_IK, 16_IK])

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%count    ", ForkJoin%UniqueProcess%count
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_count    ", 8_IK
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%Identity ", ForkJoin%UniqueProcess%Identity
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_Identity ", [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK, 7_IK, 8_IK]
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%Frequency", ForkJoin%UniqueProcess%Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_Frequency", [13_IK, 11_IK, 18_IK, 16_IK, 7_IK, 11_IK, 8_IK, 16_IK]
            write(*,"(10(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. ForkJoin%Contribution%count == 8_IK
        assertion = assertion .and. all(ForkJoin%Contribution%Identity == [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK, 7_IK, 8_IK])
        assertion = assertion .and. all(ForkJoin%Contribution%Frequency == [13_IK, 11_IK, 18_IK, 16_IK, 7_IK, 11_IK, 8_IK, 16_IK])
        assertion = assertion .and. all( abs(ForkJoin%Contribution%LogFrequency-[ 2.5649493574615367_RK &
                                                                                , 2.3978952727983707_RK &
                                                                                , 2.8903717578961645_RK &
                                                                                , 2.7725887222397811_RK &
                                                                                , 1.9459101490553132_RK &
                                                                                , 2.3978952727983707_RK &
                                                                                , 2.0794415416798357_RK &
                                                                                , 2.7725887222397811_RK ] ) < 1.e-10_RK )

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%count"
            write(*,"(10(g0,:,', '))")  ForkJoin%Contribution%count
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%Identity"
            write(*,"(10(g0,:,', '))")  ForkJoin%Contribution%Identity
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%Frequency"
            write(*,"(10(g0,:,', '))")  ForkJoin%Contribution%Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%LogFrequency"
            write(*,"(10(g0,:,', '))")  ForkJoin%Contribution%LogFrequency
            write(*,"(10(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. abs(ForkJoin%SuccessProb%current - successProb) < tolerance
        assertion = assertion .and. abs(ForkJoin%SuccessProb%effective - ForkJoin%SuccessProb%PowellMinimum%xmin(1)) < tolerance
        assertion = assertion .and. all( abs(ForkJoin%SuccessProb%PowellMinimum%xmin - [0.28663337425270718E-1_RK, 4.5593657754033101_RK]) < tolerance )
        assertion = assertion .and. all( abs(ForkJoin%Contribution%LogFrequency-[ 2.5649493574615367_RK &
                                                                                , 2.3978952727983707_RK &
                                                                                , 2.8903717578961645_RK &
                                                                                , 2.7725887222397811_RK &
                                                                                , 1.9459101490553132_RK &
                                                                                , 2.3978952727983707_RK &
                                                                                , 2.0794415416798357_RK &
                                                                                , 2.7725887222397811_RK ] ) < tolerance )

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "ForkJoin%SuccessProb%current            ", ForkJoin%SuccessProb%current
            write(*,"(10(g0,:,', '))") "ForkJoin_SuccessProb_current            ", successProb
            write(*,"(10(g0,:,', '))") "ForkJoin%SuccessProb%effective          ", ForkJoin%SuccessProb%effective
            write(*,"(10(g0,:,', '))") "ForkJoin_SuccessProb_effective          ", ForkJoin%SuccessProb%PowellMinimum%xmin(1)
            write(*,"(10(g0,:,', '))") "ForkJoin%SuccessProb%PowellMinimum%xmin ", ForkJoin%SuccessProb%PowellMinimum%xmin
            write(*,"(10(g0,:,', '))") "ForkJoin_SuccessProb_PowellMinimum_xmin ", [0.28663337425270718E-1_RK, 4.5593657754033101_RK]
            write(*,"(10(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. ForkJoin%Speedup%count == size(ForkJoinSpeedupScaling)
        assertion = assertion .and. ForkJoin%Speedup%Maximum%nproc == 79_IK
        assertion = assertion .and. abs(ForkJoin%Speedup%Maximum%value - 25.109589458344086_RK) < tolerance
        assertion = assertion .and. abs(ForkJoin%Speedup%current - 7.2045569830497680_RK) < tolerance
        assertion = assertion .and. all(abs(ForkJoin%Speedup%Scaling - ForkJoinSpeedupScaling) < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))")
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%count"
            write(*,"(10(g0,:,', '))")  ForkJoin%Speedup%count
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%current"
            write(*,"(10(g0,:,', '))")  ForkJoin%Speedup%current
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%maximum%value"
            write(*,"(10(g0,:,', '))")  ForkJoin%Speedup%maximum%value
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%maximum%nproc"
            write(*,"(10(g0,:,', '))")  ForkJoin%Speedup%maximum%nproc
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%Scaling"
            write(*,"(10(g0,:,', '))")  ForkJoin%Speedup%Scaling
            write(*,"(10(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructForkJoin_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with a valid input `processCount == 1`.
    function test_constructForkJoin_2() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin
        real(RK), parameter :: tolerance = 1.e-10_RK

        integer(IK) , parameter :: processCount = 1_IK
        real(RK)    , parameter :: ForkJoin_Speedup_Maximum_value = 1._RK
        integer(IK) , parameter :: ForkJoin_Speedup_Maximum_nproc = 1_IK
        real(RK)    , parameter :: ForkJoin_Speedup_current = 1._RK
        integer(IK) , parameter :: ForkJoin_Speedup_count = 1_IK
        real(RK)    , parameter :: ForkJoin_Speedup_Scaling(*) = [1._RK]
        integer(IK) , parameter :: ForkJoin_Contribution_count = processCount
        integer(IK) , parameter :: ForkJoin_Contribution_Identity(*) = [1_IK]
        integer(IK) , parameter :: ForkJoin_Contribution_Frequency(*) = [lenProcessID]
        real(RK)    , parameter :: ForkJoin_Contribution_LogFrequency(*) = log(real(ForkJoin_Contribution_Frequency,kind=RK))
        integer(IK) , parameter :: ForkJoin_UniqueProcess_count = 1_IK
        integer(IK) , parameter :: ForkJoin_UniqueProcess_Identity(*) = [1_IK]
        integer(IK) , parameter :: ForkJoin_UniqueProcess_Frequency(*) = [lenProcessID]
        real(RK)    , parameter :: ForkJoin_SuccessProb_current = successProb
        real(RK)    , parameter :: ForkJoin_SuccessProb_effective = successProb

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = successProb &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = .not. ForkJoin%Err%occurred
        assertion = assertion .and. ForkJoin%Speedup%Maximum%value == ForkJoin_Speedup_Maximum_value
        assertion = assertion .and. ForkJoin%Speedup%Maximum%nproc == ForkJoin_Speedup_Maximum_nproc
        assertion = assertion .and. ForkJoin%Speedup%current == ForkJoin_Speedup_current
        assertion = assertion .and. ForkJoin%Speedup%count == ForkJoin_Speedup_count
        assertion = assertion .and. all(ForkJoin%Speedup%Scaling == ForkJoin_Speedup_Scaling)
        assertion = assertion .and. ForkJoin%Contribution%count == ForkJoin_Contribution_count
        assertion = assertion .and. all(ForkJoin%Contribution%Identity == ForkJoin_Contribution_Identity)
        assertion = assertion .and. all(ForkJoin%Contribution%Frequency == ForkJoin_Contribution_Frequency)
        assertion = assertion .and. all(ForkJoin%Contribution%LogFrequency == ForkJoin_Contribution_LogFrequency)
        assertion = assertion .and. ForkJoin%UniqueProcess%count == ForkJoin_UniqueProcess_count
        assertion = assertion .and. all(ForkJoin%UniqueProcess%Identity == ForkJoin_UniqueProcess_Identity)
        assertion = assertion .and. all(ForkJoin%UniqueProcess%Frequency == ForkJoin_UniqueProcess_Frequency)
        assertion = assertion .and. abs(ForkJoin%SuccessProb%current - ForkJoin_SuccessProb_current) < tolerance
        assertion = assertion .and. abs(ForkJoin%SuccessProb%effective - ForkJoin_SuccessProb_effective) < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%Maximum%value    ", ForkJoin%Speedup%Maximum%value
            write(*,"(10(g0,:,', '))") "ForkJoin_Speedup_Maximum_value    ", ForkJoin_Speedup_Maximum_value
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%Maximum%nproc    ", ForkJoin%Speedup%Maximum%nproc
            write(*,"(10(g0,:,', '))") "ForkJoin_Speedup_Maximum_nproc    ", ForkJoin_Speedup_Maximum_nproc
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%current          ", ForkJoin%Speedup%current
            write(*,"(10(g0,:,', '))") "ForkJoin_Speedup_current          ", ForkJoin_Speedup_current
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%count            ", ForkJoin%Speedup%count
            write(*,"(10(g0,:,', '))") "ForkJoin_Speedup_count            ", ForkJoin_Speedup_count
            write(*,"(10(g0,:,', '))") "ForkJoin%Speedup%Scaling          ", ForkJoin%Speedup%Scaling
            write(*,"(10(g0,:,', '))") "ForkJoin_Speedup_Scaling          ", ForkJoin_Speedup_Scaling
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%count       ", ForkJoin%Contribution%count
            write(*,"(10(g0,:,', '))") "ForkJoin_Contribution_count       ", ForkJoin_Contribution_count
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%Identity    ", ForkJoin%Contribution%Identity
            write(*,"(10(g0,:,', '))") "ForkJoin_Contribution_Identity    ", ForkJoin_Contribution_Identity
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%Frequency   ", ForkJoin%Contribution%Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin_Contribution_Frequency   ", ForkJoin_Contribution_Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin%Contribution%LogFrequency", ForkJoin%Contribution%LogFrequency
            write(*,"(10(g0,:,', '))") "ForkJoin_Contribution_LogFrequency", ForkJoin_Contribution_LogFrequency
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%count      ", ForkJoin%UniqueProcess%count
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_count      ", ForkJoin_UniqueProcess_count
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%Identity   ", ForkJoin%UniqueProcess%Identity
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_Identity   ", ForkJoin_UniqueProcess_Identity
            write(*,"(10(g0,:,', '))") "ForkJoin%UniqueProcess%Frequency  ", ForkJoin%UniqueProcess%Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin_UniqueProcess_Frequency  ", ForkJoin_UniqueProcess_Frequency
            write(*,"(10(g0,:,', '))") "ForkJoin%SuccessProb%current      ", ForkJoin%SuccessProb%current
            write(*,"(10(g0,:,', '))") "ForkJoin_SuccessProb_current      ", ForkJoin_SuccessProb_current
            write(*,"(10(g0,:,', '))") "ForkJoin%SuccessProb%effective    ", ForkJoin%SuccessProb%effective
            write(*,"(10(g0,:,', '))") "ForkJoin_SuccessProb_effective    ", ForkJoin_SuccessProb_effective
        end if
        ! LCOV_EXCL_STOP

    end function test_constructForkJoin_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with an invalid input `processCount < 1`.
    function test_constructForkJoin_3() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = 0_IK &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = successProb &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = ForkJoin%Err%occurred

    end function test_constructForkJoin_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with an invalid input `successProb = 0`.
    function test_constructForkJoin_4() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = -0.1_RK &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = ForkJoin%Err%occurred

    end function test_constructForkJoin_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with an invalid input `successProb = 1`.
    function test_constructForkJoin_5() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = 2._RK &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = ForkJoin%Err%occurred

    end function test_constructForkJoin_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with a valid but extreme input for `successProb`.
    function test_constructForkJoin_6() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = 1.00000000001_RK &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = .not. ForkJoin%Err%occurred

    end function test_constructForkJoin_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the ForkJoin constructor with a valid but extreme input for `successProb`.
    function test_constructForkJoin_7() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical :: assertion
        type(ForkJoin_type) :: ForkJoin

        ForkJoin = ForkJoin_type( processCount = processCount &
                                , lenProcessID = lenProcessID &
                                , ProcessID = ProcessID &
                                , successProb = -0.00000000001_RK &
                                , seqSecTime = epsilon(1._RK) &
                                , parSecTime = parSecTime &
                                , comSecTime = comSecTime &
                                )
        assertion = .not. ForkJoin%Err%occurred

    end function test_constructForkJoin_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Parallelism_mod