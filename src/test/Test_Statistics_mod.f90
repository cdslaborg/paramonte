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

module Test_Statistics_mod

    use Statistics_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Statistics

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Statistics()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_doSortedKS2()
        call Test%finalize()
    end subroutine test_Statistics

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_doSortedKS2()

        use Constants_mod, only: RK, IK
        use Sort_mod, only: sortAscending

        implicit none
        integer, parameter      :: NDATA = 50_IK
        real(RK) :: StdNormRnd1(NDATA) =    [ 0.537667139546100_RK &
                                            , 1.83388501459509_RK &
                                            , -2.25884686100365_RK &
                                            , 0.862173320368121_RK &
                                            , 0.318765239858981_RK &
                                            , -1.30768829630527_RK &
                                            , -0.433592022305684_RK &
                                            , 0.342624466538650_RK &
                                            , 3.57839693972576_RK &
                                            , 2.76943702988488_RK &
                                            , -1.34988694015652_RK &
                                            , 3.03492346633185_RK &
                                            , 0.725404224946106_RK &
                                            , -0.0630548731896562_RK &
                                            , 0.714742903826096_RK &
                                            , -0.204966058299775_RK &
                                            , -0.124144348216312_RK &
                                            , 1.48969760778547_RK &
                                            , 1.40903448980048_RK &
                                            , 1.41719241342961_RK &
                                            , 0.671497133608081_RK &
                                            , -1.20748692268504_RK &
                                            , 0.717238651328839_RK &
                                            , 1.63023528916473_RK &
                                            , 0.488893770311789_RK &
                                            , 1.03469300991786_RK &
                                            , 0.726885133383238_RK &
                                            , -0.303440924786016_RK &
                                            , 0.293871467096658_RK &
                                            , -0.787282803758638_RK &
                                            , 0.888395631757642_RK &
                                            , -1.14707010696915_RK &
                                            , -1.06887045816803_RK &
                                            , -0.809498694424876_RK &
                                            , -2.94428416199490_RK &
                                            , 1.43838029281510_RK &
                                            , 0.325190539456198_RK &
                                            , -0.754928319169703_RK &
                                            , 1.37029854009523_RK &
                                            , -1.71151641885370_RK &
                                            , -0.102242446085491_RK &
                                            , -0.241447041607358_RK &
                                            , 0.319206739165502_RK &
                                            , 0.312858596637428_RK &
                                            , -0.864879917324457_RK &
                                            , -0.0300512961962686_RK &
                                            , -0.164879019209038_RK &
                                            , 0.627707287528727_RK &
                                            , 1.09326566903948_RK &
                                            , 1.10927329761440_RK &
                                            ]
        real(RK) :: StdNormRnd2(NDATA) =    [ -0.863652821988714_RK &
                                            , 0.0773590911304249_RK &
                                            , -1.21411704361541_RK &
                                            , -1.11350074148676_RK &
                                            , -0.00684932810334806_RK &
                                            , 1.53263030828475_RK &
                                            , -0.769665913753682_RK &
                                            , 0.371378812760058_RK &
                                            , -0.225584402271252_RK &
                                            , 1.11735613881447_RK &
                                            , -1.08906429505224_RK &
                                            , 0.0325574641649735_RK &
                                            , 0.552527021112224_RK &
                                            , 1.10061021788087_RK &
                                            , 1.54421189550395_RK &
                                            , 0.0859311331754255_RK &
                                            , -1.49159031063761_RK &
                                            , -0.742301837259857_RK &
                                            , -1.06158173331999_RK &
                                            , 2.35045722400204_RK &
                                            , -0.615601881466894_RK &
                                            , 0.748076783703985_RK &
                                            , -0.192418510588264_RK &
                                            , 0.888610425420721_RK &
                                            , -0.764849236567874_RK &
                                            , -1.40226896933876_RK &
                                            , -1.42237592509150_RK &
                                            , 0.488193909859941_RK &
                                            , -0.177375156618825_RK &
                                            , -0.196053487807333_RK &
                                            , 1.41931015064255_RK &
                                            , 0.291584373984183_RK &
                                            , 0.197811053464361_RK &
                                            , 1.58769908997406_RK &
                                            , -0.804465956349547_RK &
                                            , 0.696624415849607_RK &
                                            , 0.835088165072682_RK &
                                            , -0.243715140377952_RK &
                                            , 0.215670086403744_RK &
                                            , -1.16584393148205_RK &
                                            , -1.14795277889859_RK &
                                            , 0.104874716016494_RK &
                                            , 0.722254032225002_RK &
                                            , 2.58549125261624_RK &
                                            , -0.666890670701386_RK &
                                            , 0.187331024578940_RK &
                                            , -0.0824944253709554_RK &
                                            , -1.93302291785099_RK &
                                            , -0.438966153934773_RK &
                                            , -1.79467884145512_RK &
                                            ]
        real(RK), parameter     :: REFERENCE_PROB_KS = 0.056045859714425_RK
        real(RK), parameter     :: REFERENCE_STAT_KS = 0.260000000000000_RK
        real(RK)                :: difference, statKS, probKS
        type(Err_type)          :: Err

        if (Test%Image%isFirst) call Test%testing("doSortedKS2()")

        call sortAscending( np = NDATA, Point = StdNormRnd1, Err = Err )
        if (Err%occurred) then; Test%assertion = .false.; call Test%verify(); end if

        call sortAscending( np = NDATA, Point = StdNormRnd2, Err = Err )
        if (Err%occurred) then; Test%assertion = .false.; call Test%verify(); end if
        

        call doSortedKS2( np1 = NDATA &
                        , np2 = NDATA &
                        , SortedPoint1 = StdNormRnd1 &
                        , SortedPoint2 = StdNormRnd2 &
                        , statKS = statKS &
                        , probKS = probKS &
                        )

        difference = 2 * abs(REFERENCE_PROB_KS - probKS) / (REFERENCE_PROB_KS + probKS)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "REFERENCE_PROB_KS, probKS, relative-difference:", REFERENCE_PROB_KS, probKS, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < 1.e-7_RK
        call Test%verify()

        difference = 2 * abs(REFERENCE_STAT_KS - statKS) / (REFERENCE_STAT_KS + statKS)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "REFERENCE_STAT_KS, statKS, relative-difference:", REFERENCE_STAT_KS, statKS, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < 1.e-7_RK
        call Test%verify()

        !call Test%skipping()

    end subroutine test_doSortedKS2

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Statistics_mod