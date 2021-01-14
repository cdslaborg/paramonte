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

!>  \brief This module contains tests of the module [Statistics_mod](@ref statistics_mod).
!>  \author Amir Shahmoradi

module Test_Statistics_mod

    use Statistics_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Statistics

    type(Test_type) :: Test

    integer(IK) , parameter :: lenRnd = 50_IK

    real(RK) :: UnifRnd(lenRnd) =   [ 0.0759666916908419_RK &
                                    , 0.2399161535536580_RK &
                                    , 0.1233189348351660_RK &
                                    , 0.1839077882824170_RK &
                                    , 0.2399525256649030_RK &
                                    , 0.4172670690843700_RK &
                                    , 0.0496544303257421_RK &
                                    , 0.9027161099152810_RK &
                                    , 0.9447871897216460_RK &
                                    , 0.4908640924680800_RK &
                                    , 0.4892526384000190_RK &
                                    , 0.3377194098213770_RK &
                                    , 0.9000538464176620_RK &
                                    , 0.3692467811202150_RK &
                                    , 0.1112027552937870_RK &
                                    , 0.7802520683211380_RK &
                                    , 0.3897388369612530_RK &
                                    , 0.2416912859138330_RK &
                                    , 0.4039121455881150_RK &
                                    , 0.0964545251683886_RK &
                                    , 0.1319732926063350_RK &
                                    , 0.9420505907754850_RK &
                                    , 0.9561345402298020_RK &
                                    , 0.5752085950784660_RK &
                                    , 0.0597795429471558_RK &
                                    , 0.2347799133724060_RK &
                                    , 0.3531585712220710_RK &
                                    , 0.8211940401979590_RK &
                                    , 0.0154034376515551_RK &
                                    , 0.0430238016578078_RK &
                                    , 0.1689900294627040_RK &
                                    , 0.6491154749564520_RK &
                                    , 0.7317223856586700_RK &
                                    , 0.6477459631363070_RK &
                                    , 0.4509237064309450_RK &
                                    , 0.5470088922863450_RK &
                                    , 0.2963208056077730_RK &
                                    , 0.7446928070741560_RK &
                                    , 0.1889550150325450_RK &
                                    , 0.6867754333653150_RK &
                                    , 0.1835111557372700_RK &
                                    , 0.3684845964903370_RK &
                                    , 0.6256185607296900_RK &
                                    , 0.7802274351513770_RK &
                                    , 0.0811257688657853_RK &
                                    , 0.9293859709687300_RK &
                                    , 0.7757126786084020_RK &
                                    , 0.4867916324031720_RK &
                                    , 0.4358585885809190_RK &
                                    , 0.4467837494298060_RK ]

    real(RK) :: StdNormRnd1(lenRnd) =   [ 0.537667139546100_RK &
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
                                        , 1.10927329761440_RK ]

    real(RK) :: StdNormRnd2(lenRnd) =   [ -0.863652821988714_RK &
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
                                        , -1.79467884145512_RK ]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Statistics()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_doKS1_1, "test_doKS1_1")
        call Test%run(test_flatten_1, "test_flatten_1")
        call Test%run(test_getMVNDev_1, "test_getMVNDev_1")
        call Test%run(test_getMVUDev_1, "test_getMVUDev_1")
        call Test%run(test_getHist1D_1, "test_getHist1D_1")
        call Test%run(test_getHist1D_2, "test_getHist1D_2")
        call Test%run(test_getHist1D_3, "test_getHist1D_3")
        call Test%run(test_getHist2D_1, "test_getHist2D_1")
        call Test%run(test_getHist2D_2, "test_getHist2D_2")
        call Test%run(test_getHist2D_3, "test_getHist2D_3")
        call Test%run(test_getHist2D_4, "test_getHist2D_4")
        call Test%run(test_getHist2D_5, "test_getHist2D_5")
        call Test%run(test_getRandMVN_1, "test_getRandMVN_1")
        call Test%run(test_getRandMVU_1, "test_getRandMVU_1")
        call Test%run(test_getMean_2D_1, "test_getMean_2D_1")
        call Test%run(test_getMean_2D_2, "test_getMean_2D_2")
        call Test%run(test_getRandExp_1, "test_getRandExp_1")
        call Test%run(test_getNormPDF_1, "test_getNormPDF_1")
        call Test%run(test_getNormCDF_1, "test_getNormCDF_1")
        call Test%run(test_getQuantile_1, "test_getQuantile_1")
        call Test%run(test_getQuantile_2, "test_getQuantile_2")
        call Test%run(test_getRandGaus_1, "test_getRandGaus_1")
        call Test%run(test_getRandNorm_1, "test_getRandNorm_1")
        call Test%run(test_getRandLogn_1, "test_getRandLogn_1")
        call Test%run(test_getRandBeta_1, "test_getRandBeta_1")
        call Test%run(test_getSNormPDF_1, "test_getSNormPDF_1")
        call Test%run(test_doSortedKS2_1, "test_doSortedKS2_1")
        call Test%run(test_doUniformKS1_1, "test_doUniformKS1_1")
        call Test%run(test_mergeMeanCov_1, "test_mergeMeanCov_1")
        call Test%run(test_getRandGamma_1, "test_getRandGamma_1")
        call Test%run(test_getLogProbMVU_1, "test_getLogProbMVU_1")
        call Test%run(test_getSamCholFac_1, "test_getSamCholFac_1")
        call Test%run(test_getSamCovMean_1, "test_getSamCovMean_1")
        call Test%run(test_getRandCorMat_1, "test_getRandCorMat_1")
        call Test%run(test_getLogProbGeo_1, "test_getLogProbGeo_1")
        call Test%run(test_getUniformCDF_1, "test_getUniformCDF_1")
        call Test%run(test_getRandUniform_1, "test_getRandUniform_1")
        call Test%run(test_getNormData_2D_1, "test_getNormData_2D_1")
        call Test%run(test_getVariance_1D_1, "test_getVariance_1D_1")
        call Test%run(test_getVariance_1D_2, "test_getVariance_1D_2")
        call Test%run(test_getVariance_2D_1, "test_getVariance_1D_1")
        call Test%run(test_getVariance_2D_2, "test_getVariance_2D_2")
        call Test%run(test_getBetaCDF_SPR_1, "test_getBetaCDF_SPR_1")
        call Test%run(test_getBetaCDF_DPR_1, "test_getBetaCDF_DPR_1")
        call Test%run(test_getSNormCDF_SPR_1, "test_getSNormCDF_SPR_1")
        call Test%run(test_getSNormCDF_DPR_1, "test_getSNormCDF_DPR_1")
        call Test%run(test_getMahalSqSP_RK_1, "test_getMahalSqSP_RK_1")
        call Test%run(test_getMahalSqMP_RK_1, "test_getMahalSqMP_RK_1")
        call Test%run(test_getMahalSqSP_CK_1, "test_getMahalSqSP_CK_1")
        call Test%run(test_getMahalSqMP_CK_1, "test_getMahalSqMP_CK_1")
        call Test%run(test_getLogProbLognSP_1, "test_getLogProbLognSP_1")
        call Test%run(test_getLogProbLognMP_1, "test_getLogProbLognMP_1")
        call Test%run(test_isInsideEllipsoid_1, "test_isInsideEllipsoid_1")
        call Test%run(test_mergeMeanCovUpper_1, "test_mergeMeanCovUpper_1")
        call Test%run(test_getRandIntLecuyer_1, "test_getRandIntLecuyer_1")
        call Test%run(test_getRandRealLecuyer_1, "test_getRandRealLecuyer_1")
        call Test%run(test_getSamCovMeanTrans_1, "test_getSamCovMeanTrans_1")
        call Test%run(test_getLogProbMVNSP_RK_1, "test_getLogProbMVNSP_RK_1")
        call Test%run(test_getLogProbMVNMP_RK_1, "test_getLogProbMVNMP_RK_1")
        call Test%run(test_getLogProbMVNSP_CK_1, "test_getLogProbMVNSP_CK_1")
        call Test%run(test_getLogProbMVNMP_CK_1, "test_getLogProbMVNMP_CK_1")
        call Test%run(test_getLogProbNormSP_RK_1, "test_getLogProbNormSP_RK_1")
        call Test%run(test_getLogProbNormMP_RK_1, "test_getLogProbNormMP_RK_1")
        call Test%run(test_getLogProbNormSP_CK_1, "test_getLogProbNormSP_CK_1")
        call Test%run(test_getLogProbNormMP_CK_1, "test_getLogProbNormMP_CK_1")
        call Test%run(test_getLogProbGeoCyclic_1, "test_getLogProbGeoCyclic_1")
        call Test%run(test_getRandGammaIntShape_1, "test_getRandGammaIntShape_1")
        call Test%run(test_getRandExpWithInvMean_1, "test_getRandExpWithInvMean_1")
       !call Test%run(test_mergeMeanCovUpperSlow_1, "test_mergeMeanCovUpperSlow_1")
        call Test%run(test_getLogProbMixMVNSP_RK_1, "test_getLogProbMixMVNSP_RK_1")
        call Test%run(test_getLogProbMixMVNMP_RK_1, "test_getLogProbMixMVNMP_RK_1")
        call Test%run(test_getLogProbMixMVNSP_CK_1, "test_getLogProbMixMVNSP_CK_1")
        call Test%run(test_getLogProbMixMVNMP_CK_1, "test_getLogProbMixMVNMP_CK_1")
        call Test%run(test_getRandCorMatRejection_1, "test_getRandCorMatRejection_1")
        call Test%run(test_mergeMeanCovUpperDense_1, "test_mergeMeanCovUpperDense_1")
        call Test%run(test_getLogProbMixNormSP_RK_1, "test_getLogProbMixNormSP_RK_1")
        call Test%run(test_getLogProbMixNormMP_RK_1, "test_getLogProbMixNormMP_RK_1")
        call Test%run(test_getLogProbMixNormSP_CK_1, "test_getLogProbMixNormSP_CK_1")
        call Test%run(test_getLogProbMixNormMP_CK_1, "test_getLogProbMixNormMP_CK_1")
        call Test%run(test_getRandPointOnEllipsoid_1, "test_getRandPointOnEllipsoid_1")
        call Test%run(test_getSamCovUpperMeanTrans_1, "test_getSamCovUpperMeanTrans_1")
        call Test%run(test_getWeiSamCovUppMeanTrans_1, "test_getWeiSamCovUppMeanTrans_1")
        call Test%run(test_getCovMatFromCorMatUpper_1, "test_getCovMatFromCorMatUpper_1")
        call Test%run(test_getCorMatUpperFromCovMatUpper_1, "test_getCorMatUpperFromCovMatUpper_1")
        call Test%run(test_getCovMatUpperFromCorMatUpper_1, "test_getCovMatUpperFromCorMatUpper_1")
        call Test%run(test_getCovMatUpperFromCorMatLower_1, "test_getCovMatUpperFromCorMatLower_1")
        call Test%run(test_getCovMatLowerFromCorMatUpper_1, "test_getCovMatLowerFromCorMatUpper_1")
        call Test%run(test_getCovMatLowerFromCorMatLower_1, "test_getCovMatLowerFromCorMatLower_1")
        call Test%finalize()
    end subroutine test_Statistics

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMahalSqSP_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mahalSq_ref = 180._RK
        real(RK)    , parameter :: Point(nd) = [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(InvCovMat) )
        real(RK)                :: mahalSq
        real(RK)                :: difference
        mahalSq = getMahalSqSP_RK(nd = nd, MeanVec = MeanVec, InvCovMat = InvCovMat, Point = Point)
        difference = abs(mahalSq - mahalSq_ref) / mahalSq_ref
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "mahalSq_ref    ", mahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMahalSqSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMahalSqMP_RK_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK, np = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: MahalSq_ref(np) = [180._RK, 36._RK]
        real(RK)    , parameter :: Point(nd,np) = reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point))
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                , 0._RK, 2._RK, 0._RK &
                                                                , 1._RK, 0._RK, 3._RK ], shape = shape(InvCovMat) )
        real(RK)                :: MahalSq(np)
        real(RK)                :: Difference(np)
        MahalSq = getMahalSqMP_RK(nd = nd, np = np, MeanVec = MeanVec, InvCovMat = InvCovMat, Point = Point)
        Difference = abs(MahalSq - MahalSq_ref) / MahalSq_ref
        assertion = all(Difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_ref    ", MahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq        ", MahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMahalSqMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMahalSqSP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        complex(CK) , parameter :: mahalSq_ref = 180._CK
        complex(CK) , parameter :: Point(nd) = [(cmplx(i,0.,kind=RK),i=1,nd)]
        complex(CK) , parameter :: MeanVec(nd) = [(cmplx(i**2+1._RK,0.,kind=RK),i=1,nd)]
        complex(CK) , parameter :: InvCovMat(nd,nd) = cmplx(reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK ], shape = shape(InvCovMat) ), kind = RK )
        complex(CK)             :: mahalSq
        real(RK)                :: difference
        mahalSq = getMahalSqSP_CK(nd = nd, MeanVec = MeanVec, InvCovMat = InvCovMat, Point = Point)
        difference = abs(real(mahalSq - mahalSq_ref,RK)) / real(mahalSq_ref,RK)
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "mahalSq_ref    ", mahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "mahalSq        ", mahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMahalSqSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMahalSqMP_CK_1() result(assertion)
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK)             :: i
        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK, np = 2_IK
        complex(CK) , parameter :: MahalSq_ref(np) = [180._CK, 36._CK]
        complex(CK) , parameter :: Point(nd,np) = reshape([(cmplx(i,0.,kind=RK),i=1,nd*np)], shape = shape(Point))
        complex(CK) , parameter :: MeanVec(nd) = [(cmplx(i**2+1._RK,0.,kind=RK),i=1,nd)]
        complex(CK) , parameter :: InvCovMat(nd,nd) = cmplx(reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK ], shape = shape(InvCovMat) ), kind=RK )

        complex(CK)             :: MahalSq(np)
        real(RK)                :: Difference(np)
        MahalSq = getMahalSqMP_CK(nd = nd, np = np, MeanVec = MeanVec, InvCovMat = InvCovMat, Point = Point)
        Difference = abs(real(MahalSq - MahalSq_ref,RK) / real(MahalSq_ref,RK))
        assertion = all(Difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_ref    ", MahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq        ", MahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMahalSqMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormSP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: mean = 3._RK
        real(RK)    , parameter :: point = 2._RK
        real(RK)    , parameter :: inverseVariance = 1._RK / 16._RK
        real(RK)    , parameter :: logProbNorm_ref = -2.336482894324563_RK
        real(RK)    , parameter :: logSqrtInverseVariance = log(sqrt(inverseVariance))
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)                :: logProbNorm
        real(RK)                :: difference

        logProbNorm = getLogProbNormSP_RK   ( mean = mean &
                                            , inverseVariance = inverseVariance &
                                            , logSqrtInverseVariance = logSqrtInverseVariance &
                                            , point = point &
                                            )

        difference = abs( (logProbNorm - logProbNorm_ref) / logProbNorm_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormMP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: np = 2_IK
        real(RK)    , parameter :: mean = 3._RK
        real(RK)    , parameter :: point(np) = [2._RK, 3._RK]
        real(RK)    , parameter :: inverseVariance = 1._RK / 16._RK
        real(RK)    , parameter :: logProbNorm_ref(np) = [ -2.336482894324563_RK, -2.305232894324563_RK ]
        real(RK)    , parameter :: logSqrtInverseVariance = log(sqrt(inverseVariance))
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)                :: logProbNorm(np)
        real(RK)                :: difference(np)

        logProbNorm = getLogProbNormMP_RK   ( np = np &
                                            , mean = mean &
                                            , inverseVariance = inverseVariance &
                                            , logSqrtInverseVariance = logSqrtInverseVariance &
                                            , point = point &
                                            )

        difference = abs( (logProbNorm - logProbNorm_ref) / logProbNorm_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVNSP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd) = [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                , shape = shape(InvCovMat) )
        real(RK)    , parameter :: logSqrtDetInvCovMat = -0.693147180559945_RK
        real(RK)    , parameter :: logProbNorm_ref = -15.19996278017396_RK
        real(RK)                :: logProbNorm
        real(RK)                :: difference

        logProbNorm = getLogProbMVNSP_RK( nd = nd &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbNorm - logProbNorm_ref) / logProbNorm_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMVNSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVNMP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd,np) = reshape( [(real(i,RK),i=1,nd*np)], shape = shape(Point) )
        real(RK)    , parameter :: MeanVec(nd) = [(real(i**2+1._RK,RK),i=1,nd)]
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                , shape = shape(InvCovMat) )
        real(RK)    , parameter :: logSqrtDetInvCovMat = -0.693147180559945_RK
        real(RK)    , parameter :: logProbNorm_ref(np) = [ -15.19996278017396_RK, -14.44996278017396_RK ]
        real(RK)                :: logProbNorm(np)
        real(RK)                :: difference(np)

        logProbNorm = getLogProbMVNMP_RK( nd = nd &
                                        , np = np &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbNorm - logProbNorm_ref) / logProbNorm_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMVNMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormSP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        complex(CK) , parameter :: mean = 3._RK
        complex(CK) , parameter :: point = cmplx(2._RK,kind=RK)
        complex(CK) , parameter :: inverseVariance = cmplx(1._RK / 16._RK,kind=RK)
        complex(CK) , parameter :: logProbNorm_ref = cmplx(-2.336482894324563_RK,kind=RK)
        complex(CK) , parameter :: logSqrtInverseVariance = cmplx(log(sqrt(inverseVariance)),kind=RK)
        complex(CK)             :: logProbNorm
        real(RK)                :: difference

        logProbNorm = getLogProbNormSP_CK   ( mean = mean &
                                            , inverseVariance = inverseVariance &
                                            , logSqrtInverseVariance = logSqrtInverseVariance &
                                            , point = point &
                                            )

        difference = abs( real(logProbNorm - logProbNorm_ref,kind=RK) / real(logProbNorm_ref,kind=RK) )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormMP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: np = 2_IK
        complex(CK) , parameter :: mean = 3._RK
        complex(CK) , parameter :: point(np) = cmplx([2._RK, 3._RK], kind=RK)
        complex(CK) , parameter :: inverseVariance = cmplx(1._RK / 16._RK, kind=RK)
        complex(CK) , parameter :: logProbNorm_ref(np) = cmplx([ -2.336482894324563_RK, -2.305232894324563_RK ], kind=RK)
        complex(CK) , parameter :: logSqrtInverseVariance = cmplx(log(sqrt(inverseVariance)), kind=RK)
        complex(CK)             :: logProbNorm(np)
        real(RK)                :: difference(np)

        logProbNorm = getLogProbNormMP_CK   ( np = np &
                                            , mean = mean &
                                            , inverseVariance = inverseVariance &
                                            , logSqrtInverseVariance = logSqrtInverseVariance &
                                            , point = point &
                                            )

        difference = abs( (real(logProbNorm, kind=RK) - real(logProbNorm_ref, kind=RK)) / real(logProbNorm_ref, kind=RK) )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVNSP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        complex(CK) , parameter :: Point(nd) = cmplx([(real(i,RK),i=1,nd)], kind=RK)
        complex(CK) , parameter :: MeanVec(nd) = cmplx([(real(i**2+1._RK,RK),i=1,nd)], kind=RK)
        complex(CK) , parameter :: InvCovMat(nd,nd) = cmplx( reshape( [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(InvCovMat) ), kind=RK)
        complex(CK) , parameter :: logSqrtDetInvCovMat = cmplx(-0.693147180559945_RK, kind=RK)
        complex(CK) , parameter :: logProbNorm_ref = cmplx(-15.19996278017396_RK, kind=RK)
        complex(CK)             :: logProbNorm
        real(RK)                :: difference

        logProbNorm = getLogProbMVNSP_CK( nd = nd &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (real(logProbNorm, kind=RK) - real(logProbNorm_ref, kind=RK)) / real(logProbNorm_ref, kind=RK) )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMVNSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVNMP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nd = 3_IK
        complex(CK) , parameter :: Point(nd,np) = reshape( cmplx([(real(i,RK),i=1,nd*np)], kind=RK), shape = shape(Point) )
        complex(CK) , parameter :: MeanVec(nd) = cmplx([(real(i**2+1._RK,RK),i=1,nd)], kind=RK)
        complex(CK) , parameter :: InvCovMat(nd,nd) = cmplx( reshape(  [ 1.500000000000000_RK, 0.000000000000000_RK, -0.50000000000000_RK &
                                                                    , 0.000000000000000_RK, 0.500000000000000_RK, 0.000000000000000_RK &
                                                                    , -0.50000000000000_RK, 0.000000000000000_RK, 0.500000000000000_RK ] &
                                                                    , shape = shape(InvCovMat) ), kind=RK)
        complex(CK) , parameter :: logSqrtDetInvCovMat = cmplx(-0.693147180559945_RK, kind=RK)
        complex(CK) , parameter :: logProbNorm_ref(np) = cmplx([ -15.19996278017396_RK, -14.44996278017396_RK ], kind=RK)
        complex(CK)             :: logProbNorm(np)
        real(RK)                :: difference(np)

        logProbNorm = getLogProbMVNMP_CK( nd = nd &
                                        , np = np &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (real(logProbNorm, kind=RK) - real(logProbNorm_ref, kind=RK)) / real(logProbNorm_ref, kind=RK) )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMVNMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormSP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: point = 2._RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: MeanVec(nmode) = [ 0.25_RK, 0.75_RK ]
        real(RK)    , parameter :: InvCovMat(nmode) = [ 1._RK / 16._RK, 1._RK / 32._RK ]
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: LogSqrtDetInvCovMat(nmode) = log(sqrt(InvCovMat))
        real(RK)    , parameter :: logProbMixNorm_ref = 1.718832134253714_RK
        real(RK)                :: logProbMixNorm
        real(RK)                :: difference

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , LogAmplitude = LogAmplitude &
                                            , MeanVec = MeanVec &
                                            , InvCovMat = InvCovMat &
                                            , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (logProbMixNorm - logProbMixNorm_ref) / logProbMixNorm_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm_ref ", logProbMixNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm     ", logProbMixNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormMP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: point(np) = [ 2._RK, 3._RK ]
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: MeanVec(nmode) = [ 0.25_RK, 0.75_RK ]
        real(RK)    , parameter :: InvCovMat(nmode) = [ 1._RK / 16._RK, 1._RK / 32._RK ]
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: LogSqrtDetInvCovMat(nmode) = log(sqrt(InvCovMat))
        real(RK)    , parameter :: logProbMixNorm_ref(np) = [ 1.718832134253714_RK, 1.636902047052812_RK ]
        real(RK)                :: logProbMixNorm(np)
        real(RK)                :: difference(np)

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , np = np &
                                            , LogAmplitude = LogAmplitude &
                                            , MeanVec = MeanVec &
                                            , InvCovMat = InvCovMat &
                                            , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (logProbMixNorm - logProbMixNorm_ref) / logProbMixNorm_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm_ref ", logProbMixNorm_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm     ", logProbMixNorm
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNSP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd) = [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: MeanVec(nd,nmode) = reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(MeanVec))
        real(RK)    , parameter :: InvCovMat(nd,nd,nmode) = reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK &
                                                                    , 2._RK, 0._RK, 0._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 0._RK, 0._RK, 2._RK &
                                                                    ], shape = shape(InvCovMat) )
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: LogSqrtDetInvCovMat(nmode) = [-0.693147180559945_RK, -1.039720770839918_RK]
        real(RK)    , parameter :: logProbMixMVN_ref = -90.44996278017396_RK
        real(RK)                :: logProbMixMVN
        real(RK)                :: difference

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , LogAmplitude = LogAmplitude &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbMixMVN - logProbMixMVN_ref) / logProbMixMVN_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN_ref ", logProbMixMVN_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN     ", logProbMixMVN
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNMP_RK_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd,np) = reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point))
        real(RK)    , parameter :: MeanVec(nd,nmode) = reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(MeanVec))
        real(RK)    , parameter :: InvCovMat(nd,nd,nmode) = reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK &
                                                                    , 2._RK, 0._RK, 0._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 0._RK, 0._RK, 2._RK &
                                                                    ], shape = shape(InvCovMat) )
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: LogSqrtDetInvCovMat(nmode) = [-0.693147180559945_RK, -1.039720770839918_RK]
        real(RK)    , parameter :: logProbMixMVN_ref(np) = [ -90.44996278017396_RK, -18.44996278017396_RK ]
        real(RK)                :: logProbMixMVN(np)
        real(RK)                :: difference(np)

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , np = np &
                                        , LogAmplitude = LogAmplitude &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbMixMVN - logProbMixMVN_ref) / logProbMixMVN_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN_ref ", logProbMixMVN_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN     ", logProbMixMVN
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormSP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nmode = 2_IK
        complex(CK) , parameter :: point = 2._RK
        complex(CK) , parameter :: MeanVec(nmode) = cmplx([ 0.25_RK, 0.75_RK ], kind = RK)
        complex(CK) , parameter :: InvCovMat(nmode) = cmplx([ 1._RK / 16._RK, 1._RK / 32._RK ], kind = RK)
        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
        complex(CK) , parameter :: LogSqrtDetInvCovMat(nmode) = cmplx(log(sqrt(InvCovMat)), kind = RK)
        complex(CK) , parameter :: logProbMixNorm_ref = cmplx(1.718832134253714_RK, kind = RK)
        complex(CK)             :: logProbMixNorm
        real(RK)                :: difference

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , LogAmplitude = LogAmplitude &
                                            , MeanVec = MeanVec &
                                            , InvCovMat = InvCovMat &
                                            , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (real(logProbMixNorm, kind = RK) - real(logProbMixNorm_ref, kind = RK)) / real(logProbMixNorm_ref, kind = RK) )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm_ref ", real(logProbMixNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm     ", real(logProbMixNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormMP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        complex(CK) , parameter :: point(np) = cmplx([ 2._RK, 3._RK ], kind = RK)
        complex(CK) , parameter :: MeanVec(nmode) = cmplx([ 0.25_RK, 0.75_RK ], kind = RK)
        complex(CK) , parameter :: InvCovMat(nmode) = cmplx([ 1._RK / 16._RK, 1._RK / 32._RK ], kind = RK)
        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
        complex(CK) , parameter :: LogSqrtDetInvCovMat(nmode) = cmplx(log(sqrt(InvCovMat)), kind = RK)
        complex(CK) , parameter :: logProbMixNorm_ref(np) = cmplx([ 1.718832134253714_RK, 1.636902047052812_RK ], kind = RK)
        complex(CK)             :: logProbMixNorm(np)
        real(RK)                :: difference(np)

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , np = np &
                                            , LogAmplitude = LogAmplitude &
                                            , MeanVec = MeanVec &
                                            , InvCovMat = InvCovMat &
                                            , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (real(logProbMixNorm, kind = RK) - real(logProbMixNorm_ref, kind = RK)) / real(logProbMixNorm_ref, kind = RK) )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm_ref ", real(logProbMixNorm_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixNorm     ", real(logProbMixNorm, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNSP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: nmode = 2_IK
        complex(CK) , parameter :: Point(nd) = cmplx([(real(i,RK),i=1,nd)], kind = RK)
        complex(CK) , parameter :: MeanVec(nd,nmode) = cmplx(reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(MeanVec)), kind = RK)
        complex(CK) , parameter :: InvCovMat(nd,nd,nmode) = cmplx( reshape( [ 1._RK, 0._RK, 1._RK &
                                                                            , 0._RK, 2._RK, 0._RK &
                                                                            , 1._RK, 0._RK, 3._RK &
                                                                            , 2._RK, 0._RK, 0._RK &
                                                                            , 0._RK, 2._RK, 0._RK &
                                                                            , 0._RK, 0._RK, 2._RK &
                                                                            ], shape = shape(InvCovMat) ), kind = RK)
        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
        complex(CK) , parameter :: LogSqrtDetInvCovMat(nmode) = cmplx([-0.693147180559945_RK, -1.039720770839918_RK], kind = RK)
        complex(CK) , parameter :: logProbMixMVN_ref = cmplx(-90.44996278017396_RK, kind = RK)
        complex(CK)             :: logProbMixMVN
        real(RK)                :: difference

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , LogAmplitude = LogAmplitude &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (real(logProbMixMVN, kind = RK) - real(logProbMixMVN_ref, kind = RK)) / real(logProbMixMVN_ref, kind = RK) )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN_ref  ", real(logProbMixMVN_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN      ", real(logProbMixMVN, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNMP_CK_1() result(assertion)

        use Constants_mod, only: IK, RK, CK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        complex(CK) , parameter :: Point(nd,np) = cmplx(reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point)), kind = RK)
        complex(CK) , parameter :: MeanVec(nd,nmode) = cmplx(reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(MeanVec)), kind = RK)
        complex(CK) , parameter :: InvCovMat(nd,nd,nmode) = cmplx(reshape(  [ 1._RK, 0._RK, 1._RK &
                                                                            , 0._RK, 2._RK, 0._RK &
                                                                            , 1._RK, 0._RK, 3._RK &
                                                                            , 2._RK, 0._RK, 0._RK &
                                                                            , 0._RK, 2._RK, 0._RK &
                                                                            , 0._RK, 0._RK, 2._RK &
                                                                            ], shape = shape(InvCovMat) ), kind = RK)
        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
        complex(CK) , parameter :: LogSqrtDetInvCovMat(nmode) = cmplx([-0.693147180559945_RK, -1.039720770839918_RK], kind = RK)
        complex(CK) , parameter :: logProbMixMVN_ref(np) = cmplx([ -90.44996278017396_RK, -18.44996278017396_RK ], kind = RK)
        complex(CK)             :: logProbMixMVN(np)
        real(RK)                :: difference(np)

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , np = np &
                                        , LogAmplitude = LogAmplitude &
                                        , MeanVec = MeanVec &
                                        , InvCovMat = InvCovMat &
                                        , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (real(logProbMixMVN, kind = RK) - real(logProbMixMVN_ref, kind = RK)) / real(logProbMixMVN_ref, kind = RK) )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN_ref  ", real(logProbMixMVN_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMixMVN      ", real(logProbMixMVN, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMean_2D_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean_ref(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK), allocatable   :: Mean(:)
        real(RK), allocatable   :: Difference(:)
        Mean = getMean(nd,np,Point)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = Mean)

        Difference = abs( (Mean - Mean_ref) / Mean_ref )
        assertion = all(Difference < tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref   ", real(Mean_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean       ", real(Mean, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMean_2D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMean_2D_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        integer(IK) , parameter :: Weight(np) = [(i,i=0,np-1)]
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean_ref(nd) = [.4601234030687523_RK, .5228363424875015_RK, .4616067715500543_RK]
        real(RK), allocatable   :: Mean(:)
        real(RK), allocatable   :: Difference(:)
        Mean = getMean(nd, np, Point, Weight)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = Mean)

        Difference = abs( (Mean - Mean_ref) / Mean_ref )
        assertion = all(Difference < tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref   ", real(Mean_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean       ", real(Mean, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getMean_2D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getNormData_2D_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK)    , parameter :: NormData_ref(nd,np) = reshape(   [ .2566442939639110_RK &
                                                                    , -.3041688273160970_RK &
                                                                    , -.2468837997938950_RK &
                                                                    , -.4032304034245440_RK &
                                                                    , -.2388698924576700_RK &
                                                                    , .2996510435725079_RK &
                                                                    , .2454268289201190_RK &
                                                                    , -.1890219363265699E-01_RK &
                                                                    , .4264152640835700_RK &
                                                                    , -.4149557135527890_RK &
                                                                    , .1027426859628800_RK &
                                                                    , -.1422483276617770_RK &
                                                                    , .3161149940933040_RK &
                                                                    , .4591982274435449_RK &
                                                                    , -.3369341802004060_RK &
                                                                    ], shape = shape(NormData_ref) )
        real(RK), allocatable   :: NormData(:,:)
        real(RK), allocatable   :: Difference(:,:)
        NormData = getNormData(nd,np,Mean,Point)

        ! Gfortran 7.1 fails to automatically reallocate this array. This is not implemented in Gfortran 7.0.0
        if (allocated(Difference)) deallocate(Difference); allocate(Difference, mold = NormData_ref)

        Difference = abs( (transpose(NormData) - NormData_ref) / NormData_ref )
        assertion = all(Difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "NormData_ref   ", real(NormData_ref, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "NormData       ", real(transpose(NormData), kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", real(difference, kind = RK)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getNormData_2D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getVariance_1D_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(*) = [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK ]
        integer(IK) , parameter :: np = size(Point)
        real(RK)    , parameter :: mean = 0.436403417501334_RK
        real(RK)    , parameter :: variance_ref = 0.106281559054025_RK
        real(RK)                :: variance
        real(RK)                :: difference

        variance = getVariance(np,mean,Point)
        difference = abs(variance - variance_ref) / variance_ref
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "variance_ref   ", variance_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "variance       ", variance
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getVariance_1D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getVariance_1D_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        integer(IK)             :: i
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(*) = [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK ]
        integer(IK) , parameter :: np = size(Point)
        integer(IK) , parameter :: Weight(np) = [(i,i=0,np-1)]
        real(RK)    , parameter :: mean = 0.436403417501334_RK
        real(RK)    , parameter :: variance_ref = .9447699069025800E-01_RK
        real(RK)                :: variance
        real(RK)                :: difference

        variance = getVariance(np,mean,Point,Weight,sum(Weight))
        difference = abs(variance - variance_ref) / variance_ref
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "variance_ref   ", variance_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "variance       ", variance
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getVariance_1D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getVariance_2D_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK), parameter     :: Variance_ref(nd) = [.1402030784811636_RK, .9283846639096889E-01_RK, .1165828911170294_RK]
        real(RK), allocatable   :: Variance(:)
        real(RK), allocatable   :: Difference(:)
        Variance = getVariance(nd,np,Mean,Point)
        Difference = abs(Variance - Variance_ref) / Variance_ref
        assertion = all(Difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Variance_ref   ", Variance_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Variance       ", Variance
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getVariance_2D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getVariance_2D_2() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        integer(IK)             :: i
        integer(IK) , parameter :: Weight(np) = [(i,i=0,np-1)]
        real(RK), parameter     :: Variance_ref(nd) = [.1332603228384714_RK, .1036548486974186_RK, .1075836700131122_RK]
        real(RK), allocatable   :: Variance(:)
        real(RK), allocatable   :: Difference(:)
        Variance = getVariance(nd,np,Mean,Point,Weight)
        Difference = abs(Variance - Variance_ref) / Variance_ref
        assertion = all(Difference <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Variance_ref   ", Variance_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Variance       ", Variance
            write(Test%outputUnit,"(*(g0,:,', '))") "Difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getVariance_2D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSamCholFac_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )

        real(RK)    , parameter :: Mean(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
       !integer(IK)             :: i
       !integer(IK) , parameter :: Weight(np) = [(i,i=0,np-1)]
        real(RK), parameter     :: CholeskyDiago_ref(nd) = [ +0.374437015372631_RK, +0.294661191115248_RK,  +0.306127969111099_RK ]
        real(RK), parameter     :: CholeskyLower_ref(nd,nd) = transpose(reshape([ +0.140203078481164_RK, +0.029035771029083_RK, -0.031754793421423_RK &
                                                                                , +0.077545140669886_RK, +0.092838466390969_RK, -0.043469498536710_RK &
                                                                                , -0.084806768876261_RK, -0.125205309782421_RK, +0.116582891117029_RK &
                                                                                ] , shape = shape(CholeskyLower_ref) ))
        real(RK)                :: CholeskyLower_diff(nd,nd)
        real(RK)                :: CholeskyDiago_diff(nd)
        real(RK)                :: CholeskyLower(nd,nd)
        real(RK)                :: CholeskyDiago(nd)

        call getSamCholFac(nd,np,Mean,Point,CholeskyLower,CholeskyDiago)
        CholeskyLower_diff = abs(CholeskyLower - CholeskyLower_ref) / abs(CholeskyLower_ref)
        CholeskyDiago_diff = abs(CholeskyDiago - CholeskyDiago_ref) / abs(CholeskyDiago_ref)
        assertion = all(CholeskyLower_diff <= tolerance) .and. all(CholeskyDiago_diff <= tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower_ref  ", CholeskyLower_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower      ", CholeskyLower
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyLower_diff ", CholeskyLower_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiago_ref  ", CholeskyDiago_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiago      ", CholeskyDiago
            write(Test%outputUnit,"(*(g0,:,', '))") "CholeskyDiago_diff ", CholeskyDiago_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSamCholFac_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSamCovMean_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-11_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(np,nd) = transpose(reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                                    , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                                    , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                                    , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                                    , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                                    ], shape = [nd,np] ) )

        real(RK), parameter     :: Mean_ref(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK), parameter     :: CovMat_ref(nd,nd) = reshape( [ +0.140203078481164_RK, +0.029035771029083_RK, -0.031754793421423_RK &
                                                                , +0.029035771029083_RK, +0.092838466390969_RK, -0.043469498536710_RK &
                                                                , -0.031754793421423_RK, -0.043469498536710_RK, +0.116582891117029_RK &
                                                                ] , shape = shape(CovMat_ref) )
        real(RK), parameter     :: InvCovMat_ref(nd,nd) = reshape(  [ +7.831154269923546_RK, -1.757283899385404_RK, +1.477819211284223_RK &
                                                                    , -1.757283899385404_RK, 13.444000278032743_RK, +4.534128105255763_RK &
                                                                    , +1.477819211284223_RK, +4.534128105255763_RK, 10.670726269401085_RK &
                                                                    ] , shape = shape(CovMat_ref) )
        real(RK), parameter     :: MahalSq_ref(np) = [3.178088257444105_RK, 1.653804994353691_RK, 2.669296951657121_RK, 1.8980344204538842_RK, 2.6007753760912014_RK]
        real(RK), parameter     :: sqrtDetInvCovMat_ref = 29.607059382128476_RK
        real(RK)                :: CovMat(nd,nd)
        real(RK)                :: CovMat_diff(nd,nd)
        real(RK)                :: Mean_diff(nd), MahalSq_diff(np), InvCovMat_diff(nd,nd), sqrtDetInvCovMat_diff
        real(RK)                :: Mean(nd), MahalSq(np), InvCovMat(nd,nd), sqrtDetInvCovMat

        assertion = .true.

        call getSamCovMean(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
        CovMat_diff = abs(CovMat - CovMat_ref) / abs(CovMat_ref)
        assertionCurrent = all(CovMat_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat_ref  ", CovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat      ", CovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat_diff ", CovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Mean_diff = abs(Mean - Mean_ref) / abs(Mean_ref)
        assertionCurrent = all(Mean_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref  ", Mean_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean      ", Mean
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_diff ", Mean_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        MahalSq_diff = abs(MahalSq - MahalSq_ref) / abs(MahalSq_ref)
        assertionCurrent = all(MahalSq_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_ref  ", MahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq      ", MahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_diff ", MahalSq_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        InvCovMat_diff = abs(InvCovMat - InvCovMat_ref) / abs(InvCovMat_ref)
        assertionCurrent = all(InvCovMat_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat_ref  ", InvCovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat      ", InvCovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat_diff ", InvCovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        sqrtDetInvCovMat_diff = abs(sqrtDetInvCovMat - sqrtDetInvCovMat_ref) / abs(sqrtDetInvCovMat_ref)
        assertionCurrent = sqrtDetInvCovMat_diff <= tolerance
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat_ref  ", sqrtDetInvCovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat      ", sqrtDetInvCovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat_diff ", sqrtDetInvCovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSamCovMean_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSamCovMeanTrans_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-11_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )
        real(RK), parameter     :: Mean_ref(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK), parameter     :: CovMat_ref(nd,nd) = reshape( [ +0.140203078481164_RK, +0.029035771029083_RK, -0.031754793421423_RK &
                                                                , +0.029035771029083_RK, +0.092838466390969_RK, -0.043469498536710_RK &
                                                                , -0.031754793421423_RK, -0.043469498536710_RK, +0.116582891117029_RK &
                                                                ] , shape = shape(CovMat_ref) )
        real(RK), parameter     :: InvCovMat_ref(nd,nd) = reshape(  [ +7.831154269923546_RK, -1.757283899385404_RK, +1.477819211284223_RK &
                                                                    , -1.757283899385404_RK, 13.444000278032743_RK, +4.534128105255763_RK &
                                                                    , +1.477819211284223_RK, +4.534128105255763_RK, 10.670726269401085_RK &
                                                                    ] , shape = shape(CovMat_ref) )
        real(RK), parameter     :: MahalSq_ref(np) = [3.178088257444105_RK, 1.653804994353691_RK, 2.669296951657121_RK, 1.8980344204538842_RK, 2.6007753760912014_RK]
        real(RK), parameter     :: sqrtDetInvCovMat_ref = 29.607059382128476_RK
        real(RK)                :: CovMat(nd,nd)
        real(RK)                :: CovMat_diff(nd,nd)
        real(RK)                :: Mean_diff(nd), MahalSq_diff(np), InvCovMat_diff(nd,nd), sqrtDetInvCovMat_diff
        real(RK)                :: Mean(nd), MahalSq(np), InvCovMat(nd,nd), sqrtDetInvCovMat

        assertion = .true.

        call getSamCovMeanTrans(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
        CovMat_diff = abs(CovMat - CovMat_ref) / abs(CovMat_ref)
        assertionCurrent = all(CovMat_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat_ref  ", CovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat      ", CovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMat_diff ", CovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Mean_diff = abs(Mean - Mean_ref) / abs(Mean_ref)
        assertionCurrent = all(Mean_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref  ", Mean_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean      ", Mean
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_diff ", Mean_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        MahalSq_diff = abs(MahalSq - MahalSq_ref) / abs(MahalSq_ref)
        assertionCurrent = all(MahalSq_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_ref  ", MahalSq_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq      ", MahalSq
            write(Test%outputUnit,"(*(g0,:,', '))") "MahalSq_diff ", MahalSq_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        InvCovMat_diff = abs(InvCovMat - InvCovMat_ref) / abs(InvCovMat_ref)
        assertionCurrent = all(InvCovMat_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat_ref  ", InvCovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat      ", InvCovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "InvCovMat_diff ", InvCovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        sqrtDetInvCovMat_diff = abs(sqrtDetInvCovMat - sqrtDetInvCovMat_ref) / abs(sqrtDetInvCovMat_ref)
        assertionCurrent = sqrtDetInvCovMat_diff <= tolerance
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat_ref  ", sqrtDetInvCovMat_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat      ", sqrtDetInvCovMat
            write(Test%outputUnit,"(*(g0,:,', '))") "sqrtDetInvCovMat_diff ", sqrtDetInvCovMat_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSamCovMeanTrans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSamCovUpperMeanTrans_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        integer(IK)             :: i, j
        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )
        real(RK), parameter     :: Mean_ref(nd) = [0.449401794055698_RK, 0.336001673693518_RK, 0.523806784754785_RK]
        real(RK), parameter     :: CovMatUpper_ref(nd,nd) = reshape([ +0.140203078481164_RK, +0.029035771029083_RK, -0.031754793421423_RK &
                                                                    , +0.029035771029083_RK, +0.092838466390969_RK, -0.043469498536710_RK &
                                                                    , -0.031754793421423_RK, -0.043469498536710_RK, +0.116582891117029_RK &
                                                                    ] , shape = shape(CovMatUpper_ref) )
        real(RK)                :: CovMatUpper(nd,nd)
        real(RK)                :: CovMatUpper_diff(nd,nd)
        real(RK)                :: Mean_diff(nd)
        real(RK)                :: Mean(nd)

        assertion = .true.

        CovMatUpper = CovMatUpper_ref ! to avoid un-initialization error in the computation of `CovMatUpper_diff`.
        call getSamCovUpperMeanTrans(np,nd,Point,CovMatUpper,Mean)
        CovMatUpper_diff = abs(CovMatUpper - CovMatUpper_ref) / abs(CovMatUpper_ref)
        assertionCurrent = all([ ((CovMatUpper_diff(i,j) <= tolerance, i=1,j), j=1,nd) ])
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper_ref  ", ((CovMatUpper_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper      ", ((CovMatUpper(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper_diff ", ((CovMatUpper_diff(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Mean_diff = abs(Mean - Mean_ref) / abs(Mean_ref)
        assertionCurrent = all(Mean_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref  ", Mean_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean      ", Mean
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_diff ", Mean_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSamCovUpperMeanTrans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWeiSamCovUppMeanTrans_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        integer(IK)             :: i, j
        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 5_IK
        integer(IK) , parameter :: Weight(np) = [(i,i=1,np)]
        real(RK)    , parameter :: Point(nd,np) = reshape(  [ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(Point) )
        real(RK), parameter     :: Mean_ref(nd) = [.4565495333977344_RK, .4605581195561738_RK, .4823401092849646_RK]
        real(RK), parameter     :: CovMatUpper_ref(nd,nd) = reshape([ +.1256706333432408E+00_RK, +.4589762909514189E-01_RK, -.2021812246706341E-01_RK &
                                                                    , +.4589762909514189E-01_RK, +.7653806291082801E-01_RK, -.6048750626649663E-01_RK &
                                                                    , -.2021812246706341E-01_RK, -.6048750626649663E-01_RK, +.1006280226405955E+00_RK &
                                                                    ] , shape = shape(CovMatUpper_ref) )
        real(RK)                :: CovMatUpper(nd,nd)
        real(RK)                :: CovMatUpper_diff(nd,nd)
        real(RK)                :: Mean_diff(nd)
        real(RK)                :: Mean(nd)

        assertion = .true.

        call getWeiSamCovUppMeanTrans(np,sum(Weight),nd,Point,Weight,CovMatUpper,Mean)
        assertionCurrent = .true.
        do j = 1, nd
            do i = 1, j
                CovMatUpper_diff(i,j) = abs(CovMatUpper(i,j) - CovMatUpper_ref(i,j)) / abs(CovMatUpper_ref(i,j))
                assertionCurrent = assertionCurrent .and. CovMatUpper_diff(i,j) <= tolerance
            end do
        end do
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper_ref  ", ((CovMatUpper_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper      ", ((CovMatUpper(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpper_diff ", ((CovMatUpper_diff(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Mean_diff = abs(Mean - Mean_ref) / abs(Mean_ref)
        assertionCurrent = all(Mean_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_ref  ", Mean_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean      ", Mean
            write(Test%outputUnit,"(*(g0,:,', '))") "Mean_diff ", Mean_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getWeiSamCovUppMeanTrans_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_mergeMeanCov_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        integer(IK)             :: i, j
        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK

        integer(IK) , parameter :: npA = 5_IK
        real(RK)                :: MeanA(nd), CovMatA(nd,nd)
        real(RK)    , parameter :: PointA(nd,npA) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(PointA) )

        integer(IK) , parameter :: npB = 10_IK
        real(RK)                :: MeanB(nd), CovMatB(nd,nd)
        real(RK)    , parameter :: PointB(nd,npB) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            , 1.706046088019609_RK, 2.031832846377421_RK, 3.276922984960890_RK &
                                                            , 1.046171390631154_RK, 2.097131781235848_RK, 3.823457828327293_RK &
                                                            , 1.694828622975817_RK, 2.317099480060861_RK, 3.950222048838355_RK &
                                                            , 1.034446080502909_RK, 2.438744359656398_RK, 3.381558457093008_RK &
                                                            , 1.765516788149002_RK, 2.795199901137063_RK, 3.186872604554379_RK &
                                                            ], shape = shape(PointB) )


        integer(IK) , parameter :: npAB = npA + npB
        real(RK)                :: PointAB(nd,npAB)
        real(RK)                :: CovMatAB_diff(nd,nd)
        real(RK)    , parameter :: CovMatAB_refMATLAB(nd,nd) = reshape( [ 0.358269305364807_RK, 0.501078279929690_RK, 0.687067319924494_RK &
                                                                        , 0.501078279929690_RK, 1.031956780716069_RK, 1.391311858397106_RK &
                                                                        , 0.687067319924494_RK, 1.391311858397106_RK, 2.242785335243168_RK &
                                                                        ], shape = shape(CovMatAB_refMATLAB) )
        real(RK)                :: CovMatAB_ref(nd,nd)
        real(RK)                :: CovMatAB(nd,nd)
        real(RK)                :: MeanAB_diff(nd)
        real(RK)                :: MeanAB_ref(nd)
        real(RK)                :: MeanAB(nd)

        assertion = .true.

        call getSamCovMeanTrans(npA,nd,PointA,CovMatA,MeanA)
        call getSamCovMeanTrans(npB,nd,PointB,CovMatB,MeanB)

        PointAB(:,1:npA) = PointA
        PointAB(:,npA+1:npAB) = PointB

        call getSamCovMeanTrans(npAB,nd,PointAB,CovMatAB_ref,MeanAB_ref)

        !do i = 1, 10000000
        call mergeMeanCov(nd,npA,MeanA,CovMatA,npB,MeanB,CovMatB,MeanAB,CovMatAB)
        !end do

        CovMatAB_diff = abs(CovMatAB - CovMatAB_ref) / abs(CovMatAB_ref)
        assertionCurrent = all([ ((CovMatAB_diff(i,j) <= tolerance, i=1,j), j=1,nd) ])
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB_refMATLAB ", ((CovMatAB_refMATLAB(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB_ref       ", ((CovMatAB_ref(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB_diff      ", ((CovMatAB_refMATLAB(i,j)-CovMatAB_ref(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB_ref       ", ((CovMatAB_ref(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB           ", ((CovMatAB(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatAB_diff      ", ((CovMatAB_diff(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        MeanAB_diff = abs(MeanAB - MeanAB_ref) / abs(MeanAB_ref)
        assertionCurrent = all(MeanAB_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_ref  ", MeanAB_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB      ", MeanAB
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_diff ", MeanAB_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_mergeMeanCov_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_mergeMeanCovUpperSlow_1() result(assertion)
!
!        use Constants_mod, only: IK, RK
!        implicit none
!
!        integer(IK)             :: i, j
!        logical                 :: assertion, assertionCurrent
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK) , parameter :: nd = 3_IK
!
!        integer(IK) , parameter :: npA = 5_IK
!        real(RK)                :: MeanA(nd), CovMatUpperA(nd,nd)
!        real(RK)    , parameter :: PointA(nd,npA) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
!                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
!                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
!                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
!                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
!                                                            ], shape = shape(PointA) )
!
!        integer(IK) , parameter :: npB = 10_IK
!        real(RK)                :: MeanB(nd), CovMatUpperB(nd,nd)
!        real(RK)    , parameter :: PointB(nd,npB) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
!                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
!                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
!                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
!                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
!                                                            , 1.706046088019609_RK, 2.031832846377421_RK, 3.276922984960890_RK &
!                                                            , 1.046171390631154_RK, 2.097131781235848_RK, 3.823457828327293_RK &
!                                                            , 1.694828622975817_RK, 2.317099480060861_RK, 3.950222048838355_RK &
!                                                            , 1.034446080502909_RK, 2.438744359656398_RK, 3.381558457093008_RK &
!                                                            , 1.765516788149002_RK, 2.795199901137063_RK, 3.186872604554379_RK &
!                                                            ], shape = shape(PointB) )
!
!
!        integer(IK) , parameter :: npAB = npA + npB
!        real(RK)                :: PointAB(nd,npAB)
!        real(RK)                :: CovMatUpperAB_diff(nd,nd)
!        real(RK)    , parameter :: CovMatUpperAB_refMATLAB(nd,nd) = reshape([ 0.358269305364807_RK, 0.501078279929690_RK, 0.687067319924494_RK &
!                                                                            , 0.501078279929690_RK, 1.031956780716069_RK, 1.391311858397106_RK &
!                                                                            , 0.687067319924494_RK, 1.391311858397106_RK, 2.242785335243168_RK &
!                                                                            ], shape = shape(CovMatUpperAB_refMATLAB) )
!        real(RK)                :: CovMatUpperAB_ref(nd,nd)
!        real(RK)                :: CovMatUpperAB(nd,nd)
!        real(RK)                :: MeanAB_diff(nd)
!        real(RK)                :: MeanAB_ref(nd)
!        real(RK)                :: MeanAB(nd)
!
!        assertion = .true.
!
!        call getSamCovUpperMeanTrans(npA,nd,PointA,CovMatUpperA,MeanA)
!        call getSamCovUpperMeanTrans(npB,nd,PointB,CovMatUpperB,MeanB)
!
!        PointAB(:,1:npA) = PointA
!        PointAB(:,npA+1:npAB) = PointB
!
!        call getSamCovUpperMeanTrans(npAB,nd,PointAB,CovMatUpperAB_ref,MeanAB_ref)
!
!        !do i = 1, 10000000
!        call mergeMeanCovUpperSlow(nd,npA,MeanA,CovMatUpperA,npB,MeanB,CovMatUpperB,MeanAB,CovMatUpperAB)
!        !end do
!        assertionCurrent = .true.
!        do j = 1, nd
!            do i = 1, j
!                CovMatUpperAB_diff(i,j) = abs(CovMatUpperAB(i,j) - CovMatUpperAB_ref(i,j)) / abs(CovMatUpperAB_ref(i,j))
!                assertionCurrent = assertionCurrent .and. CovMatUpperAB_diff(i,j) <= tolerance
!            end do
!        end do
!        assertion = assertion .and. assertionCurrent
!
!        LCOV_EXCL_START
!        if (Test%isVerboseMode .and. .not. assertionCurrent) then
!            write(Test%outputUnit,"(*(g0,:,', '))")
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_refMATLAB ", ((CovMatUpperAB_refMATLAB(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_refMATLAB(i,j)-CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))")
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB           ", ((CovMatUpperAB(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_diff(i,j), i=1,j), j=1,nd)
!            write(Test%outputUnit,"(*(g0,:,', '))")
!        end if
!        LCOV_EXCL_STOP
!
!        MeanAB_diff = abs(MeanAB - MeanAB_ref) / abs(MeanAB_ref)
!        assertionCurrent = all(MeanAB_diff <= tolerance)
!        assertion = assertion .and. assertionCurrent
!
!        LCOV_EXCL_START
!        if (Test%isVerboseMode .and. .not. assertionCurrent) then
!            write(Test%outputUnit,"(*(g0,:,', '))")
!            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_ref  ", MeanAB_ref
!            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB      ", MeanAB
!            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_diff ", MeanAB_diff
!            write(Test%outputUnit,"(*(g0,:,', '))")
!        end if
!        LCOV_EXCL_STOP
!
!    end function test_mergeMeanCovUpperSlow_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_mergeMeanCovUpper_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        integer(IK)             :: i, j
        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK

        integer(IK) , parameter :: npA = 5_IK
        real(RK)                :: MeanA(nd), CovMatUpperA(nd,nd)
        real(RK)    , parameter :: PointA(nd,npA) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(PointA) )

        integer(IK) , parameter :: npB = 10_IK
        real(RK)                :: MeanB(nd), CovMatUpperB(nd,nd)
        real(RK)    , parameter :: PointB(nd,npB) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            , 1.706046088019609_RK, 2.031832846377421_RK, 3.276922984960890_RK &
                                                            , 1.046171390631154_RK, 2.097131781235848_RK, 3.823457828327293_RK &
                                                            , 1.694828622975817_RK, 2.317099480060861_RK, 3.950222048838355_RK &
                                                            , 1.034446080502909_RK, 2.438744359656398_RK, 3.381558457093008_RK &
                                                            , 1.765516788149002_RK, 2.795199901137063_RK, 3.186872604554379_RK &
                                                            ], shape = shape(PointB) )


        integer(IK) , parameter :: npAB = npA + npB
        real(RK)                :: PointAB(nd,npAB)
        real(RK)                :: CovMatUpperAB_diff(nd,nd)
        real(RK)    , parameter :: CovMatUpperAB_refMATLAB(nd,nd) = reshape([ 0.358269305364807_RK, 0.501078279929690_RK, 0.687067319924494_RK &
                                                                            , 0.501078279929690_RK, 1.031956780716069_RK, 1.391311858397106_RK &
                                                                            , 0.687067319924494_RK, 1.391311858397106_RK, 2.242785335243168_RK &
                                                                            ], shape = shape(CovMatUpperAB_refMATLAB) )
        real(RK)                :: CovMatUpperAB_ref(nd,nd)
        real(RK)                :: CovMatUpperAB(nd,nd)
        real(RK)                :: MeanAB_diff(nd)
        real(RK)                :: MeanAB_ref(nd)
        real(RK)                :: MeanAB(nd)

        assertion = .true.

        call getSamCovUpperMeanTrans(npA,nd,PointA,CovMatUpperA,MeanA)
        call getSamCovUpperMeanTrans(npB,nd,PointB,CovMatUpperB,MeanB)

        PointAB(:,1:npA) = PointA
        PointAB(:,npA+1:npAB) = PointB

        call getSamCovUpperMeanTrans(npAB,nd,PointAB,CovMatUpperAB_ref,MeanAB_ref)

        !do i = 1, 10000000
        call mergeMeanCovUpper(nd,npA,MeanA,CovMatUpperA,npB,MeanB,CovMatUpperB,MeanAB,CovMatUpperAB)
        !end do

        assertionCurrent = .true.
        do j = 1, nd
            do i = 1, j
                CovMatUpperAB_diff(i,j) = abs(CovMatUpperAB(i,j) - CovMatUpperAB_ref(i,j)) / abs(CovMatUpperAB_ref(i,j))
                assertionCurrent = assertionCurrent .and. CovMatUpperAB_diff(i,j) <= tolerance
            end do
        end do
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_refMATLAB ", ((CovMatUpperAB_refMATLAB(i,j), i=1,nd), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_refMATLAB(i,j)-CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB           ", ((CovMatUpperAB(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_diff(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        MeanAB_diff = abs(MeanAB - MeanAB_ref) / abs(MeanAB_ref)
        assertionCurrent = all(MeanAB_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_ref  ", MeanAB_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB      ", MeanAB
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_diff ", MeanAB_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_mergeMeanCovUpper_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_mergeMeanCovUpperDense_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none

        integer(IK)             :: i, j
        logical                 :: assertion, assertionCurrent
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        integer(IK) , parameter :: nd = 3_IK

        integer(IK) , parameter :: npA = 5_IK
        real(RK)                :: MeanA(nd), CovMatUpperA(nd,nd)
        real(RK)    , parameter :: PointA(nd,npA) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            ], shape = shape(PointA) )

        integer(IK) , parameter :: npB = 10_IK
        real(RK)                :: MeanB(nd), CovMatUpperB(nd,nd)
        real(RK)    , parameter :: PointB(nd,npB) = reshape([ 0.706046088019609_RK, 0.031832846377421_RK, 0.276922984960890_RK &
                                                            , 0.046171390631154_RK, 0.097131781235848_RK, 0.823457828327293_RK &
                                                            , 0.694828622975817_RK, 0.317099480060861_RK, 0.950222048838355_RK &
                                                            , 0.034446080502909_RK, 0.438744359656398_RK, 0.381558457093008_RK &
                                                            , 0.765516788149002_RK, 0.795199901137063_RK, 0.186872604554379_RK &
                                                            , 1.706046088019609_RK, 2.031832846377421_RK, 3.276922984960890_RK &
                                                            , 1.046171390631154_RK, 2.097131781235848_RK, 3.823457828327293_RK &
                                                            , 1.694828622975817_RK, 2.317099480060861_RK, 3.950222048838355_RK &
                                                            , 1.034446080502909_RK, 2.438744359656398_RK, 3.381558457093008_RK &
                                                            , 1.765516788149002_RK, 2.795199901137063_RK, 3.186872604554379_RK &
                                                            ], shape = shape(PointB) )


        integer(IK) , parameter :: npAB = npA + npB
        real(RK)                :: PointAB(nd,npAB)
        real(RK)                :: CovMatUpperAB_diff(nd,nd)
        real(RK)                :: MeanAB(nd), CovMatUpperAB(nd,nd)
        real(RK)    , parameter :: CovMatUpperAB_refMATLAB(nd,nd) = reshape([ 0.358269305364807_RK, 0.501078279929690_RK, 0.687067319924494_RK &
                                                                            , 0.501078279929690_RK, 1.031956780716069_RK, 1.391311858397106_RK &
                                                                            , 0.687067319924494_RK, 1.391311858397106_RK, 2.242785335243168_RK &
                                                                            ], shape = shape(CovMatUpperAB_refMATLAB) )
        real(RK)                :: CovMatUpperAB_ref(nd,nd)
        real(RK)                :: MeanAB_diff(nd)
        real(RK)                :: MeanAB_ref(nd)

        assertion = .true.

        call getSamCovUpperMeanTrans(npA,nd,PointA,CovMatUpperA,MeanA)
        call getSamCovUpperMeanTrans(npB,nd,PointB,CovMatUpperB,MeanB)

        PointAB(:,1:npA) = PointA
        PointAB(:,npA+1:npAB) = PointB

        call getSamCovUpperMeanTrans(npAB,nd,PointAB,CovMatUpperAB_ref,MeanAB_ref)

        MeanAB = MeanB
        CovMatUpperAB = CovMatUpperB
        !do i = 1, 10000000-1
        call mergeMeanCovUpperDense(nd,npA,MeanA,CovMatUpperA,npB,MeanAB,CovMatUpperAB)
        !end do

        MeanAB = MeanB
        CovMatUpperAB = CovMatUpperB
        call mergeMeanCovUpperDense(nd,npA,MeanA,CovMatUpperA,npB,MeanAB,CovMatUpperAB)

        assertionCurrent = .true.
        do j = 1, nd
            do i = 1, j
                CovMatUpperAB_diff(i,j) = abs(CovMatUpperAB(i,j) - CovMatUpperAB_ref(i,j)) / abs(CovMatUpperAB_ref(i,j))
                assertionCurrent = assertionCurrent .and. CovMatUpperAB_diff(i,j) <= tolerance
            end do
        end do
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_refMATLAB ", ((CovMatUpperAB_refMATLAB(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_refMATLAB(i,j)-CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_ref       ", ((CovMatUpperAB_ref(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB           ", ((CovMatUpperAB(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))") "CovMatUpperAB_diff      ", ((CovMatUpperAB_diff(i,j), i=1,j), j=1,nd)
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        MeanAB_diff = abs(MeanAB - MeanAB_ref) / abs(MeanAB_ref)
        assertionCurrent = all(MeanAB_diff <= tolerance)
        assertion = assertion .and. assertionCurrent

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertionCurrent) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_ref  ", MeanAB_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB      ", MeanAB
            write(Test%outputUnit,"(*(g0,:,', '))") "MeanAB_diff ", MeanAB_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_mergeMeanCovUpperDense_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandGaus_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        real(RK)                :: stdNormRnd
        stdNormRnd = getRandGaus()
        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getRandGaus_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandNorm_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        real(RK)                :: normRnd
        normRnd = getRandNorm(0._RK, 1._RK)
        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getRandNorm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandLogn_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        real(RK)                :: lognRnd
        lognRnd = getRandLogn(0._RK, 1._RK)
        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getRandLogn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getMVNDev_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)                :: X(nd)
        call getMVNDev(nd, MeanVec, CovMat, X)
        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getMVNDev_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMVUDev_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape( [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                            , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                            , shape = shape(InvCovMat) )
        real(RK)                :: X(nd), NormedPoint(nd)
        call getMVUDev(nd, MeanVec, CovMat, X)
        NormedPoint = X - MeanVec
        assertion = isInsideEllipsoid(nd, NormedPoint, InvCovMat)
    end function test_getMVUDev_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandMVN_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: CholeskyLower(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: CholeskyDiago(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)                :: X(nd)
        X = getRandMVN(nd, MeanVec, CholeskyLower, CholeskyDiago)
        assertion = .true. !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getRandMVN_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandMVU_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: CholeskyDiago(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)    , parameter :: CholeskyLower(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape( [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                            , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                            , shape = shape(InvCovMat) )
        real(RK)                :: X(nd), NormedPoint(nd)
        X = getRandMVU(nd, MeanVec, CholeskyLower, CholeskyDiago)
        NormedPoint = X - MeanVec
        assertion = isInsideEllipsoid(nd, NormedPoint, InvCovMat)
    end function test_getRandMVU_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInsideEllipsoid_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: CholeskyDiago(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)    , parameter :: CholeskyLower(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape( [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                            , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                            , shape = shape(InvCovMat) )
        real(RK)                :: X(nd), NormedPoint(nd)
        X = getRandMVU(nd, MeanVec, CholeskyLower, CholeskyDiago)
        NormedPoint = X - MeanVec
        assertion = isInsideEllipsoid(nd, NormedPoint, InvCovMat)
        NormedPoint = [-1.e2_RK, 1.e2_RK]
        assertion = assertion .and. .not. isInsideEllipsoid(nd, NormedPoint, InvCovMat)
    end function test_isInsideEllipsoid_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVU_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: logSqrtDetCovMat = log(1._RK)
        real(RK)    , parameter :: logProbMVU_ref = -1.144729885849400_RK
        real(RK)                :: difference
        real(RK)                :: logProbMVU
        logProbMVU = getLogProbMVU(nd, logSqrtDetCovMat)
        difference = abs( (logProbMVU - logProbMVU_ref) / logProbMVU_ref)
        assertion = difference < tolerance
        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMVU_ref ", logProbMVU_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbMVU     ", logProbMVU
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogProbMVU_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandPointOnEllipsoid_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: MeanVec(nd) = [ 1._RK, 2._RK ]
        real(RK)    , parameter :: CovMat(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: CholeskyDiago(nd) = [ 1._RK, 0.866025403784439_RK ]
        real(RK)    , parameter :: CholeskyLower(nd,nd) = reshape( [ 1._RK, 0.5_RK, 0.5_RK, 1._RK ], shape = shape(CovMat) )
        real(RK)    , parameter :: InvCovMat(nd,nd) = reshape( [ +1.333333333333333_RK, -0.666666666666667_RK &
                                                            , -0.666666666666667_RK, +1.333333333333333_RK ] &
                                                            , shape = shape(InvCovMat) )
        real(RK)                :: X(nd), NormedPoint(nd)
        X = getRandPointOnEllipsoid(nd,MeanVec,CholeskyLower,CholeskyDiago)
        NormedPoint = X - MeanVec
        assertion = dot_product(NormedPoint,matmul(InvCovMat,NormedPoint)) - 1._RK < tolerance
        NormedPoint = [-1.e2_RK, 1.e2_RK]
        assertion = assertion .and. .not. isInsideEllipsoid(nd, NormedPoint, InvCovMat)
        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "RandPointOnEllipsoid   ", X
            write(Test%outputUnit,"(*(g0,:,', '))") "distance from center   ", dot_product(NormedPoint,matmul(InvCovMat,NormedPoint))
            write(Test%outputUnit,"(*(g0,:,', '))") "expected distance      ", 1._RK
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getRandPointOnEllipsoid_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbLognSP_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: logMean = 2._RK
        real(RK)    , parameter :: logPoint = log(5._RK)
        real(RK)    , parameter :: inverseVariance = 1.e-2_RK
        real(RK)    , parameter :: logSqrtInverseVariance = log(sqrt(inverseVariance))
        real(RK)    , parameter :: logProbLogn_ref = -4.831724232354038_RK
        real(RK)                :: difference
        real(RK)                :: logProbLogn

        logProbLogn = getLogProbLogn(logMean = logMean, inverseVariance = inverseVariance, logSqrtInverseVariance = logSqrtInverseVariance, logPoint = logPoint)
        difference = abs( (logProbLogn - logProbLogn_ref) / logProbLogn_ref )
        assertion = difference < tolerance

        logProbLogn = getLogProbLognorm(logMean = logMean, inverseVariance = inverseVariance, logSqrtInverseVariance = logSqrtInverseVariance, logPoint = logPoint)
        difference = abs( (logProbLogn - logProbLogn_ref) / logProbLogn_ref )
        assertion = assertion .and. difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbLogn_ref    ", logProbLogn_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "logProbLogn        ", logProbLogn
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbLognSP_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbLognMP_1() result(assertion)

        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK, np = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: logMean = 2._RK
        real(RK)    , parameter :: LogPoint(np) = log([5._RK, 6._RK, 7._RK])
        real(RK)    , parameter :: inverseVariance = 1.e-2_RK
        real(RK)    , parameter :: logSqrtInverseVariance = log(sqrt(inverseVariance))
        real(RK)    , parameter :: LogProbLogn_ref(np) = [ -4.831724232354038_RK, -5.013499916020054_RK, -5.167448403813908_RK ]
        real(RK)                :: LogProbLogn(np)
        real(RK)                :: difference(np)

        LogProbLogn = getLogProbLogn(np = np, logMean = logMean, inverseVariance = inverseVariance, logSqrtInverseVariance = logSqrtInverseVariance, logPoint = LogPoint)
        difference = abs( (LogProbLogn - LogProbLogn_ref) / LogProbLogn_ref )
        assertion = all(difference < tolerance)

        LogProbLogn = getLogProbLognorm(np = np, logMean = logMean, inverseVariance = inverseVariance, logSqrtInverseVariance = logSqrtInverseVariance, logPoint = LogPoint)
        difference = abs( (LogProbLogn - LogProbLogn_ref) / LogProbLogn_ref )
        assertion = assertion .and. all(difference < tolerance)

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbLogn_ref    ", LogProbLogn_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "LogProbLogn        ", LogProbLogn
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbLognMP_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandRealLecuyer_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        integer(IK)             :: idum = 3333, i
        real(RK)                :: RandRealLecuyer(np)
        assertion = .true.
        do i = 1, np
            RandRealLecuyer(i) = getRandRealLecuyer(idum)
            assertion = assertion .and. RandRealLecuyer(i) <= 1._RK .and. RandRealLecuyer(i) >= 0._RK
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandRealLecuyer(",i,") =", RandRealLecuyer(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandRealLecuyer_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandIntLecuyer_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: lowerBound = -2, upperBound = 9
        integer(IK) , parameter :: np = 100
        integer(IK)             :: idum = 3333, i
        integer(IK)             :: RandIntLecuyer(np)
        assertion = .true.
        do i = 1, np
            RandIntLecuyer(i) = getRandIntLecuyer(lowerBound,upperBound,idum)
            assertion = assertion .and. RandIntLecuyer(i) <= upperBound .and. RandIntLecuyer(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandIntLecuyer(",i,") =", RandIntLecuyer(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandIntLecuyer_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandUniform_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: lowerBound = -2._RK, upperBound = 9._RK
        integer(IK)             :: i
        real(RK)                :: RandUniform(np)
        assertion = .true.
        do i = 1, np
            RandUniform(i) = getRandUniform(lowerBound,upperBound)
            assertion = assertion .and. RandUniform(i) <= upperBound .and. RandUniform(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandUniform(",i,") =", RandUniform(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandUniform_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandGamma_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: alpha = 2._RK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = HUGE_RK
        integer(IK)             :: i
        real(RK)                :: RandGamma(np)
        assertion = .true.
        do i = 1, np
            RandGamma(i) = getRandGamma(alpha)
            assertion = assertion .and. RandGamma(i) <= upperBound .and. RandGamma(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandGamma(",i,") =", RandGamma(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandGamma(alpha=-1._RK) < 0._RK
    end function test_getRandGamma_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandGammaIntShape_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        integer(IK) , parameter :: alpha = 2_IK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = HUGE_RK
        integer(IK)             :: i
        real(RK)                :: RandGamma(np)
        assertion = .true.
        do i = 1, np
            RandGamma(i) = getRandGammaIntShape(alpha)
            assertion = assertion .and. RandGamma(i) <= upperBound .and. RandGamma(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandGamma(",i,") =", RandGamma(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandGammaIntShape(alpha=-1_IK) < 0._RK
    end function test_getRandGammaIntShape_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandBeta_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: alpha = 2._RK, beta = 3._RK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = 1._RK
        integer(IK)             :: i
        real(RK)                :: RandBeta(np)
        assertion = .true.
        do i = 1, np
            RandBeta(i) = getRandBeta(alpha, beta)
            assertion = assertion .and. RandBeta(i) <= upperBound .and. RandBeta(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandBeta(",i,") =", RandBeta(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandBeta(alpha=-1._RK, beta=+2._RK) < 0._RK
        assertion = assertion .and. getRandBeta(alpha=+1._RK, beta=-2._RK) < 0._RK
        assertion = assertion .and. getRandBeta(alpha=-1._RK, beta=-2._RK) < 0._RK
    end function test_getRandBeta_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandExpWithInvMean_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: invMean = 0.1_RK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = HUGE_RK
        integer(IK)             :: i
        real(RK)                :: RandExp(np)
        assertion = .true.
        do i = 1, np
            RandExp(i) = getRandExpWithInvMean(invMean)
            assertion = assertion .and. RandExp(i) <= upperBound .and. RandExp(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandExp(",i,") =", RandExp(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandExpWithInvMean_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandExp_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = HUGE_RK
        integer(IK)             :: i
        real(RK)                :: RandExp(np)
        assertion = .true.
        do i = 1, np
            RandExp(i) = getRandExp()
            assertion = assertion .and. RandExp(i) <= upperBound .and. RandExp(i) >= lowerBound
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandExp(",i,") =", RandExp(i)
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandExp_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandCorMat_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        use Matrix_mod, only: isPosDef
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        integer(IK) , parameter :: np = 100_IK
        real(RK)    , parameter :: eta = 5._RK
        real(RK)    , parameter :: lowerBound = -1._RK, upperBound = 1._RK
        integer(IK)             :: i
        real(RK)                :: RandCorMat(nd,nd)
        assertion = .true.
        do i = 1, np
            RandCorMat = getRandCorMat(nd,eta)
            assertion = assertion .and. all(RandCorMat <= upperBound) .and. all(RandCorMat >= lowerBound)
            assertion = assertion .and. isPosDef(nd,RandCorMat)
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandCorMat(:,:,) =", RandCorMat
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        ! do not uncomment this as GNU Fortran 9.1 crashes on this in debug mode with error message:
        ! Fortran runtime error: Dimension 1 of array 'randcormat' has extent 0 instead of 2...
        !RandCorMat = getRandCorMat(nd = 0_IK, eta = 1._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        !RandCorMat = getRandCorMat(nd = -1_IK, eta = +2._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        !RandCorMat = getRandCorMat(nd = +1_IK, eta = +0._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        RandCorMat = getRandCorMat(nd = +2_IK, eta = +0._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        RandCorMat = getRandCorMat(nd = +2_IK, eta = -1._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        RandCorMat = getRandCorMat(nd = +2_IK, eta = -2._RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
    end function test_getRandCorMat_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandCorMatRejection_1() result(assertion)
        use Constants_mod, only: IK, RK, HUGE_RK
        use Matrix_mod, only: isPosDef
        implicit none
        logical                 :: assertion, assertionCurrent
        integer(IK) , parameter :: nd = 2_IK
        integer(IK) , parameter :: np = 100_IK
        real(RK)    , parameter :: minRho = -.3_RK, maxRho = .6_RK
        real(RK)    , parameter :: lowerBound = -1._RK, upperBound = 1._RK
        integer(IK)             :: i, j, k
        real(RK)                :: RandCorMat(nd,nd)
        assertion = .true.
        do i = 1, np
            RandCorMat = getRandCorMatRejection(nd,minRho,maxRho)
            assertion = assertion .and. all(RandCorMat <= upperBound) .and. all(RandCorMat >= lowerBound)
            do j = 1, nd
                do k = 1, nd
                    if (j==k) then
                        assertionCurrent = RandCorMat(j,k) == upperBound
                    else
                        assertionCurrent = RandCorMat(j,k) <= maxRho .and. RandCorMat(j,k) >= minRho
                    end if
                    assertion = assertion .and. assertionCurrent
                end do
            end do
            assertionCurrent = isPosDef(nd,RandCorMat)
            ! LCOV_EXCL_START
            if (Test%isVerboseMode .and. .not. assertion) then
                write(Test%outputUnit,"(*(g0,:,' '))")
                write(Test%outputUnit,"(*(g0,:,' '))") "RandCorMat(:,:,) =", RandCorMat
                write(Test%outputUnit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
            assertion = assertion .and. assertionCurrent
        end do
        ! do not uncomment this as GNU Fortran 9.1 crashes on this in debug mode with error message:
        ! Fortran runtime error: Dimension 1 of array 'randcormat' has extent 1 instead of 2...
        !RandCorMat = getRandCorMatRejection(nd = +1_IK, minRho = +.2_RK, maxRho = -.5_RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        RandCorMat = getRandCorMatRejection(nd = +2_IK, minRho = -.2_RK, maxRho = -.5_RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        RandCorMat = getRandCorMatRejection(nd = +2_IK, minRho = +.2_RK, maxRho = -.5_RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        !RandCorMat = getRandCorMatRejection(nd = +0_IK, minRho = +.2_RK, maxRho = +.5_RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
        !RandCorMat = getRandCorMatRejection(nd = -3_IK, minRho = +.2_RK, maxRho = +.5_RK); assertion = assertion .and. RandCorMat(1,1) < 0._RK
    end function test_getRandCorMatRejection_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCorMatUpperFromCovMatUpper_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMatUpper    (nd,nd) = reshape( [ 2._RK, 0.00_RK, 0.50_RK, 2._RK ], shape = shape(CovMatUpper) )
        real(RK)    , parameter :: CorMatUpper_ref(nd,nd) = reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatUpper_ref) )
        integer(IK)             :: i, j
        real(RK)                :: CorMatUpper(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CorMatUpper = getCorMatUpperFromCovMatUpper(nd, CovMatUpper)

        assertion = .true.
        do j = 1, nd
            do i = 1, j
                Difference(i,j) = abs(CorMatUpper(i,j) - CorMatUpper_ref(i,j)) / abs(CorMatUpper_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CorMatUpper_ref =", CorMatUpper_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CorMatUpper     =", CorMatUpper
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCorMatUpperFromCovMatUpper_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCovMatUpperFromCorMatUpper_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMatUpper_ref(nd,nd) = reshape( [ 2._RK, 0.00_RK, 0.50_RK, 2._RK ], shape = shape(CovMatUpper_ref) )
        real(RK)    , parameter :: CorMatUpper(nd,nd) = reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatUpper) )
        real(RK)    , parameter :: StdVec(nd) = [ 1.414213562373095_RK, 1.414213562373095_RK ]
        integer(IK)             :: i, j
        real(RK)                :: CovMatUpper(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CovMatUpper = getCovMatUpperFromCorMatUpper(nd, StdVec, CorMatUpper)

        assertion = .true.
        do j = 1, nd
            do i = 1, j
                Difference(i,j) = abs(CovMatUpper(i,j) - CovMatUpper_ref(i,j)) / abs(CovMatUpper_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        ! LCOV_EXCL_START
        if (Test%isVerboseMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatUpper_ref =", CovMatUpper_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatUpper     =", CovMatUpper
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCovMatUpperFromCorMatUpper_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCovMatUpperFromCorMatLower_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMatUpper_ref(nd,nd) = reshape( [ 2._RK, 0.00_RK, 0.50_RK, 2._RK ], shape = shape(CovMatUpper_ref) )
        real(RK)    , parameter :: CorMatLower(nd,nd) = transpose(reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatLower) ))
        real(RK)    , parameter :: StdVec(nd) = [ 1.414213562373095_RK, 1.414213562373095_RK ]
        integer(IK)             :: i, j
        real(RK)                :: CovMatUpper(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CovMatUpper = getCovMatUpperFromCorMatLower(nd, StdVec, CorMatLower)

        assertion = .true.
        do j = 1, nd
            do i = 1, j
                Difference(i,j) = abs(CovMatUpper(i,j) - CovMatUpper_ref(i,j)) / abs(CovMatUpper_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatUpper_ref =", CovMatUpper_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatUpper     =", CovMatUpper
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCovMatUpperFromCorMatLower_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCovMatLowerFromCorMatUpper_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMatLower_ref(nd,nd) = transpose( reshape( [ 2._RK, 0.00_RK, 0.50_RK, 2._RK ], shape = shape(CovMatLower_ref) ) )
        real(RK)    , parameter :: CorMatUpper(nd,nd) = reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatUpper) )
        real(RK)    , parameter :: StdVec(nd) = [ 1.414213562373095_RK, 1.414213562373095_RK ]
        integer(IK)             :: i, j
        real(RK)                :: CovMatLower(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CovMatLower = getCovMatLowerFromCorMatUpper(nd, StdVec, CorMatUpper)

        assertion = .true.
        do j = 1, nd
            do i = j, nd
                Difference(i,j) = abs(CovMatLower(i,j) - CovMatLower_ref(i,j)) / abs(CovMatLower_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatLower_ref =", CovMatLower_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatLower     =", CovMatLower
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCovMatLowerFromCorMatUpper_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCovMatLowerFromCorMatLower_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMatLower_ref(nd,nd) = transpose( reshape( [ 2._RK, 0.00_RK, 0.50_RK, 2._RK ], shape = shape(CovMatLower_ref) ) )
        real(RK)    , parameter :: CorMatLower(nd,nd) = transpose( reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatLower) ) )
        real(RK)    , parameter :: StdVec(nd) = [ 1.414213562373095_RK, 1.414213562373095_RK ]
        integer(IK)             :: i, j
        real(RK)                :: CovMatLower(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CovMatLower = getCovMatLowerFromCorMatLower(nd, StdVec, CorMatLower)

        assertion = .true.
        do j = 1, nd
            do i = j, nd
                Difference(i,j) = abs(CovMatLower(i,j) - CovMatLower_ref(i,j)) / abs(CovMatLower_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatLower_ref =", CovMatLower_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMatLower     =", CovMatLower
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCovMatLowerFromCorMatLower_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCovMatFromCorMatUpper_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: CovMat_ref(nd,nd) = reshape( [ 2._RK, 0.50_RK, 0.50_RK, 2._RK ], shape = shape(CovMat_ref) )
        real(RK)    , parameter :: CorMatUpper(nd,nd) = reshape( [ 1._RK, 0.00_RK, 0.25_RK, 1._RK ], shape = shape(CorMatUpper) )
        real(RK)    , parameter :: StdVec(nd) = [ 1.414213562373095_RK, 1.414213562373095_RK ]
        integer(IK)             :: i, j
        real(RK)                :: CovMat(nd,nd)
        real(RK)                :: Difference(nd,nd)

        CovMat = getCovMatFromCorMatUpper(nd, StdVec, CorMatUpper)

        assertion = .true.
        do j = 1, nd
            do i = 1, nd
                Difference(i,j) = abs(CovMat(i,j) - CovMat_ref(i,j)) / abs(CovMat_ref(i,j))
                assertion = assertion .and. Difference(i,j) < tolerance
            end do
        end do

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMat_ref  =", CovMat_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "CovMat      =", CovMat
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference  =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getCovMatFromCorMatUpper_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbGeo_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: numTrial = 10_IK
        integer(IK) , parameter :: SuccessStep(numTrial) = [ (i, i = 1, numTrial) ]
        real(RK)    , parameter :: successProb = 0.7_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: LogProbGeo_ref(numTrial) =   [ -.3566749439387324_RK &
                                                                , -1.560647748264668_RK &
                                                                , -2.764620552590604_RK &
                                                                , -3.968593356916540_RK &
                                                                , -5.172566161242476_RK &
                                                                , -6.376538965568412_RK &
                                                                , -7.580511769894348_RK &
                                                                , -8.784484574220283_RK &
                                                                , -9.988457378546219_RK &
                                                                , -11.19243018287215_RK ]
        real(RK)                :: LogProbGeo(numTrial)
        real(RK)                :: Difference(numTrial)

        LogProbGeo = getLogProbGeo(numTrial, SuccessStep, successProb)
        Difference = abs(LogProbGeo - LogProbGeo_ref) / abs(LogProbGeo_ref)
        assertion = all( Difference < tolerance )

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogProbGeo_ref  =", LogProbGeo_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogProbGeo      =", LogProbGeo
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference      =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbGeo_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbGeoCyclic_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: numTrial = 10_IK
        integer(IK) , parameter :: maxNumTrial = 3_IK
        integer(IK) , parameter :: SuccessStep(numTrial) = [ (i, i = 1, numTrial) ]
        real(RK)    , parameter :: successProb = 0.7_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: LogProbGeoCyclic_ref(numTrial) = [ -0.32930374714260041_RK &
                                                                    , -1.53327655146853630_RK &
                                                                    , -2.73724935579447240_RK &
                                                                    , -3.94122216012040830_RK &
                                                                    , -5.14519496444634420_RK &
                                                                    , -6.34916776877228010_RK &
                                                                    , -7.55314057309821600_RK &
                                                                    , -8.75711337742415100_RK &
                                                                    , -9.96108618175008690_RK &
                                                                    , -11.1650589860760230_RK ]
        real(RK)                :: LogProbGeoCyclic(numTrial)
        real(RK)                :: Difference(numTrial)

        LogProbGeoCyclic = getLogProbGeoCyclic(successProb, maxNumTrial, numTrial, SuccessStep)
        Difference = abs(LogProbGeoCyclic - LogProbGeoCyclic_ref) / abs(LogProbGeoCyclic_ref)
        assertion = all( Difference < tolerance )

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogProbGeoCyclic_ref =", LogProbGeoCyclic_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogProbGeoCyclic     =", LogProbGeoCyclic
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference           =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbGeoCyclic_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSNormPDF_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: snormPDF_ref = 0.004431848411938_RK
        real(RK)                :: difference
        real(RK)                :: snormPDF
        snormPDF = getSNormPDF(3._RK)
        difference = abs( (snormPDF - snormPDF_ref) / snormPDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "snormPDF_ref   ", snormPDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "snormPDF       ", snormPDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getSNormPDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSNormCDF_SPR_1() result(assertion)
        use iso_fortran_env, only: RK => real32
        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-5_RK
        real(RK)    , parameter :: snormCDF_ref = 0.998650101968370_RK
        real(RK)                :: difference
        real(RK)                :: snormCDF
        snormCDF = getSNormCDF(3._RK)
        difference = abs( (snormCDF - snormCDF_ref) / snormCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "snormCDF_ref   ", snormCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "snormCDF       ", snormCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getSNormCDF_SPR_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getSNormCDF_DPR_1() result(assertion)
        use iso_fortran_env, only: RK => real64
        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: snormCDF_ref = 0.998650101968370_RK
        real(RK)                :: difference
        real(RK)                :: snormCDF
        snormCDF = getSNormCDF(3._RK)
        difference = abs( (snormCDF - snormCDF_ref) / snormCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "snormCDF_ref   ", snormCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "snormCDF       ", snormCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getSNormCDF_DPR_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getNormPDF_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: avg = 2._RK
        real(RK)    , parameter :: std = 5._RK
        real(RK)    , parameter :: val = 10._RK
        real(RK)    , parameter :: normPDF_ref = 0.022184166935891_RK
        real(RK)                :: difference
        real(RK)                :: normPDF
        normPDF = getNormPDF(avg,std,std**2,val)
        difference = abs( (normPDF - normPDF_ref) / normPDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "normPDF_ref    ", normPDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "normPDF        ", normPDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getNormPDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getNormCDF_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: avg = 2._RK
        real(RK)    , parameter :: std = 5._RK
        real(RK)    , parameter :: val = 10._RK
        real(RK)    , parameter :: normCDF_ref = 0.945200708300442_RK
        real(RK)                :: difference
        real(RK)                :: normCDF
        normCDF = getNormCDF(avg,std,val)
        difference = abs( (normCDF - normCDF_ref) / normCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "normCDF_ref    ", normCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "normCDF        ", normCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getNormCDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBetaCDF_SPR_1() result(assertion)
        use iso_fortran_env, only: RK => real32
        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: alpha = 2._RK
        real(RK)    , parameter :: beta = 5._RK
        real(RK)    , parameter :: val = 0.5_RK
        real(RK)    , parameter :: betaCDF_ref = 0.890625000000000_RK
        real(RK)                :: difference
        real(RK)                :: betaCDF
        betaCDF = getBetaCDF(alpha,beta,val)
        difference = abs( (betaCDF - betaCDF_ref) / betaCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "betaCDF_ref    ", betaCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "betaCDF        ", betaCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
        assertion = assertion .and. getBetaCDF(alpha,beta,-.01_RK) < 0._RK .and. getBetaCDF(alpha,beta,+1.01_RK) < 0._RK
    end function test_getBetaCDF_SPR_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBetaCDF_DPR_1() result(assertion)
        use iso_fortran_env, only: RK => real64
        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: alpha = 2._RK
        real(RK)    , parameter :: beta = 5._RK
        real(RK)    , parameter :: val = 0.5_RK
        real(RK)    , parameter :: betaCDF_ref = 0.890625000000000_RK
        real(RK)                :: difference
        real(RK)                :: betaCDF
        betaCDF = getBetaCDF(alpha,beta,val)
        difference = abs( (betaCDF - betaCDF_ref) / betaCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "betaCDF_ref    ", betaCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "betaCDF        ", betaCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
        assertion = assertion .and. getBetaCDF(alpha,beta,-.01_RK) < 0._RK .and. getBetaCDF(alpha,beta,+1.01_RK) < 0._RK
    end function test_getBetaCDF_DPR_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUniformCDF_1() result(assertion)
        use Constants_mod, only: IK, RK
        implicit none
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: avg = 2._RK
        real(RK)    , parameter :: std = 5._RK
        real(RK)    , parameter :: val = 10._RK
        real(RK)    , parameter :: uniformCDF_ref = val
        real(RK)                :: difference
        real(RK)                :: uniformCDF
        uniformCDF = getUniformCDF(val)
        difference = abs( (uniformCDF - uniformCDF_ref) / uniformCDF_ref )
        assertion = difference < tolerance
        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "uniformCDF_ref ", uniformCDF_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "uniformCDF     ", uniformCDF
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getUniformCDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doKS1_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-7_RK

        real(RK), parameter     :: probKS_ref = .3763758622852317E-01_RK
        real(RK), parameter     :: statKS_ref = .1955719390701096_RK
        real(RK)                :: statKS
        real(RK)                :: probKS
        real(RK)                :: difference
        type(Err_type)          :: Err

        call doKS1  ( np = lenRnd &
                    , Point = StdNormRnd1 &
                    , getCDF = getSNormCDF_DPR &
                    , statKS = statKS &
                    , probKS = probKS &
                    , Err = Err &
                    )
        assertion = .not. Err%occurred
        if (.not. assertion) return

        difference = abs( (probKS - probKS_ref) / probKS_ref )
        assertion = difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS     :", probKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (statKS - statKS_ref) / statKS_ref )
        assertion = assertion .and. difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS     :", statKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doKS1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doUniformKS1_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: tolerance = 1.e-7_RK

        real(RK), parameter     :: probKS_ref = .1982797523608350_RK
        real(RK), parameter     :: statKS_ref = .1491359075319200_RK
        real(RK)                :: statKS
        real(RK)                :: probKS
        real(RK)                :: difference
        type(Err_type)          :: Err

        call doUniformKS1   ( np = lenRnd &
                            , Point = UnifRnd &
                            , statKS = statKS &
                            , probKS = probKS &
                            , Err = Err &
                            )
        assertion = .not. Err%occurred
        if (.not. assertion) return

        difference = abs( (probKS - probKS_ref) / probKS_ref )
        assertion = difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS     :", probKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (statKS - statKS_ref) / statKS_ref )
        assertion = assertion .and. difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS     :", statKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doUniformKS1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doSortedKS2_1() result(assertion)

        use Constants_mod, only: RK, IK
        use Sort_mod, only: sortAscending

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: lenRnd = 50_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK

        real(RK), parameter     :: probKS_ref = 0.056045859714425_RK
        real(RK), parameter     :: statKS_ref = 0.260000000000000_RK
        real(RK)                :: difference, statKS, probKS
        type(Err_type)          :: Err

        call sortAscending( np = lenRnd, Point = StdNormRnd1, Err = Err )
        assertion = .not. Err%occurred
        if (.not. assertion) return

        call sortAscending( np = lenRnd, Point = StdNormRnd2, Err = Err )
        assertion = .not. Err%occurred
        if (.not. assertion) return

        call doSortedKS2( np1 = lenRnd &
                        , np2 = lenRnd &
                        , SortedPoint1 = StdNormRnd1 &
                        , SortedPoint2 = StdNormRnd2 &
                        , statKS = statKS &
                        , probKS = probKS &
                        )

        difference = 2 * abs(probKS - probKS_ref) / (probKS_ref + probKS)
        assertion = assertion .and. difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "probKS     :", probKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = 2 * abs(statKS - statKS_ref) / (statKS_ref + statKS)
        assertion = difference < tolerance

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "statKS     :", statKS
            write(Test%outputUnit,"(*(g0,:,', '))") "difference :", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_doSortedKS2_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist1D_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nxbin) = [ 1.000000000000000_RK &
                                                        , 1.000000000000000_RK &
                                                        , 4.000000000000000_RK &
                                                        , 6.000000000000000_RK &
                                                        , 11.00000000000000_RK &
                                                        , 15.00000000000000_RK &
                                                        , 8.000000000000000_RK &
                                                        , 1.000000000000000_RK &
                                                        , 2.000000000000000_RK &
                                                        , 1.000000000000000_RK ]
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -3.068150106908867_RK &
                                                        , -2.315881996736801_RK &
                                                        , -1.563613886564735_RK &
                                                        , -.8113457763926690_RK &
                                                        , -.5907766622060301e-01_RK &
                                                        , +.6931904439514631_RK &
                                                        , +1.445458554123529_RK &
                                                        , +2.197726664295595_RK &
                                                        , +2.949994774467661_RK &
                                                        , +3.702262884639727_RK ]
        real(RK)                :: Density(nxbin)
        real(RK)                :: Xbin(nxbin)
        real(RK)                :: Difference(nxbin)

        call getHist1D  ( method = "count" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , nxbin = nxbin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Xbin = Xbin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Difference = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Difference < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Difference = abs( (Density - Density_ref) / Density_ref )
        assertion = assertion .and. all(Difference < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density    ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist1D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist1D_2() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 10_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nxbin) = [ .2000000000000000E-01_RK &
                                                        , .2000000000000000E-01_RK &
                                                        , .8000000000000000E-01_RK &
                                                        , .1200000000000000_RK &
                                                        , .2200000000000000_RK &
                                                        , .3000000000000000_RK &
                                                        , .1600000000000000_RK &
                                                        , .2000000000000000E-01_RK &
                                                        , .4000000000000000E-01_RK &
                                                        , .2000000000000000E-01_RK ]
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -3.068150106908867_RK &
                                                        , -2.315881996736801_RK &
                                                        , -1.563613886564735_RK &
                                                        , -.8113457763926690_RK &
                                                        , -.5907766622060301e-01_RK &
                                                        , +.6931904439514631_RK &
                                                        , +1.445458554123529_RK &
                                                        , +2.197726664295595_RK &
                                                        , +2.949994774467661_RK &
                                                        , +3.702262884639727_RK ]
        real(RK)                :: Density(nxbin)
        real(RK)                :: Xbin(nxbin)
        real(RK)                :: Difference(nxbin)

        call getHist1D  ( method = "pdf" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , nxbin = nxbin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Xbin = Xbin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Difference = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Difference < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Difference = abs( (Density - Density_ref) / Density_ref )
        assertion = assertion .and. all(Difference < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density    ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist1D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist1D_3() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 10_IK
        real(RK)                :: Density(nxbin)
        real(RK)                :: Xbin(nxbin)

        call getHist1D  ( method = "nonsense" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , nxbin = nxbin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Xbin = Xbin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )

    end function test_getHist1D_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist2D_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 8_IK
        integer(IK) , parameter :: nybin = 7_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nybin,nxbin) = reshape(  [ .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 3.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 3.000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 3.000000000000000_RK &
                                                                        , 5.000000000000000_RK &
                                                                        , 5.000000000000000_RK &
                                                                        , 3.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 2.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        ], shape = shape(Density_ref) )
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -2.974116593137359_RK &
                                                        , -2.033781455422276_RK &
                                                        , -1.093446317707194_RK &
                                                        , -.1531111799921113_RK &
                                                        , .7872239577229714_RK &
                                                        , 1.727559095438054_RK &
                                                        , 2.667894233153136_RK &
                                                        , 3.608229370868219_RK ]
        real(RK)    , parameter :: Ybin_ref(nybin) =    [ -2.038843334246188_RK &
                                                        , -1.250484167036584_RK &
                                                        , -.4621249998269794_RK &
                                                        , .3262341673826251_RK &
                                                        , 1.114593334592229_RK &
                                                        , 1.902952501801833_RK &
                                                        , 2.691311669011438_RK ]
        real(RK)                :: Density(nybin,nxbin)
        real(RK)                :: Xbin(nxbin),Xbin_diff(nxbin)
        real(RK)                :: Ybin(nybin),Ybin_diff(nybin)
        real(RK)                :: Density_diff(nybin,nxbin)

        call getHist2D  ( histType = "count" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , ymin = minval(StdNormRnd2) - 0.5_RK &
                        , ymax = maxval(StdNormRnd2) + 0.5_RK &
                        , nxbin = nxbin &
                        , nybin = nybin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Y = StdNormRnd2 &
                        , Xbin = Xbin &
                        , Ybin = Ybin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Xbin_diff = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Xbin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_diff  ", Xbin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Ybin_diff = abs( (Ybin - Ybin_ref) / Ybin_ref )
        assertion = all(Ybin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_ref   ", Ybin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin       ", Ybin
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_diff  ", Ybin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Density_diff = abs(Density - Density_ref)
        assertion = assertion .and. all(Density_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref    ", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density        ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_diff   ", Density_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist2D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist2D_2() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 8_IK
        integer(IK) , parameter :: nybin = 7_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nybin,nxbin) = reshape(  [ .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .6000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .6000000000000000E-01_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .6000000000000000E-01_RK &
                                                                        , .1000000000000000_RK &
                                                                        , .1000000000000000_RK &
                                                                        , .6000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .4000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        ], shape = shape(Density_ref) )
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -2.974116593137359_RK &
                                                        , -2.033781455422276_RK &
                                                        , -1.093446317707194_RK &
                                                        , -.1531111799921113_RK &
                                                        , .7872239577229714_RK &
                                                        , 1.727559095438054_RK &
                                                        , 2.667894233153136_RK &
                                                        , 3.608229370868219_RK ]
        real(RK)    , parameter :: Ybin_ref(nybin) =    [ -2.038843334246188_RK &
                                                        , -1.250484167036584_RK &
                                                        , -.4621249998269794_RK &
                                                        , .3262341673826251_RK &
                                                        , 1.114593334592229_RK &
                                                        , 1.902952501801833_RK &
                                                        , 2.691311669011438_RK ]
        real(RK)                :: Density(nybin,nxbin)
        real(RK)                :: Xbin(nxbin),Xbin_diff(nxbin)
        real(RK)                :: Ybin(nybin),Ybin_diff(nybin)
        real(RK)                :: Density_diff(nybin,nxbin)

        call getHist2D  ( histType = "pdf" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , ymin = minval(StdNormRnd2) - 0.5_RK &
                        , ymax = maxval(StdNormRnd2) + 0.5_RK &
                        , nxbin = nxbin &
                        , nybin = nybin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Y = StdNormRnd2 &
                        , Xbin = Xbin &
                        , Ybin = Ybin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Xbin_diff = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Xbin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_diff  ", Xbin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Ybin_diff = abs( (Ybin - Ybin_ref) / Ybin_ref )
        assertion = all(Ybin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_ref   ", Ybin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin       ", Ybin
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_diff  ", Ybin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Density_diff = abs(Density - Density_ref)
        assertion = assertion .and. all(Density_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref    ", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density        ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_diff   ", Density_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist2D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist2D_3() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 8_IK
        integer(IK) , parameter :: nybin = 7_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nybin,nxbin) = reshape(  [ .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2222222222222222_RK &
                                                                        , .2222222222222222_RK &
                                                                        , .3333333333333333_RK &
                                                                        , .1111111111111111_RK &
                                                                        , .1111111111111111_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1818181818181818_RK &
                                                                        , .2727272727272727_RK &
                                                                        , .1818181818181818_RK &
                                                                        , .1818181818181818_RK &
                                                                        , .9090909090909091E-01_RK &
                                                                        , .9090909090909091E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1764705882352941_RK &
                                                                        , .2941176470588235_RK &
                                                                        , .2941176470588235_RK &
                                                                        , .1764705882352941_RK &
                                                                        , .5882352941176471E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1428571428571428_RK &
                                                                        , .2857142857142857_RK &
                                                                        , .2857142857142857_RK &
                                                                        , .1428571428571428_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1428571428571428_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , 1.000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        ], shape = shape(Density_ref) )
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -2.974116593137359_RK &
                                                        , -2.033781455422276_RK &
                                                        , -1.093446317707194_RK &
                                                        , -.1531111799921113_RK &
                                                        , .7872239577229714_RK &
                                                        , 1.727559095438054_RK &
                                                        , 2.667894233153136_RK &
                                                        , 3.608229370868219_RK ]
        real(RK)    , parameter :: Ybin_ref(nybin) =    [ -2.038843334246188_RK &
                                                        , -1.250484167036584_RK &
                                                        , -.4621249998269794_RK &
                                                        , .3262341673826251_RK &
                                                        , 1.114593334592229_RK &
                                                        , 1.902952501801833_RK &
                                                        , 2.691311669011438_RK ]
        real(RK)                :: Density(nybin,nxbin)
        real(RK)                :: Xbin(nxbin),Xbin_diff(nxbin)
        real(RK)                :: Ybin(nybin),Ybin_diff(nybin)
        real(RK)                :: Density_diff(nybin,nxbin)

        call getHist2D  ( histType = "pdf(y|x)" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , ymin = minval(StdNormRnd2) - 0.5_RK &
                        , ymax = maxval(StdNormRnd2) + 0.5_RK &
                        , nxbin = nxbin &
                        , nybin = nybin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Y = StdNormRnd2 &
                        , Xbin = Xbin &
                        , Ybin = Ybin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Xbin_diff = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Xbin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_diff  ", Xbin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Ybin_diff = abs( (Ybin - Ybin_ref) / Ybin_ref )
        assertion = all(Ybin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_ref   ", Ybin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin       ", Ybin
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_diff  ", Ybin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Density_diff = abs(Density - Density_ref)
        assertion = assertion .and. all(Density_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref    ", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density        ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_diff   ", Density_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist2D_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist2D_4() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 8_IK
        integer(IK) , parameter :: nybin = 7_IK
        real(RK)    , parameter :: tolerance = 1.e-7_RK
        real(RK)    , parameter :: Density_ref(nybin,nxbin) = reshape(  [ .000000000000000_RK &
                                                                        , .1000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .7692307692307693E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000_RK &
                                                                        , .1538461538461539_RK &
                                                                        , .2307692307692308_RK &
                                                                        , .1428571428571428_RK &
                                                                        , .3333333333333333_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .2000000000000000_RK &
                                                                        , .2307692307692308_RK &
                                                                        , .1538461538461539_RK &
                                                                        , .2857142857142857_RK &
                                                                        , .3333333333333333_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .3000000000000000_RK &
                                                                        , .3846153846153846_RK &
                                                                        , .3846153846153846_RK &
                                                                        , .4285714285714285_RK &
                                                                        , .3333333333333333_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .1000000000000000_RK &
                                                                        , .1538461538461539_RK &
                                                                        , .1538461538461539_RK &
                                                                        , .1428571428571428_RK &
                                                                        , .000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .7692307692307693E-01_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .5000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        , .000000000000000_RK &
                                                                        ], shape = shape(Density_ref) )
        real(RK)    , parameter :: Xbin_ref(nxbin) =    [ -2.974116593137359_RK &
                                                        , -2.033781455422276_RK &
                                                        , -1.093446317707194_RK &
                                                        , -.1531111799921113_RK &
                                                        , .7872239577229714_RK &
                                                        , 1.727559095438054_RK &
                                                        , 2.667894233153136_RK &
                                                        , 3.608229370868219_RK ]
        real(RK)    , parameter :: Ybin_ref(nybin) =    [ -2.038843334246188_RK &
                                                        , -1.250484167036584_RK &
                                                        , -.4621249998269794_RK &
                                                        , .3262341673826251_RK &
                                                        , 1.114593334592229_RK &
                                                        , 1.902952501801833_RK &
                                                        , 2.691311669011438_RK ]
        real(RK)                :: Density(nybin,nxbin)
        real(RK)                :: Xbin(nxbin),Xbin_diff(nxbin)
        real(RK)                :: Ybin(nybin),Ybin_diff(nybin)
        real(RK)                :: Density_diff(nybin,nxbin)

        call getHist2D  ( histType = "pdf(x|y)" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , ymin = minval(StdNormRnd2) - 0.5_RK &
                        , ymax = maxval(StdNormRnd2) + 0.5_RK &
                        , nxbin = nxbin &
                        , nybin = nybin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Y = StdNormRnd2 &
                        , Xbin = Xbin &
                        , Ybin = Ybin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )
        assertion = .not. assertion
        if (.not. assertion) return

        Xbin_diff = abs( (Xbin - Xbin_ref) / Xbin_ref )
        assertion = all(Xbin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_ref   ", Xbin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin       ", Xbin
            write(Test%outputUnit,"(*(g0,:,', '))") "Xbin_diff  ", Xbin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Ybin_diff = abs( (Ybin - Ybin_ref) / Ybin_ref )
        assertion = all(Ybin_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_ref   ", Ybin_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin       ", Ybin
            write(Test%outputUnit,"(*(g0,:,', '))") "Ybin_diff  ", Ybin_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Density_diff = abs(Density - Density_ref)
        assertion = assertion .and. all(Density_diff < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_ref    ", Density_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Density        ", Density
            write(Test%outputUnit,"(*(g0,:,', '))") "Density_diff   ", Density_diff
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getHist2D_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getHist2D_5() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nxbin = 8_IK
        integer(IK) , parameter :: nybin = 7_IK
        real(RK)                :: Density(nybin,nxbin)
        real(RK)                :: Xbin(nxbin)
        real(RK)                :: Ybin(nybin)

        call getHist2D  ( histType = "nonsense" &
                        , xmin = minval(StdNormRnd1) - 0.5_RK &
                        , xmax = maxval(StdNormRnd1) + 0.5_RK &
                        , ymin = minval(StdNormRnd2) - 0.5_RK &
                        , ymax = maxval(StdNormRnd2) + 0.5_RK &
                        , nxbin = nxbin &
                        , nybin = nybin &
                        , np = lenRnd &
                        , X = StdNormRnd1 &
                        , Y = StdNormRnd2 &
                        , Xbin = Xbin &
                        , Ybin = Ybin &
                        , Density = Density &
                        , errorOccurred = assertion &
                        )

    end function test_getHist2D_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getQuantile_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nq = 9_IK
        integer(IK) , parameter :: Weight(lenRnd) = [ 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 0_IK, 1_IK, 2_IK &
                                                    , 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK &
                                                    , 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK &
                                                    , 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK ]
        real(RK)    , parameter :: SortedQuantileProbability(nq) = [0._RK, 0.05_RK, 0.1_RK, 0.25_RK, 0.5_RK, 0.75_RK, 0.9_RK, 0.95_RK, 1._RK]
        real(RK)    , parameter :: Quantile_ref(nq) =   [ -2.944284161994900_RK &
                                                        , -2.258846861003650_RK &
                                                        , -1.711516418853700_RK &
                                                        , -.3034409247860160_RK &
                                                        , .3251905394561980_RK &
                                                        , 1.034693009917860_RK &
                                                        , 1.489697607785470_RK &
                                                        , 1.630235289164730_RK &
                                                        , 3.578396939725760_RK ]
        real(RK)                :: Quantile(nq)
        real(RK)                :: Difference(nq)

        Quantile = getQuantile  ( np = lenRnd &
                                , nq = nq &
                                , SortedQuantileProbability = SortedQuantileProbability &
                                , Point = StdNormRnd1 &
                                , Weight = Weight &
                                , sumWeight = sum(Weight) &
                                )

        Difference = (Quantile - Quantile_ref)
        assertion = all(Difference==0._RK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Quantile_ref   ", Quantile_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Quantile       ", Quantile
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getQuantile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getQuantile_2() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        integer(IK) , parameter :: nq = 9_IK
        real(RK)    , parameter :: SortedQuantileProbability(nq) = [0._RK, 0.05_RK, 0.1_RK, 0.25_RK, 0.5_RK, 0.75_RK, 0.9_RK, 0.95_RK, 1._RK]
        real(RK)    , parameter :: Quantile_ref(nq) =   [ -2.944284161994900_RK &
                                                        , -1.711516418853700_RK &
                                                        , -1.307688296305270_RK &
                                                        , -.4335920223056840_RK &
                                                        , .3192067391655020_RK &
                                                        , 1.034693009917860_RK &
                                                        , 1.489697607785470_RK &
                                                        , 2.769437029884880_RK &
                                                        , 3.578396939725760_RK ]
        real(RK)                :: Quantile(nq)
        real(RK)                :: Difference(nq)

        Quantile = getQuantile  ( np = lenRnd &
                                , nq = nq &
                                , SortedQuantileProbability = SortedQuantileProbability &
                                , Point = StdNormRnd1 &
                                )

        Difference = (Quantile - Quantile_ref)
        assertion = all(Difference==0._RK)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Quantile_ref   ", Quantile_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "Quantile       ", Quantile
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getQuantile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether [flatten_2D()](@ref statistics_mod::flatten_2d) can successfully flatten an input weighted 2D array.
    function test_flatten_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        integer(IK)             :: i
        logical                 :: assertion
        integer(IK) , parameter :: nd = 2_IK
        integer(IK) , parameter :: np = 4_IK
        integer(IK) , parameter :: Weight(np) = [-1_IK, 2_IK, 0_IK, 1_IK]
        integer(IK) , parameter :: sumWeight = 3_IK
        real(RK)    , parameter :: tolerance = 1.e-14_RK
        real(RK)    , parameter :: Point(nd,np) = reshape( [( real(i,RK),i=1,nd*np )], shape = shape(Point) )
        real(RK)    , parameter :: FlattenedPoint_ref(nd,3) = reshape([ 3._RK, 4._RK, 3._RK, 4._RK, 7._RK, 8._RK, 7._RK, 8._RK ], shape = shape(FlattenedPoint_ref) )
        real(RK), allocatable   :: FlattenedPoint(:,:)
        real(RK), allocatable   :: Difference(:,:)

        FlattenedPoint = flatten_2D(nd, np, Point, Weight)

        Difference = abs(FlattenedPoint - FlattenedPoint_ref)
        assertion = all(shape(FlattenedPoint)==[nd,sumWeight])
        assertion = assertion .and. all(Difference < tolerance)

        if (Test%isVerboseMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "FlattenedPoint_ref ", FlattenedPoint_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "FlattenedPoint     ", FlattenedPoint
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", Difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_flatten_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Statistics_mod