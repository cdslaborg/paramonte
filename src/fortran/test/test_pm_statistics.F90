!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains tests of the statistics modules.
!>
!>  \final
!>
!>  \author
!>  Amir Shahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_statistics

    use pm_statistics
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

   !integer(IK) , parameter :: lenRnd = 50_IK
   !real(RK) :: UnifRnd(lenRnd) =   [ 0.0759666916908419_RK &
   !                                , 0.2399161535536580_RK &
   !                                , 0.1233189348351660_RK &
   !                                , 0.1839077882824170_RK &
   !                                , 0.2399525256649030_RK &
   !                                , 0.4172670690843700_RK &
   !                                , 0.0496544303257421_RK &
   !                                , 0.9027161099152810_RK &
   !                                , 0.9447871897216460_RK &
   !                                , 0.4908640924680800_RK &
   !                                , 0.4892526384000190_RK &
   !                                , 0.3377194098213770_RK &
   !                                , 0.9000538464176620_RK &
   !                                , 0.3692467811202150_RK &
   !                                , 0.1112027552937870_RK &
   !                                , 0.7802520683211380_RK &
   !                                , 0.3897388369612530_RK &
   !                                , 0.2416912859138330_RK &
   !                                , 0.4039121455881150_RK &
   !                                , 0.0964545251683886_RK &
   !                                , 0.1319732926063350_RK &
   !                                , 0.9420505907754850_RK &
   !                                , 0.9561345402298020_RK &
   !                                , 0.5752085950784660_RK &
   !                                , 0.0597795429471558_RK &
   !                                , 0.2347799133724060_RK &
   !                                , 0.3531585712220710_RK &
   !                                , 0.8211940401979590_RK &
   !                                , 0.0154034376515551_RK &
   !                                , 0.0430238016578078_RK &
   !                                , 0.1689900294627040_RK &
   !                                , 0.6491154749564520_RK &
   !                                , 0.7317223856586700_RK &
   !                                , 0.6477459631363070_RK &
   !                                , 0.4509237064309450_RK &
   !                                , 0.5470088922863450_RK &
   !                                , 0.2963208056077730_RK &
   !                                , 0.7446928070741560_RK &
   !                                , 0.1889550150325450_RK &
   !                                , 0.6867754333653150_RK &
   !                                , 0.1835111557372700_RK &
   !                                , 0.3684845964903370_RK &
   !                                , 0.6256185607296900_RK &
   !                                , 0.7802274351513770_RK &
   !                                , 0.0811257688657853_RK &
   !                                , 0.9293859709687300_RK &
   !                                , 0.7757126786084020_RK &
   !                                , 0.4867916324031720_RK &
   !                                , 0.4358585885809190_RK &
   !                                , 0.4467837494298060_RK ]

   !real(RK) :: StdNormRnd1(lenRnd) =   [ 0.537667139546100_RK &
   !                                    , 1.83388501459509_RK &
   !                                    , -2.25884686100365_RK &
   !                                    , 0.862173320368121_RK &
   !                                    , 0.318765239858981_RK &
   !                                    , -1.30768829630527_RK &
   !                                    , -0.433592022305684_RK &
   !                                    , 0.342624466538650_RK &
   !                                    , 3.57839693972576_RK &
   !                                    , 2.76943702988488_RK &
   !                                    , -1.34988694015652_RK &
   !                                    , 3.03492346633185_RK &
   !                                    , 0.725404224946106_RK &
   !                                    , -0.0630548731896562_RK &
   !                                    , 0.714742903826096_RK &
   !                                    , -0.204966058299775_RK &
   !                                    , -0.124144348216312_RK &
   !                                    , 1.48969760778547_RK &
   !                                    , 1.40903448980048_RK &
   !                                    , 1.41719241342961_RK &
   !                                    , 0.671497133608081_RK &
   !                                    , -1.20748692268504_RK &
   !                                    , 0.717238651328839_RK &
   !                                    , 1.63023528916473_RK &
   !                                    , 0.488893770311789_RK &
   !                                    , 1.03469300991786_RK &
   !                                    , 0.726885133383238_RK &
   !                                    , -0.303440924786016_RK &
   !                                    , 0.293871467096658_RK &
   !                                    , -0.787282803758638_RK &
   !                                    , 0.888395631757642_RK &
   !                                    , -1.14707010696915_RK &
   !                                    , -1.06887045816803_RK &
   !                                    , -0.809498694424876_RK &
   !                                    , -2.94428416199490_RK &
   !                                    , 1.43838029281510_RK &
   !                                    , 0.325190539456198_RK &
   !                                    , -0.754928319169703_RK &
   !                                    , 1.37029854009523_RK &
   !                                    , -1.71151641885370_RK &
   !                                    , -0.102242446085491_RK &
   !                                    , -0.241447041607358_RK &
   !                                    , 0.319206739165502_RK &
   !                                    , 0.312858596637428_RK &
   !                                    , -0.864879917324457_RK &
   !                                    , -0.0300512961962686_RK &
   !                                    , -0.164879019209038_RK &
   !                                    , 0.627707287528727_RK &
   !                                    , 1.09326566903948_RK &
   !                                    , 1.10927329761440_RK ]
   !
   !real(RK) :: StdNormRnd2(lenRnd) =   [ -0.863652821988714_RK &
   !                                    , 0.0773590911304249_RK &
   !                                    , -1.21411704361541_RK &
   !                                    , -1.11350074148676_RK &
   !                                    , -0.00684932810334806_RK &
   !                                    , 1.53263030828475_RK &
   !                                    , -0.769665913753682_RK &
   !                                    , 0.371378812760058_RK &
   !                                    , -0.225584402271252_RK &
   !                                    , 1.11735613881447_RK &
   !                                    , -1.08906429505224_RK &
   !                                    , 0.0325574641649735_RK &
   !                                    , 0.552527021112224_RK &
   !                                    , 1.10061021788087_RK &
   !                                    , 1.54421189550395_RK &
   !                                    , 0.0859311331754255_RK &
   !                                    , -1.49159031063761_RK &
   !                                    , -0.742301837259857_RK &
   !                                    , -1.06158173331999_RK &
   !                                    , 2.35045722400204_RK &
   !                                    , -0.615601881466894_RK &
   !                                    , 0.748076783703985_RK &
   !                                    , -0.192418510588264_RK &
   !                                    , 0.888610425420721_RK &
   !                                    , -0.764849236567874_RK &
   !                                    , -1.40226896933876_RK &
   !                                    , -1.42237592509150_RK &
   !                                    , 0.488193909859941_RK &
   !                                    , -0.177375156618825_RK &
   !                                    , -0.196053487807333_RK &
   !                                    , 1.41931015064255_RK &
   !                                    , 0.291584373984183_RK &
   !                                    , 0.197811053464361_RK &
   !                                    , 1.58769908997406_RK &
   !                                    , -0.804465956349547_RK &
   !                                    , 0.696624415849607_RK &
   !                                    , 0.835088165072682_RK &
   !                                    , -0.243715140377952_RK &
   !                                    , 0.215670086403744_RK &
   !                                    , -1.16584393148205_RK &
   !                                    , -1.14795277889859_RK &
   !                                    , 0.104874716016494_RK &
   !                                    , 0.722254032225002_RK &
   !                                    , 2.58549125261624_RK &
   !                                    , -0.666890670701386_RK &
   !                                    , 0.187331024578940_RK &
   !                                    , -0.0824944253709554_RK &
   !                                    , -1.93302291785099_RK &
   !                                    , -0.438966153934773_RK &
   !                                    , -1.79467884145512_RK ]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_getRandLogn_1, SK_"test_getRandLogn_1")
        call test%run(test_getRandBeta_1, SK_"test_getRandBeta_1")
        call test%run(test_getRandGamma_1, SK_"test_getRandGamma_1")
        call test%run(test_getLogProbMVU_1, SK_"test_getLogProbMVU_1")
        call test%run(test_getRandCorMat_1, SK_"test_getRandCorMat_1")
        call test%run(test_getLogProbGeo_1, SK_"test_getLogProbGeo_1")
        call test%run(test_getUniformCDF_1, SK_"test_getUniformCDF_1")
        call test%run(test_getUnifRandorm_1, SK_"test_getUnifRandorm_1")
        call test%run(test_getBetaCDF_RK1_1, SK_"test_getBetaCDF_RK1_1")
        call test%run(test_getBetaCDF_RK2_1, SK_"test_getBetaCDF_RK2_1")
        call test%run(test_getLogProbLognSP_1, SK_"test_getLogProbLognSP_1")
        call test%run(test_getLogProbLognMP_1, SK_"test_getLogProbLognMP_1")
        call test%run(test_getRandIntLecuyer_1, SK_"test_getRandIntLecuyer_1")
        call test%run(test_getRandRealLecuyer_1, SK_"test_getRandRealLecuyer_1")
        call test%run(test_getLogProbNormSP_RK_1, SK_"test_getLogProbNormSP_RK_1")
        call test%run(test_getLogProbNormMP_RK_1, SK_"test_getLogProbNormMP_RK_1")
       !call test%run(test_getLogProbNormSP_CK_1, SK_"test_getLogProbNormSP_CK_1")
       !call test%run(test_getLogProbNormMP_CK_1, SK_"test_getLogProbNormMP_CK_1")
        call test%run(test_getLogProbGeoCyclic_1, SK_"test_getLogProbGeoCyclic_1")
        call test%run(test_getRandGammaIntShape_1, SK_"test_getRandGammaIntShape_1")
        call test%run(test_getLogProbMixMVNSP_RK_1, SK_"test_getLogProbMixMVNSP_RK_1")
        call test%run(test_getLogProbMixMVNMP_RK_1, SK_"test_getLogProbMixMVNMP_RK_1")
       !call test%run(test_getLogProbMixMVNSP_CK_1, SK_"test_getLogProbMixMVNSP_CK_1")
       !call test%run(test_getLogProbMixMVNMP_CK_1, SK_"test_getLogProbMixMVNMP_CK_1")
        call test%run(test_getRandCorMatRejection_1, SK_"test_getRandCorMatRejection_1")
        call test%run(test_getLogProbMixNormSP_RK_1, SK_"test_getLogProbMixNormSP_RK_1")
        call test%run(test_getLogProbMixNormMP_RK_1, SK_"test_getLogProbMixNormMP_RK_1")
       !call test%run(test_getLogProbMixNormSP_CK_1, SK_"test_getLogProbMixNormSP_CK_1")
       !call test%run(test_getLogProbMixNormMP_CK_1, SK_"test_getLogProbMixNormMP_CK_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormSP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormMP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm_ref    ", logProbNorm_ref
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm        ", logProbNorm
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormSP_CK_1() result(assertion)

        use pm_kind, only: IK, RK, CK
        implicit none

        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormSP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbNormMP_CK_1() result(assertion)

        use pm_kind, only: IK, RK, CK
        implicit none

        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm_ref    ", real(logProbNorm_ref, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbNorm        ", real(logProbNorm, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbNormMP_CK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormSP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: point = 2._RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mean(nmode) = [ 0.25_RK, 0.75_RK ]
        real(RK)    , parameter :: invCov(nmode) = [ 1._RK / 16._RK, 1._RK / 32._RK ]
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: logSqrtDetInvCovMat(nmode) = log(sqrt(invCov))
        real(RK)    , parameter :: logProbMixNorm_ref = 1.718832134253714_RK
        real(RK)                :: logProbMixNorm
        real(RK)                :: difference

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , LogAmplitude = LogAmplitude &
                                            , mean = mean &
                                            , invCov = invCov &
                                            , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (logProbMixNorm - logProbMixNorm_ref) / logProbMixNorm_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm_ref ", logProbMixNorm_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm     ", logProbMixNorm
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixNormMP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: point(np) = [ 2._RK, 3._RK ]
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: mean(nmode) = [ 0.25_RK, 0.75_RK ]
        real(RK)    , parameter :: invCov(nmode) = [ 1._RK / 16._RK, 1._RK / 32._RK ]
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: logSqrtDetInvCovMat(nmode) = log(sqrt(invCov))
        real(RK)    , parameter :: logProbMixNorm_ref(np) = [ 1.718832134253714_RK, 1.636902047052812_RK ]
        real(RK)                :: logProbMixNorm(np)
        real(RK)                :: difference(np)

        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
                                            , np = np &
                                            , LogAmplitude = LogAmplitude &
                                            , mean = mean &
                                            , invCov = invCov &
                                            , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                            , point = point &
                                            )

        difference = abs( (logProbMixNorm - logProbMixNorm_ref) / logProbMixNorm_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm_ref ", logProbMixNorm_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm     ", logProbMixNorm
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixNormMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNSP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd) = [(real(i,RK),i=1,nd)]
        real(RK)    , parameter :: mean(nd,nmode) = reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(mean))
        real(RK)    , parameter :: invCov(nd,nd,nmode) = reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK &
                                                                    , 2._RK, 0._RK, 0._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 0._RK, 0._RK, 2._RK &
                                                                    ], shape = shape(invCov) )
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: logSqrtDetInvCovMat(nmode) = [-0.693147180559945_RK, -1.039720770839918_RK]
        real(RK)    , parameter :: logProbMixMVN_ref = -90.44996278017396_RK
        real(RK)                :: logProbMixMVN
        real(RK)                :: difference

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , LogAmplitude = LogAmplitude &
                                        , mean = mean &
                                        , invCov = invCov &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbMixMVN - logProbMixMVN_ref) / logProbMixMVN_ref )
        assertion = difference <= tolerance

        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN_ref ", logProbMixMVN_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN     ", logProbMixMVN
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNSP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMixMVNMP_RK_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none

        logical(LK)             :: assertion
        integer(IK)             :: i
        integer(IK) , parameter :: nd = 3_IK
        integer(IK) , parameter :: np = 2_IK
        integer(IK) , parameter :: nmode = 2_IK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: Point(nd,np) = reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point))
        real(RK)    , parameter :: mean(nd,nmode) = reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(mean))
        real(RK)    , parameter :: invCov(nd,nd,nmode) = reshape([ 1._RK, 0._RK, 1._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 1._RK, 0._RK, 3._RK &
                                                                    , 2._RK, 0._RK, 0._RK &
                                                                    , 0._RK, 2._RK, 0._RK &
                                                                    , 0._RK, 0._RK, 2._RK &
                                                                    ], shape = shape(invCov) )
        real(RK)    , parameter :: LogAmplitude(nmode) = [ 3._RK, 4._RK ]
        real(RK)    , parameter :: logSqrtDetInvCovMat(nmode) = [-0.693147180559945_RK, -1.039720770839918_RK]
        real(RK)    , parameter :: logProbMixMVN_ref(np) = [ -90.44996278017396_RK, -18.44996278017396_RK ]
        real(RK)                :: logProbMixMVN(np)
        real(RK)                :: difference(np)

        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
                                        , nd = nd &
                                        , np = np &
                                        , LogAmplitude = LogAmplitude &
                                        , mean = mean &
                                        , invCov = invCov &
                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
                                        , point = point &
                                        )

        difference = abs( (logProbMixMVN - logProbMixMVN_ref) / logProbMixMVN_ref )
        assertion = all(difference <= tolerance)

        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN_ref ", logProbMixMVN_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN     ", logProbMixMVN
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbMixMVNMP_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getLogProbMixNormSP_CK_1() result(assertion)
!
!        use pm_kind, only: IK, RK, CK
!        implicit none
!
!        logical(LK)             :: assertion
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK) , parameter :: nmode = 2_IK
!        complex(CK) , parameter :: point = 2._RK
!        complex(CK) , parameter :: mean(nmode) = cmplx([ 0.25_RK, 0.75_RK ], kind = RK)
!        complex(CK) , parameter :: invCov(nmode) = cmplx([ 1._RK / 16._RK, 1._RK / 32._RK ], kind = RK)
!        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
!        complex(CK) , parameter :: logSqrtDetInvCovMat(nmode) = cmplx(log(sqrt(invCov)), kind = RK)
!        complex(CK) , parameter :: logProbMixNorm_ref = cmplx(1.718832134253714_RK, kind = RK)
!        complex(CK)             :: logProbMixNorm
!        real(RK)                :: difference
!
!        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
!                                            , LogAmplitude = LogAmplitude &
!                                            , mean = mean &
!                                            , invCov = invCov &
!                                            , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
!                                            , point = point &
!                                            )
!
!        difference = abs( (real(logProbMixNorm, kind = RK) - real(logProbMixNorm_ref, kind = RK)) / real(logProbMixNorm_ref, kind = RK) )
!        assertion = difference <= tolerance
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm_ref ", real(logProbMixNorm_ref, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm     ", real(logProbMixNorm, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getLogProbMixNormSP_CK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getLogProbMixNormMP_CK_1() result(assertion)
!
!        use pm_kind, only: IK, RK, CK
!        implicit none
!
!        logical(LK)             :: assertion
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK) , parameter :: np = 2_IK
!        integer(IK) , parameter :: nmode = 2_IK
!        complex(CK) , parameter :: point(np) = cmplx([ 2._RK, 3._RK ], kind = RK)
!        complex(CK) , parameter :: mean(nmode) = cmplx([ 0.25_RK, 0.75_RK ], kind = RK)
!        complex(CK) , parameter :: invCov(nmode) = cmplx([ 1._RK / 16._RK, 1._RK / 32._RK ], kind = RK)
!        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
!        complex(CK) , parameter :: logSqrtDetInvCovMat(nmode) = cmplx(log(sqrt(invCov)), kind = RK)
!        complex(CK) , parameter :: logProbMixNorm_ref(np) = cmplx([ 1.718832134253714_RK, 1.636902047052812_RK ], kind = RK)
!        complex(CK)             :: logProbMixNorm(np)
!        real(RK)                :: difference(np)
!
!        logProbMixNorm = getLogProbMixNorm  ( nmode = nmode &
!                                            , np = np &
!                                            , LogAmplitude = LogAmplitude &
!                                            , mean = mean &
!                                            , invCov = invCov &
!                                            , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
!                                            , point = point &
!                                            )
!
!        difference = abs( (real(logProbMixNorm, kind = RK) - real(logProbMixNorm_ref, kind = RK)) / real(logProbMixNorm_ref, kind = RK) )
!        assertion = all(difference <= tolerance)
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm_ref ", real(logProbMixNorm_ref, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixNorm     ", real(logProbMixNorm, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getLogProbMixNormMP_CK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getLogProbMixMVNSP_CK_1() result(assertion)
!
!        use pm_kind, only: IK, RK, CK
!        implicit none
!
!        logical(LK)             :: assertion
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK)             :: i
!        integer(IK) , parameter :: nd = 3_IK
!        integer(IK) , parameter :: nmode = 2_IK
!        complex(CK) , parameter :: Point(nd) = cmplx([(real(i,RK),i=1,nd)], kind = RK)
!        complex(CK) , parameter :: mean(nd,nmode) = cmplx(reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(mean)), kind = RK)
!        complex(CK) , parameter :: invCov(nd,nd,nmode) = cmplx( reshape( [ 1._RK, 0._RK, 1._RK &
!                                                                            , 0._RK, 2._RK, 0._RK &
!                                                                            , 1._RK, 0._RK, 3._RK &
!                                                                            , 2._RK, 0._RK, 0._RK &
!                                                                            , 0._RK, 2._RK, 0._RK &
!                                                                            , 0._RK, 0._RK, 2._RK &
!                                                                            ], shape = shape(invCov) ), kind = RK)
!        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
!        complex(CK) , parameter :: logSqrtDetInvCovMat(nmode) = cmplx([-0.693147180559945_RK, -1.039720770839918_RK], kind = RK)
!        complex(CK) , parameter :: logProbMixMVN_ref = cmplx(-90.44996278017396_RK, kind = RK)
!        complex(CK)             :: logProbMixMVN
!        real(RK)                :: difference
!
!        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
!                                        , nd = nd &
!                                        , LogAmplitude = LogAmplitude &
!                                        , mean = mean &
!                                        , invCov = invCov &
!                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
!                                        , point = point &
!                                        )
!
!        difference = abs( (real(logProbMixMVN, kind = RK) - real(logProbMixMVN_ref, kind = RK)) / real(logProbMixMVN_ref, kind = RK) )
!        assertion = difference <= tolerance
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN_ref  ", real(logProbMixMVN_ref, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN      ", real(logProbMixMVN, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getLogProbMixMVNSP_CK_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getLogProbMixMVNMP_CK_1() result(assertion)
!
!        use pm_kind, only: IK, RK, CK
!        implicit none
!
!        logical(LK)             :: assertion
!        integer(IK)             :: i
!        real(RK)    , parameter :: tolerance = 1.e-12_RK
!        integer(IK) , parameter :: nd = 3_IK
!        integer(IK) , parameter :: np = 2_IK
!        integer(IK) , parameter :: nmode = 2_IK
!        complex(CK) , parameter :: Point(nd,np) = cmplx(reshape([(real(i,RK),i=1,nd*np)], shape = shape(Point)), kind = RK)
!        complex(CK) , parameter :: mean(nd,nmode) = cmplx(reshape([(real(i**2+1._RK,RK),i=1,nd*nmode)], shape = shape(mean)), kind = RK)
!        complex(CK) , parameter :: invCov(nd,nd,nmode) = cmplx(reshape(  [ 1._RK, 0._RK, 1._RK &
!                                                                            , 0._RK, 2._RK, 0._RK &
!                                                                            , 1._RK, 0._RK, 3._RK &
!                                                                            , 2._RK, 0._RK, 0._RK &
!                                                                            , 0._RK, 2._RK, 0._RK &
!                                                                            , 0._RK, 0._RK, 2._RK &
!                                                                            ], shape = shape(invCov) ), kind = RK)
!        complex(CK) , parameter :: LogAmplitude(nmode) = cmplx([ 3._RK, 4._RK ], kind = RK)
!        complex(CK) , parameter :: logSqrtDetInvCovMat(nmode) = cmplx([-0.693147180559945_RK, -1.039720770839918_RK], kind = RK)
!        complex(CK) , parameter :: logProbMixMVN_ref(np) = cmplx([ -90.44996278017396_RK, -18.44996278017396_RK ], kind = RK)
!        complex(CK)             :: logProbMixMVN(np)
!        real(RK)                :: difference(np)
!
!        logProbMixMVN = getlogProbMixMVN( nmode = nmode &
!                                        , nd = nd &
!                                        , np = np &
!                                        , LogAmplitude = LogAmplitude &
!                                        , mean = mean &
!                                        , invCov = invCov &
!                                        , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
!                                        , point = point &
!                                        )
!
!        difference = abs( (real(logProbMixMVN, kind = RK) - real(logProbMixMVN_ref, kind = RK)) / real(logProbMixMVN_ref, kind = RK) )
!        assertion = all(difference <= tolerance)
!
!        ! LCOV_EXCL_START
!        if (test%traceable .and. .not. assertion) then
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN_ref  ", real(logProbMixMVN_ref, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "logProbMixMVN      ", real(logProbMixMVN, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))") "Difference         ", real(difference, kind = RK)
!            write(test%disp%unit,"(*(g0,:,', '))")
!        end if
!        ! LCOV_EXCL_STOP
!
!    end function test_getLogProbMixMVNMP_CK_1
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! \todo
    ! What is the best method of testing for randomness?
    function test_getRandLogn_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        real(RK)                :: lognRnd
        lognRnd = getRandLogn(0._RK, 1._RK)
        assertion = .true._LK !< There is really no easy testing of randomness. For now, we trust the giants.
    end function test_getRandLogn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbMVU_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMVU_ref ", logProbMVU_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbMVU     ", logProbMVU
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getLogProbMVU_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbLognSP_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "logProbLogn_ref    ", logProbLogn_ref
            write(test%disp%unit,"(*(g0,:,', '))") "logProbLogn        ", logProbLogn
            write(test%disp%unit,"(*(g0,:,', '))") "difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbLognSP_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbLognMP_1() result(assertion)

        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbLogn_ref    ", LogProbLogn_ref
            write(test%disp%unit,"(*(g0,:,', '))") "LogProbLogn        ", LogProbLogn
            write(test%disp%unit,"(*(g0,:,', '))") "difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbLognMP_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandRealLecuyer_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 100
        integer(IK)             :: idum = 3333, i
        real(RK)                :: RandRealLecuyer(np)
        assertion = .true._LK
        do i = 1, np
            RandRealLecuyer(i) = getRandRealLecuyer(idum)
            assertion = assertion .and. RandRealLecuyer(i) <= 1._RK .and. RandRealLecuyer(i) >= 0._RK
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandRealLecuyer(",i,") =", RandRealLecuyer(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandRealLecuyer_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandIntLecuyer_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: lowerBound = -2, upperBound = 9
        integer(IK) , parameter :: np = 100
        integer(IK)             :: idum = 3333, i
        integer(IK)             :: RandIntLecuyer(np)
        assertion = .true._LK
        do i = 1, np
            RandIntLecuyer(i) = getRandIntLecuyer(lowerBound,upperBound,idum)
            assertion = assertion .and. RandIntLecuyer(i) <= upperBound .and. RandIntLecuyer(i) >= lowerBound
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandIntLecuyer(",i,") =", RandIntLecuyer(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getRandIntLecuyer_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUnifRandorm_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: lowerBound = -2._RK, upperBound = 9._RK
        integer(IK)             :: i
        real(RK)                :: RandUniform(np)
        assertion = .true._LK
        do i = 1, np
            RandUniform(i) = getUnifRandorm(lowerBound,upperBound)
            assertion = assertion .and. RandUniform(i) <= upperBound .and. RandUniform(i) >= lowerBound
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandUniform(",i,") =", RandUniform(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
    end function test_getUnifRandorm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandGamma_1() result(assertion)
        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: alpha = 2._RK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = huge(0._RK)
        integer(IK)             :: i
        real(RK)                :: RandGamma(np)
        assertion = .true._LK
        do i = 1, np
            RandGamma(i) = getRandGamma(alpha)
            assertion = assertion .and. RandGamma(i) <= upperBound .and. RandGamma(i) >= lowerBound
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandGamma(",i,") =", RandGamma(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandGamma(alpha=-1._RK) < 0._RK
    end function test_getRandGamma_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandGammaIntShape_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 100
        integer(IK) , parameter :: alpha = 2_IK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = huge(0._RK)
        integer(IK)             :: i
        real(RK)                :: RandGamma(np)
        assertion = .true._LK
        do i = 1, np
            RandGamma(i) = getRandGammaIntShape(alpha)
            assertion = assertion .and. RandGamma(i) <= upperBound .and. RandGamma(i) >= lowerBound
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandGamma(",i,") =", RandGamma(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandGammaIntShape(alpha=-1_IK) < 0._RK
    end function test_getRandGammaIntShape_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandBeta_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: np = 100
        real(RK)    , parameter :: alpha = 2._RK, beta = 3._RK
        real(RK)    , parameter :: lowerBound = 0._RK, upperBound = 1._RK
        integer(IK)             :: i
        real(RK)                :: RandBeta(np)
        assertion = .true._LK
        do i = 1, np
            RandBeta(i) = getRandBeta(alpha, beta)
            assertion = assertion .and. RandBeta(i) <= upperBound .and. RandBeta(i) >= lowerBound
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandBeta(",i,") =", RandBeta(i)
                write(test%disp%unit,"(*(g0,:,' '))")
            end if
            ! LCOV_EXCL_STOP
        end do
        assertion = assertion .and. getRandBeta(alpha=-1._RK, beta=+2._RK) < 0._RK
        assertion = assertion .and. getRandBeta(alpha=+1._RK, beta=-2._RK) < 0._RK
        assertion = assertion .and. getRandBeta(alpha=-1._RK, beta=-2._RK) < 0._RK
    end function test_getRandBeta_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getRandCorMat_1() result(assertion)
        use pm_kind, only: IK, RK
        use pm_matrixDet, only: isPosDef
        implicit none
        logical(LK)             :: assertion
        integer(IK) , parameter :: nd = 2_IK
        integer(IK) , parameter :: np = 100_IK
        real(RK)    , parameter :: eta = 5._RK
        real(RK)    , parameter :: lowerBound = -1._RK, upperBound = 1._RK
        integer(IK)             :: i
        real(RK)                :: RandCorMat(nd,nd)
        assertion = .true._LK
        do i = 1, np
            RandCorMat = getRandCorMat(nd,eta)
            assertion = assertion .and. all(RandCorMat <= upperBound) .and. all(RandCorMat >= lowerBound)
            assertion = assertion .and. isPosDef(nd,RandCorMat)
            ! LCOV_EXCL_START
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandCorMat(:,:,) =", RandCorMat
                write(test%disp%unit,"(*(g0,:,' '))")
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
        use pm_kind, only: IK, RK
        use pm_matrixDet, only: isPosDef
        implicit none
        logical(LK)             :: assertion, assertionCurrent
        integer(IK) , parameter :: nd = 2_IK
        integer(IK) , parameter :: np = 100_IK
        real(RK)    , parameter :: minRho = -.3_RK, maxRho = .6_RK
        real(RK)    , parameter :: lowerBound = -1._RK, upperBound = 1._RK
        integer(IK)             :: i, j, k
        real(RK)                :: RandCorMat(nd,nd)
        assertion = .true._LK
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
            if (test%traceable .and. .not. assertion) then
                write(test%disp%unit,"(*(g0,:,' '))")
                write(test%disp%unit,"(*(g0,:,' '))") "RandCorMat(:,:,) =", RandCorMat
                write(test%disp%unit,"(*(g0,:,' '))")
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

    function test_getLogProbGeo_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "LogProbGeo_ref  =", LogProbGeo_ref
            write(test%disp%unit,"(*(g0,:,' '))") "LogProbGeo      =", LogProbGeo
            write(test%disp%unit,"(*(g0,:,' '))") "Difference      =", Difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbGeo_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogProbGeoCyclic_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "LogProbGeoCyclic_ref =", LogProbGeoCyclic_ref
            write(test%disp%unit,"(*(g0,:,' '))") "LogProbGeoCyclic     =", LogProbGeoCyclic
            write(test%disp%unit,"(*(g0,:,' '))") "Difference           =", Difference
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogProbGeoCyclic_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBetaCDF_RK1_1() result(assertion)
        use pm_kind, only: RK => RK1
        use pm_kind, only: IK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "betaCDF_ref    ", betaCDF_ref
            write(test%disp%unit,"(*(g0,:,', '))") "betaCDF        ", betaCDF
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
        assertion = assertion .and. getBetaCDF(alpha,beta,-.01_RK) < 0._RK .and. getBetaCDF(alpha,beta,+1.01_RK) < 0._RK
    end function test_getBetaCDF_RK1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBetaCDF_RK2_1() result(assertion)
        use pm_kind, only: RK => RK2
        use pm_kind, only: IK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "betaCDF_ref    ", betaCDF_ref
            write(test%disp%unit,"(*(g0,:,', '))") "betaCDF        ", betaCDF
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if
        assertion = assertion .and. getBetaCDF(alpha,beta,-.01_RK) < 0._RK .and. getBetaCDF(alpha,beta,+1.01_RK) < 0._RK
    end function test_getBetaCDF_RK2_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUniformCDF_1() result(assertion)
        use pm_kind, only: IK, RK
        implicit none
        logical(LK)             :: assertion
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "uniformCDF_ref ", uniformCDF_ref
            write(test%disp%unit,"(*(g0,:,', '))") "uniformCDF     ", uniformCDF
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getUniformCDF_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_statistics