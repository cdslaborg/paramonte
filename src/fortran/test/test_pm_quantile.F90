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

!>  \brief This module contains tests of the module [pm_quantile](@ref pm_quantile).
!>  \author Amir Shahmoradi

module test_pm_quantile

    use pm_quantile
    use pm_err, only: err_type
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

    integer(IK) , parameter :: lenRnd = 50_IK

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
   !
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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer, parameter  :: ndata = 50
   !integer(IK), parameter :: DataUnsorted_IK(ndata)= &
   !                                                [ 1201_IK &
   !                                                , 1187_IK &
   !                                                , 1188_IK &
   !                                                , 1193_IK &
   !                                                , 1177_IK &
   !                                                , 1153_IK &
   !                                                , 1134_IK &
   !                                                , 1146_IK &
   !                                                , 1172_IK &
   !                                                , 1181_IK &
   !                                                , 1197_IK &
   !                                                , 1172_IK &
   !                                                , 1141_IK &
   !                                                , 1216_IK &
   !                                                , 1158_IK &
   !                                                , 1174_IK &
   !                                                , 1189_IK &
   !                                                , 1211_IK &
   !                                                , 1157_IK &
   !                                                , 1184_IK &
   !                                                , 1177_IK &
   !                                                , 1157_IK &
   !                                                , 1191_IK &
   !                                                , 1176_IK &
   !                                                , 1196_IK &
   !                                                , 1150_IK &
   !                                                , 1185_IK &
   !                                                , 1190_IK &
   !                                                , 1172_IK &
   !                                                , 1161_IK &
   !                                                , 1179_IK &
   !                                                , 1189_IK &
   !                                                , 1136_IK &
   !                                                , 1148_IK &
   !                                                , 1176_IK &
   !                                                , 1142_IK &
   !                                                , 1146_IK &
   !                                                , 1202_IK &
   !                                                , 1156_IK &
   !                                                , 1170_IK &
   !                                                , 1146_IK &
   !                                                , 1170_IK &
   !                                                , 1169_IK &
   !                                                , 1211_IK &
   !                                                , 1168_IK &
   !                                                , 1189_IK &
   !                                                , 1170_IK &
   !                                                , 1162_IK &
   !                                                , 1167_IK &
   !                                                , 1180_IK ]
    real(RK), parameter :: DataUnsorted_RK(ndata)   = &
                                                    [ 5.28935260000000_RK &
                                                    , 5.50145870000000_RK &
                                                    , 5.89022390000000_RK &
                                                    , 5.06549460000000_RK &
                                                    , 5.62128260000000_RK &
                                                    , 4.49246930000000_RK &
                                                    , 3.54559920000000_RK &
                                                    , 4.17171310000000_RK &
                                                    , 5.34432780000000_RK &
                                                    , 4.30855910000000_RK &
                                                    , 6.12466330000000_RK &
                                                    , 4.45103540000000_RK &
                                                    , 4.08259680000000_RK &
                                                    , 7.64761290000000_RK &
                                                    , 6.53095480000000_RK &
                                                    , 6.07550490000000_RK &
                                                    , 7.32100850000000_RK &
                                                    , 5.82501650000000_RK &
                                                    , 4.19347540000000_RK &
                                                    , 4.89687790000000_RK &
                                                    , 5.61290890000000_RK &
                                                    , 5.70994940000000_RK &
                                                    , 5.00047920000000_RK &
                                                    , 5.47741520000000_RK &
                                                    , 4.99151560000000_RK &
                                                    , 5.08172850000000_RK &
                                                    , 5.98773500000000_RK &
                                                    , 6.97849360000000_RK &
                                                    , 6.91612860000000_RK &
                                                    , 4.90595890000000_RK &
                                                    , 5.71852950000000_RK &
                                                    , 4.12146660000000_RK &
                                                    , 5.51241440000000_RK &
                                                    , 5.26293780000000_RK &
                                                    , 5.14932990000000_RK &
                                                    , 4.14738170000000_RK &
                                                    , 5.55786790000000_RK &
                                                    , 7.08800450000000_RK &
                                                    , 6.08987380000000_RK &
                                                    , 4.73697940000000_RK &
                                                    , 3.80934450000000_RK &
                                                    , 6.03942270000000_RK &
                                                    , 5.96600840000000_RK &
                                                    , 6.06674510000000_RK &
                                                    , 5.84361600000000_RK &
                                                    , 6.19013970000000_RK &
                                                    , 4.43891700000000_RK &
                                                    , 4.45833300000000_RK &
                                                    , 5.47659170000000_RK &
                                                    , 4.65761920000000_RK ]

   !integer(IK) , parameter :: DataUnsorted2_IK(ndata) = DataUnsorted_IK
   !real(RK)    , parameter :: DataUnsorted2_RK(ndata) = DataUnsorted_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_getQuantile_1, SK_"test_getQuantile_1")
        call test%run(test_getQuantile_2, SK_"test_getQuantile_2")
        call test%run(test_getMedian_RK_1, SK_"test_getMedian_RK_1")
        call test%run(test_getMedian_RK_2, SK_"test_getMedian_RK_2")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getQuantile_1() result(assertion)

        use pm_kind, only: RK, IK
        use pm_arrayVerbose, only: getVerbose

        implicit none

        logical(LK)                 :: assertion
        integer(IK) , parameter     :: nq = 9_IK
        integer(IK) , parameter     :: Weight(lenRnd) = [ 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 0_IK, 1_IK, 2_IK &
                                                        , 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK &
                                                        , 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK &
                                                        , 0_IK, 1_IK, 2_IK, 2_IK, 1_IK, 0_IK, 0_IK, 1_IK ]
        real(RK)    , parameter     :: SortedQuantileProbability(nq) = [0._RK, 0.05_RK, 0.1_RK, 0.25_RK, 0.5_RK, 0.75_RK, 0.9_RK, 0.95_RK, 1._RK]
        real(RK)                    :: Quantile_ref(nq)
        real(RK)                    :: Quantile(nq)
        real(RK)                    :: Difference(nq)

        ! First, getVerbose the input `Point` vector to compute the reference quantiles.

        call getQuantile( SortedQuantileProbability = SortedQuantileProbability & ! LCOV_EXCL_LINE
                        , Point = getVerbose(StdNormRnd1, Weight) & ! LCOV_EXCL_LINE
                        , Quantile = Quantile_ref & ! LCOV_EXCL_LINE
                        )
        if (assertion) return
        assertion = .not. assertion

        ! Now computed the quantiles using the input weighted `Point`.

        call getQuantile( SortedQuantileProbability = SortedQuantileProbability & ! LCOV_EXCL_LINE
                        , Point = StdNormRnd1 & ! LCOV_EXCL_LINE
                        , Weight = Weight & ! LCOV_EXCL_LINE
                        , sumWeight = sum(Weight) & ! LCOV_EXCL_LINE
                        , Quantile = Quantile & ! LCOV_EXCL_LINE
                        )
        if (assertion) return
        assertion = .not. assertion

        Difference = (Quantile - Quantile_ref)
        assertion = all(Difference == 0._RK)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "Quantile_ref   ", Quantile_ref
            write(test%disp%unit,"(*(g0,:,', '))") "Quantile       ", Quantile
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", Difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getQuantile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getQuantile_2() result(assertion)

        use pm_kind, only: RK, IK

        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: nq = 9_IK
        real(RK)    , parameter :: SortedQuantileProbability(nq) = [0._RK, 0.05_RK, 0.1_RK, 0.25_RK, 0.5_RK, 0.75_RK, 0.9_RK, 0.95_RK, 1._RK]
        real(RK)    , parameter :: Quantile_ref(nq) =   [ -2.944284161994900_RK &
                                                        , -1.711516418853700_RK &
                                                        , -1.307688296305270_RK &
                                                        , -.4335920223056840_RK &
                                                        , +.3192067391655020_RK &
                                                        , +1.034693009917860_RK &
                                                        , +1.489697607785470_RK &
                                                        , +2.769437029884880_RK &
                                                        , +3.578396939725760_RK ]
        real(RK)                :: Quantile(nq)
        real(RK)                :: Difference(nq)

        call getQuantile( SortedQuantileProbability = SortedQuantileProbability & ! LCOV_EXCL_LINE
                        , Point = StdNormRnd1 & ! LCOV_EXCL_LINE
                        , Quantile = Quantile & ! LCOV_EXCL_LINE
                        )

        Difference = (Quantile - Quantile_ref)
        assertion = all(Difference == 0._RK)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "Quantile_ref   ", Quantile_ref
            write(test%disp%unit,"(*(g0,:,', '))") "Quantile       ", Quantile
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", Difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getQuantile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMedian_RK_1() result(assertion)

        use pm_err, only: err_type
        implicit none

        logical(LK)                 :: assertion
        real(RK)    , parameter     :: median_ref = 5.477003450000000_RK
        real(RK)                    :: median

        call getMedian(lenArray = ndata, Array = DataUnsorted_RK, median = median)

        assertion = median == median_ref

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "median_ref  =", median_ref
            write(test%disp%unit,"(*(g0,:,' '))") "median      =", median
            write(test%disp%unit,"(*(g0,:,' '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getMedian_RK_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getMedian_RK_2() result(assertion)

        use pm_err, only: err_type
        implicit none

        logical(LK)                 :: assertion
        real(RK)    , parameter     :: median_ref = 5.477415200000000_RK
        real(RK)                    :: median

        call getMedian_RK(lenArray = ndata-1, Array = DataUnsorted_RK(1:ndata-1), median = median)

        assertion = median == median_ref

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") "median_ref  =", median_ref
            write(test%disp%unit,"(*(g0,:,' '))") "median      =", median
            write(test%disp%unit,"(*(g0,:,' '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_getMedian_RK_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_quantile