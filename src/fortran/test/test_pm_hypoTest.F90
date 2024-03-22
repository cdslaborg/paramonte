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

!>  \brief This module contains tests of the module [pm_statest](@ref pm_statest).
!>  \author Amir Shahmoradi

module test_pm_statest

    use pm_statest
    use pm_test, only: test_type, LK
    use pm_kind, only: LK

    implicit none

    private
    public :: setTest
    type(test_type) :: test

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

    integer(IK) , parameter :: lenDistSortedDiff = 200_IK
    real(RK)    , parameter :: DistSortedDiff(*) =  [ 0.00137390044027030_RK & ! LCOV_EXCL_LINE
                                                    , 0.00106568311889532_RK & ! LCOV_EXCL_LINE
                                                    , 0.00261528537896705_RK & ! LCOV_EXCL_LINE
                                                    , 0.00418324458049812_RK & ! LCOV_EXCL_LINE
                                                    , 0.00173791427946435_RK & ! LCOV_EXCL_LINE
                                                    , 0.00683744981995738_RK & ! LCOV_EXCL_LINE
                                                    , 0.00122509129524828_RK & ! LCOV_EXCL_LINE
                                                    , 0.000963897105374811_RK & ! LCOV_EXCL_LINE
                                                    , 0.00659617936766466_RK & ! LCOV_EXCL_LINE
                                                    , 0.00121099678988135_RK & ! LCOV_EXCL_LINE
                                                    , 0.00152627239164915_RK & ! LCOV_EXCL_LINE
                                                    , 5.59784103162375e-05_RK & ! LCOV_EXCL_LINE
                                                    , 0.00396276701851184_RK & ! LCOV_EXCL_LINE
                                                    , 0.00585926454186925_RK & ! LCOV_EXCL_LINE
                                                    , 0.00183541725616831_RK & ! LCOV_EXCL_LINE
                                                    , 0.00213464221769843_RK & ! LCOV_EXCL_LINE
                                                    , 0.00456380640221332_RK & ! LCOV_EXCL_LINE
                                                    , 0.000300222645855053_RK & ! LCOV_EXCL_LINE
                                                    , 0.00113346207247378_RK & ! LCOV_EXCL_LINE
                                                    , 0.00403412924770930_RK & ! LCOV_EXCL_LINE
                                                    , 0.00876326613184120_RK & ! LCOV_EXCL_LINE
                                                    , 0.00240591898034481_RK & ! LCOV_EXCL_LINE
                                                    , 0.00111791498744240_RK & ! LCOV_EXCL_LINE
                                                    , 0.000314496869014858_RK & ! LCOV_EXCL_LINE
                                                    , 0.00166278866932723_RK & ! LCOV_EXCL_LINE
                                                    , 0.00123356868877644_RK & ! LCOV_EXCL_LINE
                                                    , 0.00136174700279745_RK & ! LCOV_EXCL_LINE
                                                    , 0.00180033865772067_RK & ! LCOV_EXCL_LINE
                                                    , 0.000188696505842079_RK & ! LCOV_EXCL_LINE
                                                    , 0.00447901518516258_RK & ! LCOV_EXCL_LINE
                                                    , 0.000649229339662938_RK & ! LCOV_EXCL_LINE
                                                    , 0.00665965285198666_RK & ! LCOV_EXCL_LINE
                                                    , 0.000639392345076595_RK & ! LCOV_EXCL_LINE
                                                    , 0.00439231287788411_RK & ! LCOV_EXCL_LINE
                                                    , 0.00514580427373124_RK & ! LCOV_EXCL_LINE
                                                    , 0.00325615977507276_RK & ! LCOV_EXCL_LINE
                                                    , 0.00262444209665191_RK & ! LCOV_EXCL_LINE
                                                    , 0.00287172623926713_RK & ! LCOV_EXCL_LINE
                                                    , 0.00726524941509710_RK & ! LCOV_EXCL_LINE
                                                    , 0.000679801673850622_RK & ! LCOV_EXCL_LINE
                                                    , 0.000771848959060684_RK & ! LCOV_EXCL_LINE
                                                    , 0.00236687990071305_RK & ! LCOV_EXCL_LINE
                                                    , 0.00279738279147568_RK & ! LCOV_EXCL_LINE
                                                    , 0.00134507052319155_RK & ! LCOV_EXCL_LINE
                                                    , 0.00192351678728353_RK & ! LCOV_EXCL_LINE
                                                    , 0.000403579255582320_RK & ! LCOV_EXCL_LINE
                                                    , 0.00466075625864260_RK & ! LCOV_EXCL_LINE
                                                    , 0.00128562661588161_RK & ! LCOV_EXCL_LINE
                                                    , 0.00240934864593989_RK & ! LCOV_EXCL_LINE
                                                    , 0.00379834817045643_RK & ! LCOV_EXCL_LINE
                                                    , 0.00626106613651900_RK & ! LCOV_EXCL_LINE
                                                    , 0.00155470502430244_RK & ! LCOV_EXCL_LINE
                                                    , 0.00386451831278178_RK & ! LCOV_EXCL_LINE
                                                    , 0.00633652247623751_RK & ! LCOV_EXCL_LINE
                                                    , 0.00212716547361247_RK & ! LCOV_EXCL_LINE
                                                    , 0.000117519040001346_RK & ! LCOV_EXCL_LINE
                                                    , 0.00268124387647972_RK & ! LCOV_EXCL_LINE
                                                    , 0.00330190345409753_RK & ! LCOV_EXCL_LINE
                                                    , 0.00181556281964235_RK & ! LCOV_EXCL_LINE
                                                    , 0.00680782445764760_RK & ! LCOV_EXCL_LINE
                                                    , 0.00344385667182678_RK & ! LCOV_EXCL_LINE
                                                    , 0.00341876847109601_RK & ! LCOV_EXCL_LINE
                                                    , 0.00156480476836285_RK & ! LCOV_EXCL_LINE
                                                    , 0.00392916582362890_RK & ! LCOV_EXCL_LINE
                                                    , 0.000536610168473173_RK & ! LCOV_EXCL_LINE
                                                    , 0.00308379345606979_RK & ! LCOV_EXCL_LINE
                                                    , 0.00382371993615627_RK & ! LCOV_EXCL_LINE
                                                    , 0.00143530941667058_RK & ! LCOV_EXCL_LINE
                                                    , 0.000273761654861704_RK & ! LCOV_EXCL_LINE
                                                    , 0.000218939488190961_RK & ! LCOV_EXCL_LINE
                                                    , 0.00395298813907463_RK & ! LCOV_EXCL_LINE
                                                    , 0.00461683348784969_RK & ! LCOV_EXCL_LINE
                                                    , 0.000587765974902510_RK & ! LCOV_EXCL_LINE
                                                    , 0.00130109023546354_RK & ! LCOV_EXCL_LINE
                                                    , 0.000446086978761473_RK & ! LCOV_EXCL_LINE
                                                    , 0.00102849200817856_RK & ! LCOV_EXCL_LINE
                                                    , 0.00475349511066359_RK & ! LCOV_EXCL_LINE
                                                    , 0.00177004768590394_RK & ! LCOV_EXCL_LINE
                                                    , 0.00118757227387534_RK & ! LCOV_EXCL_LINE
                                                    , 0.00511260207997721_RK & ! LCOV_EXCL_LINE
                                                    , 0.000966468749441840_RK & ! LCOV_EXCL_LINE
                                                    , 0.00153506873793674_RK & ! LCOV_EXCL_LINE
                                                    , 0.000163947408605591_RK & ! LCOV_EXCL_LINE
                                                    , 0.00159865940856585_RK & ! LCOV_EXCL_LINE
                                                    , 0.00319324233871277_RK & ! LCOV_EXCL_LINE
                                                    , 0.00182318600045628_RK & ! LCOV_EXCL_LINE
                                                    , 0.00270039749021500_RK & ! LCOV_EXCL_LINE
                                                    , 0.00233393067978271_RK & ! LCOV_EXCL_LINE
                                                    , 0.000304246250464657_RK & ! LCOV_EXCL_LINE
                                                    , 0.00107643657879242_RK & ! LCOV_EXCL_LINE
                                                    , 0.00149167058888822_RK & ! LCOV_EXCL_LINE
                                                    , 3.22029433108551e-05_RK & ! LCOV_EXCL_LINE
                                                    , 0.000507552164931036_RK & ! LCOV_EXCL_LINE
                                                    , 0.00284976448963692_RK & ! LCOV_EXCL_LINE
                                                    , 0.00261297396275373_RK & ! LCOV_EXCL_LINE
                                                    , 0.00115230561975765_RK & ! LCOV_EXCL_LINE
                                                    , 0.000589761776301434_RK & ! LCOV_EXCL_LINE
                                                    , 0.00254655116972891_RK & ! LCOV_EXCL_LINE
                                                    , 0.000736372303162147_RK & ! LCOV_EXCL_LINE
                                                    , 0.00298938200986687_RK & ! LCOV_EXCL_LINE
                                                    , 0.00115908892632211_RK & ! LCOV_EXCL_LINE
                                                    , 0.0118304957984586_RK & ! LCOV_EXCL_LINE
                                                    , 0.00427409203243478_RK & ! LCOV_EXCL_LINE
                                                    , 0.00408018550920808_RK & ! LCOV_EXCL_LINE
                                                    , 0.00137769679358790_RK & ! LCOV_EXCL_LINE
                                                    , 0.00172460242999972_RK & ! LCOV_EXCL_LINE
                                                    , 0.000397268483832813_RK & ! LCOV_EXCL_LINE
                                                    , 0.00640511077756611_RK & ! LCOV_EXCL_LINE
                                                    , 0.00160078461814450_RK & ! LCOV_EXCL_LINE
                                                    , 0.00116751842007934_RK & ! LCOV_EXCL_LINE
                                                    , 0.00689243594329803_RK & ! LCOV_EXCL_LINE
                                                    , 5.96420573087952e-05_RK & ! LCOV_EXCL_LINE
                                                    , 0.00120353557763264_RK & ! LCOV_EXCL_LINE
                                                    , 0.00542771429675204_RK & ! LCOV_EXCL_LINE
                                                    , 0.00610071142686630_RK & ! LCOV_EXCL_LINE
                                                    , 0.00213356732282444_RK & ! LCOV_EXCL_LINE
                                                    , 0.000281572449424172_RK & ! LCOV_EXCL_LINE
                                                    , 0.00129330702548425_RK & ! LCOV_EXCL_LINE
                                                    , 0.000169205982398446_RK & ! LCOV_EXCL_LINE
                                                    , 0.000101693454083507_RK & ! LCOV_EXCL_LINE
                                                    , 0.000396590512815487_RK & ! LCOV_EXCL_LINE
                                                    , 0.000271783190621933_RK & ! LCOV_EXCL_LINE
                                                    , 0.00207176434429135_RK & ! LCOV_EXCL_LINE
                                                    , 0.00195525703558364_RK & ! LCOV_EXCL_LINE
                                                    , 0.00133407286142151_RK & ! LCOV_EXCL_LINE
                                                    , 0.00146450204323612_RK & ! LCOV_EXCL_LINE
                                                    , 0.00172854780161391_RK & ! LCOV_EXCL_LINE
                                                    , 0.00175309009682001_RK & ! LCOV_EXCL_LINE
                                                    , 0.00591372011950542_RK & ! LCOV_EXCL_LINE
                                                    , 0.00383246845630059_RK & ! LCOV_EXCL_LINE
                                                    , 0.00478224831434315_RK & ! LCOV_EXCL_LINE
                                                    , 0.00140155544663756_RK & ! LCOV_EXCL_LINE
                                                    , 0.00311278264609216_RK & ! LCOV_EXCL_LINE
                                                    , 0.00117361727622267_RK & ! LCOV_EXCL_LINE
                                                    , 0.00197417012646872_RK & ! LCOV_EXCL_LINE
                                                    , 0.00124210536927416_RK & ! LCOV_EXCL_LINE
                                                    , 0.000362752683200629_RK & ! LCOV_EXCL_LINE
                                                    , 0.000939106199745576_RK & ! LCOV_EXCL_LINE
                                                    , 0.000794249402029323_RK & ! LCOV_EXCL_LINE
                                                    , 0.000913827979191928_RK & ! LCOV_EXCL_LINE
                                                    , 0.00720437268872554_RK & ! LCOV_EXCL_LINE
                                                    , 0.000205371175825086_RK & ! LCOV_EXCL_LINE
                                                    , 0.00715773811743015_RK & ! LCOV_EXCL_LINE
                                                    , 0.00251452853672052_RK & ! LCOV_EXCL_LINE
                                                    , 0.000861445629684599_RK & ! LCOV_EXCL_LINE
                                                    , 0.00418790413483983_RK & ! LCOV_EXCL_LINE
                                                    , 0.000370461989028015_RK & ! LCOV_EXCL_LINE
                                                    , 0.00125679370708720_RK & ! LCOV_EXCL_LINE
                                                    , 0.00349303438590243_RK & ! LCOV_EXCL_LINE
                                                    , 0.000724203659035916_RK & ! LCOV_EXCL_LINE
                                                    , 0.000491862014471267_RK & ! LCOV_EXCL_LINE
                                                    , 0.00287196585168448_RK & ! LCOV_EXCL_LINE
                                                    , 0.00185611717084277_RK & ! LCOV_EXCL_LINE
                                                    , 0.000150069773691919_RK & ! LCOV_EXCL_LINE
                                                    , 0.000848534416059033_RK & ! LCOV_EXCL_LINE
                                                    , 0.00123442270407526_RK & ! LCOV_EXCL_LINE
                                                    , 0.00103148223612315_RK & ! LCOV_EXCL_LINE
                                                    , 0.00132400045208614_RK & ! LCOV_EXCL_LINE
                                                    , 0.00435456041987981_RK & ! LCOV_EXCL_LINE
                                                    , 0.00117222896204028_RK & ! LCOV_EXCL_LINE
                                                    , 0.000873409715677065_RK & ! LCOV_EXCL_LINE
                                                    , 0.00512459761807316_RK & ! LCOV_EXCL_LINE
                                                    , 0.000423255105326148_RK & ! LCOV_EXCL_LINE
                                                    , 0.00801979475421000_RK & ! LCOV_EXCL_LINE
                                                    , 0.00452332834787750_RK & ! LCOV_EXCL_LINE
                                                    , 0.000662717679231761_RK & ! LCOV_EXCL_LINE
                                                    , 0.00612674945114589_RK & ! LCOV_EXCL_LINE
                                                    , 0.000241023672676199_RK & ! LCOV_EXCL_LINE
                                                    , 0.00241368458759861_RK & ! LCOV_EXCL_LINE
                                                    , 0.00721401234545482_RK & ! LCOV_EXCL_LINE
                                                    , 0.000617822886748609_RK & ! LCOV_EXCL_LINE
                                                    , 4.20361472183162e-05_RK & ! LCOV_EXCL_LINE
                                                    , 0.000176033863938718_RK & ! LCOV_EXCL_LINE
                                                    , 0.00261140736592547_RK & ! LCOV_EXCL_LINE
                                                    , 0.00161441323061740_RK & ! LCOV_EXCL_LINE
                                                    , 0.00119563944997603_RK & ! LCOV_EXCL_LINE
                                                    , 0.000161186670205815_RK & ! LCOV_EXCL_LINE
                                                    , 0.00186893825670897_RK & ! LCOV_EXCL_LINE
                                                    , 0.0156568954429827_RK & ! LCOV_EXCL_LINE
                                                    , 0.00333188871092183_RK & ! LCOV_EXCL_LINE
                                                    , 0.00460140509714901_RK & ! LCOV_EXCL_LINE
                                                    , 0.00106350049875970_RK & ! LCOV_EXCL_LINE
                                                    , 0.000278484374323207_RK & ! LCOV_EXCL_LINE
                                                    , 0.000109843896551887_RK & ! LCOV_EXCL_LINE
                                                    , 0.00517794656024040_RK & ! LCOV_EXCL_LINE
                                                    , 0.000764568679501587_RK & ! LCOV_EXCL_LINE
                                                    , 0.00106947825226744_RK & ! LCOV_EXCL_LINE
                                                    , 0.00475394751147906_RK & ! LCOV_EXCL_LINE
                                                    , 0.00700310227105183_RK & ! LCOV_EXCL_LINE
                                                    , 0.00280107346803815_RK & ! LCOV_EXCL_LINE
                                                    , 0.00594747756623304_RK & ! LCOV_EXCL_LINE
                                                    , 0.00268586958531547_RK & ! LCOV_EXCL_LINE
                                                    , 0.00120417364527792_RK & ! LCOV_EXCL_LINE
                                                    , 0.00344542024642325_RK & ! LCOV_EXCL_LINE
                                                    , 0.00217032516670990_RK & ! LCOV_EXCL_LINE
                                                    , 0.00170145556378476_RK & ! LCOV_EXCL_LINE
                                                    , 0.00159968663729304_RK & ! LCOV_EXCL_LINE
                                                    , 0.00153233985478085_RK & ! LCOV_EXCL_LINE
                                                    , 0.00219750724001799_RK & ! LCOV_EXCL_LINE
                                                    , 0.00230064745992042_RK & ! LCOV_EXCL_LINE
                                                    ]
    integer(IK) , parameter :: SampleSize(*) =  [ 1_IK &
                                                , 101_IK &
                                                , 102_IK &
                                                , 103_IK &
                                                , 104_IK &
                                                , 105_IK &
                                                , 106_IK &
                                                , 107_IK &
                                                , 108_IK &
                                                , 109_IK &
                                                , 110_IK &
                                                ]
    integer(IK) , parameter :: lenSampleSize = size(SampleSize,kind=IK)
    real(RK)    , parameter :: ProbKS_REF(*) =  [ 1._RK &
                                                , 0.482474176998994_RK &
                                                , 0.370817196700337_RK &
                                                , 0.493146633387423_RK &
                                                , 0.702120402082310_RK &
                                                , 0.671010642014021_RK &
                                                , 0.713612126474551_RK &
                                                , 0.683487750289677_RK &
                                                , 0.559321682773637_RK &
                                                , 0.530113153332594_RK &
                                                , 0.571577075807556_RK &
                                                ] ! The KS two-sample test probabilities.
    real(RK)    , parameter :: StatKS_REF(*) =  [ 0._RK &
                                                , 0.162352941176471_RK &
                                                , 0.176470588235294_RK &
                                                , 0.159502262443439_RK &
                                                , 0.134615384615385_RK &
                                                , 0.137518142235123_RK &
                                                , 0.132075471698113_RK &
                                                , 0.134870719776380_RK &
                                                , 0.148148148148148_RK &
                                                , 0.150841750841751_RK &
                                                , 0.145454545454545_RK &
                                                ] ! The KS two-sample test statistics.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test_Uniformity()
        return
        call test%run(test_doKS1_1, SK_"test_doKS1_1")
        call test%run(test_doSortedKS2_1, SK_"test_doSortedKS2_1")
        call test%run(test_doUniformKS1_1, SK_"test_doUniformKS1_1")
        call test%run(test_getCumDenComKS_1, SK_"test_getCumDenComKS_1")
        call test%run(test_getSampleDisparity_1, SK_"test_getSampleDisparity_1")
        call test%run(test_getSampleDisparity_2, SK_"test_getSampleDisparity_2")
       !call test%run(test_performance_1, SK_"test_performance_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doKS1_1() result(assertion)

        use pm_distNorm, only: getNormCDF
        use pm_kind, only: RK, IK

        implicit none

        logical(LK)             :: assertion
        real(RK)    , parameter :: TOLERANCE = 1.e-7_RK

        real(RK), parameter     :: probKS_ref = .3763758622852317E-01_RK
        real(RK), parameter     :: statKS_ref = .1955719390701096_RK
        real(RK)                :: statKS
        real(RK)                :: probKS
        real(RK)                :: difference

        call doKS1  ( np = lenRnd &
                    , Point = StdNormRnd1 &
                    , getCDF = getCDF &
                    , statKS = statKS &
                    , probKS = probKS &
                    )

        difference = abs( (probKS - probKS_ref) / probKS_ref )
        assertion = difference < TOLERANCE

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "probKS     :", probKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        difference = abs( (statKS - statKS_ref) / statKS_ref )
        assertion = assertion .and. difference < TOLERANCE
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "statKS     :", statKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    contains

        PURE function getCDF(x) result(cdf)
            implicit none
            real(RK)    , intent(in)    :: x
            real(RK)                    :: cdf
            cdf = getNormCDF(x)
        end function

    end function test_doKS1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCumDenComKS_1() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        integer(IK)                 :: i
        logical(LK)                 :: assertion
        real(RK)    , parameter     :: TOLERANCE = 1.e-11_RK
        real(RK)    , parameter     :: Score(*) = [(0.01_RK * real(i,RK), i = -1, 100)]
        real(RK)                    :: Difference(size(Score))
        real(RK)                    :: ProbKS_ref(size(Score))
        real(RK)                    :: ProbKS(size(Score))
        do i = 1, size(Score)
            ProbKS_ref(i) = getProbKS(Score(i))
            ProbKS(i) = getCumDenComKS(Score(i))
            Difference(i) = abs(ProbKS(i) - ProbKS_ref(i))
        end do
        assertion = all(Difference < TOLERANCE)
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "ProbKS_ref         :", ProbKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "ProbKS             :", ProbKS
            write(test%disp%unit,"(*(g0,:,', '))") "Difference         :", Difference
            write(test%disp%unit,"(*(g0,:,', '))") "maxval(Difference) :", maxval(Difference)
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if
    end function test_getCumDenComKS_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Compare performance of [getCumDenComKS()](@ref pm_statest)` with [pm_statest::getProbKS()](@ref pm_statest::getProbKS).
    !> Results indicated that [pm_statest::getProbKS()](@ref pm_statest::getProbKS) is approximately 2-3 times faster.
    function test_performance_1() result(assertion)
        use pm_kind, only: RK, IK
        use pm_timer, only: timerSYS_type
        use pm_err, only: err_type
        implicit none
        integer(IK)                 :: i
        integer(IK) , parameter     :: NSCORE = 10**6_IK
        logical(LK)                 :: assertion
        real(RK)                    :: Score(NSCORE)
        real(RK)                    :: ProbKS(NSCORE), dummy
        type(timerSYS_type)            :: timer
        timer = timerSYS_type()
        assertion = .true.
        do i = -1, NSCORE - 2
            Score(i+2) = 0.01_RK * real(i,RK)
        end do
        timer%start = timer%time()
        do i = 1, NSCORE
            ProbKS(i) = getCumDenComKS(Score(i))
        end do
        timer%delta = timer%time(since = timer%start)
        write(test%disp%unit,"(*(g0,:,', '))")
        write(test%disp%unit,"(*(g0,:,', '))") "time(getCumDenComKS):", timer%delta
        dummy = timer%delta
        timer%start = timer%time()
        do i = 1, NSCORE
            ProbKS(i) = getProbKS(Score(i))
        end do
        timer%delta = timer%time(since = timer%start)
        write(test%disp%unit,"(*(g0,:,', '))") "time(getProbKS):", timer%delta
        write(test%disp%unit,"(*(g0,:,', '))") "ratio(getProbKS/getCumDenComKS):", timer%delta / dummy
        write(test%disp%unit,"(*(g0,:,', '))")
    end function test_performance_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doUniformKS1_1() result(assertion)

        use pm_kind, only: RK, IK

        implicit none

        logical(LK)             :: assertion
        real(RK)    , parameter :: TOLERANCE = 1.e-7_RK

        real(RK), parameter     :: probKS_ref = .1982797523608350_RK
        real(RK), parameter     :: statKS_ref = .1491359075319200_RK
        real(RK)                :: statKS
        real(RK)                :: probKS
        real(RK)                :: difference

        call doUniformKS1   ( np = lenRnd &
                            , Point = UnifRnd &
                            , statKS = statKS &
                            , probKS = probKS &
                            )

        difference = abs( (probKS - probKS_ref) / probKS_ref )
        assertion = difference < TOLERANCE
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "probKS     :", probKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        difference = abs( (statKS - statKS_ref) / statKS_ref )
        assertion = assertion .and. difference < TOLERANCE
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "statKS     :", statKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_doUniformKS1_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_doSortedKS2_1() result(assertion)

        use pm_kind, only: RK, IK
        use pm_arraySort, only: setSorted

        implicit none

        logical(LK)             :: assertion
        integer(IK) , parameter :: lenRnd = 50_IK
        real(RK)    , parameter :: TOLERANCE = 1.e-12_RK

        real(RK), parameter     :: probKS_ref = 0.056045859714425_RK
        real(RK), parameter     :: statKS_ref = 0.260000000000000_RK
        real(RK)                :: difference, probKS, statKS, lambda

        assertion = .true._LK

        call setSorted(StdNormRnd1)
        call setSorted(StdNormRnd2)

        call doSortedKS2( sortedSample1 = StdNormRnd1 &
                        , SortedSample2 = StdNormRnd2 &
                        , probKS = probKS &
                        , statKS = statKS &
                        , lambda = lambda &
                        )

        difference = 2 * abs(probKS - probKS_ref) / (probKS_ref + probKS)
        assertion = assertion .and. difference < TOLERANCE

        !write(*,*) "statKS, diff", statKS, abs(statKS - statKS_ref)
        !write(*,*) "probKS, diff", probKS, abs(probKS - probKS_ref)
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "probKS_ref :", probKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "probKS     :", probKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))") "TOLERANCE  :", TOLERANCE
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion)

        difference = 2 * abs(statKS - statKS_ref) / (statKS_ref + statKS)
        assertion = assertion .and. difference < TOLERANCE
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "statKS_ref :", statKS_ref
            write(test%disp%unit,"(*(g0,:,', '))") "statKS     :", statKS
            write(test%disp%unit,"(*(g0,:,', '))") "difference :", difference
            write(test%disp%unit,"(*(g0,:,', '))") "TOLERANCE  :", TOLERANCE
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

    end function test_doSortedKS2_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Any value of `SampleSize` smaller than one or larger than the size of the input `Point` should lead to an error.
    function test_getSampleDisparity_1() result(assertion)

        use pm_kind, only: RK, IK

        implicit none

        logical(LK)                 :: assertion
        real(RK)                    :: ProbKS(2)
        real(RK)                    :: StatKS(2)
        real(RK)                    :: Lambda(2)

        assertion = .true._LK
        call setSampleDisparityAuto ( SortedSample = DistSortedDiff & ! LCOV_EXCL_LINE
                                    , SampleSize = [3, 4] & ! [0_IK, 1_IK] & ! LCOV_EXCL_LINE
                                    , StatKS = StatKS & ! LCOV_EXCL_LINE
                                    , ProbKS = ProbKS & ! LCOV_EXCL_LINE
                                    , Lambda = Lambda & ! LCOV_EXCL_LINE
                                    )
        call test%assert(assertion)

        call setSampleDisparityAuto ( SortedSample = DistSortedDiff & ! LCOV_EXCL_LINE
                                    , SampleSize = [ 1, lenDistSortedDiff] & ! [1_IK, lenDistSortedDiff + 1_IK] & ! LCOV_EXCL_LINE
                                    , StatKS = StatKS & ! LCOV_EXCL_LINE
                                    , ProbKS = ProbKS & ! LCOV_EXCL_LINE
                                    , Lambda = Lambda & ! LCOV_EXCL_LINE
                                    )
        call test%assert(assertion)

    end function test_getSampleDisparity_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Any value of `SampleSize` smaller than one or larger than the size of the input `Point` should lead to an error.
    function test_getSampleDisparity_2() result(assertion)

        use pm_kind, only: RK, IK

        implicit none

        logical(LK)                 :: assertion
        real(RK)    , allocatable   :: Difference(:)
        real(RK)    , parameter     :: TOLERANCE = 1.e-11_RK
        real(RK)                    :: ProbKS(lenSampleSize)
        real(RK)                    :: StatKS(lenSampleSize)
        real(RK)                    :: Lambda(lenSampleSize)

        assertion = .true.
        call setSampleDisparityAuto ( SortedSample = DistSortedDiff & ! LCOV_EXCL_LINE
                                    , SampleSize = SampleSize & ! LCOV_EXCL_LINE
                                    , ProbKS = ProbKS & ! LCOV_EXCL_LINE
                                    , StatKS = StatKS & ! LCOV_EXCL_LINE
                                    , Lambda = Lambda & ! LCOV_EXCL_LINE
                                    )
        call test%assert(assertion)
        Difference = 2 * abs(ProbKS - ProbKS_REF) / (ProbKS_REF + ProbKS)
        assertion = assertion .and. all(Difference < TOLERANCE)
        call test%assert(assertion)

        !write(*,*) "ProbKS, diff", ProbKS, abs(ProbKS - ProbKS_REF)
        !write(*,*) "StatKS, diff", StatKS, abs(StatKS - StatKS_REF)
        !write(*,*) "Lambda, diff", Lambda, abs(Lambda - StatKS_REF)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "ProbKS_REF:", ProbKS_REF
            write(test%disp%unit,"(*(g0,:,', '))") "ProbKS     :", ProbKS
            write(test%disp%unit,"(*(g0,:,', '))") "Difference :", Difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        Difference = abs(StatKS - StatKS_REF) !* 2 / (StatKS_REF + StatKS)
        assertion = assertion .and. all(Difference < TOLERANCE)
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "StatKS_REF :", StatKS_REF
            write(test%disp%unit,"(*(g0,:,', '))") "StatKS     :", StatKS
            write(test%disp%unit,"(*(g0,:,', '))") "Difference :", Difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getSampleDisparity_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Uniformity()

        use pm_kind, only: IK, LK, RK, RKS
        use pm_sysPath, only: getDirCurrent
        use pm_val2str, only: getStr
        use pm_matrixInit, only: getMatInit
        use pm_matrixInit, only: uppLowDia
        use pm_matrixDet, only: setMatDetSqrtLog, uppDia, lowDia, transHerm
        use pm_distUnifEll, only: setUnifEllRand
        use pm_statest, only: setSampleDisparityAuto
        !use pm_statest, only: doSortedKS2
        use pm_arrayRange, only: getRange
        use pm_mathCumSum, only: setCumSum
        use pm_distanceEuclid, only: getDisEuclid, euclidsq
        use pm_knn, only: setKnnSorted
        use pm_knn, only: setDisSortedExpDiff
        use pm_distanceEuclid, only: setDisMatEuclid
        use pm_distanceEuclid, only: rdpack, uppLowDia, euclidu
        use pm_ellipsoid, only: getLogVolUnitBall
        use pm_timer, only: timerSYS_type

        type(timerSYS_type)             :: timer
        integer(IK)                     :: ipnt, idim, iref
        integer(IK)                     :: npntMinusOne
        integer(IK)                     :: npntHalf
        integer(IK)                     :: npnt
        integer(IK)                     :: ndim
        real(RK)                        :: minDistanceSq
        real(RK)                        :: sqrtDetInvCovMat
        real(RK)                        :: rho
        real(RK)        , allocatable   :: refDistMat(:,:)
       !real(RK)        , allocatable   :: PairDistMat(:,:)
       !real(RK)        , allocatable   :: PairDistVec(:)
        real(RK)        , allocatable   :: Reference(:)
        real(RK)        , allocatable   :: Center(:)
        real(RK)        , allocatable   :: Volume(:)
        real(RK)        , allocatable   :: VolumeWeighted(:)
        real(RK)        , allocatable   :: disSortedExpDiff(:)
        real(RK)        , allocatable   :: sample(:,:), SampleDum(:,:) ! ndim, npnt
        real(RK)        , allocatable   :: chocov(:,:)
        real(RK)        , allocatable   :: DistanceSq(:)
        real(RK)        , allocatable   :: ProbKS(:), StatKS(:), LambdaKS(:)
        real(RK)        , allocatable   :: CumSumProbKS(:)
        integer(IK)     , allocatable   :: ProxyIndex(:,:)
        integer(IK)     , allocatable   :: IndexNN(:)
        integer(IK)     , allocatable   :: IndexKS(:)
        integer(IK)                     :: fileUnit
        integer(IK)                     :: counter
        integer(IK)                     :: info
        character(*, SK), parameter     :: CWD = SK_"/mnt/d/Dropbox/Projects/20220601_nsfCareer/codes/logz"
        character(:, SK), allocatable   :: filePath

        rho = 0.9_RK
        ndim = 2_IK
        npnt = 2000_IK
        Center = [(0._RK, idim = 1, ndim)]

        ! Get the Covariance.

        chocov = getMatInit([ndim, ndim + 1_IK], uppLowDia, rho, rho, 1._RK, doff = 1_IK)

        ! Get the Cholesky lower factor.

        call setMatDetSqrtLog(chocov(:, 2 : ndim + 1), uppDia, sqrtDetInvCovMat, info, chocov(:, 1 : ndim), transHerm)
        if (info /= 0_IK) error stop "setMatChol() failed."

        ! Normalize the volume to 1.

        sqrtDetInvCovMat = exp(-sqrtDetInvCovMat / ndim)
        do idim = 1, ndim
            chocov(idim + 1 : ndim, idim) = chocov(idim + 1 : ndim, idim) * sqrtDetInvCovMat! / exp(getLogVolUnitBall(ndim)/ndim)
        end do

        ! Get the Cholesky lower factor.

        allocate(sample(ndim, npnt))
        call setUnifEllRand(sample, chocov, lowDia)

        ! Find the closest point to the center.

        minDistanceSq = +huge(minDistanceSq)
        allocate(DistanceSq(npnt))
        do ipnt = 1_IK, npnt
            DistanceSq(ipnt) = getDisEuclid(center, sample(1:ndim, ipnt), euclidsq)
            if (minDistanceSq > DistanceSq(ipnt)) then
                minDistanceSq = DistanceSq(ipnt)
                iref = ipnt
            end if
        end do
        Reference = sample(1:ndim, iref)

        ! Create a donut. \warning this may not be consistent with the above center point.

        if (.false.) then
            allocate(SampleDum, mold = sample)
            counter = 0_IK
            !DistanceSq = sum(sample**2, dim = 1)
            !call setDisEuclid(distancesq, sample, spread(center, 2_ik, npnt))
            do ipnt = 1_IK, npnt
                !DistanceSq(ipnt) = getDisEuclid([(0._RK, idim = 1, ndim)], sample(:, ipnt), euclidsq)
                if (DistanceSq(ipnt) < 0.4_RK**2 .or. DistanceSq(ipnt) > 0.6_RK**2) then
                    counter = counter + 1_IK
                    SampleDum(1:ndim, counter) = sample(1:ndim, ipnt)
                end if
            end do
            npnt = counter
            sample = SampleDum(1:ndim,1:counter)
            deallocate(SampleDum)
        end if

        ! Define the sample sizes for which the KS test must be run.

        npntMinusOne = npnt - 1
        npntHalf = npntMinusOne / 2
        allocate(disSortedExpDiff(0:npnt-1), IndexNN(npnt))

        ! Write the points to the output file

        filePath = CWD//SK_"/sample.txt"
        open(newunit = fileUnit, file = filePath, status = "replace")
        !print *, filePath
        do ipnt = 1, npnt
            !write(fileUnit, "("//getStr(ndim)//"(g0,:,','))") sample
            write(fileUnit, "(*(g0,:,','))") sample(:, ipnt)
        end do
        close(fileUnit)

        ! Sort the sample from the closest to the farthest with respect to the best reference point.

        !print *, maxval(Reference)
        !Center = [(0._RK, idim = 1, ndim)]
        call setDisSortedExpDiff(sample, Center, disSortedExpDiff(0:npnt-1), IndexNN)
        sample(1:ndim, 1:npnt) = sample(1:ndim, IndexNN)

        ! Compute the pairwise distance matrix.

        !!allocate(PairDistVec((npnt**2 - npnt) / 2))
        !!call setPairDistSqVec(PairDistVec, sample)
        !!PairDistVec = sqrt(PairDistVec)
        !allocate(PairDistMat(npnt, npnt))
        !call setPairDistSqMat(PairDistMat, sample)
        !PairDistMat = sqrt(PairDistMat)

        ! Compute the pairwise distance differences matrix.

        allocate(refDistMat(npnt, npnt), ProxyIndex(npnt, npnt))
        call setDisMatEuclid(refDistMat, rdpack, uppLowDia, sample, euclidu)
        call setKnnSorted(refDistMat, ProxyIndex)

        ! Estimate the volume.

        call doCrossKS()
        call doAutoKS()

    contains

        subroutine doCrossKS()

            use pm_statest, only: doKS2
            integer(IK) :: counter
            integer(IK) :: lenIndexKS
            integer(IK) :: maxIndexKS
            integer(IK) :: oldIndexKS
            integer(IK) :: index, knn
            integer(IK) :: skip, jpnt, kpnt
            integer(IK) :: jpntIndexInRefDistMat
            real(RK)    :: VolSortedDiff(npnt, npnt)
            real(RK)    :: VolSortedDiffAlt((npntMinusOne**2 - npntMinusOne))
            real(RK)    :: lambdaKS, statKS
            real(RK)    :: margin
            real(RK)    , allocatable :: ProbKS(:), Volume(:), VolumeWeighted(:), CumSumProbKS(:)

            ! Compute the distance difference matrix.

            do ipnt = 1_IK, npnt
                VolSortedDiff(1:npnt, ipnt) = refDistMat(ProxyIndex(1:npnt, ipnt), ipnt)**ndim
                do jpnt = npnt, 2_IK, -1_IK
                    VolSortedDiff(jpnt, ipnt) = VolSortedDiff(jpnt, ipnt) - VolSortedDiff(jpnt - 1, ipnt)
                end do
            end do
            !VolSortedDiff(1:npnt, 1) = disSortedExpDiff

            skip = 1_IK
            IndexKS = [(idim, idim = 8_IK, npnt, skip)]
            maxIndexKS = maxval(IndexKS, dim = 1)
            lenIndexKS = size(IndexKS,1,IK)
            allocate( ProbKS(maxIndexKS) &
                    , Volume(maxIndexKS) &
                    , VolumeWeighted(maxIndexKS) &
                    )

            ! Get the KS tests

            timer%start = timer%time()

            ! First find all bounded neighbors within the specified index.

            loopOverReference: do ipnt = 1_IK, 1!npnt

                oldIndexKS = 1_IK
                do index = 1_IK, lenIndexKS
                    knn = IndexKS(index)
                    !print *, ipnt, knn
                    counter = 0_IK
                    do jpnt = 1_IK, knn - 1_IK
                        ! find the farthest neighbor for subcircle jpnt within reference and knn.
                        ! The last point on the border has no neighbor within the super-circle, so do not loop over it.
                        ! The jpnt'th neighbor with respect to the ipnt'th reference point has
                        ! distance refDistMat(ProxyIndex(jpnt, ipnt), ipnt) from the reference.
                        ! ProxyIndex(jpnt, ipnt) is the column of jpnt'th neighbor of ipnt reference in refDistMat.
                        jpntIndexInRefDistMat = ProxyIndex(jpnt, ipnt) ! Note that `refDistMat` is symmetric.
                        margin = refDistMat(ProxyIndex(knn, ipnt), ipnt) - refDistMat(jpntIndexInRefDistMat, ipnt) ! distance to border.
                        !print *, margin
                        !if (ProxyIndex(knn, ipnt) < jpntIndexInRefDistMat) error stop "nothing else matters."
                        if (margin <= 0._RK) error stop "negative margin = "//getStr(margin)
                        ! Loop over the neighbors of jpnt'th reference within the neighborhood of ipnt.
                        ! the first point is the reference jpnt, so ignore it.
                        ! But we do not know which point is the last neighbor of jpnt that is closer than `margin` to the jpnt'th reference.
                        ! So loop over all points in the column. This is effectively a linear search,
                        ! but can be done much faster using binary search.
                        do kpnt = 2_IK, npnt
                            if (refDistMat(ProxyIndex(kpnt, jpntIndexInRefDistMat), jpntIndexInRefDistMat) > margin) exit
                            counter = counter + 1_IK
                            VolSortedDiffAlt(counter) = VolSortedDiff(kpnt, jpnt)
                        end do
                        !print *, counter, refDistMat(ProxyIndex(knn, ipnt), ipnt), refDistMat(ProxyIndex(jpnt, ipnt), ipnt), margin
                    end do
                    if (counter == 0_IK) error stop "counter == 0_IK"
                    call doKS2(VolSortedDiff(1:knn, ipnt), VolSortedDiffAlt(1:counter), ProbKS(knn), statKS, lambdaKS)
                    ProbKS(oldIndexKS : knn - 1) = ProbKS(knn)
                    oldIndexKS = knn + 1_IK
                end do
                write(*,*) "time = ", timer%time(since = timer%start)

                ! Write the KS probabilities to the output file

                filePath = CWD//SK_"/CrossProbKS.txt"
                open(newunit = fileUnit, file = filePath, status = "replace")
                do index = 1, maxIndexKS
                    write(fileUnit, "(*(g0,:,','))") index, ProbKS(index)
                end do
                close(fileUnit)

                ! Compute the mean volume

                call setCumSum(Volume, VolSortedDiff(1:maxIndexKS,1))
                Volume = npnt * Volume / getRange(1_IK, maxIndexKS)

                allocate(CumSumProbKS, mold = ProbKS)
                call setCumSum(CumSumProbKS, ProbKS)
                call setCumSum(VolumeWeighted, VolSortedDiff(1:maxIndexKS,1) * ProbKS)
                VolumeWeighted = npnt * VolumeWeighted / CumSumProbKS

                filePath = CWD//SK_"/CrossVolumeKS.txt"
                open(newunit = fileUnit, file = filePath, status = "replace")
                do index = 1, maxIndexKS
                    write(fileUnit, "(*(g0,:,','))") Volume(index), VolumeWeighted(index)
                end do
                close(fileUnit)

            end do loopOverReference

        end subroutine

        subroutine doAutoKS()

            allocate(Volume(npntHalf), VolumeWeighted(npntHalf), ProbKS(npntHalf), StatKS(npntHalf), LambdaKS(npntHalf))
            IndexKS = [(idim, idim = 1_IK, npntHalf, 1_IK)]

            ! Get the KS tests

            timer%start = timer%time()
            iref = 1_IK
            disSortedExpDiff(0:npntMinusOne) = refDistMat(ProxyIndex(:,iref), iref)**ndim
            disSortedExpDiff(1:npntMinusOne) = disSortedExpDiff(1:npntMinusOne) - disSortedExpDiff(0:npntMinusOne-1)
            call setSampleDisparityAuto(disSortedExpDiff(1:npntMinusOne), IndexKS, ProbKS, StatKS, LambdaKS)
            write(*,*) "time = ", timer%time(since = timer%start)

            ! Write the KS probabilities to the output file

            filePath = CWD//SK_"/AutoProbKS.txt"
            open(newunit = fileUnit, file = filePath, status = "replace")
            do ipnt = 1, size(IndexKS, 1, IK)
                write(fileUnit, "(*(g0,:,','))") 2*IndexKS(ipnt), ProbKS(ipnt)
            end do
            close(fileUnit)

            ! Compute the mean volume

            call setCumSum(Volume, disSortedExpDiff(1:npntHalf))
            Volume = npntMinusOne * Volume / getRange(1_IK, size(Volume, kind = IK))

            allocate(CumSumProbKS, mold = ProbKS)
            call setCumSum(CumSumProbKS, ProbKS)
            call setCumSum(VolumeWeighted, disSortedExpDiff(1:npntHalf) * ProbKS)
            VolumeWeighted = npntMinusOne * VolumeWeighted / CumSumProbKS

            filePath = CWD//SK_"/AutoVolumeKS.txt"
            open(newunit = fileUnit, file = filePath, status = "replace")
            do ipnt = 1, size(Volume, kind = IK)
                write(fileUnit, "(*(g0,:,','))") Volume(ipnt), VolumeWeighted(ipnt)
            end do
            close(fileUnit)

        end subroutine

    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_statest