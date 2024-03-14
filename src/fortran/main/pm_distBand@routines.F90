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
!>  This file contains procedure implementations of [pm_distBand](@ref pm_distBand).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distBand) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
    use pm_arrayUnique, only: getUnique
    use pm_arraySort, only: isAscending
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_quadPack, only: getQuadErr, GK21, weps
    use pm_distGamma, only: setGammaCDF

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBandEpeak_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBandEpeak_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBandEpeak_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBandEpeak_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBandEpeak_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBandEpeak_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBandEpeak_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBandEbreak_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBandEbreak_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBandEbreak_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBandEbreak_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBandEbreak_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBandEbreak_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBandEbreak_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBandZeta_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBandZeta_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBandZeta_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBandZeta_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBandZeta_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBandZeta_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBandZeta_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBandUDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Any_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBandUDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBandUDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBandUDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBandUDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBandUDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Any_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBandUDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBandUCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Any_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandUCDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandUCDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandUCDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandUCDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandUCDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Any_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBandUCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBandMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandMeanDef_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandMeanDef_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandMeanDef_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandMeanDef_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandMeanDef_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandMeanNew_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandMeanNew_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandMeanNew_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandMeanNew_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandMeanNew_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBandMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBandPhoton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FromEnergy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OldB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandPhotonFromEnergyOldB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandPhotonFromEnergyOldB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandPhotonFromEnergyOldB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandPhotonFromEnergyOldB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandPhotonFromEnergyOldB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OldB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NewB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandPhotonFromEnergyNewB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandPhotonFromEnergyNewB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandPhotonFromEnergyNewB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandPhotonFromEnergyNewB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandPhotonFromEnergyNewB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NewB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FromEnergy_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FromPhoton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NewB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandPhotonFromPhotonNewB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandPhotonFromPhotonNewB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandPhotonFromPhotonNewB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandPhotonFromPhotonNewB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandPhotonFromPhotonNewB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NewB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FromPhoton_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBandPhoton_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBandEnergy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FromPhoton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OldB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandEnergyFromPhotonOldB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandEnergyFromPhotonOldB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandEnergyFromPhotonOldB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandEnergyFromPhotonOldB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandEnergyFromPhotonOldB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OldB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NewB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandEnergyFromPhotonNewB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandEnergyFromPhotonNewB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandEnergyFromPhotonNewB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandEnergyFromPhotonNewB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandEnergyFromPhotonNewB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NewB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FromPhoton_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FromEnergy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NewB_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBandEnergyFromEnergyNewB_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBandEnergyFromEnergyNewB_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBandEnergyFromEnergyNewB_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBandEnergyFromEnergyNewB_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBandEnergyFromEnergyNewB_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBand@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NewB_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FromEnergy_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBandEnergy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines ! LCOV_EXCL_LINE