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
!>  This file contains procedure implementations of [pm_ellipsoid](@ref pm_ellipsoid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_ellipsoid) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_matrixCopy, only: setMatCopy, rdpack, uppDia
    use pm_mathFactorial, only: getLogFactorial
    use pm_matrixDet, only: getMatDetSqrtLog

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVolUnitBall_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Iter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getVolUnitBallIter_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getVolUnitBallIter_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getVolUnitBallIter_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getVolUnitBallIter_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getVolUnitBallIter_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Iter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVolUnitBall_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVolUnitBall_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Iter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setVolUnitBallIter_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setVolUnitBallIter_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setVolUnitBallIter_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setVolUnitBallIter_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setVolUnitBallIter_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Iter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVolUnitBall_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogVolUnitBall_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Gamm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogVolUnitBallGamm_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogVolUnitBallGamm_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogVolUnitBallGamm_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogVolUnitBallGamm_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogVolUnitBallGamm_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Gamm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogVolUnitBall_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogVolUnitBall_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Gamm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogVolUnitBallGamm_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogVolUnitBallGamm_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogVolUnitBallGamm_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogVolUnitBallGamm_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogVolUnitBallGamm_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Gamm_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Iter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogVolUnitBallIter_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogVolUnitBallIter_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogVolUnitBallIter_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogVolUnitBallIter_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogVolUnitBallIter_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Iter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogVolUnitBall_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogVolEll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogVolEll_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogVolEll_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogVolEll_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogVolEll_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogVolEll_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogVolEll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCountMemberEll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sph_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCountMemberSphOrg_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCountMemberSphOrg_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCountMemberSphOrg_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCountMemberSphOrg_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCountMemberSphOrg_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCountMemberSphCen_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCountMemberSphCen_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCountMemberSphCen_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCountMemberSphCen_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCountMemberSphCen_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sph_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCountMemberEllOrg_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCountMemberEllOrg_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCountMemberEllOrg_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCountMemberEllOrg_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCountMemberEllOrg_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCountMemberEllCen_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCountMemberEllCen_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCountMemberEllCen_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCountMemberEllCen_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCountMemberEllCen_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ell_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCountMemberEll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isMemberEll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sph_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMemberSphOrg_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMemberSphOrg_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMemberSphOrg_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMemberSphOrg_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMemberSphOrg_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMemberSphCen_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMemberSphCen_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMemberSphCen_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMemberSphCen_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMemberSphCen_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sph_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMemberEllOrg_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMemberEllOrg_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMemberEllOrg_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMemberEllOrg_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMemberEllOrg_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Cen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isMemberEllCen_RK5
        use pm_kind, only: RKG => RK5
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isMemberEllCen_RK4
        use pm_kind, only: RKG => RK4
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isMemberEllCen_RK3
        use pm_kind, only: RKG => RK3
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isMemberEllCen_RK2
        use pm_kind, only: RKG => RK2
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isMemberEllCen_RK1
        use pm_kind, only: RKG => RK1
#include "pm_ellipsoid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Cen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ell_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isMemberEll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
