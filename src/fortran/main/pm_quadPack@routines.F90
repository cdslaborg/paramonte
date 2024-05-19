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
!>  This file contains procedure implementations of [pm_quadPack](@ref pm_quadPack).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_quadPack) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_err, only: setAsserted
    use pm_arrayUnique, only: getUnique
    use pm_arraySort, only: isAscending
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayMerge, only: setMerged
    use pm_val2str, only: getStr

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define wcauchy_typer_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure wcauchy_typer_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure wcauchy_typer_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure wcauchy_typer_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure wcauchy_typer_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure wcauchy_typer_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef wcauchy_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define wsin_typer_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure wsin_typer_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure wsin_typer_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure wsin_typer_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure wsin_typer_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure wsin_typer_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef wsin_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define wcos_typer_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure wcos_typer_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure wcos_typer_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure wcos_typer_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure wcos_typer_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure wcos_typer_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef wcos_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNodeWeightGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNodeWeightGKFixed_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNodeWeightGKFixed_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNodeWeightGKFixed_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNodeWeightGKFixed_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNodeWeightGKFixed_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Alloc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNodeWeightGKAlloc_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNodeWeightGKAlloc_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNodeWeightGKAlloc_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNodeWeightGKAlloc_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNodeWeightGKAlloc_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Alloc_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNodeWeightGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK15_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK15_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK15_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK15_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK15_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK15_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK15_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK15_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK15_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK15_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK15_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK15_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK15_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK15_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK15_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK15_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK15_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK15_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK15_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK15_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK21_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK21_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK21_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK21_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK21_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK21_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK21_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK21_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK21_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK21_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK21_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK21_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK21_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK21_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK21_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK21_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK21_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK21_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK21_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK21_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK31_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK31_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK31_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK31_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK31_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK31_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK31_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK31_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK31_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK31_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK31_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK31_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK31_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK31_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK31_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK31_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK31_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK31_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK31_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK31_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK31_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK31_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK41_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK41_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK41_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK41_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK41_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK41_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK41_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK41_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK41_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK41_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK41_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK41_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK41_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK41_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK41_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK41_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK41_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK41_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK41_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK41_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK41_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK41_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK51_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK51_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK51_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK51_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK51_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK51_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK51_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK51_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK51_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK51_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK51_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK51_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK51_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK51_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK51_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK51_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK51_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK51_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK51_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK51_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK51_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK51_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK61_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK61_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK61_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK61_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK61_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK61_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK61_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK61_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK61_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK61_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK61_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK61_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK61_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK61_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK61_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK61_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGK61_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGK61_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGK61_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGK61_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGK61_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK61_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadGK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GKXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGKXX_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGKXX_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGKXX_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGKXX_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGKXX_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGKXX_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGKXX_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGKXX_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGKXX_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGKXX_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGKXX_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGKXX_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGKXX_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGKXX_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGKXX_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getQuadGKXX_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getQuadGKXX_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getQuadGKXX_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getQuadGKXX_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getQuadGKXX_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GKXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadGK_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSeqLimEps_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setSeqLimEps_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setSeqLimEps_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setSeqLimEps_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setSeqLimEps_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setSeqLimEps_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSeqLimEps_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setErrSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setErrSorted_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setErrSorted_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setErrSorted_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setErrSorted_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setErrSorted_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setErrSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChebExpan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setChebExpan_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setChebExpan_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setChebExpan_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setChebExpan_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setChebExpan_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setChebExpan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedQuad_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedQuadQAGD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedQuadQAGD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedQuadQAGD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedQuadQAGD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedQuadQAGD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedQuadQAGS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedQuadQAGS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedQuadQAGS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedQuadQAGS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedQuadQAGS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedQuadQAGP_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedQuadQAGP_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedQuadQAGP_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedQuadQAGP_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedQuadQAGP_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAWC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedQuadQAWC_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedQuadQAWC_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedQuadQAWC_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedQuadQAWC_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedQuadQAWC_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAWC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isFailedQuad_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadErr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK15_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK15_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK15_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK15_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK15_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK15_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK15_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK15_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK15_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK15_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK15_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK15_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK15_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK15_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK15_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK15_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK15_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK15_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK15_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK15_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK21_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK21_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK21_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK21_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK21_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK21_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK21_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK21_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK21_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK21_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK21_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK21_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK21_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK21_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK21_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK21_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK21_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK21_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK21_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK21_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK31_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK31_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK31_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK31_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK31_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK31_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK31_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK31_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK31_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK31_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK31_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK31_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK31_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK31_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK31_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK31_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK31_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK31_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK31_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK31_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK31_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK31_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK41_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK41_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK41_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK41_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK41_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK41_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK41_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK41_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK41_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK41_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK41_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK41_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK41_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK41_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK41_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK41_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK41_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK41_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK41_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK41_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK41_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK41_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK51_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK51_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK51_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK51_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK51_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK51_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK51_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK51_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK51_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK51_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK51_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK51_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK51_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK51_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK51_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK51_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK51_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK51_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK51_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK51_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK51_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK51_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK61_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK61_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK61_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK61_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK61_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK61_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK61_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK61_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK61_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK61_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK61_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK61_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK61_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK61_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK61_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK61_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GK61_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GK61_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GK61_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GK61_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GK61_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK61_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GKXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GKXX_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GKXX_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GKXX_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GKXX_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GKXX_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GKXX_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GKXX_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GKXX_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GKXX_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GKXX_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GKXX_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GKXX_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GKXX_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GKXX_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GKXX_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGD_GKXX_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGD_GKXX_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGD_GKXX_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGD_GKXX_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGD_GKXX_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GKXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadErr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadErr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK15_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK15_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK15_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK15_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK15_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK15_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK15_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK15_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK15_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK15_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK15_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK15_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK15_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK15_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK15_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK15_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK15_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK15_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK15_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK15_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK21_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK21_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK21_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK21_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK21_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK21_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK21_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK21_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK21_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK21_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK21_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK21_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK21_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK21_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK21_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK21_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK21_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK21_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK21_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK21_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK31_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK31_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK31_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK31_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK31_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK31_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK31_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK31_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK31_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK31_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK31_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK31_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK31_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK31_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK31_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK31_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK31_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK31_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK31_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK31_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK31_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK31_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK41_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK41_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK41_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK41_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK41_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK41_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK41_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK41_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK41_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK41_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK41_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK41_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK41_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK41_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK41_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK41_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK41_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK41_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK41_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK41_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK41_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK41_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK51_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK51_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK51_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK51_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK51_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK51_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK51_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK51_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK51_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK51_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK51_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK51_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK51_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK51_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK51_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK51_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK51_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK51_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK51_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK51_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK51_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK51_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK61_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK61_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK61_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK61_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK61_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK61_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK61_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK61_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK61_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK61_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK61_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK61_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK61_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK61_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK61_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK61_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GK61_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GK61_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GK61_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GK61_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GK61_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK61_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GKXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GKXX_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GKXX_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GKXX_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GKXX_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GKXX_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GKXX_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GKXX_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GKXX_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GKXX_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GKXX_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GKXX_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GKXX_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GKXX_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GKXX_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GKXX_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGS_GKXX_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGS_GKXX_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGS_GKXX_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGS_GKXX_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGS_GKXX_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GKXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadErr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadErr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAGP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK15_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK15_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK15_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK15_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK15_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK15_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK15_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK15_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK15_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK15_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK15_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK15_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK15_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK15_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK15_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK15_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK15_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK15_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK15_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK15_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK21_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK21_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK21_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK21_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK21_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK21_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK21_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK21_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK21_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK21_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK21_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK21_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK21_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK21_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK21_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK21_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK21_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK21_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK21_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK21_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK31_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK31_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK31_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK31_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK31_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK31_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK31_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK31_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK31_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK31_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK31_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK31_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK31_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK31_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK31_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK31_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK31_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK31_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK31_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK31_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK31_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK31_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK41_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK41_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK41_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK41_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK41_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK41_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK41_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK41_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK41_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK41_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK41_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK41_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK41_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK41_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK41_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK41_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK41_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK41_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK41_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK41_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK41_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK41_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK51_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK51_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK51_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK51_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK51_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK51_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK51_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK51_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK51_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK51_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK51_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK51_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK51_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK51_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK51_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK51_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK51_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK51_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK51_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK51_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK51_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK51_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK61_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK61_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK61_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK61_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK61_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK61_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK61_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK61_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK61_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK61_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK61_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK61_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK61_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK61_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK61_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK61_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GK61_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GK61_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GK61_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GK61_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GK61_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK61_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GKXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GKXX_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GKXX_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GKXX_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GKXX_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GKXX_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GKXX_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GKXX_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GKXX_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GKXX_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GKXX_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GKXX_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GKXX_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GKXX_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GKXX_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GKXX_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAGP_GKXX_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAGP_GKXX_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAGP_GKXX_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAGP_GKXX_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAGP_GKXX_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GKXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAGP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadErr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getQuadErr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define QAWC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK15_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK15_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK15_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK15_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK15_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK15_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK15_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK15_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK15_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK15_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK15_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK15_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK15_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK15_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK15_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK15_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK15_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK15_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK15_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK15_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK21_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK21_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK21_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK21_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK21_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK21_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK21_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK21_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK21_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK21_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK21_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK21_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK21_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK21_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK21_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK21_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK21_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK21_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK21_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK21_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK31_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK31_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK31_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK31_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK31_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK31_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK31_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK31_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK31_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK31_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK31_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK31_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK31_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK31_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK31_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK31_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK31_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK31_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK31_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK31_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK31_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK31_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK41_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK41_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK41_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK41_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK41_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK41_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK41_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK41_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK41_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK41_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK41_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK41_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK41_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK41_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK41_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK41_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK41_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK41_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK41_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK41_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK41_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK41_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK51_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK51_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK51_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK51_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK51_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK51_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK51_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK51_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK51_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK51_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK51_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK51_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK51_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK51_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK51_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK51_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK51_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK51_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK51_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK51_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK51_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK51_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK61_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK61_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK61_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK61_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK61_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK61_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK61_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK61_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK61_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK61_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK61_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK61_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK61_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK61_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK61_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK61_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GK61_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GK61_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GK61_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GK61_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GK61_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK61_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GKXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GKXX_FF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GKXX_FF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GKXX_FF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GKXX_FF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GKXX_FF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GKXX_FI_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GKXX_FI_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GKXX_FI_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GKXX_FI_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GKXX_FI_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GKXX_IF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GKXX_IF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GKXX_IF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GKXX_IF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GKXX_IF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure QAWC_GKXX_II_RK5
        use pm_kind, only: RKG => RK5
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure QAWC_GKXX_II_RK4
        use pm_kind, only: RKG => RK4
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure QAWC_GKXX_II_RK3
        use pm_kind, only: RKG => RK3
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure QAWC_GKXX_II_RK2
        use pm_kind, only: RKG => RK2
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure QAWC_GKXX_II_RK1
        use pm_kind, only: RKG => RK1
#include "pm_quadPack@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GKXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef QAWC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getQuadErr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines