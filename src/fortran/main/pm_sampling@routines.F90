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
!>  This file contains the implementations of the generic interfaces in [pm_sampling]@(ref pm_sampling).
!>
!>  \test
!>  [test_pm_sampling](@ref test_pm_sampling)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampling) routines

    !   \bug
    !   Avoid Intel ifort bug for too many `use` statements in a submodule by placing them all in the submodule header.
#if CHECK_ENABLED
    use pm_err, only: setAsserted, getFine
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

   !use pm_str, only: UNDEFINED
    use pm_sampleCor, only: getCor
    use pm_statest, only: getProbKS
    use pm_sampleMean, only: setMean
    use pm_distanceKolm, only: ascending
    use pm_distanceKolm, only: getDisKolm
    use pm_arrayVerbose, only: getVerbose
    use pm_matrixCopy, only: setMatCopy, rdpack, transHerm
    use pm_sampleCov, only: setCovMean, setCov, lowDia, uppDia
    use pm_distGeomCyclic, only: isFailedGeomCyclicFit
    use pm_arrayComplement, only: getComplementRange
    use pm_parallelism, only: setForkJoinScaling
    use pm_arrayUnique, only: setUnique
    use pm_arraySort, only: setSorted
    use pm_kind, only: modelr_type
    use pm_container, only: css_type
    use pm_io, only: getErrTableRead, trans
    use, intrinsic :: iso_fortran_env, only: output_unit
    use pm_arrayFill, only: getFilled
    use pm_mathCumSum, only: setCumSum
    use pm_distUnif, only: getUnifRand
    use pm_arrayInit, only: setCoreHalo
    use pm_arrayCenter, only: setCentered
    use pm_arrayReplace, only: getReplaced
    use pm_io, only: INDENT, display_type
    use pm_distUnif, only: setUnifRand
    use pm_distUnif, only: getUnifRandState
    use pm_distUnif, only: getUnifRandStateSize
    use pm_mathNumSys, only: getCountDigit
    use pm_matrixInit, only: getMatInit, setMatInit, uppLowDia
    use pm_parallelism, only: setImageCount, isFailedImage
    use pm_sysPath, only: LENPATH => MAX_LEN_FILE_PATH
    use pm_strASCII, only: getStrLower, setStrLower
    use pm_sysShell, only: isFailedGetShellWidth
    use pm_option, only: getOption
    use pm_val2str, only: getStr
    use pm_err, only: err_type
    use pm_except, only: setNAN
    use pm_except, only: setInfPos
    use pm_arraySplit, only: setSplit
    use pm_arrayResize, only: setResized
    use pm_arrayRebind, only: setRebound
    use pm_matrixChol, only: setMatChol
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_dateTime, only: getDateTime
   !use pm_except, only: isInfNeg, setInfNeg
    use pm_arrayRefine, only: setRefined
    use pm_sampleWeight, only: setReweight
    use pm_sampleWeight, only: getReweight
    use pm_mathLogSumExp, only: getLogSumExp
    use pm_sampleQuan, only: getQuan, neimean!, piwilin ! piwilin is problematic, because some states can be repeated, although highly unlikely.
    use pm_io, only: getFileUnit
    use pm_err, only: getLine

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getErrSampling_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDRAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrParaDRAM_RK5
#define pm_sampling_kernel pm_sampling_kernel_dram_RK5
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrParaDRAM_RK4
#define pm_sampling_kernel pm_sampling_kernel_dram_RK4
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrParaDRAM_RK3
#define pm_sampling_kernel pm_sampling_kernel_dram_RK3
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrParaDRAM_RK2
#define pm_sampling_kernel pm_sampling_kernel_dram_RK2
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrParaDRAM_RK1
#define pm_sampling_kernel pm_sampling_kernel_dram_RK1
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaDRAM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDISE_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrParaDISE_RK5
#define pm_sampling_kernel pm_sampling_kernel_dram_RK5
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrParaDISE_RK4
#define pm_sampling_kernel pm_sampling_kernel_dram_RK4
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrParaDISE_RK3
#define pm_sampling_kernel pm_sampling_kernel_dram_RK3
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrParaDISE_RK2
#define pm_sampling_kernel pm_sampling_kernel_dram_RK2
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrParaDISE_RK1
#define pm_sampling_kernel pm_sampling_kernel_dram_RK1
#include "pm_sampling@routines.inc.F90"
#undef pm_sampling_kernel
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ParaDISE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getErrSampling_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine killMeAlreadyCMake1_RK5(); use pm_sampling_kernel_dise_RK5, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK4(); use pm_sampling_kernel_dise_RK4, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK3(); use pm_sampling_kernel_dise_RK3, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK2(); use pm_sampling_kernel_dise_RK2, only: RKG; end subroutine
    subroutine killMeAlreadyCMake1_RK1(); use pm_sampling_kernel_dise_RK1, only: RKG; end subroutine

    subroutine killMeAlreadyCMake2_RK5(); use pm_sampling_kernel_dram_RK5, only: RKG; end subroutine
    subroutine killMeAlreadyCMake2_RK4(); use pm_sampling_kernel_dram_RK4, only: RKG; end subroutine
    subroutine killMeAlreadyCMake2_RK3(); use pm_sampling_kernel_dram_RK3, only: RKG; end subroutine
    subroutine killMeAlreadyCMake2_RK2(); use pm_sampling_kernel_dram_RK2, only: RKG; end subroutine
    subroutine killMeAlreadyCMake2_RK1(); use pm_sampling_kernel_dram_RK1, only: RKG; end subroutine

    !subroutine killMeAlreadyCMake3_RK5(); use pm_sampling_kernel_nest_RK5, only: RKG; end subroutine
    !subroutine killMeAlreadyCMake3_RK4(); use pm_sampling_kernel_nest_RK4, only: RKG; end subroutine
    !subroutine killMeAlreadyCMake3_RK3(); use pm_sampling_kernel_nest_RK3, only: RKG; end subroutine
    !subroutine killMeAlreadyCMake3_RK2(); use pm_sampling_kernel_nest_RK2, only: RKG; end subroutine
    !subroutine killMeAlreadyCMake3_RK1(); use pm_sampling_kernel_nest_RK1, only: RKG; end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines