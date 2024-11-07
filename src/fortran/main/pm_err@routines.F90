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
!>  This file contains procedure implementations of [pm_err](@ref pm_err).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if MEXPRINT_ENABLED
#include "fintrf.h"
#endif

submodule (pm_err) routines ! LCOV_EXCL_LINE

    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure err_typer
        if (present(occurred)) err%occurred = occurred
        if (present(stat)) err%stat = stat
        if (present(msg)) then
            err%msg = msg
        else
            err%msg = repeat(SK_"A", 255)
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define mark_typer_ENABLED 1
    module procedure mark_typer
#include "pm_err@routines.inc.F90"
    end procedure
#undef mark_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define note_typer_ENABLED 1
    module procedure note_typer
#include "pm_err@routines.inc.F90"
    end procedure
#undef note_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define warn_typer_ENABLED 1
    module procedure warn_typer
#include "pm_err@routines.inc.F90"
    end procedure
#undef warn_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define stop_typer_ENABLED 1
    module procedure stop_typer
#include "pm_err@routines.inc.F90"
    end procedure
#undef stop_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFine_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getFine_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getFine_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getFine_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getFine_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getFine_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFine_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFile_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getFile_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getFile_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getFile_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getFile_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getFile_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFile_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLine_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getLine_IK5
        use pm_kind, only: IKG => IK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getLine_IK4
        use pm_kind, only: IKG => IK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getLine_IK3
        use pm_kind, only: IKG => IK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getLine_IK2
        use pm_kind, only: IKG => IK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getLine_IK1
        use pm_kind, only: IKG => IK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLine_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setAsserted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setAsserted_LK5
        use pm_kind, only: LKG => LK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setAsserted_LK4
        use pm_kind, only: LKG => LK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setAsserted_LK3
        use pm_kind, only: LKG => LK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setAsserted_LK2
        use pm_kind, only: LKG => LK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setAsserted_LK1
        use pm_kind, only: LKG => LK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setAsserted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMarked_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Static_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMarkedStatic_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMarkedStatic_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMarkedStatic_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMarkedStatic_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMarkedStatic_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Static_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Method_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMarkedMethod_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMarkedMethod_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMarkedMethod_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMarkedMethod_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMarkedMethod_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Method_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMarked_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNoted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Static_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setNotedStatic_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setNotedStatic_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setNotedStatic_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setNotedStatic_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setNotedStatic_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Static_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Method_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setNotedMethod_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setNotedMethod_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setNotedMethod_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setNotedMethod_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setNotedMethod_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Method_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNoted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setWarned_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Static_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setWarnedStatic_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setWarnedStatic_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setWarnedStatic_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setWarnedStatic_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setWarnedStatic_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Static_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Method_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setWarnedMethod_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setWarnedMethod_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setWarnedMethod_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setWarnedMethod_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setWarnedMethod_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Method_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setWarned_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setAborted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Static_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setAbortedStatic_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setAbortedStatic_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setAbortedStatic_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setAbortedStatic_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setAbortedStatic_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Static_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Method_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setAbortedMethod_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setAbortedMethod_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setAbortedMethod_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setAbortedMethod_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setAbortedMethod_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_err@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Method_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setAborted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
