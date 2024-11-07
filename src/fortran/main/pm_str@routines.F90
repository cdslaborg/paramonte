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
!>  This file contains procedure implementations of [pm_str](@ref pm_str).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_str) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define alleq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure alleq_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure alleq_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure alleq_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure alleq_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure alleq_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef alleq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isEndedWith_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isEndedWith_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isEndedWith_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isEndedWith_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isEndedWith_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isEndedWith_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure isEndedWith_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isEndedWith_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isEndedWith_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isEndedWith_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isEndedWith_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1
    module procedure isEndedWith_BSSK
        use pm_kind, only: SKG => SK
#include "pm_str@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isEndedWith_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinLoc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NoMask_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinLocNoMask_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinLocNoMask_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinLocNoMask_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinLocNoMask_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinLocNoMask_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NoMask_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Masked_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinLocMasked_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinLocMasked_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinLocMasked_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinLocMasked_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinLocMasked_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Masked_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMinLoc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMaxLoc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NoMask_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMaxLocNoMask_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMaxLocNoMask_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMaxLocNoMask_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMaxLocNoMask_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMaxLocNoMask_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NoMask_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Masked_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMaxLocMasked_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMaxLocMasked_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMaxLocMasked_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMaxLocMasked_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMaxLocMasked_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Masked_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMaxLoc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinVal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NoMask_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinValNoMask_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinValNoMask_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinValNoMask_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinValNoMask_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinValNoMask_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NoMask_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Masked_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinValMasked_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinValMasked_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinValMasked_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinValMasked_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinValMasked_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Masked_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMinVal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMaxVal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NoMask_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMaxValNoMask_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMaxValNoMask_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMaxValNoMask_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMaxValNoMask_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMaxValNoMask_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NoMask_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Masked_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMaxValMasked_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMaxValMasked_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMaxValMasked_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMaxValMasked_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMaxValMasked_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Masked_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMaxVal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCharSeq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getCharSeq_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getCharSeq_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getCharSeq_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getCharSeq_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getCharSeq_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCharSeq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCharVec_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getCharVec_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getCharVec_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getCharVec_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getCharVec_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getCharVec_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCharVec_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getTrimmedTZ_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getTrimmedTZ_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getTrimmedTZ_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getTrimmedTZ_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getTrimmedTZ_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getTrimmedTZ_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getTrimmedTZ_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStrWrapped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStrWrapped_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStrWrapped_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStrWrapped_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStrWrapped_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStrWrapped_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getStrWrapped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLenIndent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getLenIndent_SK5
        use pm_kind, only: SKG => SK5
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getLenIndent_SK4
        use pm_kind, only: SKG => SK4
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getLenIndent_SK3
        use pm_kind, only: SKG => SK3
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getLenIndent_SK2
        use pm_kind, only: SKG => SK2
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getLenIndent_SK1
        use pm_kind, only: SKG => SK1
#include "pm_str@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLenIndent_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#define setStrWrapped_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!#if 0
!#define List_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getStrWrappedList_SK5
!        use pm_kind, only: SKG => SK5
!#include "pm_str@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getStrWrappedList_SK4
!        use pm_kind, only: SKG => SK4
!#include "pm_str@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getStrWrappedList_SK3
!        use pm_kind, only: SKG => SK3
!#include "pm_str@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getStrWrappedList_SK2
!        use pm_kind, only: SKG => SK2
!#include "pm_str@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getStrWrappedList_SK1
!        use pm_kind, only: SKG => SK1
!#include "pm_str@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef List_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines