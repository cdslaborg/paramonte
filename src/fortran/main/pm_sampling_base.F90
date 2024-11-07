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
!>  This module contains the **internal** classes and procedures for setting up the attributes of the ParaMonte library samplers and optimizers.<br>
!>
!>  \details
!>  For more information, see the description of the attributes within the bodies of their constructors in this module.<br>
!>  Alternatively, a description of these simulation specifications is always printed out the in `_report.txt` files of each ParaMonte exploration simulation.
!>
!>  \note
!>  The contents of this module are not meant to be used by the end users of the ParaMonte library.<br>
!>
!>  \devnote
!>  CMake is totally incapable of understanding the dependency between the modules declared in this file and the included files.<br>
!>  As such, any change to the interfaces in the included implementation file requires
!>  manual change in this file to trigger interface module file generation by CMake.<br>
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin<br>

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampling_base_RK5
#if RK5_ENABLED
    use pm_kind, only: RKG => RK5
#define pm_sampling_scio pm_sampling_scio_RK5
#include "pm_sampling_base.imp.F90"
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_base_RK4
#if RK4_ENABLED
    use pm_kind, only: RKG => RK4
#define pm_sampling_scio pm_sampling_scio_RK4
#include "pm_sampling_base.imp.F90"
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_base_RK3
#if RK3_ENABLED
    use pm_kind, only: RKG => RK3
#define pm_sampling_scio pm_sampling_scio_RK3
#include "pm_sampling_base.imp.F90"
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_base_RK2
#if RK2_ENABLED
    use pm_kind, only: RKG => RK2
#define pm_sampling_scio pm_sampling_scio_RK2
#include "pm_sampling_base.imp.F90"
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_base_RK1
#if RK1_ENABLED
    use pm_kind, only: RKG => RK1
#define pm_sampling_scio pm_sampling_scio_RK1
#include "pm_sampling_base.imp.F90"
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

!#if RK5_ENABLED
!module pm_sampling_base_RK5
!    use pm_kind, only: RKG => RK5
!#define pm_sampling_scio pm_sampling_scio_RK5
!#include "pm_sampling_base.imp.F90"
!#undef pm_sampling_scio
!end module
!#else
!module pm_sampling_base_RK5
!    use pm_kind, only: RKG => RK
!end module
!#endif
!
!#if RK4_ENABLED
!module pm_sampling_base_RK4
!    use pm_kind, only: RKG => RK4
!#define pm_sampling_scio pm_sampling_scio_RK4
!#include "pm_sampling_base.imp.F90"
!#undef pm_sampling_scio
!end module
!#else
!module pm_sampling_base_RK4
!    use pm_kind, only: RKG => RK
!end module
!#endif
!
!#if RK3_ENABLED
!module pm_sampling_base_RK3
!    use pm_kind, only: RKG => RK3
!#define pm_sampling_scio pm_sampling_scio_RK3
!#include "pm_sampling_base.imp.F90"
!#undef pm_sampling_scio
!end module
!#else
!module pm_sampling_base_RK3
!    use pm_kind, only: RKG => RK
!end module
!#endif
!
!#if RK2_ENABLED
!module pm_sampling_base_RK2
!    use pm_kind, only: RKG => RK2
!#define pm_sampling_scio pm_sampling_scio_RK2
!#include "pm_sampling_base.imp.F90"
!#undef pm_sampling_scio
!end module
!#else
!module pm_sampling_base_RK2
!    use pm_kind, only: RKG => RK
!end module
!#endif
!
!#if RK1_ENABLED
!module pm_sampling_base_RK1
!    use pm_kind, only: RKG => RK1
!#define pm_sampling_scio pm_sampling_scio_RK1
!#include "pm_sampling_base.imp.F90"
!#undef pm_sampling_scio
!end module
!#else
!module pm_sampling_base_RK1
!    use pm_kind, only: RKG => RK
!end module
!#endif
