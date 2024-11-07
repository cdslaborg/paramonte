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
!>  This module contains the **internal** classes and procedures for setting up the attributes of the ParaMonte library sampler kernels.<br>
!>
!>  \note
!>  The contents of this module are not meant to be used by the end users of the ParaMonte library.<br>
!>
!>  \devnote
!>  The madness seen here with module-level generics is due to the lack of support for PDTs in \gfortran{13.1} and older versions.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin<br>

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDISE_ENABLED 1

module pm_sampling_kernel_dise_RK5
#if RK5_ENABLED
    use pm_kind, only: RKG => RK5
#define pm_sampling_proposal pm_sampling_proposal_dise_RK5
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dise_RK4
#if RK4_ENABLED
    use pm_kind, only: RKG => RK4
#define pm_sampling_proposal pm_sampling_proposal_dise_RK4
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dise_RK3
#if RK3_ENABLED
    use pm_kind, only: RKG => RK3
#define pm_sampling_proposal pm_sampling_proposal_dise_RK3
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dise_RK2
#if RK2_ENABLED
    use pm_kind, only: RKG => RK2
#define pm_sampling_proposal pm_sampling_proposal_dise_RK2
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dise_RK1
#if RK1_ENABLED
    use pm_kind, only: RKG => RK1
#define pm_sampling_proposal pm_sampling_proposal_dise_RK1
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

#undef ParaDISE_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaDRAM_ENABLED 1

module pm_sampling_kernel_dram_RK5
#if RK5_ENABLED
    use pm_kind, only: RKG => RK5
#define pm_sampling_proposal pm_sampling_proposal_dram_RK5
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dram_RK4
#if RK4_ENABLED
    use pm_kind, only: RKG => RK4
#define pm_sampling_proposal pm_sampling_proposal_dram_RK4
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dram_RK3
#if RK3_ENABLED
    use pm_kind, only: RKG => RK3
#define pm_sampling_proposal pm_sampling_proposal_dram_RK3
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dram_RK2
#if RK2_ENABLED
    use pm_kind, only: RKG => RK2
#define pm_sampling_proposal pm_sampling_proposal_dram_RK2
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_kernel_dram_RK1
#if RK1_ENABLED
    use pm_kind, only: RKG => RK1
#define pm_sampling_proposal pm_sampling_proposal_dram_RK1
#include "pm_sampling_kernel.imp.F90"
#undef pm_sampling_proposal
#else
    use pm_kind, only: RKG => RK
#endif
end module

#undef ParaDRAM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ParaNest_ENABLED 1
!
!module pm_sampling_kernel_nest_RK5
!#if RK5_ENABLED
!    use pm_kind, only: RKG => RK5
!#define pm_sampling_proposal pm_sampling_proposal_dram_RK5
!#include "pm_sampling_kernel.imp.F90"
!#undef pm_sampling_proposal
!!#else
!    use pm_kind, only: RKG => RK
!#endif
!end module
!
!module pm_sampling_kernel_nest_RK4
!#if RK4_ENABLED
!    use pm_kind, only: RKG => RK4
!#define pm_sampling_proposal pm_sampling_proposal_dram_RK4
!#include "pm_sampling_kernel.imp.F90"
!#undef pm_sampling_proposal
!!#else
!    use pm_kind, only: RKG => RK
!#endif
!end module
!
!module pm_sampling_kernel_nest_RK3
!#if RK3_ENABLED
!    use pm_kind, only: RKG => RK3
!#define pm_sampling_proposal pm_sampling_proposal_dram_RK3
!#include "pm_sampling_kernel.imp.F90"
!#undef pm_sampling_proposal
!!#else
!    use pm_kind, only: RKG => RK
!#endif
!end module
!
!module pm_sampling_kernel_nest_RK2
!#if RK2_ENABLED
!    use pm_kind, only: RKG => RK2
!#define pm_sampling_proposal pm_sampling_proposal_dram_RK2
!#include "pm_sampling_kernel.imp.F90"
!#undef pm_sampling_proposal
!!#else
!    use pm_kind, only: RKG => RK
!#endif
!end module
!
!module pm_sampling_kernel_nest_RK1
!#if RK1_ENABLED
!    use pm_kind, only: RKG => RK1
!#define pm_sampling_proposal pm_sampling_proposal_dram_RK1
!#include "pm_sampling_kernel.imp.F90"
!#undef pm_sampling_proposal
!!#else
!    use pm_kind, only: RKG => RK
!#endif
!end module

#undef ParaNest_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
