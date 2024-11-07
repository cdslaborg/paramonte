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
!>  This module contains the **internal** classes and procedures for setting up the attributes of the ParaMonte library DRAM-MCMC samplers.<br>
!>
!>  \details
!>  For more information, see the description of the attributes within the bodies of their constructors in this module.<br>
!>  Alternatively, a description of these simulation specifications is always printed out the in `_report.txt` files of each ParaMonte DRAM-MCMC simulation.
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

module pm_sampling_dram_RK5
#if RK5_ENABLED
    use pm_kind, only: RKG => RK5
#define pm_sampling_scio pm_sampling_scio_RK5
#define pm_sampling_mcmc pm_sampling_mcmc_RK5
#include "pm_sampling_dram.imp.F90"
#undef pm_sampling_mcmc
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_dram_RK4
#if RK4_ENABLED
    use pm_kind, only: RKG => RK4
#define pm_sampling_scio pm_sampling_scio_RK4
#define pm_sampling_mcmc pm_sampling_mcmc_RK4
#include "pm_sampling_dram.imp.F90"
#undef pm_sampling_mcmc
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_dram_RK3
#if RK3_ENABLED
    use pm_kind, only: RKG => RK3
#define pm_sampling_scio pm_sampling_scio_RK3
#define pm_sampling_mcmc pm_sampling_mcmc_RK3
#include "pm_sampling_dram.imp.F90"
#undef pm_sampling_mcmc
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_dram_RK2
#if RK2_ENABLED
    use pm_kind, only: RKG => RK2
#define pm_sampling_scio pm_sampling_scio_RK2
#define pm_sampling_mcmc pm_sampling_mcmc_RK2
#include "pm_sampling_dram.imp.F90"
#undef pm_sampling_mcmc
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module

module pm_sampling_dram_RK1
#if RK1_ENABLED
    use pm_kind, only: RKG => RK1
#define pm_sampling_scio pm_sampling_scio_RK1
#define pm_sampling_mcmc pm_sampling_mcmc_RK1
#include "pm_sampling_dram.imp.F90"
#undef pm_sampling_mcmc
#undef pm_sampling_scio
#else
    use pm_kind, only: RKG => RK
#endif
end module
