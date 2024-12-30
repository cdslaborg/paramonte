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
!>  This module contains the **internal** classes and procedures for setting up the attributes of the ParaMonte library MCMC samplers.<br>
!>
!>  \details
!>  For more information, see the description of the attributes within the bodies of their constructors in this module.<br>
!>  Alternatively, a description of these simulation specifications is always printed out the in `_report.txt` files of each ParaMonte MCMC simulation.
!>
!>  \note
!>  The contents of this module are not meant to be used by the end users of the ParaMonte library.<br>
!>
!>  \todo
!>  \phigh
!>  The simulations specifications for a globular proposal start domain must be added.<br>
!>  Currently, only a cubical domain is supported, which can be difficult to deal with when the density function domain is globular.<br>
!>  This can be fixed by adding additional specifications:<br>
!>  <ol>
!>      <li>    `proposalStartDomainBallAvg`
!>      <li>    `proposalStartDomainBallCor`
!>      <li>    `proposalStartDomainBallCov`
!>      <li>    `proposalStartDomainBallStd`
!>  </ol>
!>
!>  \devnote
!>  The madness seen here with module-level generics is due to the lack of support for PDTs in \gfortran{13.1} and older versions.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin<br>

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampling_mcmc_RK5
    use pm_kind, only: RKG => RK5
#if RK5_ENABLED
#define pm_sampling_scio pm_sampling_scio_RK5
#define pm_sampling_base pm_sampling_base_RK5
#include "pm_sampling_mcmc.imp.F90"
#undef pm_sampling_base
#undef pm_sampling_scio
#endif
end module

module pm_sampling_mcmc_RK4
    use pm_kind, only: RKG => RK4
#if RK4_ENABLED
#define pm_sampling_scio pm_sampling_scio_RK4
#define pm_sampling_base pm_sampling_base_RK4
#include "pm_sampling_mcmc.imp.F90"
#undef pm_sampling_base
#undef pm_sampling_scio
#endif
end module

module pm_sampling_mcmc_RK3
    use pm_kind, only: RKG => RK3
#if RK3_ENABLED
#define pm_sampling_scio pm_sampling_scio_RK3
#define pm_sampling_base pm_sampling_base_RK3
#include "pm_sampling_mcmc.imp.F90"
#undef pm_sampling_base
#undef pm_sampling_scio
#endif
end module

module pm_sampling_mcmc_RK2
    use pm_kind, only: RKG => RK2
#if RK2_ENABLED
#define pm_sampling_scio pm_sampling_scio_RK2
#define pm_sampling_base pm_sampling_base_RK2
#include "pm_sampling_mcmc.imp.F90"
#undef pm_sampling_base
#undef pm_sampling_scio
#endif
end module

module pm_sampling_mcmc_RK1
    use pm_kind, only: RKG => RK1
#if RK1_ENABLED
#define pm_sampling_scio pm_sampling_scio_RK1
#define pm_sampling_base pm_sampling_base_RK1
#include "pm_sampling_mcmc.imp.F90"
#undef pm_sampling_base
#undef pm_sampling_scio
#endif
end module
