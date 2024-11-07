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
!>  This module contains classes and procedures for computing the Bhattacharyya statistical distance between two probability distributions.
!>
!>  \details
!>  The Bhattacharyya distance is a quantity which represents a notion of similarity between two probability distributions.<br>
!>  It is closely related to the **Bhattacharyya coefficient**, which is a measure of the amount of overlap between two statistical samples or populations.<br>
!>  The Bhattacharyya distance is **not a metric**, despite being named a *distance*, since **it does not obey the triangle inequality**.<br>
!>
!>  History
!>  -------
!>
!>  Both the Bhattacharyya distance and the Bhattacharyya coefficient are named after Anil Kumar Bhattacharyya,
!>  a statistician who worked in the 1930s at the Indian Statistical Institute.<br>
!>  He developed the method to measure the distance between two non-normal distributions and illustrated this with
!>  the classical multinomial populations, this work despite being submitted for publication in 1941, appeared almost five years later in Sankhya.<br>
!>  Consequently, Professor Bhattacharyya started working toward developing a distance metric for probability distributions that are absolutely
!>  continuous with respect to the Lebesgue measure and published his progress in 1942, at Proceedings of the Indian Science Congress
!>  and the final work has appeared in 1943 in the Bulletin of the Calcutta Mathematical Society.<br>
!>
!>  Definition using Probability Theory
!>  -----------------------------------
!>
!>  For probability distributions \f$P\f$ and \f$Q\f$ on the same domain \f$\mathcal{X}\f$, the Bhattacharyya distance is defined as,
!>  \f{equation}{
!>      D_{B}(P, Q) = -\ln\left(\up{BC}(P,Q)\right) ~,
!>  \f}
!>
!>  where
!>  \f{equation}{
!>      \up{BC}(P, Q) = \sum_{x \in {\mathcal{X}}}{\sqrt {P(x)Q(x)}} ~,
!>  \f}
!>
!>  is the Bhattacharyya coefficient for discrete probability distributions.<br>
!>  For continuous probability distributions, \f$P(dx) = p(x)dx\f$ and \f$Q(dx) = q(x)dx\f$ where \f$p(x)\f$ and \f$q(x)\f$ are the probability density functions, the Bhattacharyya coefficient is defined as,
!>  \f{equation}{
!>      \up{BC}(P, Q) = \int_{\mathcal{X}}{\sqrt{p(x)q(x)}}\,dx ~.
!>  \f}
!>
!>  Definition using Measure Theory
!>  -------------------------------
!>
!>  More generally, given two probability measures \f$P, Q\f$ on a measurable space \f$(\mathcal{X}, \mathcal{B})\f$, let \f$\lambda\f$ be a sigma finite measure such that \f$P\f$ and \f$Q\f$
!>  are absolutely continuous with respect to \f$\lambda\f$, that is, such that \f$P(dx) = p(x) \lambda(dx)\f$, and \f$Q(dx) = q(x)\lambda(dx)\f$ for probability density functions \f$p, q\f$
!>  with respect to \f$\lambda\f$ defined \f$\lambda\f$-almost everywhere.<br>
!>  Such a measure, even such a probability measure, always exists, for example, \f$\lambda = \frac{1}{2}(P + Q)\f$.<br>
!>  Then the Bhattacharyya measure on \f$({\mathcal {X}},{\mathcal {B}})\f$ is defined by,
!>  \f{equation}{
!>      \up{bc}(dx|P,Q)={\sqrt {p(x)q(x)}}\,\lambda (dx)={\sqrt {{\frac {P(dx)}{\lambda (dx)}}(x){\frac {Q(dx)}{\lambda (dx)}}(x)}}\lambda (dx) ~.
!>  \f}
!>
!>  The definition does not depend on the measure \f$\lambda\f$, for if we choose a measure \f$\mu\f$ such that \f$\lambda\f$ and another measure choice \f$\lambda'\f$
!>  are absolutely continuous, i.e., \f$\lambda = l(x)\mu\f$ and \f$\lambda '=l'(x)\mu\f$, then,
!>  \f{equation}{
!>      P(dx) = p(x)\lambda (dx) = p'(x)\lambda '(dx)=p(x)l(x)\mu (dx)=p'(x)l'(x)\mu (dx) ~,
!>  \f}
!>
!>  and similarly for \f$Q\f$.<br>
!>  We then have,
!>  \f{equation}{
!>      \up{bc}(dx | P, Q) = {\sqrt{p(x)q(x)}} \, \lambda (dx) = {\sqrt {p(x)q(x)}} \, l(x) \mu(x) =
!>      {\sqrt{p(x)l(x)q(x) \, l(x)}} \mu(dx) = {\sqrt{p'(x)l'(x)q'(x)l'(x)}} \, \mu(dx) =
!>      {\sqrt{p'(x)q'(x)}} \, \lambda'(dx) ~.
!>  \f}
!>
!>  Then define the Bhattacharyya coefficient as,
!>  \f{equation}{
!>      \up{BC}(P,Q) = \int_{\mathcal{X}} \up{bc}(dx|P, Q) = \int_{\mathcal{X}}{\sqrt {p(x)q(x)}}\,\lambda (dx) ~.
!>  \f}
!>
!>  By the above, the quantity \f$\up{BC}(P,Q)\f$ does not depend on \f$\lambda\f$, and by the Cauchy inequality \f$0\leq \up{BC}(P,Q)\leq 1\f$.<br>
!>  In particular if \f$P(dx) = p(x)Q(dx)\f$ is absolutely continuous w.r.t. to \f$Q\f$ with Radon Nikodym derivative \f$p(x) = {\frac {P(dx)}{Q(dx)}}(x)\f$, then,
!>  \f{equation}{
!>      \up{BC}(P,Q) = \int_{\mathcal{X}}{\sqrt {p(x)}}Q(dx) = \int_{\mathcal{X}} {\sqrt {\frac{P(dx)}{Q(dx)}}} Q(dx) = E_{Q} \left[{\sqrt{\frac{P(dx)}{Q(dx)}}}\right] ~.
!>  \f}
!>
!>  Properties of the Bhattacharyya Distance
!>  ----------------------------------------
!>
!>  <ol>
!>      <li>    The conditions \f$0\leq \up{BC}\leq 1\f$ and \f$0\leq D_{B}\leq \infty\f$ hold for the Bhattacharyya coefficient and distance, respectively.
!>      <li>    The Bhattacharyya distance \f$D_{B}\f$ does not obey the triangle inequality, though the [Hellinger distance](@ref pm_distanceHellinger) \f$\sqrt{1 - \up{BC}(p,q)}\f$ does.<br>
!>  </ol>
!>
!>  Connection with Total Variation Distance (TVD)
!>  ----------------------------------------------
!>
!>  The Bhattacharyya distance \f$H(P, Q)\f$ and the total variation distance (or statistical distance) \f$\delta(P,Q)\f$ are related as follows,
!>  \f{equation}{
!>      H^{2}(P, Q)\leq \delta(P, Q)\leq {\sqrt{2}}H(P, Q) ~.
!>  \f}
!>
!>  These inequalities follow immediately from the inequalities between the 1-norm and the 2-norm.
!>
!>  \see
!>  [pm_distanceBhat](@ref pm_distanceBhat)<br>
!>  [pm_distanceEuclid](@ref pm_distanceEuclid)<br>
!>  [pm_distanceHellinger](@ref pm_distanceHellinger)<br>
!>  [pm_distanceKolm](@ref pm_distanceKolm)<br>
!>  [pm_distanceMahal](@ref pm_distanceMahal)<br>
!>
!>  \test
!>  [test_pm_distanceBhat](@ref test_pm_distanceBhat)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distanceBhat

    use pm_kind, only: SK, IK, LK
    use pm_mathConst, only: ninf_type
    use pm_mathConst, only: pinf_type
    use pm_except, only: getInfNeg
    use pm_except, only: getInfPos

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distanceBhat"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Bhattacharyya distance of two univariate (discrete or continuous) distributions.
    !>
    !>  \details
    !>  See [pm_distanceBhat](@ref pm_distanceBhat) for the mathematical definition of the Bhattacharyya distance.
    !>
    !>  \param[in]  p       :   The vector of the same size as the non-zero size of the input argument `q`, of the same type and kind as the output `bhat`,
    !>                          representing the Probability Mass Function (PMF) of the first distribution in the computation of the squared Bhattacharyya distance.<br>
    !>                          Alternatively, the input `p` can be a non-elemental function that takes a scalar of the same type and kind as the output `bhat`
    !>                          and returns the output `pdf` of the same type and kind as `x`, representing the Probability Density Function (PDF)
    !>                          of the first distribution in the computation of the squared Bhattacharyya distance.<br>
    !>                          The following illustrates the generic interface of `p`,
    !>                          \code{.F90}
    !>                              function p(x) result(pdf)
    !>                                  real(RKG), intent(in) :: x
    !>                                  real(RKG) :: pdf
    !>                              end function
    !>                          \endcode
    !>                          where `RKG` is the kind type parameter of the output `bhat`.<br>
    !>                          The input arguments `p` and `q` must be of the same type and kind (and size if they are vectors).<br>
    !>  \param[in]  q       :   The vector of the same size as the non-zero size of the input argument `p`, of the same type and kind as the output `bhat`,
    !>                          representing the Probability Mass Function (PMF) of the second distribution in the computation of the squared Bhattacharyya distance.<br>
    !>                          Alternatively, the input `p` can be a non-elemental function that takes a scalar of the same type and kind as the output `bhat`
    !>                          and returns the output `pdf` of the same type and kind as `x`, representing the Probability Density Function (PDF)
    !>                          of the second distribution in the computation of the squared Bhattacharyya distance.<br>
    !>                          The following illustrates the generic interface of `p`,
    !>                          \code{.F90}
    !>                              function p(x) result(pdf)
    !>                                  real(RKG), intent(in) :: x
    !>                                  real(RKG) :: pdf
    !>                              end function
    !>                          \endcode
    !>                          where `RKG` is the kind type parameter of the output `bhat`.<br>
    !>                          The input arguments `p` and `q` must be of the same type and kind (and size if they are vectors).<br>
    !>  \param[in]  lb      :   The input scalar of type `real` of the same kind as `integral`, representing the lower limit of integration.<br>
    !>                          Set `lb = -huge(lb)` or to the IEEE-compliant negative infinity (`lb = `[getInfNeg(lb)](@ref pm_except::getInfNeg))
    !>                          to imply \f$-\infty\f$ as the lower bound of integration.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**, default = [getInfNeg(lb)](@ref pm_except::getInfNeg). It can be present **only if** the input arguments `p` and `q` are procedures.)
    !>  \param[in]  ub      :   The input scalar of type `real` of the same kind as `integral`, representing the upper limit of integration.<br>
    !>                          Set `ub = huge(ub)` or to the IEEE-compliant positive infinity (`ub = `[getInfPos(lb)](@ref pm_except::getInfPos))
    !>                          to imply \f$+\infty\f$ as the upper bound of integration.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**, default = [getInfPos(lb)](@ref pm_except::getInfPos). It can be present **only if** the input arguments `p` and `q` are procedures.)
    !>  \param[out] failed  :   The output scalar of type `logical` of default kind \LK, that is set to `.true.` <b>if and only if</b> the integration Bhattacharyya **fails to converge within the tolerances**.<br>
    !>                          Otherwise, it is set `.false.` if the integration succeeds with no errors.<br>
    !>                          See also the corresponding output argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          See also the description of the output argument `err` of [getQuadErr](@ref pm_quadPack::getQuadErr) for information on the kinds of integration failures that can happen.<br>
    !>                          (**optional**. If missing and the integration fails, the procedure will halt the program by calling `error stop`.)<br>
    !>  \param[out] msg     :   The output scalar argument of type `character` of default kind \SK of arbitrary length type parameter that is set to a diagnostic message if the integration fails to converge.<br>
    !>                          A length type parameter of `127` is sufficient to capture all error messages.<br>
    !>                          If `msg` has shorter length parameter, the output message will be trimmed from the end, otherwise padded with blanks as necessary.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**. If missing, no diagnostic message will be returned.)<br>
    !>
    !>  \return
    !>  `bhat`            :   The output scalar of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          representing the Bhattacharyya distance of the discrete or continuous distributions represented by `p` and `q`.<br>
    !>
    !>  \interface{getDisBhat}
    !>  \code{.F90}
    !>
    !>      use pm_distanceBhat, only: getDisBhat
    !>
    !>      ! discrete distributions.
    !>
    !>      bhat = getDisBhat(p(:), q(:))
    !>
    !>      ! continuous univariate distributions.
    !>
    !>      bhat = getDisBhat(p, q, lb = lb, ub = ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `sum(p) == 1` must hold for the corresponding input arguments.<br>
    !>  The condition `sum(q) == 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= p) .and. all(p <= 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= q) .and. all(q <= 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(p) == size(q)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input arguments `p` and `q` are procedures.
    !>
    !>  \example{getDisBhat}
    !>  \include{lineno} example/pm_distanceBhat/getDisBhat/main.F90
    !>  \compilef{getDisBhat}
    !>  \output{getDisBhat}
    !>  \include{lineno} example/pm_distanceBhat/getDisBhat/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceBhat](@ref test_pm_distanceBhat)<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The runtime checks for the complex input `invCov` must be implemented.<br>
    !>
    !>  \final{getDisBhat}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! PMF

    interface getDisBhat

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisBhatPMF_RK5(p, q) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPMF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: bhat
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisBhatPMF_RK4(p, q) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPMF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: bhat
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisBhatPMF_RK3(p, q) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPMF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: bhat
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisBhatPMF_RK2(p, q) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPMF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: bhat
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisBhatPMF_RK1(p, q) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPMF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: bhat
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Proc

    interface getDisBhat

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisBhatPRC_RK5(p, q, lb, ub, failed, msg) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: bhat
    end function
#endif

#if RK4_ENABLED
    module function getDisBhatPRC_RK4(p, q, lb, ub, failed, msg) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: bhat
    end function
#endif

#if RK3_ENABLED
    module function getDisBhatPRC_RK3(p, q, lb, ub, failed, msg) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: bhat
    end function
#endif

#if RK2_ENABLED
    module function getDisBhatPRC_RK2(p, q, lb, ub, failed, msg) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: bhat
    end function
#endif

#if RK1_ENABLED
    module function getDisBhatPRC_RK1(p, q, lb, ub, failed, msg) result(bhat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: bhat
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisBhatPRC_FF_RK5(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FF_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisBhatPRC_FF_RK4(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FF_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisBhatPRC_FF_RK3(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FF_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisBhatPRC_FF_RK2(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FF_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisBhatPRC_FF_RK1(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FF_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisBhatPRC_FI_RK5(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FI_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisBhatPRC_FI_RK4(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FI_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisBhatPRC_FI_RK3(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FI_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisBhatPRC_FI_RK2(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FI_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisBhatPRC_FI_RK1(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_FI_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisBhatPRC_IF_RK5(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_IF_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisBhatPRC_IF_RK4(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_IF_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisBhatPRC_IF_RK3(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_IF_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisBhatPRC_IF_RK2(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_IF_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisBhatPRC_IF_RK1(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_IF_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisBhatPRC_II_RK5(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_II_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisBhatPRC_II_RK4(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_II_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisBhatPRC_II_RK3(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_II_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisBhatPRC_II_RK2(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_II_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisBhatPRC_II_RK1(p, q, lb, ub) result(bhat)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisBhatPRC_II_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: bhat
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distanceBhat ! LCOV_EXCL_LINE