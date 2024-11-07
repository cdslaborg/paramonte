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
!>  This module contains the types, classes, and procedures relevant to weights of random samples.
!>
!>  \details
!>  Within the world of statistics, particularly, traditional frequentist statistics,
!>  there is absolute confusion on what types of weights exist and what their implications are.<br>
!>
!>  The most relevant definitions, which are also currently used in the ParaMonte library, are the following,
!>  <ol>
!>      <li>    [Frequency weights](@ref pm_sampleWeight::fweight_type), also known as **frequency**, **count**, or **repeat** weights,
!>              represent the number of duplications of each observation in a sample.<br>
!>              Therefore, the **frequency weights** are expected to be **integers** or **whole numbers**.<br>
!>      <li>    [Reliability weights](@ref pm_sampleWeight::fweight_type) represent a measure of the reliability of each observation in a sample.<br>
!>              Reliability weights are frequently normalized, but not necessarily.<br>
!>  </ol>
!>
!>  There are also other definitions of weight which special use cases, for example:<br>
!>  <ol>
!>              <li>    [Analytic or precision weights](@ref pm_sampleWeight::aweight_type), also known as **analytic weight**, **inverse variance weight**, or **regression weight**,
!>                      are weights that are inversely proportional to the variance of an observation.<br>
!>                      Typically, the observations represent averages and the weights are the number of elements that gave rise to the average.<br>
!>                      These weights are calculated by taking the *inverse of the sampling fraction*. See also<br>
!>                      <ol>
!>                          <li>    [Inverse-variance weighting](https://en.wikipedia.org/wiki/Inverse-variance_weighting).
!>                          <li>    Gatz, Donald F., and Luther Smith (1995), *The standard error of a weighted mean concentration—I. Bootstrapping vs other methods*.
!>                      </ol>
!>  </ol>
!>
!>  \note
!>  A weighted sample has an **effective sample size (ESS)** that is smaller than the actual number of observations in the sample,
!>  \f{equation}{
!>      \mathrm{ESS} = \frac{ \big( \sum_{i = 1}^{N} w_i \big)^2 }{ \sum_{i = 1}^{N} w_i^2 } ~,
!>  \f}
!>  If the weights \f$w_i\f$ are normalized such that \f$\sum_{i = 1}^{N} w_i = N\f$, then,
!>  \f{equation}{
!>      \mathrm{ESS} = \frac{ N }{ 1 + \mathrm{VAR}(w) } ~,
!>  \f}
!>  where \f$\mathrm{VAR}(w)\f$ represents the variance of the weights.<br>
!>  Therefore, the effective sample size is impacted by both the sample size \f$N\f$ and the amount of variability in weights.<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  [Stata documentation](https://www.stata.com/help.cgi?weight)<br>
!>  [Wikipedia](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Frequency_weights)<br>
!>  [Types of weights](https://stats.oarc.ucla.edu/other/mult-pkg/faq/what-types-of-weights-do-sas-stata-and-spss-support/)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleWeight

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleWeight"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different sample weights (e.g., fweight, aweight, rweight, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{weight_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: weight_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the **aweight** type of sample weights.<br>
    !>
    !>  \details
    !>  Precision (or Analytic) weights are weights that are inversely proportional to the variance of an observation.<br>
    !>  The variance of the \f$i\f$th observation is assumed to be \f$\sigma^2/w_i\f$.<br>
    !>  Typically, the observations represent averages and the weights are the number of elements that gave rise to the average.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [aweight](@ref pm_sampleWeight::aweight) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{aweight_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(weight_type) :: aweight_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [aweight_type](@ref pm_sampleWeight::aweight_type) that is exclusively used
    !>  to signify the Symmetric class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Precision (or Analytic) weights are weights that are inversely proportional to the variance of an observation.<br>
    !>  The variance of the \f$i\f$th observation is assumed to be \f$\sigma^2/w_i\f$.<br>
    !>  Typically, the observations represent averages and the weights are the number of elements that gave rise to the average.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{aweight}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(aweight_type), parameter :: aweight = aweight_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: aweight
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the **fweight** type of sample weights.<br>
    !>
    !>  \details
    !>  Frequency weights, represent the number of duplications of each observation in the sample.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [fweight](@ref pm_sampleWeight::fweight) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{fweight_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(weight_type) :: fweight_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [fweight_type](@ref pm_sampleWeight::fweight_type) that is exclusively used
    !>  to signify the Symmetric class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Frequency weights, represent the number of duplications of each observation in the sample.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{fweight}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(fweight_type), parameter :: fweight = fweight_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: fweight
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the **rweight** type of sample weights.<br>
    !>
    !>  \details
    !>  Frequency weights, represent the number of duplications of each observation in the sample.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [rweight](@ref pm_sampleWeight::rweight) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{rweight_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(weight_type) :: rweight_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [rweight_type](@ref pm_sampleWeight::rweight_type) that is exclusively used
    !>  to signify the Symmetric class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [aweight](@ref pm_sampleWeight::aweight)<br>
    !>  [fweight](@ref pm_sampleWeight::fweight)<br>
    !>  [rweight](@ref pm_sampleWeight::rweight)<br>
    !>  [weight_type](@ref pm_sampleWeight::weight_type)<br>
    !>  [aweight_type](@ref pm_sampleWeight::aweight_type)<br>
    !>  [fweight_type](@ref pm_sampleWeight::fweight_type)<br>
    !>  [rweight_type](@ref pm_sampleWeight::rweight_type)<br>
    !>
    !>  \final{rweight}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(rweight_type), parameter :: rweight = rweight_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: rweight
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a reweighting of the input `weight` vector, such that the sequence represented by output `reweight`
    !>  is the refined (skipped by the amount `skip`) version of the sequence represented by the input `weight`.<br>
    !>
    !>  \brief
    !>  The following examples illustrate the specific functioning of this generic interface:<br>
    !>  \code{.sh}
    !>               skip: 1
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 5, 0, 1, 3, 1
    !>               skip: 2
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 3, 0, 0, 2, 0
    !>               skip: 3
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 2, 0, 0, 1, 1
    !>  \endcode
    !>
    !>  \param[in]  weight  :   The input `contiguous` array of shape `(:)` of
    !>                          <ol>
    !>                              <li>    type `integer` of default kind \IK,
    !>                              <li>    type `integer` of kind \IKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the weights that are to be reweighted (refined).<br>
    !>                          Non-positive values are allowed, but ignored by definition and reset to zero in the output `reweight`.<br>
    !>  \param[in]  skip    :   The input scalar of the same type and kind as the input `weight`,
    !>                          representing the size of the skip in the sequence whose `weight` is given as input.<br>
    !>                          If the input `weight` is of type `integer` of default kind \IK, 
    !>                          then `skip` can be also of type `real` of kind \RKALL.<br>
    !>
    !>  \return
    !>  `reweight`          :   The output object of the same type, kind, rank, shape, and size as the input `weight`
    !>                          containing the refined weights by the input amount `skip`.<br>
    !>                          Any element of the input `weight` that does not contribute
    !>                          to the output refined weight is set to `0` on output.<br>
    !>
    !>  \interface{getReweight}
    !>  \code{.F90}
    !>
    !>      use pm_sampleWeight, only: getReweight
    !>
    !>      reweight(1 : size(weight)) = getReweight(weight, skip)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < skip` must hold.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  The term `reweight` here stands for **refined weight**.<br>
    !>
    !>  \see
    !>  [getReweight](@ref pm_sampleWeight::getReweight)<br>
    !>  [setReweight](@ref pm_sampleWeight::setReweight)<br>
    !>  [getRefined](@ref pm_arrayRefine::getRefined)<br>
    !>  [setRefined](@ref pm_arrayRefine::setRefined)<br>
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [setCompact](@ref pm_arrayCompact::setCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>
    !>  \example{getReweight}
    !>  \include{lineno} example/pm_sampleWeight/getReweight/main.F90
    !>  \compilef{getReweight}
    !>  \output{getReweight}
    !>  \include{lineno} example/pm_sampleWeight/getReweight/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleWeight](@ref test_pm_sampleWeight)
    !>
    !>  \final{getReweight}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface getReweight

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReweight_IK_IK5(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)    , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getReweight_IK_IK4(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)    , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getReweight_IK_IK3(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)    , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getReweight_IK_IK2(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)    , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getReweight_IK_IK1(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)    , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReweight_RK_RK5(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_RK_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        real(RKG)                                           :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getReweight_RK_RK4(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_RK_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        real(RKG)                                           :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getReweight_RK_RK3(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_RK_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        real(RKG)                                           :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getReweight_RK_RK2(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_RK_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        real(RKG)                                           :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getReweight_RK_RK1(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_RK_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        real(RKG)                                           :: reweight(size(weight, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReweight_IK_RK5(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        integer(IK)                                         :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getReweight_IK_RK4(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        integer(IK)                                         :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getReweight_IK_RK3(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        integer(IK)                                         :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getReweight_IK_RK2(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        integer(IK)                                         :: reweight(size(weight, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getReweight_IK_RK1(weight, skip) result(reweight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReweight_IK_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(in)    , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
        integer(IK)                                         :: reweight(size(weight, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a reweighting of the input `weight` vector, such that the sequence represented by output `reweight`
    !>  is the refined (skipped by the amount `skip`) version of the sequence represented by the input `weight`.<br>
    !>
    !>  \brief
    !>  The following examples illustrate the specific functioning of this generic interface:<br>
    !>  \code{.sh}
    !>               skip: 1
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 5, 0, 1, 3, 1
    !>               skip: 2
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 3, 0, 0, 2, 0
    !>               skip: 3
    !>             Weight: 5, 0, 1, 3, 1
    !>      RefinedWeight: 2, 0, 0, 1, 1
    !>  \endcode
    !>
    !>  \param[inout]   weight  :   The input/output `contiguous` array of shape `(:)` of
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK,
    !>                                  <li>    type `integer` of kind \IKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the weights that are to be reweighted (refined).<br>
    !>                              On output, each element of `weight` is reset to the refined value (or zero if the element does not contribute to the refinement).<br>
    !>                              Non-positive values are allowed, but ignored by definition and reset to zero in the output `reweight`.<br>
    !>  \param[in]      skip    :   The input scalar of the same type and kind as the input `weight`,
    !>                              representing the size of the skip in the sequence whose `weight` is given as input.<br>
    !>                              If the input `weight` is of type `integer` of default kind \IK, 
    !>                              then `skip` can be also of type `real` of kind \RKALL.<br>
    !>
    !>  \interface{setReweight}
    !>  \code{.F90}
    !>
    !>      use pm_sampleWeight, only: setReweight
    !>
    !>      reweight(1 : size(weight)) = setReweight(weight, skip)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < skip` must hold.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  The term `reweight` here stands for **refined weight**.<br>
    !>
    !>  \see
    !>  [getReweight](@ref pm_sampleWeight::getReweight)<br>
    !>  [setReweight](@ref pm_sampleWeight::setReweight)<br>
    !>  [getRefined](@ref pm_arrayRefine::getRefined)<br>
    !>  [setRefined](@ref pm_arrayRefine::setRefined)<br>
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [setCompact](@ref pm_arrayCompact::setCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>
    !>  \example{setReweight}
    !>  \include{lineno} example/pm_sampleWeight/setReweight/main.F90
    !>  \compilef{setReweight}
    !>  \output{setReweight}
    !>  \include{lineno} example/pm_sampleWeight/setReweight/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleWeight](@ref test_pm_sampleWeight)
    !>
    !>  \final{setReweight}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface setReweight

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReweight_IK_IK5(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(inout) , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReweight_IK_IK4(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(inout) , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReweight_IK_IK3(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(inout) , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReweight_IK_IK2(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(inout) , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setReweight_IK_IK1(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(inout) , contiguous    :: weight(:)
        integer(IKG)        , intent(in)                    :: skip
        integer(IKG)                                        :: reweight(size(weight, 1, IK))
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReweight_RK_RK5(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_RK_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReweight_RK_RK4(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_RK_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReweight_RK_RK3(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_RK_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReweight_RK_RK2(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_RK_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReweight_RK_RK1(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_RK_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReweight_IK_RK5(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReweight_IK_RK4(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReweight_IK_RK3(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReweight_IK_RK2(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setReweight_IK_RK1(weight, skip)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReweight_IK_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(inout) , contiguous    :: weight(:)
        real(RKG)           , intent(in)                    :: skip
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleWeight ! LCOV_EXCL_LINE