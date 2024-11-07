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
!>  This module contains procedures and generic interfaces for computing the Ziggurat set for for pseudo-random number sampling.
!>
!>  \details
!>  The ziggurat algorithm is an algorithm for pseudo-random number sampling.<br>
!>  Belonging to the class of rejection sampling algorithms, it relies on an underlying source of uniformly-distributed random numbers,
!>  typically from a pseudo-random number generator, as well as precomputed tables.<br>
!>  The algorithm is used to generate values from a monotonically decreasing probability distribution.<br>
!>  It can also be applied to symmetric unimodal distributions, such as the normal distribution,
!>  by choosing a value from one half of the distribution and then randomly choosing which half
!>  the value is considered to have been drawn from.<br>
!>  It was developed by George Marsaglia and others in the 1960s.<br>
!>
!>  \see
!>  [getZigNorm](@ref pm_distNorm::getZigNorm)<br>
!>
!>  \test
!>  [test_pm_ziggurat](@ref test_pm_ziggurat)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_ziggurat

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_ziggurat"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a Ziggurat set for the specified distribution
    !>  that can be subsequently used for random number generation from the distribution.
    !>
    !>  \details
    !>  The following figure illustrates the concept of the layers and rectangle corners.<br>
    !>  \image html pm_ziggurat.png width=700
    !>  The shared regions are the Ziggurat rectangles formed and the thin violet curve represents
    !>  the underlying density function returned by the user-specified function `getFunc()`.<br>
    !>  All layers have equal areas.<br>
    !>
    !>  \param[in]  nlay        :   The input scalar `integer` of default kind \IK, representing the number of
    !>                              Ziggurat rectangles that will be used to approximate the target distribution.<br>
    !>                              One of the most popular values for `nlay` is `256`.<br>
    !>  \param      getFunc     :   The `external` user-specified function that takes a single input **scalar**
    !>                              argument `x` of the same type and kind as the input argument `abserr`.<br>
    !>                              On output, `getFunc` must return a scalar of the same type and kind as the input `x`,
    !>                              containing the value of the target distribution at the specified input `x`.<br>
    !>                              The following illustrates the generic interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(KIND)  , intent(in)    :: x
    !>                                      real(KIND)                  :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `KIND` represents the `real` kind of the input argument `x`, which must be that of the input argument `abserr`.<br>
    !>                              The density function represented by `getFunc()` must be monotonically decreasing to \f$+\infty\f$, starting from \f$x = 0\f$.<br>
    !>  \param      getFuncInv  :   The `external` user-specified function that takes a single input **scalar**
    !>                              argument `func` of the same type and kind as the input argument `abserr`.<br>
    !>                              On output, `getFuncInv` must return a scalar of the same type and kind as the input `func`,
    !>                              containing the location `x` within the domain of the distribution such that `func = getFunc(x)`.<br>
    !>                              In other words, `getFuncInv` is the inverse of `getFunc`.<br>
    !>                              The following illustrates the generic interface of `getFuncInv`,
    !>                              \code{.F90}
    !>                                  function getFuncInv(func) result(x)
    !>                                      real(KIND)  , intent(in)    :: func
    !>                                      real(KIND)                  :: x
    !>                                  end function
    !>                              \endcode
    !>                              where `KIND` represents the `real` kind of the input argument `func`, which must be that of the input argument `abserr`.<br>
    !>  \param      getZigArea  :   The `external` user-specified function that takes a single input **scalar**
    !>                              argument `r` of the same type and kind as the input argument `abserr`.<br>
    !>                              On output, `getZigArea` must return a scalar of the same type and kind as the input `r`,
    !>                              containing the area of individual Ziggurat rectangles `area`.<br>
    !>                              This area is computed as,
    !>                              \f{equation}{
    !>                                  \ms{area} = r \times \ms{getFunc}(r) + \int_r^{+\infty} \ms{getFunc}(x) ~dx ~,
    !>                              \f}
    !>                              where \f$r\f$ is (frequently) the last point in the Ziggurat set (beyond which there is no more rectangles defined, but only a distribution tail).<br>
    !>                              The following illustrates the generic interface of `getFuncInv`,
    !>                              \code{.F90}
    !>                                  function getZigArea(r) result(area)
    !>                                      real(KIND)  , intent(in)    :: r
    !>                                      real(KIND)                  :: area
    !>                                  end function
    !>                              \endcode
    !>                              where `KIND` represents the `real` kind of the input argument `r`, which must be that of the input argument `abserr`.<br>
    !>  \param[out]     abserr  :   The output scalar of type `real` of kind \RKALL.<br>
    !>                              On output it contains an estimate of the error in computing the areas of individual rectangles in the Ziggurat set.<br>
    !>                              The returned value of `abserr` must be typically zero or near `epsilon(abserr)`.<br>
    !>                              A significant deviation strongly indicates lack of convergence to the correct Ziggurat set, in which case, the output set is likely unreliable.<br>
    !>                              As a general rule of thumb, the value of `abserr` should be typically less than `sqrt(epsilon(real(0, kind(abserr))))`.<br>
    !>  \param[in]      abstol  :   The input scalar of same type and kind as the output `abserr`, representing the absolute tolerance used as the stopping criterion of the search.<br>
    !>                              The search iterations continue until the search interval becomes smaller than `abstol` in absolute units.<br>
    !>                              Care must be taken for specifying a reasonable value for `abstol` (see the warnings below).<br>
    !>                              If no suitable value for `abstol` is known a priori, try `abstol = epsilon(0._RKG)**.8 * (abs(lb) + abs(ub))`
    !>                              where `RKG` refers to the kind of the output argument `root`.<br>
    !>                              See also the corresponding argument of [getRoot](@ref pm_mathRoot::getRoot).<br>
    !>                              (**optional** The default value is set by [getRoot](@ref pm_mathRoot::getRoot).)
    !>
    !>  \return
    !>  `zig`                   :   The output array of shape `(1 : 2, 0 : nlay)` of the same type and kind as `abserr`, containing the Ziggurat set.<br>
    !>                              <ol>
    !>                                  <li>    The subset `zig(1, 1 : nlay)` contains the rightmost corners of the rectangles corresponding to
    !>                                          the `nlay` Ziggurat layers, \f$\{x_{i} ~:~ i = 1, \ms{nlay}\}\f$ as depicted in the figure above.<br>
    !>                                  <li>    The subset `zig(2, 1 : nlay)` contains the density function values returned by the user-specified `getFunc()`
    !>                                          corresponding to the rightmost corners of the rectangles stored in the subset `zig(1, 1 : nlay)`.<br>
    !>                                  <li>    The element `zig(1, 0)` contains the fictitious point \f$x_0\f$ representing the rightmost corner
    !>                                          of a hypothetical rectangle of the same area as the lowermost Ziggurat layer (which includes the tail).<br>
    !>                                          The element `zig(2, 0)` contains the corresponding function value at this fictitious point, which is assumed to be `0`.<br>
    !>                                          The zeroth subset is important for random number generation purposes using the Ziggurat algorithm.<br>
    !>                              </ol>
    !>
    !>  \interface{getZig}
    !>  \code{.F90}
    !>
    !>      use pm_ziggurat, only: getZig
    !>
    !>      zig(1 : 2, 0 : nlay) = getZig(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol = abstol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < nlay` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < abstol` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getZig](@ref pm_ziggurat::getZig)<br>
    !>  [getZigNorm](@ref pm_distNorm::getZigNorm)<br>
    !>
    !   \example{getZig}
    !   \include{lineno} example/pm_ziggurat/getZig/main.F90
    !   \compile{getZig}
    !   \output{getZig}
    !   \include{lineno} example/pm_ziggurat/getZig/main.out.F90
    !   \postproc{getZig}
    !   \include{lineno} example/pm_ziggurat/getZig/main.py
    !   \vis{getZig}
    !   \image html pm_ziggurat/getZig/getZig.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_ziggurat](@ref test_pm_ziggurat)<br>
    !>
    !>  \todo
    !>  \pmed
    !>  Examples should be added to the documentation of this generic interface in future.<br>
    !>
    !>  \final{getZig}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getZig

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getZig_RK5(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZig_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                            :: getFunc, getFuncInv, getZigArea
        integer(IK)         , intent(in)                :: nlay
        real(RKG)           , intent(out)               :: abserr
        real(RKG)           , intent(in)    , optional  :: abstol
        real(RKG)                                       :: zig(2, 0 : nlay)
    end function
#endif

#if RK4_ENABLED
    impure module function getZig_RK4(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZig_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                            :: getFunc, getFuncInv, getZigArea
        integer(IK)         , intent(in)                :: nlay
        real(RKG)           , intent(out)               :: abserr
        real(RKG)           , intent(in)    , optional  :: abstol
        real(RKG)                                       :: zig(2, 0 : nlay)
    end function
#endif

#if RK3_ENABLED
    impure module function getZig_RK3(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZig_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                            :: getFunc, getFuncInv, getZigArea
        integer(IK)         , intent(in)                :: nlay
        real(RKG)           , intent(out)               :: abserr
        real(RKG)           , intent(in)    , optional  :: abstol
        real(RKG)                                       :: zig(2, 0 : nlay)
    end function
#endif

#if RK2_ENABLED
    impure module function getZig_RK2(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZig_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                            :: getFunc, getFuncInv, getZigArea
        integer(IK)         , intent(in)                :: nlay
        real(RKG)           , intent(out)               :: abserr
        real(RKG)           , intent(in)    , optional  :: abstol
        real(RKG)                                       :: zig(2, 0 : nlay)
    end function
#endif

#if RK1_ENABLED
    impure module function getZig_RK1(nlay, getFunc, getFuncInv, getZigArea, abserr, abstol) result(zig)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getZig_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                            :: getFunc, getFuncInv, getZigArea
        integer(IK)         , intent(in)                :: nlay
        real(RKG)           , intent(out)               :: abserr
        real(RKG)           , intent(in)    , optional  :: abstol
        real(RKG)                                       :: zig(2, 0 : nlay)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Generate and return an array of coordinates of the lower rightmost corners of rectangles corresponding to an input Ziggurat set,
!    !>  including the highest point of the target density function.
!    !>
!    !>  \brief
!    !>  This is merely a convenience generic interface.<br>
!    !>  The functionality of this interface can be identically generated via the expression
!    !>  \code{.F90}
!    !>      crd = transpose(reshape(zig(1:), [size(zig) / 2, 2]))
!    !>  \endcode
!    !>  assuming the starting index of `zig` is `0`.<br>
!    !>
!    !>  \param[in]  zig :   The input vector of size `(0 : nlay * 2)` of type `real` of kind \RKALL,
!    !>                      that is the direct output of [getZig](@ref pm_ziggurat::getZig) or similar generic interfaces.<br>
!    !>
!    !>  \return
!    !>  `crd`           :   The output vector of shape `(2, nlay)` of the same type and kind as `zig`,
!    !>                      containing the Ziggurat rectangle lower-corner coordinates.<br>
!    !>
!    !>  \interface{getZigCRD}
!    !>  \code{.F90}
!    !>
!    !>      use pm_ziggurat, only: getZigCRD
!    !>
!    !>      crd(1:2, 1:size(zig)/2) = getZigCRD(zig)
!    !>
!    !>  \endcode
!    !>
!    !>  \warning
!    !>  The condition `mod(size(zig), 2) == 1` must hold for the corresponding input arguments.<br>
!    !>  \vericon
!    !>
!    !>  \impure
!    !>
!    !>  \see
!    !>  [getZig](@ref pm_ziggurat::getZig)<br>
!    !>  [getZigNorm](@ref pm_distNorm::getZigNorm)<br>
!    !>
!    !   \example{getZigCRD}
!    !   \include{lineno} example/pm_ziggurat/getZigCRD/main.F90
!    !   \compile{getZigCRD}
!    !   \output{getZigCRD}
!    !   \include{lineno} example/pm_ziggurat/getZigCRD/main.out.F90
!    !   \postproc{getZigCRD}
!    !   \include{lineno} example/pm_ziggurat/getZigCRD/main.py
!    !   \vis{getZigCRD}
!    !   \image html pm_ziggurat/getZigCRD/getZigCRD.RK.png width=700
!    !>
!    !>  \test
!    !>  [test_pm_ziggurat](@ref test_pm_ziggurat)<br>
!    !>
!    !>  \todo
!    !>  \pmed
!    !>  Examples should be added to the documentation of this generic interface in future.<br>
!    !>
!    !>  \final{getZigCRD}
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>
!    interface getZigCRD
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    impure module function getZigCRD_RK5(zig) result(crd)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigCRD_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)                        :: zig(0:)
!        real(RKG)                                       :: crd(2, size(zig, 1, IK) / 2_IK)
!    end function
!#endif
!
!#if RK4_ENABLED
!    impure module function getZigCRD_RK4(zig) result(crd)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigCRD_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)                        :: zig(0:)
!        real(RKG)                                       :: crd(2, size(zig, 1, IK) / 2_IK)
!    end function
!#endif
!
!#if RK3_ENABLED
!    impure module function getZigCRD_RK3(zig) result(crd)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigCRD_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in)                        :: zig(0:)
!        real(RKG)                                       :: crd(2, size(zig, 1, IK) / 2_IK)
!    end function
!#endif
!
!#if RK2_ENABLED
!    impure module function getZigCRD_RK2(zig) result(crd)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigCRD_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in)                        :: zig(0:)
!        real(RKG)                                       :: crd(2, size(zig, 1, IK) / 2_IK)
!    end function
!#endif
!
!#if RK1_ENABLED
!    impure module function getZigCRD_RK1(zig) result(crd)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getZigCRD_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in)                        :: zig(0:)
!        real(RKG)                                       :: crd(2, size(zig, 1, IK) / 2_IK)
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_ziggurat ! LCOV_EXCL_LINE