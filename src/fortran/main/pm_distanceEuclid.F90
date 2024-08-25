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
!>  This module contains procedures and generic interfaces for computing the Euclidean norm of a single point (with respect to origin or a given reference)
!>  or the pairwise Euclidean distances (squared) of a collection of points with respect to another set of reference points,
!>  optionally **without undue overflow or underflow**.
!>
!>  \details
!>  In mathematics, a **norm** is a function from a real or complex vector space to
!>  the non-negative real numbers that behaves in certain ways like the distance from the origin.<br>
!>  It commutes with scaling, obeys a form of the triangle inequality, and is zero only at the origin.<br>
!>  In particular, the **Euclidean distance** of a vector from the origin is a norm, called the **Euclidean norm**, or **2-norm**.<br>
!>  It is commonly defined as the square root of the inner product of a vector with itself.<br>
!>
!>  <b>Usage Directions</b><br>
!>  <ol>
!>      <li>    Use [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid) functional generic interface
!>              to compute distance of a (set of) point(s) from the origin or a (set of) reference(s), optionally without undue overflow/underflow.<br>
!>      <li>    Use [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) functional generic interface or
!>              the equivalent faster [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) procedural generic
!>              interface to compute pairwise (squared) Euclidean distances (i.e., the distance matrix) of a set of points
!>              (with or without undue overflow/underflow).<br>
!>  </ol>
!>
!>  \see
!>  [pm_distanceMahal](@ref pm_distanceMahal)<br>
!>  [pm_distanceHellinger](@ref pm_distanceHellinger)<br>
!>  [pm_distanceManhattan](@ref pm_distanceManhattan)<br>
!>  [pm_distanceMinkowski](@ref pm_distanceMinkowski)<br>
!>
!>  \devnote
!>  While it is possible to merge the generic interfaces for computing the Euclidean, Minkowski, Manhattan,
!>  and other distances into a single module, the decision was made not to do so, as each of these metrics require
!>  certain parameters that are distinct from other methods.<br>
!>  In such a case, there is no point in merging disparate procedures under a single generic interface name.<br>
!>
!>  \test
!>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)<br>
!>
!>  \todo
!>  \phigh
!>  The connection between the different packing schemes of the distance matrix and
!>  the packing methods in [pm_matrixPack](@ref pm_matrixPack) must be further clarified.<br>
!>
!>  \todo
!>  \pmed
!>  A comparison and benchmark with faster less numerically-stable computational methods for pairwise distances might be informative here.<br>
!>  See also a relevant discussion in [stackexchange](https://scicomp.stackexchange.com/questions/30360/fast-and-numerically-stable-pairwise-distance-algorithms).<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2022, 2:38 AM, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distanceEuclid

    use pm_kind, only: SK, IK
    use pm_array, only: nothing, nothing_type
    use pm_matrixPack, only: lfpack, lfpack_type
    use pm_matrixPack, only: rdpack, rdpack_type
    use pm_matrixSubset, only: uppLow, uppLow_type
    use pm_matrixSubset, only: uppLowDia, uppLowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distanceEuclid"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to
    !>  request safe method of computing Euclidean distances (without undue overflow or underflow).<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate the procedures within
    !>  the generic interfaces of module [pm_distanceEuclid](@ref pm_distanceEuclid) generic interfaces.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users can use the sole object parameter instance of this derived type [euclid](@ref pm_distanceEuclid::euclid)
    !>  to request a robust method of computing the Euclidean distances when calling [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)
    !>  and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclid_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: euclid_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [euclid_type](@ref pm_distanceEuclid::euclid_type) that is
    !>  exclusively used to request safe method of computing Euclidean distances (without undue overflow or underflow).<br>
    !>
    !>  \details
    !>  This scalar `parameter` object is exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  See [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) for example usage.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclid}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(euclid_type), parameter :: euclid = euclid_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: euclid
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to
    !>  request unsafe method of computing Euclidean distances (with the possibility of undue overflow or underflow).<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users can use the sole object parameter instance of this derived type [euclidu](@ref pm_distanceEuclid::euclidu)
    !>  to request a non-robust method of computing the Euclidean distances when calling [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)
    !>  and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidu_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: euclidu_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [euclidu_type](@ref pm_distanceEuclid::euclidu_type)that is
    !>  exclusively used to request unsafe method of computing Euclidean distances (with the possibility of undue overflow or underflow).<br>
    !>
    !>  \details
    !>  This scalar `parameter` object is exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  See [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) for example usage.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidu}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(euclidu_type), parameter :: euclidu = euclidu_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: euclidu
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to
    !>  request computing Euclidean volumes (with the possible risk of overflow or underflow).<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users can use the sole object parameter instance of this derived type [euclidv](@ref pm_distanceEuclid::euclidv)
    !>  to request computing the Euclidean squared-distances when calling [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)
    !>  and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidv_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: euclidv_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [euclidv_type](@ref pm_distanceEuclid::euclidv_type)that is
    !>  exclusively used to request computing Euclidean <b>squared</b>-distances (without undue overflow or underflow).<br>
    !>
    !>  \details
    !>  This scalar `parameter` object is exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  See [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) for example usage.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidv}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(euclidv_type), parameter :: euclidv = euclidv_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: euclidv
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to
    !>  request computing Euclidean squared-distances (with the possible risk of overflow or underflow).<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users can use the sole object parameter instance of this derived type [euclidsq](@ref pm_distanceEuclid::euclidsq)
    !>  to request computing the Euclidean squared-distances when calling [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)
    !>  and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidsq_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: euclidsq_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)that is
    !>  exclusively used to request computing Euclidean <b>squared</b>-distances (without undue overflow or underflow).<br>
    !>
    !>  \details
    !>  This scalar `parameter` object is exclusively used to differentiate the procedures within the
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) generic interfaces.<br>
    !>  See [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) and [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid) for example usage.<br>
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidv](@ref pm_distanceEuclid::euclidv)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidv_type](@ref pm_distanceEuclid::euclidv_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \final{euclidsq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(euclidsq_type), parameter :: euclidsq = euclidsq_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: euclidsq
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (squared) Euclidean distance of a (set of) point(s) in `ndim`-dimensions from a reference point (possibly origin),
    !>  optionally robustly without underflow or overflow.<br>
    !>
    !>  \param[in]      x           :   The input scalar (or array of the same rank as other array-like arguments) of the same type and kind as the output `distance`
    !>                                  containing the `x` component of a 3D vector whose Euclidean norm is to be computed.<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `point` and `ref` are missing and `y` and `z` are present.)
    !>  \param[in]      y           :   The input scalar (or array of the same rank as other array-like arguments), of the same type and kind as `x`,
    !>                                  containing the `y` component of a 2D or 3D vector whose Euclidean norm is to be computed.<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `point` and `ref` are missing and `x` an `z` are present.)
    !>  \param[in]      z           :   The input scalar (or array of the same rank as other array-like arguments), of the same type and kind as `x`,
    !>                                  containing the `z` component of a 3D vector whose Euclidean norm is to be computed.<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `point` and `ref` are missing and `x` and `y` are present.)
    !>  \param[in]      point       :   The input `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:npnt)` of the same type and kind as the output `distance`,
    !>                                  containing a (set of `npnt`) point(s) in the `ndim`-dimensional Euclidean space whose distances with respect to the input reference point `ref` must be returned.<br>
    !>                                  (**optional**. It must be present **if and only if** the input arguments `x`, `y`, `z` are missing.)
    !>  \param[in]      ref         :   The input `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:nref)` of the same type and kind as `point`,
    !>                                  containing the (set of `nref`) reference point(s) from which the distance(s) of `point` must be computed.<br>
    !>                                  (**optional**, default = `[(0., i = 1, size(point, 1))]`. It can be present **only if** the input argument `point` is also present.)
    !>  \param[in]      method      :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [euclid](@ref pm_distanceEuclid::euclid) or an object of type [euclid_type](@ref pm_distanceEuclid::euclid_type),
    !>                                              implying that all distance calculations must be done without undue numerical overflow.<br>
    !>                                              This option is computationally the most expensive method.<br>
    !>                                      <li>    the constant [euclidu](@ref pm_distanceEuclid::euclidu) or an object of type [euclidu_type](@ref pm_distanceEuclid::euclidu_type),
    !>                                              implying that all distance calculations must be **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally faster than the [euclid](@ref pm_distanceEuclid::euclid) method.<br>
    !>                                      <li>    the constant [euclidv](@ref pm_distanceEuclid::euclidv) or an object of type [euclidv_type](@ref pm_distanceEuclid::euclidv_type),
    !>                                              implying that the volumes corresponding to all Euclidean distances must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to computing the distances because it avoid costly `sqrt()` operations and runtime overflow 
    !>                                      <li>    the constant [euclidsq](@ref pm_distanceEuclid::euclidsq) or an object of type [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)
    !>                                              implying that the **squared** values of all distance calculations must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to computing the distances because it avoid costly `sqrt()` operations and runtime overflow checks.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [euclid](@ref pm_distanceEuclid::euclid))
    !>
    !>  \return
    !>  `distance`                  :   The output object of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the requested Euclidean (squared) distance(s).<br>
    !>                                  The rank and shape of the output `distance` follows that of the interfaces illustrated below.<br>
    !>
    !>  \interface{getDisEuclid}
    !>  \code{.F90}
    !>
    !>      use pm_distanceEuclid, only: getDisEuclid
    !>
    !>      ! distance with respect to origin.
    !>
    !>      distance = getDisEuclid(x, y, z) ! elemental
    !>      distance = getDisEuclid(point(1:ndim), method = method)
    !>      distance(1:npnt) = getDisEuclid(point(1:ndim, 1:npnt), method = method)
    !>
    !>      ! distance with respect to custom reference.
    !>
    !>      distance = getDisEuclid(point(1:ndim), ref(1:ndim), method = method)
    !>      distance(1:nref) = getDisEuclid(point(1:ndim), ref(1:ndim, 1:nref), method = method)
    !>      distance(1:npnt) = getDisEuclid(point(1:ndim, 1:npnt), ref(1:ndim), method = method)
    !>      distance(1:npnt, 1:nref) = getDisEuclid(point(1:ndim, 1:npnt), ref(1:ndim, 1:nref), method = method)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(point, 1) == size(ref, 1)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>  The `elemental` attribute does not apply to interfaces that take `point` as an input argument.<br>
    !>
    !>  \note
    !>  The Fortran standard provides the intrinsic procedure `norm2()` for computing the Euclidean norm of a vector.<br>
    !>  However, the standard does not enforce robustness of the intrinsic procedure with respect to possible underflows or overflows.<br>
    !>  The procedures of this module ensure robustness of the distance computations.<br>
    !>  This will inevitably lead to worse runtime performance compared to the compiler
    !>  implementations of the intrinsic routine that do not respect robustness.<br>
    !>  Use the routines of this module in place of the Fortran intrinsics
    !>  if you believe there is a possibility of under/over-flow.
    !>
    !>  \note
    !>  The procedures of this module can be used for a robust computation of `abs(x)` when `x` is a large `complex` value.<br>
    !>  In such a case, calling [getDisEuclid([x%re, x%im])](@ref pm_distanceEuclid::getDisEuclid) would be equivalent to `abs(x)` intrinsic operation.<br>
    !>  However, note that the Fortran standard already offers a better intrinsic alternative to the routines of this procedure for this task,
    !>  namely `hypot()` which is robust against overflow and underflow.<br>
    !>
    !>  \note
    !>  This generic interface intentionally does not have explicit procedures for 2D Euclidean
    !>  distance `(x, y)` because the Fortran intrinsic procedure `hypot()` already serves the purpose.<br>
    !>
    !>  \lapack{3.11}
    !>  `dlapy3`
    !>
    !>  \see
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [setDisEuclid](@ref pm_distanceEuclid::setDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>  Intrinsic Fortran procedure `hypot(x, y)` (robust)<br>
    !>  Intrinsic Fortran procedure `norm2(x(:))` (unsafe)<br>
    !>
    !>  \lapack{3.11}
    !>  `dlapy3`
    !>
    !>  \example{getDisEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/getDisEuclid/main.F90
    !>  \compilef{getDisEuclid}
    !>  \output{getDisEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/getDisEuclid/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)
    !>
    !>  \todo
    !>  \pmed
    !>  A benchmark comparison with the equivalent compiler implementations would be informative.
    !>
    !>  \final{getDisEuclid}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! x, y, z

    interface getDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getDisEuclidXYZ_RK5(x, y, z) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisEuclidXYZ_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: x, y, z
        real(RKG)                                   :: distance
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getDisEuclidXYZ_RK4(x, y, z) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisEuclidXYZ_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: x, y, z
        real(RKG)                                   :: distance
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getDisEuclidXYZ_RK3(x, y, z) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisEuclidXYZ_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: x, y, z
        real(RKG)                                   :: distance
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getDisEuclidXYZ_RK2(x, y, z) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisEuclidXYZ_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: x, y, z
        real(RKG)                                   :: distance
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getDisEuclidXYZ_RK1(x, y, z) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisEuclidXYZ_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: x, y, z
        real(RKG)                                   :: distance
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! vector, matrix

    interface getDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D1_XX_RK5(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D1_XX_RK4(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D1_XX_RK3(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D1_XX_RK2(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D1_XX_RK1(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D2_XX_RK5(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D2_XX_RK4(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D2_XX_RK3(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D2_XX_RK2(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D2_XX_RK1(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D1_D1_RK5(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D1_D1_RK4(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D1_D1_RK3(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D1_D1_RK2(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D1_D1_RK1(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D1_D2_RK5(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance(size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D1_D2_RK4(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance(size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D1_D2_RK3(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance(size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D1_D2_RK2(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance(size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D1_D2_RK1(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D1_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)                                           :: distance(size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D2_D1_RK5(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D2_D1_RK4(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D2_D1_RK3(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D2_D1_RK2(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D2_D1_RK1(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDE_D2_D2_RK5(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDE_D2_D2_RK4(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDE_D2_D2_RK3(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDE_D2_D2_RK2(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDE_D2_D2_RK1(point, ref, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDE_D2_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(ref, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (squared) Euclidean distance of a (set of) point(s) in `ndim`-dimensions from a reference point (possibly origin),
    !>  optionally robustly without underflow or overflow.<br>
    !>
    !>  \param[out]     distance    :   The output object of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the requested Euclidean (squared) distance(s).<br>
    !>                                  The rank and shape of the output `distance` follows that of the interfaces illustrated below.<br>
    !>  \param[in]      point       :   The input `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:npnt)` of the same type and kind as the output `distance`,
    !>                                  containing a (set of `npnt`) point(s) in the `ndim`-dimensional Euclidean space whose distances with respect to the input reference point `ref` must be returned.<br>
    !>  \param[in]      ref         :   The input `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:nref)` of the same type and kind as `point`,
    !>                                  containing the (set of `nref`) reference point(s) from which the distance(s) of `point` must be computed.<br>
    !>                                  (**optional**, default = `[(0., i = 1, size(point, 1))]`.)
    !>  \param[in]      method      :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [euclid](@ref pm_distanceEuclid::euclid) or an object of type [euclid_type](@ref pm_distanceEuclid::euclid_type),
    !>                                              implying that all distance calculations must be done without undue numerical overflow.<br>
    !>                                              This option is computationally the most expensive method.<br>
    !>                                      <li>    the constant [euclidu](@ref pm_distanceEuclid::euclidu) or an object of type [euclidu_type](@ref pm_distanceEuclid::euclidu_type),
    !>                                              implying that all distance calculations must be **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally faster than the [euclid](@ref pm_distanceEuclid::euclid) method.<br>
    !>                                      <li>    the constant [euclidv](@ref pm_distanceEuclid::euclidv) or an object of type [euclidv_type](@ref pm_distanceEuclid::euclidv_type),
    !>                                              implying that the volumes corresponding to all Euclidean distances must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to computing the distances because it avoid costly `sqrt()` operations and runtime overflow 
    !>                                      <li>    the constant [euclidsq](@ref pm_distanceEuclid::euclidsq) or an object of type [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)
    !>                                              implying that the **squared** values of all distance calculations must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to computing the distances because it avoid costly `sqrt()` operations and runtime overflow checks.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [euclid](@ref pm_distanceEuclid::euclid))
    !>
    !>  \interface{setDisEuclid}
    !>  \code{.F90}
    !>
    !>      use pm_distanceEuclid, only: setDisEuclid
    !>
    !>      ! distance with respect to origin.
    !>
    !>      call setDisEuclid(distance, point(1:ndim), method)
    !>      call setDisEuclid(distance(1:npnt), point(1:ndim, 1:npnt), method)
    !>
    !>      ! distance with respect to custom reference.
    !>
    !>      call setDisEuclid(distance, point(1:ndim), ref(1:ndim), method)
    !>      call setDisEuclid(distance(1:nref), point(1:ndim), ref(1:ndim, 1:nref), method)
    !>      call setDisEuclid(distance(1:npnt), point(1:ndim, 1:npnt), ref(1:ndim), method)
    !>      call setDisEuclid(distance(1:npnt, 1:nref), point(1:npnt), ref(1:nref), method)
    !>      call setDisEuclid(distance(1:npnt, 1:nref), point(1:ndim, 1:npnt), ref(1:ndim, 1:nref), method)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The shapes of `distance`, `point`, and `ref` must be consistent as in the above interface at all times.<br>
    !>  The condition `size(point, 1) == size(ref, 1) .or. rank(distance) == 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `dlapy3`
    !>
    !>  \see
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [setDisEuclid](@ref pm_distanceEuclid::setDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>  Intrinsic Fortran procedure `hypot(x, y)` (robust)<br>
    !>  Intrinsic Fortran procedure `norm2(x(:))` (unsafe)<br>
    !>
    !>  \example{setDisEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/setDisEuclid/main.F90
    !>  \compilef{setDisEuclid}
    !>  \output{setDisEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/setDisEuclid/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setDE_MEQ_D1_D1_D2_RK5()
    !>         || ||| || || || |||
    !>         || ||| || || || The type and kind parameters.
    !>         || ||| || || The dimension of `ref` array: D1 => vector, D2 => matrix, XX => `ref` missing.
    !>         || ||| || The dimension of `point` array: D1 => vector, D2 => matrix.
    !>         || ||| The dimension of the output `distance`: D0 => scalar, D1 => vector, D2 => matrix.
    !>         || The Method of Euclidean distance computation: MED => euclid_type (default/safe), MEU => eulidu_type (unsafe), MEV => eulidv_type (volume), MEQ => eulidsq_type (squared)
    !>         The abbreviation for DisEuclid to shorten procedure names.
    !>  \endcode
    !>
    !>  \todo
    !>  \pmed
    !>  A benchmark comparison with the equivalent compiler implementations would be informative.
    !>
    !>  \final{setDisEuclid}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! euclid

    interface setDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D0_D1_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D0_D1_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D0_D1_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D0_D1_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D0_D1_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D1_D2_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D1_D2_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D1_D2_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D1_D2_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D1_D2_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D0_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D0_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D0_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D0_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D0_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D0_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D1_D1_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D1_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D1_D1_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D1_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D1_D1_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D1_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D1_D1_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D1_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D1_D1_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D1_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D1_D2_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D1_D2_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D1_D2_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D1_D2_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D1_D2_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D1_D2_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D2_D2_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D2_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D2_D2_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D2_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D2_D2_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D2_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D2_D2_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D2_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D2_D2_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D2_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MED_D2_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MED_D2_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MED_D2_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MED_D2_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MED_D2_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MED_D2_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! euclidu

    interface setDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D0_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D0_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D1_D1_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D1_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D1_D1_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D1_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D1_D1_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D1_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D1_D1_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D1_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D1_D1_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D1_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D1_D2_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D1_D2_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D2_D2_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D2_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D2_D2_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D2_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D2_D2_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D2_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D2_D2_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D2_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D2_D2_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D2_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEU_D2_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEU_D2_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEU_D2_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEU_D2_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEU_D2_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEU_D2_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! euclidv

    interface setDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D0_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D0_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D1_D1_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D1_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D1_D1_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D1_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D1_D1_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D1_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D1_D1_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D1_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D1_D1_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D1_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D1_D2_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D1_D2_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D2_D2_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D2_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D2_D2_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D2_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D2_D2_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D2_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D2_D2_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D2_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D2_D2_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D2_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEV_D2_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEV_D2_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEV_D2_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEV_D2_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEV_D2_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEV_D2_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidv_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! euclidsq

    interface setDisEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_XX_RK5(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_XX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_XX_RK4(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_XX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_XX_RK3(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_XX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_XX_RK2(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_XX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_XX_RK1(distance, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_XX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D0_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D0_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)                   :: distance
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D1_D1_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D1_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D1_D1_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D1_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D1_D1_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D1_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D1_D1_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D1_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D1_D1_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D1_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D1_D2_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D1_D2_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D2_D2_D2_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D2_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D2_D2_D2_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D2_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D2_D2_D2_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D2_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D2_D2_D2_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D2_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D2_D2_D2_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D2_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:,:)
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDE_MEQ_D2_D1_D1_RK5(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDE_MEQ_D2_D1_D1_RK4(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDE_MEQ_D2_D1_D1_RK3(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDE_MEQ_D2_D1_D1_RK2(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDE_MEQ_D2_D1_D1_RK1(distance, point, ref, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDE_MEQ_D2_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: ref(:)
        real(RKG)           , intent(in)    , contiguous    :: point(:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the full or a subset of the Euclidean (squared) distance matrix of the input set of `npnt` points in `ndim` dimensions.
    !>
    !>  \param[in]      pack        :   The input scalar that can be:
    !>                                  <ol>
    !>                                      <li>    the constant [rdpack](@ref pm_matrixPack::rdpack) or an object of type [rdpack_type](@ref pm_matrixPack::rdpack_type),
    !>                                              implying the use of Rectangular Default Packing format for the output matrix.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar that can be:
    !>                                  <ol>
    !>                                      <li>    the constant [uppLowDia](@ref pm_matrixSubset::uppLowDia) or an object of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                                              indicating that the output `distance` must contain the full distance matrix of shape `(1:npnt, 1:npnt)` including the zero-valued diagonals.<br>
    !>                                      <li>    the constant [uppLow](@ref pm_matrixSubset::uppLow) or an object of type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                                              indicating that the output `distance` must exclude the zero-valued diagonals from the distance matrix yielding a distance matrix of shape `(1:npnt - 1, 1:npnt)`.<br>
    !>                                              **Motivation:** The zero-valued diagonal elements of the distance matrix are are frequently troubling for subsequent vector operations on the output distance matrix.<br>
    !>                                              Such vector operations include but are not limited to finding the extrema of distances, for example, the nearest and farthest neighbors.<br>
    !>                                              This `subset` value offers a fast convenient method of excluding self-distance values from the output distance matrix
    !>                                              such that each column `(1:npnt-1 , i)` of the distance matrix contains only the distances of `point(1:ndim, i)` with all other `npnt - 1` points in `point`.<br>
    !>                                              For example, finding the nearest neighbor of the points using the output distance matrix would be as simple as `minval(distance, 1)`.<br>
    !>                                              Finding the actual index of the point that is the nearest neighbor to each point would be slightly more involved as a two-step process:<br>
    !>                                              \code{.F90}
    !>                                                  nn1loc(1 : npnt) = minloc(distance(1 : npnt - 1, 1 : npnt), 1)
    !>                                                  nn1loc = merge(nn1loc, nn1loc + 1, getRange(1, npnt) <= nn1loc)
    !>                                              \endcode
    !>                                              where `nn1loc` is the vector of indices of the first nearest neighbors such that `point(:,nn1loc(i))` is the nearest neighbor to `point(:,i)`.<br>
    !>                                  </ol>
    !>  \param[in]      point       :   The input `contiguous` matrix of shape `(1:ndim, 1:npnt)` of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing `npnt` points in the `ndim`-dimensional Euclidean space
    !>                                  whose distances with respect to each other must be computed and returned.<br>
    !>  \param[in]      method      :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [euclid](@ref pm_distanceEuclid::euclid) or an object of type [euclid_type](@ref pm_distanceEuclid::euclid_type),
    !>                                              implying that all distance calculations must be done without undue numerical overflow.<br>
    !>                                              This option is computationally the most expensive method.<br>
    !>                                      <li>    the constant [euclidu](@ref pm_distanceEuclid::euclidu) or an object of type [euclidu_type](@ref pm_distanceEuclid::euclidu_type),
    !>                                              implying that all distance calculations must be **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally faster than the [euclid](@ref pm_distanceEuclid::euclid) method.<br>
    !>                                      <li>    the constant [euclidsq](@ref pm_distanceEuclid::euclidsq) or an object of type [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)
    !>                                              implying that the **squared** values of all distance calculations must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to constructing the distance matrix because it avoid costly `sqrt()` operations and runtime overflow checks.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [euclid](@ref pm_distanceEuclid::euclid))
    !>
    !>  \return
    !>  `distance`                  :   The output `contiguous` array of rank `2` of the same type and kind as the input argument `point`.<br>
    !>                                  On output, it contains the requested `subset` of the (squared) distance matrix in the specified packing format `pack`.<br>
    !>                                  Any element of `distance` that is not included in the specified `subset` will remain intact, if any such element exists.<br>
    !>
    !>  \interface{getDisMatEuclid}
    !>  \code{.F90}
    !>
    !>      use pm_distanceEuclid, only: getDisMatEuclid
    !>
    !>      distance(1:npnt, 1:npnt) = getDisMatEuclid(pack, subset, point(1:ndim, 1:npnt), method) ! subset = uppLowDia, pack = rdpack
    !>      distance(1:npnt-1, 1:npnt) = getDisMatEuclid(pack, subset, point(1:ndim, 1:npnt), method) ! subset = uppLow, pack = rdpack
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(point, 1) == size(point, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition `shape(distance) == [size(point, 1), size(point, 1)] .or. .not. same_type_as(subset, uppLowDia)` must hold for the corresponding input arguments.<br>
    !>  The condition `shape(distance) == [size(point, 1) - 1, size(point, 1)] .or. .not. same_type_as(subset, uppLow)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  The input arguments `pack, subset` appear first for a good reason:
    !>  To allow the possibility of adding of similarly-named arguments for the input `point` matrix.
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [setDisEuclid](@ref pm_distanceEuclid::setDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \example{getDisMatEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/getDisMatEuclid/main.F90
    !>  \compilef{getDisMatEuclid}
    !>  \output{getDisMatEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/getDisMatEuclid/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)
    !>
    !>  \todo
    !>  \phigh
    !>  This generic interface must be extended to allow other packing and subsets of the output distance matrix.
    !>
    !>  \final{getDisMatEuclid}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! rdpack

    interface getDisMatEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDME_RDP_FUL_RK5(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_FUL_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDME_RDP_FUL_RK4(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_FUL_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDME_RDP_FUL_RK3(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_FUL_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDME_RDP_FUL_RK2(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_FUL_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDME_RDP_FUL_RK1(point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_FUL_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDME_RDP_ULD_RK5(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK4_ENABLED
    PURE module function getDME_RDP_ULD_RK4(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK3_ENABLED
    PURE module function getDME_RDP_ULD_RK3(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK2_ENABLED
    PURE module function getDME_RDP_ULD_RK2(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK1_ENABLED
    PURE module function getDME_RDP_ULD_RK1(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK), size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDME_RDP_ULX_RK5(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK) - 1, size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK4_ENABLED
    PURE module function getDME_RDP_ULX_RK4(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK) - 1, size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK3_ENABLED
    PURE module function getDME_RDP_ULX_RK3(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK) - 1, size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK2_ENABLED
    PURE module function getDME_RDP_ULX_RK2(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK) - 1, size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK1_ENABLED
    PURE module function getDME_RDP_ULX_RK1(pack, subset, point, method) result(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDME_RDP_ULX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)                                           :: distance(size(point, 2, IK) - 1, size(point, 2, IK))
        class(*)            , intent(in)    , optional      :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the full or a subset of the Euclidean (squared) distance matrix of the input set of `npnt` points in `ndim` dimensions.
    !>
    !>  \param[inout]   distance    :   The input/output `contiguous` array of rank `2` of the same type and kind as the input argument `point`.<br>
    !>                                  On output, it contains the requested `subset` of the (squared) distance matrix in the specified packing format `pack`.<br>
    !>                                  Any element of `distance` that is not included in the specified `subset` will remain intact, if any such element exists.<br>
    !>  \param[in]      pack        :   The input scalar that can be:
    !>                                  <ol>
    !>                                      <li>    the constant [rdpack](@ref pm_matrixPack::rdpack) or an object of type [rdpack_type](@ref pm_matrixPack::rdpack_type),
    !>                                              implying the use of Rectangular Default Packing format for the output matrix.<br>
    !                                       <li>    the constant [lfpack](@ref pm_matrixPack::lfpack) or an object of type [lfpack_type](@ref pm_matrixSubset::lfpack_type),
    !                                               implying the use of Linear Full Packing format for the output matrix.<br>
    !                                               This means that the output matrix must be contiguous vector of appropriate size.<br>
    !                                               On output, `distance` will be a **dense** ([linear contiguous packed](@ref pm_matrixPack::lfpack_type))
    !                                               vector containing the pairwise distances of the points in `point` from itself,
    !                                               corresponding to lower triangle of symmetric square matrix of pairwise distances of `point` from itself.<br>
    !                                               By definition, the diagonal elements of the square distance matrix are zeros and not included in the output dense vector.<br>
    !                                               The following figure illustrates the storage layout for the dense vector format compared to the corresponding symmetric square distance matrix.<br>
    !                                               \image html pm_distanceEuclid@dense.png width=500
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar that can be:
    !>                                  <ol>
    !>                                      <li>    the constant [uppLowDia](@ref pm_matrixSubset::uppLowDia) or an object of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                                              indicating that the output `distance` must contain the full distance matrix of shape `(1:npnt, 1:npnt)` including the zero-valued diagonals.<br>
    !>                                      <li>    the constant [uppLow](@ref pm_matrixSubset::uppLow) or an object of type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                                              indicating that the output `distance` must exclude the zero-valued diagonals from the distance matrix yielding a distance matrix of shape `(1:npnt - 1, 1:npnt)`.<br>
    !>                                              **Motivation:** The zero-valued diagonal elements of the distance matrix are are frequently troubling for subsequent vector operations on the output distance matrix.<br>
    !>                                              Such vector operations include but are not limited to finding the extrema of distances, for example, the nearest and farthest neighbors.<br>
    !>                                              This `subset` value offers a fast convenient method of excluding self-distance values from the output distance matrix
    !>                                              such that each column `(1:npnt-1 , i)` of the distance matrix contains only the distances of `point(1:ndim, i)` with all other `npnt - 1` points in `point`.<br>
    !>                                              For example, finding the nearest neighbor of the points using the output distance matrix would be as simple as `minval(distance, 1)`.<br>
    !>                                              Finding the actual index of the point that is the nearest neighbor to each point would be slightly more involved as a two-step process:<br>
    !>                                              \code{.F90}
    !>                                                  nn1loc(1 : npnt) = minloc(distance(1 : npnt - 1, 1 : npnt), 1)
    !>                                                  nn1loc = merge(nn1loc, nn1loc + 1, getRange(1, npnt) <= nn1loc)
    !>                                              \endcode
    !>                                              where `nn1loc` is the vector of indices of the first nearest neighbors such that `point(:,nn1loc(i))` is the nearest neighbor to `point(:,i)`.<br>
    !>  \cond excluded
    !                                       <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type),
    !                                               implying that only the upper-diagonal subset of the distance matrix must be returned.<br>
    !                                       <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type),
    !                                               implying that only the lower-diagonal subset of the distance matrix must be returned.<br>
    !                                       <li>    the constant [upp](@ref pm_matrixSubset::upp) or an object of type [upp_type](@ref pm_matrixSubset::upp_type),
    !                                               implying that only the upper-diagonal subset of the distance matrix must be returned.<br>
    !                                       <li>    the constant [low](@ref pm_matrixSubset::low) or an object of type [low_type](@ref pm_matrixSubset::low_type),
    !                                               implying that only the lower-diagonal subset of the distance matrix must be returned.<br>
    !>  \endcond excluded
    !>                                  </ol>
    !>  \param[in]      point       :   The input `contiguous` matrix of shape `(1:ndim, 1:npnt)` of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing `npnt` points in the `ndim`-dimensional Euclidean space
    !>                                  whose distances with respect to each other must be computed and returned.<br>
    !>  \param[in]      method      :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [euclid](@ref pm_distanceEuclid::euclid) or an object of type [euclid_type](@ref pm_distanceEuclid::euclid_type),
    !>                                              implying that all distance calculations must be done without undue numerical overflow.<br>
    !>                                              This option is computationally the most expensive method.<br>
    !>                                      <li>    the constant [euclidu](@ref pm_distanceEuclid::euclidu) or an object of type [euclidu_type](@ref pm_distanceEuclid::euclidu_type),
    !>                                              implying that all distance calculations must be **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally faster than the [euclid](@ref pm_distanceEuclid::euclid) method.<br>
    !>                                      <li>    the constant [euclidsq](@ref pm_distanceEuclid::euclidsq) or an object of type [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)
    !>                                              implying that the **squared** values of all distance calculations must be returned **without** runtime checks for numerical overflow.<br>
    !>                                              This option is computationally the fastest approach to constructing the distance matrix because it avoid costly `sqrt()` operations and runtime overflow checks.<br>
    !>                                  </ol>
    !>
    !>  \interface{setDisMatEuclid}
    !>  \code{.F90}
    !>
    !>      use pm_distanceEuclid, only: setDisMatEuclid
    !>
    !>      call setDisMatEuclid(distance(1:npnt, 1:npnt), pack, subset, point(1:ndim, 1:npnt), method) ! subset /= uppLow
    !>      call setDisMatEuclid(distance(1:npnt-1, 1:npnt), pack, subset, point(1:ndim, 1:npnt), method) ! subset = uppLow
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(point, 1) == size(point, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition `shape(distance) == [size(point, 1), size(point, 1)] .or. .not. same_type_as(subset, uppLowDia)` must hold for the corresponding input arguments.<br>
    !>  The condition `shape(distance) == [size(point, 1) - 1, size(point, 1)] .or. .not. same_type_as(subset, uppLow)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [euclid](@ref pm_distanceEuclid::euclid)<br>
    !>  [euclidu](@ref pm_distanceEuclid::euclidu)<br>
    !>  [euclidsq](@ref pm_distanceEuclid::euclidsq)<br>
    !>  [euclid_type](@ref pm_distanceEuclid::euclid_type)<br>
    !>  [euclidu_type](@ref pm_distanceEuclid::euclidu_type)<br>
    !>  [euclidsq_type](@ref pm_distanceEuclid::euclidsq_type)<br>
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [setDisEuclid](@ref pm_distanceEuclid::setDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \example{setDisMatEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/setDisMatEuclid/main.F90
    !>  \compilef{setDisMatEuclid}
    !>  \output{setDisMatEuclid}
    !>  \include{lineno} example/pm_distanceEuclid/setDisMatEuclid/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)
    !>
    !>  \todo
    !>  \phigh
    !>  This generic interface must be extended to allow other packing and subsets of the output distance matrix.
    !>
    !>  \final{setDisMatEuclid}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! rdpack, euclid

    interface setDisMatEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MED_RDP_ULD_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MED_RDP_ULD_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MED_RDP_ULD_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MED_RDP_ULD_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MED_RDP_ULD_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MED_RDP_ULX_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MED_RDP_ULX_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MED_RDP_ULX_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MED_RDP_ULX_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MED_RDP_ULX_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MED_RDP_ULX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclid_type)   , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! rdpack, euclidu

    interface setDisMatEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULD_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULD_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULD_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULD_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULD_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULX_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULX_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULX_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULX_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MEU_RDP_ULX_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEU_RDP_ULX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidu_type)  , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! rdpack, euclidsq

    interface setDisMatEuclid

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULD_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULD_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULD_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULD_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULD_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLowDia_type), intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULX_RK5(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULX_RK4(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULX_RK3(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULX_RK2(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDME_MEQ_RDP_ULX_RK1(distance, pack, subset, point, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDME_MEQ_RDP_ULX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: point(:,:)
        real(RKG)           , intent(out)   , contiguous    :: distance(:,:)
        type(euclidsq_type) , intent(in)                    :: method
        type(uppLow_type)   , intent(in)                    :: subset
        type(rdpack_type)   , intent(in)                    :: pack
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distanceEuclid
