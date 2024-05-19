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
!>  This module contains classes and procedures for scaling (i.e., multiplying)
!>  univariate or multivariate samples by arbitrary amounts along specific directions.
!>
!>  \devnote
!>  While it is tempting to add generic interfaces for automatic scaling by the inverse of the standard deviation of the sample (in the absence of arbitrary `amount` argument),
!>  such interfaces were not added to this module for the following reasons:<br>
!>  <ol>
!>      <li>    Why should scaling by inverse of the standard deviation be the default behavior?
!>      <li>    Even though scaling by inverse of the standard deviation is popular, its implementation as the default normalization in 
!>              the generic interfaces of this module requires inclusion of sample `weight` and variance `correction` arguments,
!>              thus significantly complicating the interfaces of this module with little gain.<br>
!>  </ol>
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
!>
!>  \test
!>  [test_pm_sampleScale](@ref test_pm_sampleScale)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleScale

    use pm_kind, only: SK, IK, LK
    use pm_matrixTrans, only: transHerm, transHerm_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleScale"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type whose instances are meant to signify a sample scaling
    !>  by an amount equal to the inverse of the sample standard deviation or an equivalent measure.<br>
    !>
    !>  \details
    !>  For example usage, see the relevant interfaces that use instances of this derived type.<br>
    !>
    !>  \see
    !>  [stdscale](@ref pm_sampleScale::stdscale)<br>
    !>
    !>  \test
    !>  [test_pm_sampleScale](@ref test_pm_sampleScale)
    !>
    !>  \final{stdscale_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX
    type :: stdscale_type; end type
    type(stdscale_type), parameter :: stdscale = stdscale_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: stdscale
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is scaled by the specified input `amount` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `amount` is the square-root of the inverse of the variance of the (zero-mean) sample, then the returned sample will have a variance of one.<br>
    !>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample to be scaled.<br>
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK,
    !>                              whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  amount      :   The input scalar or the `contiguous` vector of,
    !>                              <ol>
    !>                                  <li>    the same type and kind as `sample`,
    !>                                  <li>    type `real` of the same kind as that of `sample` if `sample` is of type `complex`,
    !>                              </ol>
    !>                              representing the amount by which the input sample must be scaled.<br>
    !>                              <ol>
    !>                                  <li>    If the input `rank(sample) = 1`, then `amount` must be a scalar.<br>
    !>                                  <li>    If the input `rank(sample) = 2`, then `amount` must be a vector of size `ndim`.<br>
    !>                              </ol>
    !>                              If the sample is to be scaled to have a unit variance, then the input `amount`
    !>                              corresponds to the square-root of the inverse of the variance of the sample that is returned by procedures
    !>                              collectively represented by the generic interface [getVar](@ref pm_sampleVar::getVar).<br>
    !>                              Note that the size of the input `amount` must be consistent with the size of the input `sample`:<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(amount) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(amount) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>  \param[in]  operation   :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [transHerm](@ref pm_matrixTrans::transHerm)
    !>                                          implying that a Hermitian transpose of the input sample must be returned.<br>
    !>                                          In other words, the actions `getScaled(sample, dim, transHerm)` and `transpose(conjg(getScaled(sample, dim)))` are equivalent.<br>
    !>                                          Specifying [transHerm](@ref pm_matrixTrans::transHerm) for `source` of type other than `complex` is identical to
    !>                                          specifying [transSymm](@ref pm_matrixTrans::transSymm) for `source` of type other than `complex`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `.false.`. It **can** be present **only if** the condition `rank(sample) == 2` holds.)
    !>
    !>  \return
    !>  `sampleScaled`          :   The output object of the same type and kind and rank as `sample`, containing the scaled sample.<br>
    !>                              <ol>
    !>                                  <li>    If the input `sample` is a vector, then `sampleScaled` has the same shape and size as that of `sample`.<br>
    !>                                  <li>    If the input `sample` is a matrix of shape `(nrow, ncol)`, then
    !>                                          <ol>
    !>                                              <li>    `sampleScaled` has the shape `(nrow, ncol)` if `operation` is missing.
    !>                                              <li>    `sampleScaled` has the shape `(ncol, nrow)` if `operation = transHerm`.
    !>                                          </ol>
    !>                              </ol>
    !>
    !>  \interface{getScaled}
    !>  \code{.F90}
    !>
    !>      use pm_sampleScale, only: getScaled
    !>
    !>      sampleScaled(:) = getScaled(sample(:), amount)
    !>      sampleScaled(:,:) = getScaled(sample(:,:), dim, amount(:))
    !>      sampleScaled(:,:) = getScaled(sample(:,:), dim, amount(:), operation)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(amount) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getScaled](@ref pm_sampleScale::getScaled)<br>
    !>  [setScaled](@ref pm_sampleScale::setScaled)<br>
    !>
    !>  \example{getScaled}
    !>  \include{lineno} example/pm_sampleScale/getScaled/main.F90
    !>  \compilef{getScaled}
    !>  \output{getScaled}
    !>  \include{lineno} example/pm_sampleScale/getScaled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleScale](@ref test_pm_sampleScale)
    !>
    !>  \todo
    !>  \pvlow
    !>  The functionality of this interface can be expanded to include scaling of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \final{getScaled}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 2:48 AM, August 22, 2021, Dallas, TX

    ! RK ARK

    interface getScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getScaled_ARK_ONO_D1_RK5(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample(:)
        real(RKG)           , intent(in)                :: amount
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getScaled_ARK_ONO_D1_RK4(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample(:)
        real(RKG)           , intent(in)                :: amount
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getScaled_ARK_ONO_D1_RK3(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample(:)
        real(RKG)           , intent(in)                :: amount
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getScaled_ARK_ONO_D1_RK2(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample(:)
        real(RKG)           , intent(in)                :: amount
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getScaled_ARK_ONO_D1_RK1(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample(:)
        real(RKG)           , intent(in)                :: amount
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getScaled_ARK_ONO_D2_RK5(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getScaled_ARK_ONO_D2_RK4(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getScaled_ARK_ONO_D2_RK3(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getScaled_ARK_ONO_D2_RK2(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getScaled_ARK_ONO_D2_RK1(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKG)                                       :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getScaled_ARK_OTH_D2_RK5(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKG)                                       :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getScaled_ARK_OTH_D2_RK4(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKG)                                       :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getScaled_ARK_OTH_D2_RK3(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKG)                                       :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getScaled_ARK_OTH_D2_RK2(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKG)                                       :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getScaled_ARK_OTH_D2_RK1(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: sample(:,:)
        real(RKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKG)                                       :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getScaled

    ! CK ARK

    interface getScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ACK_ONO_D1_CK5(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        complex(CKG)        , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ACK_ONO_D1_CK4(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        complex(CKG)        , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ACK_ONO_D1_CK3(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        complex(CKG)        , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ACK_ONO_D1_CK2(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        complex(CKG)        , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ACK_ONO_D1_CK1(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        complex(CKG)        , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ACK_ONO_D2_CK5(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ACK_ONO_D2_CK4(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ACK_ONO_D2_CK3(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ACK_ONO_D2_CK2(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ACK_ONO_D2_CK1(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_ONO_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ACK_OTH_D2_CK5(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_OTH_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ACK_OTH_D2_CK4(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_OTH_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ACK_OTH_D2_CK3(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_OTH_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ACK_OTH_D2_CK2(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_OTH_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ACK_OTH_D2_CK1(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ACK_OTH_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        complex(CKG)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getScaled

    ! CK ACK

    interface getScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ARK_ONO_D1_CK5(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        real(CKG)           , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ARK_ONO_D1_CK4(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        real(CKG)           , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ARK_ONO_D1_CK3(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        real(CKG)           , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ARK_ONO_D1_CK2(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        real(CKG)           , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ARK_ONO_D1_CK1(sample, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:)
        real(CKG)           , intent(in)                :: amount
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ARK_ONO_D2_CK5(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ARK_ONO_D2_CK4(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ARK_ONO_D2_CK3(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ARK_ONO_D2_CK2(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ARK_ONO_D2_CK1(sample, dim, amount) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_ONO_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKG)                                    :: sampleScaled(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getScaled_ARK_OTH_D2_CK5(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getScaled_ARK_OTH_D2_CK4(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getScaled_ARK_OTH_D2_CK3(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getScaled_ARK_OTH_D2_CK2(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getScaled_ARK_OTH_D2_CK1(sample, dim, amount, operation) result(sampleScaled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getScaled_ARK_OTH_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: sample(:,:)
        real(CKG)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKG)                                    :: sampleScaled(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getScaled

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is scaled by the specified input `amount` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `amount` is the square-root of the inverse of the variance of the (zero-mean) sample, then the returned sample will have a variance of one.<br>
    !>
    !>  \param[inout]   sample      :   The input/output `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the sample to be scaled.<br>
    !>  \param[in]      dim         :   The input scalar of type `integer` of default kind \IK,
    !>                                  whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]      amount      :   The input scalar or the `contiguous` vector of,
    !>                                  <ol>
    !>                                      <li>    the same type and kind as `sample`,
    !>                                      <li>    type `real` of the same kind as that of `sample` if `sample` is of type `complex`,
    !>                                  </ol>
    !>                                  representing the amount by which the input sample must be scaled.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `rank(sample) = 1`, then `amount` must be a scalar.<br>
    !>                                      <li>    If the input `rank(sample) = 2`, then `amount` must be a vector of size `ndim`.<br>
    !>                                  </ol>
    !>                                  If the sample is to be scaled to have a unit variance, then the input `amount`
    !>                                  corresponds to the square-root of the inverse of the variance of the sample that is returned by procedures
    !>                                  collectively represented by the generic interface [getVar](@ref pm_sampleVar::getVar).<br>
    !>                                  Note that the size of the input `amount` must be consistent with the size of the input `sample`:<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `dim = 1` then, `size(amount) == size(sample, 2) == ndim` must hold.<br>
    !>                                      <li>    If the input argument `dim = 2` then, `size(amount) == size(sample, 1) == ndim` must hold.<br>
    !>                                  </ol>
    !>
    !>  \interface{setScaled}
    !>  \code{.F90}
    !>
    !>      use pm_sampleScale, only: setScaled
    !>
    !>      call setScaled(sample(:), amount)
    !>
    !>      call setScaled(sample(:,:), dim, amount(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(amount) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures represented by this generic interface have the same functionality as [getScaled](@ref pm_sampleScale::getScaled).<br>
    !>  However, unlike [getScaled](@ref pm_sampleScale::getScaled) which are functions and return the scaled data as a new object,
    !>  [setScaled](@ref pm_sampleScale::setScaled) subroutines translate the input data **in place**, making it potentially more efficient and less
    !>  memory-consuming than [getScaled](@ref pm_sampleScale::getScaled) functions.<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getScaled](@ref pm_sampleScale::getScaled)<br>
    !>  [setScaled](@ref pm_sampleScale::setScaled)<br>
    !>
    !>  \example{setScaled}
    !>  \include{lineno} example/pm_sampleScale/setScaled/main.F90
    !>  \compilef{setScaled}
    !>  \output{setScaled}
    !>  \include{lineno} example/pm_sampleScale/setScaled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleScale](@ref test_pm_sampleScale)
    !>
    !>  \todo
    !>  \pmed
    !>  The functionality of this interface can be expanded to include scaling of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \final{setScaled}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 00:01 AM, August 25, 2021, Dallas, TX

    ! RK ARK

    interface setScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setScaled_ARK_D1_RK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample(:)
        real(RKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setScaled_ARK_D1_RK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: sample(:)
        real(RKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setScaled_ARK_D1_RK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: sample(:)
        real(RKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setScaled_ARK_D1_RK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: sample(:)
        real(RKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setScaled_ARK_D1_RK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: sample(:)
        real(RKG)           , intent(in)                    :: amount
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setScaled_ARK_D2_RK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(in)                    :: dim
        real(RKG)           , intent(in)    , contiguous    :: amount(:)
        real(RKG)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setScaled_ARK_D2_RK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(in)                    :: dim
        real(RKG)           , intent(in)    , contiguous    :: amount(:)
        real(RKG)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setScaled_ARK_D2_RK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(in)                    :: dim
        real(RKG)           , intent(in)    , contiguous    :: amount(:)
        real(RKG)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setScaled_ARK_D2_RK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(in)                    :: dim
        real(RKG)           , intent(in)    , contiguous    :: amount(:)
        real(RKG)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setScaled_ARK_D2_RK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(in)                    :: dim
        real(RKG)           , intent(in)    , contiguous    :: amount(:)
        real(RKG)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setScaled

    ! CK ARK

    interface setScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setScaled_ARK_D1_CK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_CK5
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        real(CKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setScaled_ARK_D1_CK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        real(CKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setScaled_ARK_D1_CK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        real(CKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setScaled_ARK_D1_CK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        real(CKG)           , intent(in)                    :: amount
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setScaled_ARK_D1_CK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        real(CKG)           , intent(in)                    :: amount
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setScaled_ARK_D2_CK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(in)                    :: dim
        real(CKG)           , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setScaled_ARK_D2_CK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(in)                    :: dim
        real(CKG)           , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setScaled_ARK_D2_CK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(in)                    :: dim
        real(CKG)           , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setScaled_ARK_D2_CK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(in)                    :: dim
        real(CKG)           , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setScaled_ARK_D2_CK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ARK_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(in)                    :: dim
        real(CKG)           , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setScaled

    ! CK ACK

    interface setScaled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setScaled_ACK_D1_CK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D1_CK5
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        complex(CKG)        , intent(in)                    :: amount
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setScaled_ACK_D1_CK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        complex(CKG)        , intent(in)                    :: amount
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setScaled_ACK_D1_CK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        complex(CKG)        , intent(in)                    :: amount
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setScaled_ACK_D1_CK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        complex(CKG)        , intent(in)                    :: amount
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setScaled_ACK_D1_CK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: sample(:)
        complex(CKG)        , intent(in)                    :: amount
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setScaled_ACK_D2_CK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(in)                    :: dim
        complex(CKG)        , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setScaled_ACK_D2_CK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(in)                    :: dim
        complex(CKG)        , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setScaled_ACK_D2_CK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(in)                    :: dim
        complex(CKG)        , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setScaled_ACK_D2_CK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(in)                    :: dim
        complex(CKG)        , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setScaled_ACK_D2_CK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaled_ACK_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(in)                    :: dim
        complex(CKG)        , intent(in)    , contiguous    :: amount(:)
        complex(CKG)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setScaled

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleScale ! LCOV_EXCL_LINE