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
!>  This module contains classes and procedures for shifting
!>  univariate or multivariate samples by arbitrary amounts along specific directions.
!>
!>  \devnote
!>  While it is tempting to add generic interfaces for automatic shifting by the negative of the mean of the sample (in the absence of arbitrary `amount` argument),
!>  such interfaces were not added to this module for the following reasons:<br>
!>  <ol>
!>      <li>    Why should shifting by the negative of the mean be the default behavior?
!>      <li>    Even though shifting by the negative of the mean is popular, its implementation as the default normalization in
!>              the generic interfaces of this module requires inclusion of sample `weight` arguments,
!>              thus significantly complicating the interfaces of this module with little gain.<br>
!>  </ol>
!>
!>  \see
!>  [pm_sampling](pm_sampling)<br>
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
!>  [test_pm_sampleShift](@ref test_pm_sampleShift)
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleShift

    use pm_kind, only: SK, IK, LK
    use pm_matrixTrans, only: transHerm, transHerm_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleShift"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type whose instances are meant to signify a sample shifting
    !>  by an amount equal to the negative of the sample mean.<br>
    !>
    !>  \details
    !>  For example usage, see the relevant interfaces that use instances of this derived type.<br>
    !>
    !>  \see
    !>  [meanshift](@ref pm_sampleShift::meanshift)
    !>
    !>  \test
    !>  [test_pm_sampleShift](@ref test_pm_sampleShift)
    !>
    !>  \finmain{meanshift_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX
    type :: meanshift_type; end type
    type(meanshift_type), parameter :: meanshift = meanshift_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: meanshift
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is shifted by the specified input `amount` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `amount` is the negative of the mean of the sample, then the returned sample will have a mean of zero.<br>
    !>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample to be shifted.<br>
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK,
    !>                              whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  amount      :   The input scalar or the `contiguous` vector of the same type and kind as `sample`,
    !>                              representing the amount by which the input sample must be shifted.<br>
    !>                              <ol>
    !>                                  <li>    If the input `rank(sample) = 1`, then `amount` must be a scalar.<br>
    !>                                  <li>    If the input `rank(sample) = 2`, then `amount` must be a vector of size `ndim`.<br>
    !>                              </ol>
    !>                              If the sample is to be shifted toward the origin to have a mean of zero, then the input `amount`
    !>                              corresponds to the negative of the current mean of the sample that is returned by procedures
    !>                              collectively represented by the generic interface [getMean](@ref pm_sampleMean::getMean).<br>
    !>                              Note that the size of the input `amount` must be consistent with the size of the input `sample`:<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(amount) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(amount) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>  \param[in]  operation   :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [transHerm](@ref pm_matrixTrans::transHerm)
    !>                                          implying that a Hermitian transpose of the input sample must be returned.<br>
    !>                                          In other words, the actions `getShifted(sample, dim, transHerm)` and `transpose(conjg(getShifted(sample, dim)))` are equivalent.<br>
    !>                                          Specifying [transHerm](@ref pm_matrixTrans::transHerm) for `source` of type other than `complex` is identical to
    !>                                          specifying [transSymm](@ref pm_matrixTrans::transSymm) for `source` of type other than `complex`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `.false.`. It **can** be present **only if** the condition `rank(sample) == 2` holds.)
    !>
    !>  \return
    !>  `sampleShifted`          :   The output object of the same type and kind and rank as `sample`, containing the shifted sample.<br>
    !>                              <ol>
    !>                                  <li>    If the input `sample` is a vector, then `sampleShifted` has the same shape and size as that of `sample`.<br>
    !>                                  <li>    If the input `sample` is a matrix of shape `(nrow, ncol)`, then
    !>                                          <ol>
    !>                                              <li>    `sampleShifted` has the shape `(nrow, ncol)` if `operation` is missing.
    !>                                              <li>    `sampleShifted` has the shape `(ncol, nrow)` if `operation = transHerm`.
    !>                                          </ol>
    !>                              </ol>
    !>
    !>  \interface{getShifted}
    !>  \code{.F90}
    !>
    !>      use pm_sampleShift, only: getShifted
    !>
    !>      sampleShifted(:) = getShifted(sample(:), amount)
    !>      sampleShifted(:,:) = getShifted(sample(:,:), dim, amount(:))
    !>      sampleShifted(:,:) = getShifted(sample(:,:), dim, amount(:), operation)
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
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>
    !>  \example{getShifted}
    !>  \include{lineno} example/pm_sampleShift/getShifted/main.F90
    !>  \compilef{getShifted}
    !>  \output{getShifted}
    !>  \include{lineno} example/pm_sampleShift/getShifted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleShift](@ref test_pm_sampleShift)
    !>
    !>  \todo
    !>  \pvlow
    !>  The functionality of this interface can be expanded to include shifting of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \finmain{getShifted}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 2:48 AM, August 22, 2021, Dallas, TX

    ! DIM

    interface getShifted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedDIM_ONO_D1_CK5(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedDIM_ONO_D1_CK4(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedDIM_ONO_D1_CK3(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedDIM_ONO_D1_CK2(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedDIM_ONO_D1_CK1(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedDIM_ONO_D1_RK5(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedDIM_ONO_D1_RK4(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedDIM_ONO_D1_RK3(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedDIM_ONO_D1_RK2(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedDIM_ONO_D1_RK1(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedDIM_ONO_D2_CK5(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedDIM_ONO_D2_CK4(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedDIM_ONO_D2_CK3(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedDIM_ONO_D2_CK2(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedDIM_ONO_D2_CK1(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedDIM_ONO_D2_RK5(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedDIM_ONO_D2_RK4(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedDIM_ONO_D2_RK3(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedDIM_ONO_D2_RK2(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedDIM_ONO_D2_RK1(sample, dim, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_ONO_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedDIM_OTH_D2_CK5(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedDIM_OTH_D2_CK4(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedDIM_OTH_D2_CK3(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedDIM_OTH_D2_CK2(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedDIM_OTH_D2_CK1(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedDIM_OTH_D2_RK5(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedDIM_OTH_D2_RK4(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedDIM_OTH_D2_RK3(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedDIM_OTH_D2_RK2(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedDIM_OTH_D2_RK1(sample, dim, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedDIM_OTH_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: amount(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getShifted

    ! ALL

    interface getShifted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedALL_ONO_D1_CK5(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedALL_ONO_D1_CK4(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedALL_ONO_D1_CK3(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedALL_ONO_D1_CK2(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedALL_ONO_D1_CK1(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedALL_ONO_D1_RK5(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedALL_ONO_D1_RK4(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedALL_ONO_D1_RK3(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedALL_ONO_D1_RK2(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedALL_ONO_D1_RK1(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedALL_ONO_D2_CK5(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedALL_ONO_D2_CK4(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedALL_ONO_D2_CK3(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedALL_ONO_D2_CK2(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedALL_ONO_D2_CK1(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        complex(CKC)                                    :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedALL_ONO_D2_RK5(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedALL_ONO_D2_RK4(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedALL_ONO_D2_RK3(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedALL_ONO_D2_RK2(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedALL_ONO_D2_RK1(sample, amount) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_ONO_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        real(RKC)                                       :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getShiftedALL_OTH_D2_CK5(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getShiftedALL_OTH_D2_CK4(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getShiftedALL_OTH_D2_CK3(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getShiftedALL_OTH_D2_CK2(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getShiftedALL_OTH_D2_CK1(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getShiftedALL_OTH_D2_RK5(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getShiftedALL_OTH_D2_RK4(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getShiftedALL_OTH_D2_RK3(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getShiftedALL_OTH_D2_RK2(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getShiftedALL_OTH_D2_RK1(sample, amount, operation) result(sampleShifted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShiftedALL_OTH_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in)                :: amount
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleShifted(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getShifted

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is shifted by the specified input `amount` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `amount` is the negative of the mean of the sample, then the returned sample will have a mean of zero.<br>
    !>
    !>  \param[inout]   sample      :   The input/output `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the sample to be shifted.<br>
    !>  \param[in]      dim         :   The input scalar of type `integer` of default kind \IK,
    !>                                  whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]      amount      :   The input scalar or the `contiguous` vector of the same type and kind as `sample`,
    !>                                  representing the amount by which the input sample must be shifted.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `rank(sample) = 1`, then `amount` must be a scalar.<br>
    !>                                      <li>    If the input `rank(sample) = 2`, then `amount` must be a vector of size `ndim`.<br>
    !>                                  </ol>
    !>                                  If the sample is to be shifted toward the origin to have a mean of zero, then the input `amount`
    !>                                  corresponds to the negative of the current mean of the sample that is returned by procedures
    !>                                  collectively represented by the generic interface [getMean](@ref pm_sampleMean::getMean).<br>
    !>                                  Note that the size of the input `amount` must be consistent with the size of the input `sample`:<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `dim = 1` then, `size(amount) == size(sample, 2) == ndim` must hold.<br>
    !>                                      <li>    If the input argument `dim = 2` then, `size(amount) == size(sample, 1) == ndim` must hold.<br>
    !>                                  </ol>
    !>
    !>  \interface{setShifted}
    !>  \code{.F90}
    !>
    !>      use pm_sampleShift, only: setShifted
    !>
    !>      call setShifted(sample(:), amount)
    !>
    !>      call setShifted(sample(:,:), dim, amount(:))
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
    !>  The procedures represented by this generic interface have the same functionality as [getShifted](@ref pm_sampleShift::getShifted).<br>
    !>  However, unlike [getShifted](@ref pm_sampleShift::getShifted) which are functions and return the shifted data as a new object,
    !>  [setShifted](@ref pm_sampleShift::setShifted) subroutines translate the input data **in place**, making it potentially more efficient and less
    !>  memory-consuming than [getShifted](@ref pm_sampleShift::getShifted) functions.<br>
    !>
    !>  \see
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>
    !>  \example{setShifted}
    !>  \include{lineno} example/pm_sampleShift/setShifted/main.F90
    !>  \compilef{setShifted}
    !>  \output{setShifted}
    !>  \include{lineno} example/pm_sampleShift/setShifted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleShift](@ref test_pm_sampleShift)
    !>
    !>  \todo
    !>  \pmed
    !>  The functionality of this interface can be expanded to include shifting of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \finmain{setShifted}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 00:01 AM, August 25, 2021, Dallas, TX

    ! DIM

    interface setShifted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setShiftedDIM_D1_CK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_CK5
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setShiftedDIM_D1_CK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setShiftedDIM_D1_CK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setShiftedDIM_D1_CK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setShiftedDIM_D1_CK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setShiftedDIM_D1_RK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_RK5
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setShiftedDIM_D1_RK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setShiftedDIM_D1_RK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setShiftedDIM_D1_RK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setShiftedDIM_D1_RK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setShiftedDIM_D2_CK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: amount(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setShiftedDIM_D2_CK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: amount(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setShiftedDIM_D2_CK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: amount(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setShiftedDIM_D2_CK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: amount(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setShiftedDIM_D2_CK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: amount(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setShiftedDIM_D2_RK5(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: amount(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setShiftedDIM_D2_RK4(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: amount(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setShiftedDIM_D2_RK3(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: amount(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setShiftedDIM_D2_RK2(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: amount(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setShiftedDIM_D2_RK1(sample, dim, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedDIM_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: amount(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setShifted

    ! ALL

    interface setShifted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setShiftedALL_D1_CK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_CK5
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setShiftedALL_D1_CK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setShiftedALL_D1_CK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setShiftedALL_D1_CK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setShiftedALL_D1_CK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setShiftedALL_D1_RK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_RK5
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setShiftedALL_D1_RK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setShiftedALL_D1_RK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setShiftedALL_D1_RK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setShiftedALL_D1_RK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setShiftedALL_D2_CK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setShiftedALL_D2_CK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setShiftedALL_D2_CK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setShiftedALL_D2_CK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setShiftedALL_D2_CK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: amount
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setShiftedALL_D2_RK5(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setShiftedALL_D2_RK4(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setShiftedALL_D2_RK3(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setShiftedALL_D2_RK2(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setShiftedALL_D2_RK1(sample, amount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setShiftedALL_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: amount
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setShifted

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleShift ! LCOV_EXCL_LINE