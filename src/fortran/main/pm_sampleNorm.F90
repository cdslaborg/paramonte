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
!>  This module contains classes and procedures for normalizing
!>  univariate or multivariate samples by arbitrary amounts along specific directions.
!>
!>  \details
!>  Normalization can have a wide variety of meanings in science.<br>
!>  In this module, it refers to the creation of a **shifted** and **scaled** version of a sample,
!>  where the intention is that these normalized values allow the comparison of corresponding normalized values
!>  for different datasets in a way that eliminates the effects of certain gross influences.<br>
!>
!>  The procedures of this module facilitate the computation of the following popular sample normalizations, among others:<br>
!>
!>  Standard score
!>  --------------
!>
!>  The standard score (**z-score**) is the number of standard deviations by which the value of a raw score (i.e., an observed value or data point)
!>  is above or below the mean value of what is being observed or measured.<br>
!>  Raw scores above the mean have positive standard scores, while those below the mean have negative standard scores.<br>
!>
!>  If the population mean and population standard deviation are known, a raw score x is converted into a standard score by,
!>  \f{equation}{
!>      z = {x - \mu \over \sigma} ~,
!>  \f}
!>  where:<br>
!>  <ol>
!>      <li>    \f$\mu\f$ is the mean of the population,
!>      <li>    \f$\sigma\f$ is the standard deviation of the population.
!>  </ol>
!>  When the population mean and the population standard deviation are unknown,
!>  the standard score may be estimated by using the sample mean and sample standard deviation as estimates of the population values.<br>
!>  In these cases, the z-score is given by,
!>  \f{equation}{
!>      z = {x - {\hat\mu} \over \hat\sigma} ~,
!>  \f}
!>  where:<br>
!>  <ol>
!>      <li>    \f$\hat\mu\f$ is the mean of the sample,
!>      <li>    \f$\hat\sigma\f$ is the standard deviation of the sample.
!>  </ol>
!>
!>  Rescaling (min-max normalization)
!>  ---------------------------------
!>
!>  Also known as **min-max scaling** or **min-max normalization**, it consists of rescaling the range of features to scale the range in \f$[0, 1]\f$ or \f$[−1, 1]\f$.<br>
!>  Selecting the target range depends on the nature of the data.<br>
!>  The general formula for a min-max of \f$[0, 1]\f$ is given as:<br>
!>  \f{equation}{
!>      \tilde x = \frac {x - {\text{min}}(x)}{{\text{max}}(x)-{\text{min}}(x)} ~,
!>  \f}
!>  where \f$x\f$ is an original value and \f$\tilde x\f$ is the normalized value.<br>
!>  For example, suppose that we have the students weight data, and the students weights span [160 pounds, 200 pounds].<br>
!>  To rescale this data, we first subtract \f$160\f$ from each student weight and divide the result by \f$40\f$ (the difference between the maximum and minimum weights).<br>
!>
!>  To rescale a range between an arbitrary set of values \f$[a, b]\f$, the formula becomes:<br>
!>  \f{equation}{
!>      \tilde x = a + {\frac {(x-{\text{min}}(x))(b-a)}{{\text{max}}(x)-{\text{min}}(x)}} ~,
!>  \f}
!>  where \f$a, b\f$ are the min-max values.<br>
!>
!>  Mean normalization
!>  ------------------
!>
!>  \f{equation}{
!>      \tilde x = {\frac {x-{\bar {x}}}{{\text{max}}(x)-{\text{min}}(x)}} ~,
!>  \f}
!>  where \f$x\f$ is an original value and \f$\tilde x\f$ is the normalized value and \f${\bar{x}} = {\text{average}}(x)\f$ is the mean of that feature vector.<br>
!>  There is another form of the means normalization which divides by the standard deviation which is also called standardization.<br>
!>
!>  \devnote
!>  While it is tempting to add generic interfaces for automatic **standard normalization** of the sample (in the absence of arbitrary `shift and `scale` arguments),
!>  such interfaces were not added to this module for the following reasons:<br>
!>  <ol>
!>      <li>    Why should standard normalization be the default behavior?
!>      <li>    Even though standard normalization is popular, its implementation as the default normalization in
!>              the generic interfaces of this module requires inclusion of sample `weight` and variance `correction` arguments,
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
!>  [Normalization](https://en.wikipedia.org/wiki/Normalization_(statistics))<br>
!>
!>  \test
!>  [test_pm_sampleNorm](@ref test_pm_sampleNorm)
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleNorm

    use pm_kind, only: SK, IK, LK
    use pm_matrixTrans, only: transHerm, transHerm_type
    use pm_sampleShift, only: meanshift, meanshift_type
    use pm_sampleScale, only: stdscale, stdscale_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleNorm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type whose instances are meant to signify a sample shifting by an amount equal to the negative of the sample mean
    !>  and scaling the result by an amount equal to the inverse of the sample standard deviation or an equivalent measure.<br>
    !>
    !>  \details
    !>  For example usage, see the relevant interfaces that use instances of this derived type.<br>
    !>
    !>  \see
    !>  [zscore](@ref pm_sampleNorm::zscore)<br>
    !>
    !>  \test
    !>  [test_pm_sampleNorm](@ref test_pm_sampleNorm)
    !>
    !>  \finmain{zscore_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX
    type :: zscore_type; end type
    type(zscore_type), parameter :: zscore = zscore_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is normalized by the specified input `shift` and `scale` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `shift` is the negative of the mean of the sample and `scale` is the square-root of the inverse of the variance of the sample,
    !>  then the returned sample will be a standard score (i.e., z-score).<br>
    !>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample to be normalized.<br>
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK,
    !>                              whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  shift       :   The input scalar or the `contiguous` vector of the same type and kind as `sample`,
    !>                              representing the amount by which the input sample must be shifted.<br>
    !>                              <ol>
    !>                                  <li>    If the input `rank(sample) = 1`, then `shift` must be a scalar.<br>
    !>                                  <li>    If the input `rank(sample) = 2`, then `shift` must be a vector of size `ndim`.<br>
    !>                              </ol>
    !>                              If the sample is to be shifted toward the origin to have a mean of zero, then the input `shift`
    !>                              corresponds to the negative of the current mean of the sample that is returned by procedures
    !>                              collectively represented by the generic interface [getMean](@ref pm_sampleMean::getMean).<br>
    !>                              Note that the size of the input `shift` must be consistent with the size of the input `sample`:<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(shift) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(shift) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>  \param[in]  scale       :   The input scalar or the `contiguous` vector of,
    !>                              <ol>
    !>                                  <li>    the same type and kind as `sample`,
    !>                                  <li>    type `real` of the same kind as that of `sample` if `sample` is of type `complex`,
    !>                              </ol>
    !>                              representing the amount by which the input sample must be scaled.<br>
    !>                              <ol>
    !>                                  <li>    If the input `rank(sample) = 1`, then `scale` must be a scalar.<br>
    !>                                  <li>    If the input `rank(sample) = 2`, then `scale` must be a vector of size `ndim`.<br>
    !>                              </ol>
    !>                              If the sample is to be scaled to have a unit variance, then the input `scale`
    !>                              corresponds to the square-root of the inverse of the variance of the sample that is returned by procedures
    !>                              collectively represented by the generic interface [getVar](@ref pm_sampleVar::getVar).<br>
    !>                              Note that the size of the input `scale` must be consistent with the size of the input `sample`:<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(scale) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(scale) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>  \param[in]  operation   :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [transHerm](@ref pm_matrixTrans::transHerm)
    !>                                          implying that a Hermitian transpose of the input sample must be returned.<br>
    !>                                          In other words, the actions `getNormed(sample, dim, transHerm)` and `transpose(conjg(getNormed(sample, dim)))` are equivalent.<br>
    !>                                          Specifying [transHerm](@ref pm_matrixTrans::transHerm) for `source` of type other than `complex` is identical to
    !>                                          specifying [transSymm](@ref pm_matrixTrans::transSymm) for `source` of type other than `complex`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `.false.`. It **can** be present **only if** the condition `rank(sample) == 2` holds.)
    !>
    !>  \return
    !>  `sampleNormed`          :   The output object of the same type and kind and rank as `sample`, containing the normalized sample.<br>
    !>                              <ol>
    !>                                  <li>    If the input `sample` is a vector, then `sampleNormed` has the same shape and size as that of `sample`.<br>
    !>                                  <li>    If the input `sample` is a matrix of shape `(nrow, ncol)`, then
    !>                                          <ol>
    !>                                              <li>    `sampleNormed` has the shape `(nrow, ncol)` if `operation` is missing.
    !>                                              <li>    `sampleNormed` has the shape `(ncol, nrow)` if `operation = transHerm`.
    !>                                          </ol>
    !>                              </ol>
    !>
    !>  \interface{getNormed}
    !>  \code{.F90}
    !>
    !>      use pm_sampleNorm, only: getNormed
    !>
    !>      sampleNormed(:) = getNormed(sample(:), shift, scale)
    !>      sampleNormed(:,:) = getNormed(sample(:,:), dim, shift(:), scale(:))
    !>      sampleNormed(:,:) = getNormed(sample(:,:), dim, shift(:), scale(:), operation)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(scale) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  The condition `size(shift) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getScaled](@ref pm_sampleScale::getScaled)<br>
    !>  [setScaled](@ref pm_sampleScale::setScaled)<br>
    !>  [getNormed](@ref pm_sampleNorm::getNormed)<br>
    !>  [setNormed](@ref pm_sampleNorm::setNormed)<br>
    !>
    !>  \example{getNormed}
    !>  \include{lineno} example/pm_sampleNorm/getNormed/main.F90
    !>  \compilef{getNormed}
    !>  \output{getNormed}
    !>  \include{lineno} example/pm_sampleNorm/getNormed/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleNorm](@ref test_pm_sampleNorm)
    !>
    !>  \todo
    !>  \pvlow
    !>  The functionality of this interface can be expanded to include normalizing of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \finmain{getNormed}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 2:48 AM, August 22, 2021, Dallas, TX

    ! CK ACK

    interface getNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ACK_ONO_D1_CK5(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        complex(CKC)        , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ACK_ONO_D1_CK4(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        complex(CKC)        , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ACK_ONO_D1_CK3(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        complex(CKC)        , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ACK_ONO_D1_CK2(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        complex(CKC)        , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ACK_ONO_D1_CK1(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        complex(CKC)        , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ACK_ONO_D2_CK5(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ACK_ONO_D2_CK4(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ACK_ONO_D2_CK3(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ACK_ONO_D2_CK2(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ACK_ONO_D2_CK1(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_ONO_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ACK_OTH_D2_CK5(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_OTH_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ACK_OTH_D2_CK4(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_OTH_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ACK_OTH_D2_CK3(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_OTH_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ACK_OTH_D2_CK2(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_OTH_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ACK_OTH_D2_CK1(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ACK_OTH_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        complex(CKC)        , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormed

    ! CK ARK

    interface getNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ARK_ONO_D1_CK5(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        real(CKC)           , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ARK_ONO_D1_CK4(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        real(CKC)           , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ARK_ONO_D1_CK3(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        real(CKC)           , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ARK_ONO_D1_CK2(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        real(CKC)           , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ARK_ONO_D1_CK1(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:)
        complex(CKC)        , intent(in)                :: shift
        real(CKC)           , intent(in)                :: scale
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ARK_ONO_D2_CK5(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ARK_ONO_D2_CK4(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ARK_ONO_D2_CK3(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ARK_ONO_D2_CK2(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ARK_ONO_D2_CK1(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        complex(CKC)                                    :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getNormed_ARK_OTH_D2_CK5(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getNormed_ARK_OTH_D2_CK4(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getNormed_ARK_OTH_D2_CK3(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getNormed_ARK_OTH_D2_CK2(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getNormed_ARK_OTH_D2_CK1(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in), contiguous    :: sample(:,:)
        complex(CKC)        , intent(in), contiguous    :: shift(:)
        real(CKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        complex(CKC)                                    :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormed

    ! RK ARK

    interface getNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormed_ARK_ONO_D1_RK5(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: shift
        real(RKC)           , intent(in)                :: scale
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormed_ARK_ONO_D1_RK4(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: shift
        real(RKC)           , intent(in)                :: scale
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormed_ARK_ONO_D1_RK3(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: shift
        real(RKC)           , intent(in)                :: scale
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormed_ARK_ONO_D1_RK2(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: shift
        real(RKC)           , intent(in)                :: scale
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormed_ARK_ONO_D1_RK1(sample, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:)
        real(RKC)           , intent(in)                :: shift
        real(RKC)           , intent(in)                :: scale
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormed_ARK_ONO_D2_RK5(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormed_ARK_ONO_D2_RK4(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormed_ARK_ONO_D2_RK3(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormed_ARK_ONO_D2_RK2(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormed_ARK_ONO_D2_RK1(sample, dim, shift, scale) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_ONO_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        real(RKC)                                       :: sampleNormed(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormed_ARK_OTH_D2_RK5(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormed_ARK_OTH_D2_RK4(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormed_ARK_OTH_D2_RK3(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormed_ARK_OTH_D2_RK2(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormed_ARK_OTH_D2_RK1(sample, dim, shift, scale, operation) result(sampleNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormed_ARK_OTH_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in), contiguous    :: sample(:,:)
        real(RKC)           , intent(in), contiguous    :: shift(:)
        real(RKC)           , intent(in), contiguous    :: scale(:)
        integer(IK)         , intent(in)                :: dim
        type(transHerm_type), intent(in)                :: operation
        real(RKC)                                       :: sampleNormed(size(sample, 2, IK), size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getNormed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a sample of shape `(nsam)`, or `(ndim, nsam)` or `(nsam, ndim)` that is normalized by the specified input `shift` and `scale` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `shift` is the negative of the mean of the sample and `scale` is the square-root of the inverse of the variance of the sample,
    !>  then the returned sample will be a standard score (i.e., z-score).<br>
    !>
    !>  \param[inout]   sample      :   The input/output `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the sample to be normalized.<br>
    !>  \param[in]      dim         :   The input scalar of type `integer` of default kind \IK,
    !>                                  whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]      shift       :   The input scalar or the `contiguous` vector of the same type and kind as `sample`,
    !>                                  representing the amount by which the input sample must be shifted.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `rank(sample) = 1`, then `shift` must be a scalar.<br>
    !>                                      <li>    If the input `rank(sample) = 2`, then `shift` must be a vector of size `ndim`.<br>
    !>                                  </ol>
    !>                                  If the sample is to be shifted toward the origin to have a mean of zero, then the input `shift`
    !>                                  corresponds to the negative of the current mean of the sample that is returned by procedures
    !>                                  collectively represented by the generic interface [getMean](@ref pm_sampleMean::getMean).<br>
    !>                                  Note that the size of the input `shift` must be consistent with the size of the input `sample`:<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `dim = 1` then, `size(shift) == size(sample, 2) == ndim` must hold.<br>
    !>                                      <li>    If the input argument `dim = 2` then, `size(shift) == size(sample, 1) == ndim` must hold.<br>
    !>                                  </ol>
    !>  \param[in]      scale       :   The input scalar or the `contiguous` vector of,
    !>                                  <ol>
    !>                                      <li>    the same type and kind as `sample`,
    !>                                      <li>    type `real` of the same kind as that of `sample` if `sample` is of type `complex`,
    !>                                  </ol>
    !>                                  representing the amount by which the input sample must be scaled.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `rank(sample) = 1`, then `scale` must be a scalar.<br>
    !>                                      <li>    If the input `rank(sample) = 2`, then `scale` must be a vector of size `ndim`.<br>
    !>                                  </ol>
    !>                                  If the sample is to be scaled to have a unit variance, then the input `scale`
    !>                                  corresponds to the square-root of the inverse of the variance of the sample that is returned by procedures
    !>                                  collectively represented by the generic interface [getVar](@ref pm_sampleVar::getVar).<br>
    !>                                  Note that the size of the input `scale` must be consistent with the size of the input `sample`:<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `dim = 1` then, `size(scale) == size(sample, 2) == ndim` must hold.<br>
    !>                                      <li>    If the input argument `dim = 2` then, `size(scale) == size(sample, 1) == ndim` must hold.<br>
    !>                                  </ol>
    !>
    !>  \interface{setNormed}
    !>  \code{.F90}
    !>
    !>      use pm_sampleNorm, only: setNormed
    !>
    !>      call setNormed(sample(:), shift, scale)
    !>
    !>      call setNormed(sample(:,:), dim, shift(:), scale(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(scale) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  The condition `size(shift) == size(sample, 3 - dim) .or. rank(sample) /= 2` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures represented by this generic interface have the same functionality as [getNormed](@ref pm_sampleNorm::getNormed).<br>
    !>  However, unlike [getNormed](@ref pm_sampleNorm::getNormed) which are functions and return the normalized data as a new object,
    !>  [setNormed](@ref pm_sampleNorm::setNormed) subroutines translate the input data **in place**, making it potentially more efficient and less
    !>  memory-consuming than [getNormed](@ref pm_sampleNorm::getNormed) functions.<br>
    !>
    !>  \see
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getNormed](@ref pm_sampleNorm::getNormed)<br>
    !>  [setNormed](@ref pm_sampleNorm::setNormed)<br>
    !>
    !>  \example{setNormed}
    !>  \include{lineno} example/pm_sampleNorm/setNormed/main.F90
    !>  \compilef{setNormed}
    !>  \output{setNormed}
    !>  \include{lineno} example/pm_sampleNorm/setNormed/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleNorm](@ref test_pm_sampleNorm)
    !>
    !>  \todo
    !>  \pmed
    !>  The functionality of this interface can be expanded to include normalizing of higher dimensional input `sample`
    !>  and whole `sample` input arrays of arbitrary shape, although the latter is trivial using the Fortran array syntax.<br>
    !>
    !>  \finmain{setNormed}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 00:01 AM, August 25, 2021, Dallas, TX

    ! CK ACK

    interface setNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setNormed_ACK_D1_CK5(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D1_CK5
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        complex(CKC)        , intent(in)                    :: scale
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setNormed_ACK_D1_CK4(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        complex(CKC)        , intent(in)                    :: scale
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setNormed_ACK_D1_CK3(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        complex(CKC)        , intent(in)                    :: scale
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setNormed_ACK_D1_CK2(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        complex(CKC)        , intent(in)                    :: scale
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setNormed_ACK_D1_CK1(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        complex(CKC)        , intent(in)                    :: scale
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setNormed_ACK_D2_CK5(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setNormed_ACK_D2_CK4(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setNormed_ACK_D2_CK3(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setNormed_ACK_D2_CK2(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setNormed_ACK_D2_CK1(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ACK_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)         , intent(in)                    :: dim
        complex(CKC)        , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormed

    ! CK ARK

    interface setNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setNormed_ARK_D1_CK5(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_CK5
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        real(CKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setNormed_ARK_D1_CK4(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        real(CKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setNormed_ARK_D1_CK3(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        real(CKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setNormed_ARK_D1_CK2(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        real(CKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setNormed_ARK_D1_CK1(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(inout) , contiguous    :: sample(:)
        complex(CKC)        , intent(in)                    :: shift
        real(CKC)           , intent(in)                    :: scale
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setNormed_ARK_D2_CK5(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)         , intent(in)                    :: dim
        real(CKC)           , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setNormed_ARK_D2_CK4(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)         , intent(in)                    :: dim
        real(CKC)           , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setNormed_ARK_D2_CK3(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)         , intent(in)                    :: dim
        real(CKC)           , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setNormed_ARK_D2_CK2(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)         , intent(in)                    :: dim
        real(CKC)           , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setNormed_ARK_D2_CK1(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)         , intent(in)                    :: dim
        real(CKC)           , intent(in)    , contiguous    :: scale(:)
        complex(CKC)        , intent(in)    , contiguous    :: shift(:)
        complex(CKC)        , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormed

    ! RK ARK

    interface setNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormed_ARK_D1_RK5(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_RK5
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
        real(RKC)           , intent(in)                    :: shift
        real(RKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormed_ARK_D1_RK4(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
        real(RKC)           , intent(in)                    :: shift
        real(RKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormed_ARK_D1_RK3(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
        real(RKC)           , intent(in)                    :: shift
        real(RKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormed_ARK_D1_RK2(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
        real(RKC)           , intent(in)                    :: shift
        real(RKC)           , intent(in)                    :: scale
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormed_ARK_D1_RK1(sample, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(inout) , contiguous    :: sample(:)
        real(RKC)           , intent(in)                    :: shift
        real(RKC)           , intent(in)                    :: scale
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setNormed_ARK_D2_RK5(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: shift(:)
        real(RKC)           , intent(in)    , contiguous    :: scale(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setNormed_ARK_D2_RK4(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: shift(:)
        real(RKC)           , intent(in)    , contiguous    :: scale(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setNormed_ARK_D2_RK3(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: shift(:)
        real(RKC)           , intent(in)    , contiguous    :: scale(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setNormed_ARK_D2_RK2(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: shift(:)
        real(RKC)           , intent(in)    , contiguous    :: scale(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setNormed_ARK_D2_RK1(sample, dim, shift, scale)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNormed_ARK_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)         , intent(in)                    :: dim
        real(RKC)           , intent(in)    , contiguous    :: shift(:)
        real(RKC)           , intent(in)    , contiguous    :: scale(:)
        real(RKC)           , intent(inout) , contiguous    :: sample(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setNormed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleNorm ! LCOV_EXCL_LINE