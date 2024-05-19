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
!>  This module contains classes and procedures for affine transformation of multivariate samples.
!>
!>  \details
!>
!>  In Euclidean geometry, an **affine transformation** or **affinity** (from the Latin, **affinis**, *connected with*)
!>  is a geometric transformation that preserves **lines** and **parallelism**, but not necessarily Euclidean distances and angles.<br>
!>
!>  Definition
!>  ----------
!>
!>  If \f$X\f$ is the point set of an affine space, then every affine transformation on \f$X\f$
!>  can be represented as the composition of a linear transformation on \f$X\f$ and a translation of \f$X\f$.<br>
!>
!>  An affine map is the composition of two functions: a translation and a linear map.<br>
!>  Ordinary vector algebra uses matrix multiplication to represent linear maps, and vector addition to represent translations.<br>
!>  Formally, in the finite-dimensional case, if the linear map is represented as a multiplication by an invertible matrix \f$A\f$ and
!>  the translation as the addition of a vector \f$\mathbf{b}\f$, an affine map \f$f\f$ acting on a vector \f$\mathbf{x}\f$ can be represented as,
!>  \f{equation}{
!>      \mathbf{y} = f(\mathbf{x}) = A\mathbf{x} + \mathbf{b} ~.
!>  \f}
!>
!>  Augmented matrix representation
!>  -------------------------------
!>
!>  Using an augmented matrix and an augmented vector, it is possible to represent both the translation and the linear map using a single matrix multiplication.<br>
!>  The technique requires that all vectors be augmented with a **1** at the end, and all matrices be augmented with an extra row of zeros at the bottom,
!>  an extra column—the translation vector—to the right, and a **1** in the lower right corner.<br>
!>  If \f$A\f$ is a matrix,
!>  \f{equation}
!>      \begin{bmatrix}
!>          \mathbf{y} \\
!>          1
!>      \end{bmatrix}
!>      =
!>      \left[
!>          \begin{array}{ccc|c}
!>              &A&&\mathbf{b} \\
!>              0&\cdots&0&1
!>          \end{array}
!>      \right]
!>      \begin{bmatrix}
!>          \mathbf{x} \\
!>          1
!>      \end{bmatrix}
!>  \f}
!>  is equivalent to the following
!>  \f{equation}
!>      \mathbf{y} = A\mathbf{x} + \mathbf{b} ~.
!>  \f}
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
!>  [affine transformation](https://en.wikipedia.org/wiki/Affine_transformation)<br>
!>
!>  \test
!>  [test_pm_sampleAffinity](@ref test_pm_sampleAffinity)
!>
!>  \final
!>
!>  \author
!>  Amir Shahmoradi, Monday 4:00 AM, August 23, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleAffinity

    use pm_kind, only: SK, IK, LK
    use pm_matrixClass, only: genrecmat, genrecmat_type
    use pm_matrixClass, only: upperDiag, upperDiag_type
    use pm_matrixClass, only: lowerDiag, lowerDiag_type
    use pm_matrixClass, only: upperUnit, upperUnit_type
    use pm_matrixClass, only: lowerUnit, lowerUnit_type
    use pm_matrixClass, only: upperZero, upperZero_type
    use pm_matrixClass, only: lowerZero, lowerZero_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleAffinity"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an affine-transformation of the input `sample` of shape `(1:ndim)` or `(1:ndim, 1:nsam)` or `(1:nsam, 1:ndim)`
    !>  based on the specified values for the translation `tlate` and tranformations `tform` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `tlate` is the negative of the mean of the sample, then the returned sample will have a mean of zero.<br>
    !>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample to be affine-transformed.<br>
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK,
    !>                              whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(1:nsam, 1:ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(1:ndim, 1:nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  tform       :   The input `contiguous` matrix of the same type and kind as `sample` of shape `(1:ndim, 1:ndim)`,
    !>                              containing the transformation matrix of the affine-transformation.<br>
    !>  \param[in]  class       :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The scalar constant [genrecmat](@ref pm_matrixClass::genrecmat) or equivalently an object of type [genrecmat_type](@ref pm_matrixClass::genrecmat_type)
    !>                                          signifying the **general rectangular matrix** as the class of the input transformation matrix.<br>
    !>                                          All elements of the input matrix will be used in constructing the output sample.<br>
    !>                                  <li>    The scalar constant [upperDiag](@ref pm_matrixClass::upperDiag) or equivalently an object of type [upperDiag_type](@ref pm_matrixClass::upperDiag_type)
    !>                                          signifying the **upper-diagonal matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the upper-diagonal elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The lower triangle elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [lowerDiag](@ref pm_matrixClass::lowerDiag) or equivalently an object of type [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)
    !>                                          signifying the **lower-diagonal matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the lower-diagonal elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The upper triangle elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [upperUnit](@ref pm_matrixClass::upperUnit) or equivalently an object of type [upperUnit_type](@ref pm_matrixClass::upperUnit_type)
    !>                                          signifying the **unit-triangular matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the upper elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The lower triangle elements of `tform` are assumed to be zero.<br>
    !>                                          The diagonal elements of `tform` are assumed to be one.<br>
    !>                                          The diagonal elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [lowerUnit](@ref pm_matrixClass::lowerUnit) or equivalently an object of type [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)
    !>                                          signifying the **unit-triangular matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the lower elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The upper triangle elements of `tform` are assumed to be zero.<br>
    !>                                          The diagonal elements of `tform` are assumed to be zero.<br>
    !>                              </ol>
    !>                              (**optional**, default = [genrecmat](@ref pm_matrixClass::genrecmat).)
    !>  \param[in]  tlate       :   The input `contiguous` vector of the same type and kind as `sample` of shape `(1:ndim)`,
    !>                              containing the amount by which the input sample must be translated (i.e., shifted).<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(tlate) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(tlate) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>                              (**optional**, default = `[(0, i = 1, ndim)]`.)
    !>
    !>  \return
    !>  `affinity`              :   The output object of the same type and kind and rank and shape as `sample`, containing the affine-transformed sample.<br>
    !>
    !>  \interface{getAffinity}
    !>  \code{.F90}
    !>
    !>      use pm_sampleAffinity, only: getAffinity
    !>
    !>      affinity(1:ndim) = getAffinity(sample(1:ndim), tform(1:ndim, 1:ndim), class = class, tlate = tlate(1:ndim))
    !>      affinity(@shape(sample)) = getAffinity(sample(:,:), dim, tform(1:size(sample, 3 - dim), 1:size(sample, 3 - dim)), class = class, tlate = tlate(1:size(sample, 3 - dim)))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(tlate) == size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(tform) == size(sample, 3 - dim))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [setAffinity](@ref pm_sampleAffinity::setAffinity)<br>
    !>  [setAffinity](@ref pm_sampleAffinity::setAffinity)<br>
    !>
    !>  \example{getAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/getAffinity/main.F90
    !>  \compilef{getAffinity}
    !>  \output{getAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/getAffinity/main.out.F90
    !>  \postproc{getAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/getAffinity/main.py
    !>  \vis{getAffinity}
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.circle.sample.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.circle.affinity.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.circle.affinInv.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.square.sample.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.square.affinity.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.warp.square.affinInv.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.circle.sample.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.circle.affinity.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.circle.affinInv.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.square.sample.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.square.affinity.png width=700
    !>  \image html pm_sampleAffinity/getAffinity/getAffinity.rotation.square.affinInv.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleAffinity](@ref test_pm_sampleAffinity)
    !>
    !>  \final{getAffinity}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 2:48 AM, August 22, 2021, Dallas, TX

    interface getAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getAffinity_D1_RK5(sample, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: sample(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getAffinity_D1_RK4(sample, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: sample(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getAffinity_D1_RK3(sample, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: sample(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getAffinity_D1_RK2(sample, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: sample(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getAffinity_D1_RK1(sample, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: sample(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getAffinity_D2_RK5(sample, dim, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(in)                            :: dim
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)           , intent(in), contiguous                :: sample(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getAffinity_D2_RK4(sample, dim, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(in)                            :: dim
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)           , intent(in), contiguous                :: sample(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getAffinity_D2_RK3(sample, dim, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(in)                            :: dim
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)           , intent(in), contiguous                :: sample(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getAffinity_D2_RK2(sample, dim, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(in)                            :: dim
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)           , intent(in), contiguous                :: sample(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getAffinity_D2_RK1(sample, dim, tform, class, tlate) result(affinity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAffinity_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(in)                            :: dim
        class(*)            , intent(in)            , optional      :: class
        real(RKG)           , intent(in), contiguous, optional      :: tlate(:)
        real(RKG)           , intent(in), contiguous                :: tform(:,:)
        real(RKG)           , intent(in), contiguous                :: sample(:,:)
        real(RKG)                                                   :: affinity(size(sample, 1, IK), size(sample, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getAffinity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return an affine-transformation of the input `sample` of shape `(1:ndim)` or `(1:ndim, 1:nsam)` or `(1:nsam, 1:ndim)`
    !>  based on the specified values for the translation `tlate` and tranformations `tform` along the specified axis `dim`.<br>
    !>
    !>  \details
    !>  Here, `ndim` stands for the number of dimensions (data attributes) of the input `sample` and `nsam` represents the number of data points in the `sample`.<br>
    !>  If the input `tlate` is the negative of the mean of the sample, then the returned sample will have a mean of zero.<br>
    !>
    !>  \param[out] affinity    :   The output object of the same type and kind and rank and shape as `sample`, containing the affine-transformed sample.<br>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample to be affine-transformed.<br>
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK,
    !>                              whose value represents the dimension of the input `sample` containing different `nsam` observations:<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(1:nsam, 1:ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(1:ndim, 1:nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input arguments the condition `rank(sample) > 1` holds.)
    !>  \param[in]  tform       :   The input `contiguous` matrix of the same type and kind as `sample` of shape `(1:ndim, 1:ndim)`,
    !>                              containing the transformation matrix of the affine-transformation.<br>
    !>  \param[in]  class       :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The scalar constant [genrecmat](@ref pm_matrixClass::genrecmat) or equivalently an object of type [genrecmat_type](@ref pm_matrixClass::genrecmat_type)
    !>                                          signifying the **general rectangular matrix** as the class of the input transformation matrix.<br>
    !>                                          All elements of the input matrix will be used in constructing the output sample.<br>
    !>                                  <li>    The scalar constant [upperDiag](@ref pm_matrixClass::upperDiag) or equivalently an object of type [upperDiag_type](@ref pm_matrixClass::upperDiag_type)
    !>                                          signifying the **upper-diagonal matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the upper-diagonal elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The lower triangle elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [lowerDiag](@ref pm_matrixClass::lowerDiag) or equivalently an object of type [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)
    !>                                          signifying the **lower-diagonal matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the lower-diagonal elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The upper triangle elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [upperUnit](@ref pm_matrixClass::upperUnit) or equivalently an object of type [upperUnit_type](@ref pm_matrixClass::upperUnit_type)
    !>                                          signifying the **unit-triangular matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the upper elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The lower triangle elements of `tform` are assumed to be zero.<br>
    !>                                          The diagonal elements of `tform` are assumed to be one.<br>
    !>                                          The diagonal elements of `tform` are assumed to be zero.<br>
    !>                                  <li>    The scalar constant [lowerUnit](@ref pm_matrixClass::lowerUnit) or equivalently an object of type [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)
    !>                                          signifying the **unit-triangular matrix** as the class of the input transformation matrix.<br>
    !>                                          Only the lower elements of `tform` will be used in constructing the output sample.<br>
    !>                                          The upper triangle elements of `tform` are assumed to be zero.<br>
    !>                                          The diagonal elements of `tform` are assumed to be zero.<br>
    !>                              </ol>
    !>  \param[in]  tlate       :   The input `contiguous` vector of the same type and kind as `sample` of shape `(1:ndim)`,
    !>                              containing the amount by which the input sample must be translated (i.e., shifted).<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `dim = 1` then, `size(tlate) == size(sample, 2) == ndim` must hold.<br>
    !>                                  <li>    If the input argument `dim = 2` then, `size(tlate) == size(sample, 1) == ndim` must hold.<br>
    !>                              </ol>
    !>                              (**optional**, default = `[(0, i = 1, ndim)]`.)
    !>
    !>  \interface{setAffinity}
    !>  \code{.F90}
    !>
    !>      use pm_sampleAffinity, only: setAffinity
    !>
    !>      call setAffinity(affinity(1:ndim), sample(1:ndim), tform(1:ndim, 1:ndim), class)
    !>      call setAffinity(affinity(1:ndim), sample(1:ndim), tform(1:ndim, 1:ndim), class, tlate(1:ndim))
    !>
    !>      call setAffinity(affinity(@shape(sample)), sample(:,:), dim, tform(1:size(sample, 3 - dim), 1:size(sample, 3 - dim)), class)
    !>      call setAffinity(affinity(@shape(sample)), sample(:,:), dim, tform(1:size(sample, 3 - dim), 1:size(sample, 3 - dim)), class, tlate(1:size(sample, 3 - dim)))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(tlate) == size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(tform) == size(sample, 3 - dim))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [setAffinity](@ref pm_sampleAffinity::setAffinity)<br>
    !>  [setAffinity](@ref pm_sampleAffinity::setAffinity)<br>
    !>
    !>  \example{setAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/setAffinity/main.F90
    !>  \compilef{setAffinity}
    !>  \output{setAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/setAffinity/main.out.F90
    !>  \postproc{setAffinity}
    !>  \include{lineno} example/pm_sampleAffinity/setAffinity/main.py
    !>  \vis{setAffinity}
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.circle.sample.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.circle.affinity.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.circle.affinInv.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.square.sample.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.square.affinity.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.warp.square.affinInv.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.circle.sample.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.circle.affinity.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.circle.affinInv.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.square.sample.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.square.affinity.png width=700
    !>  \image html pm_sampleAffinity/setAffinity/setAffinity.rotation.square.affinInv.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleAffinity](@ref test_pm_sampleAffinity)
    !>
    !>  \todo
    !>  \phigh
    !>  The performance of the algorithm for the case `dim = 1` can be
    !>  further improved by looping over the sample in the innermost layer.<br>
    !>
    !>  \final{setAffinity}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 2:45 AM, August 19, 2021, Dallas, TX<br>
    !>  \FatemehBagheri, Wednesday 00:01 AM, August 25, 2021, Dallas, TX<br>

    ! CGR_ATL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D1_RK5(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D1_RK4(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D1_RK3(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D1_RK2(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D1_RK1(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D2_RK5(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D2_RK4(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D2_RK3(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D2_RK2(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CGR_ATL_D2_RK1(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_ATL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CUD_ATL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D1_RK5(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D1_RK4(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D1_RK3(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D1_RK2(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D1_RK1(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D2_RK5(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D2_RK4(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D2_RK3(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D2_RK2(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUD_ATL_D2_RK1(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_ATL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CLD_ATL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D1_RK5(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D1_RK4(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D1_RK3(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D1_RK2(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D1_RK1(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D2_RK5(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D2_RK4(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D2_RK3(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D2_RK2(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLD_ATL_D2_RK1(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_ATL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CUU_ATL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D1_RK5(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D1_RK4(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D1_RK3(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D1_RK2(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D1_RK1(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D2_RK5(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D2_RK4(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D2_RK3(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D2_RK2(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUU_ATL_D2_RK1(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_ATL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CLU_ATL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D1_RK5(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D1_RK4(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D1_RK3(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D1_RK2(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D1_RK1(affinity, sample, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D2_RK5(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D2_RK4(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D2_RK3(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D2_RK2(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLU_ATL_D2_RK1(affinity, sample, dim, tform, class, tlate)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_ATL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tlate(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CGR_DTL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D1_RK5(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D1_RK4(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D1_RK3(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D1_RK2(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D1_RK1(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D2_RK5(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D2_RK4(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D2_RK3(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D2_RK2(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CGR_DTL_D2_RK1(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CGR_DTL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(genrecmat_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CUD_DTL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D1_RK5(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D1_RK4(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D1_RK3(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D1_RK2(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D1_RK1(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D2_RK5(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D2_RK4(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D2_RK3(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D2_RK2(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUD_DTL_D2_RK1(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUD_DTL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(upperDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CLD_DTL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D1_RK5(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D1_RK4(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D1_RK3(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D1_RK2(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D1_RK1(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D2_RK5(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D2_RK4(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D2_RK3(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D2_RK2(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLD_DTL_D2_RK1(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLD_DTL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(lowerDiag_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CUU_DTL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D1_RK5(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D1_RK4(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D1_RK3(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D1_RK2(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D1_RK1(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D2_RK5(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D2_RK4(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D2_RK3(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D2_RK2(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CUU_DTL_D2_RK1(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CUU_DTL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(upperUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

    ! CLU_DTL translation arbitrary.

    interface setAffinity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D1_RK5(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D1_RK5
#endif
        use pm_kind, only: RKG => RK4
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D1_RK4(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D1_RK3(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D1_RK2(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D1_RK1(affinity, sample, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D2_RK5(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D2_RK4(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D2_RK3(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D2_RK2(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setAffinity_CLU_DTL_D2_RK1(affinity, sample, dim, tform, class)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAffinity_CLU_DTL_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                        :: dim
        type(lowerUnit_type)    , intent(in)                        :: class
        real(RKG)               , intent(out)   , contiguous        :: affinity(:,:)
        real(RKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(RKG)               , intent(in)    , contiguous        :: tform(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setAffinity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleAffinity ! LCOV_EXCL_LINE