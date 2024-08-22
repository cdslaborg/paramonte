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
!>  This module contains classes and procedures for computing the Mahalanobis statistical distance.
!>
!>  \details
!>  The Mahalanobis distance of an observation \f$\vec{x} = (x_1, x_2, x_3, \ldots, x_N)^\mathsf{H}\f$ from
!>  a set of observations represented by a Multivariate Normal (**MVN**) distribution in \f$N\f$ dimensions
!>  with \f$(\bu{\mu}, \bu{\Sigma})\f$ as its mean vector and covariance matrix is defined as,
!>  \f{equation}{
!>      \large
!>      D_M( \vec{x} ) = \sqrt{
!>          (\vec{x} - \bu{\mu})^\mathsf{H} ~ \bu{\Sigma}^{-1} (\vec{x} - \bu{\mu})
!>      }~,
!>
!>  \f}
!>  where \f$^{H}\f$ stands for the Hermitian transpose.<br>
!>  When the Covariance of the MVN distribution is the Identity matrix,
!>  the Mahalanobis distance simply becomes the Euclidean norm.
!>
!>  \benchmarks
!>
!>  \benchmark{row_vs_col_major, The runtime performance of [getDisMahalSq](@ref pm_distanceMahal::getDisMahalSq) vs. [setDisMahalSq](@ref pm_distanceMahal::setDisMahalSq)}
!>  \include{lineno} benchmark/pm_distanceMahal/row_vs_col_major/main.F90
!>  \compilefb{row_vs_col_major}
!>  \postprocb{row_vs_col_major}
!>  \include{lineno} benchmark/pm_distanceMahal/row_vs_col_major/main.py
!>  \visb{row_vs_col_major}
!>  \image html benchmark/pm_distanceMahal/row_vs_col_major/benchmark.row_vs_col_major.runtime.png width=1000
!>  \image html benchmark/pm_distanceMahal/row_vs_col_major/benchmark.row_vs_col_major.runtime.ratio.png width=1000
!>  \moralb{row_vs_col_major}
!>      -#  Fortran is a column-major language, meaning that matrix elements are stored column-wise in computer memory.<br>
!>      -#  As such, matrix multiplication format that respects column-major order of Fortran,
!>          is significantly faster than the row-major matrix multiplication.<br>
!>      -#  This is particularly relevant when one matrix is symmetric square and the other is a vector,
!>          which is the case with the procedures of the generic interface [getDisMahalSq](@ref pm_distanceMahal::getDisMahalSq).<br>
!>
!>  \benchmark{looping_vs_intrinsic, The runtime performance of [getDisMahalSq](@ref pm_distanceMahal::getDisMahalSq) vs. [setDisMahalSq](@ref pm_distanceMahal::setDisMahalSq)}
!>  \include{lineno} benchmark/pm_distanceMahal/looping_vs_intrinsic/main.F90
!>  \compilefb{looping_vs_intrinsic}
!>  \postprocb{looping_vs_intrinsic}
!>  \include{lineno} benchmark/pm_distanceMahal/looping_vs_intrinsic/main.py
!>  \visb{looping_vs_intrinsic}
!>  \image html benchmark/pm_distanceMahal/looping_vs_intrinsic/benchmark.looping_vs_intrinsic.runtime.png width=1000
!>  \image html benchmark/pm_distanceMahal/looping_vs_intrinsic/benchmark.looping_vs_intrinsic.runtime.ratio.png width=1000
!>  \moralb{looping_vs_intrinsic}
!>      -#  The procedures in this benchmark compute the Mahalanobis distance using two different implementations.<br>
!>          <ol>
!>              <li>    The procedure named `loop_and_dotp` computes the distance via looping and Fortran intrinsic `dot_product()`.<br>
!>                      This approach avoids temporary array creations.<br>
!>              <li>    The procedure named `dotp_matmul` uses the all-intrinsic expression `dot_product(vec, matmul(vec, mat))`
!>                      to compute the distance.
!>          </ol>
!>      -#  Based on the benchmark results, it appears that the looping version offers a faster implementation.<br>
!>      -#  Additionally, the specification of the slice of the matrix in the dot product
!>          of the looping approach significantly improves the performance.<br>
!>
!>  \test
!>  [test_pm_distanceMahal](@ref test_pm_distanceMahal)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distanceMahal

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distanceMahal"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the square of the Mahalanobis distance of a (set of `npnt`) point(s) from a single
    !>  (or a set of `nsam` independent) sample(s) characterized by a (set of) Multivariate Normal (**MVN**) distribution(s) in `ndim` dimensions.
    !>
    !>  \details
    !>  See [pm_distanceMahal](@ref pm_distanceMahal) for the mathematical definition of the Mahalanobis distance.
    !>
    !>  \param[in]  point   :   The input scalar, or vector of size `ndim`, or array of shape `(ndim, npnt)` of either <br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          representing the `ndim`-dimensional coordinates of (`npnt`) point(s)
    !>                          whose distance(s) from the specified MVN distribution(s) should be computed.
    !>  \param[in]  invCov  :   The input scalar or matrix of shape `(ndim, ndim)` or
    !>                          array of shape `(ndim, ndim, nsam)` of the same type and kind as the input `point`,
    !>                          representing the inverse covariance matrix(es) of the MVN distribution(s).<br>
    !>                          (**optional**, default = `1` <b>if and only if</b> `point` is a scalar, otherwise, the Identity matrix of shape `(ndim, ndim)`).
    !>  \param[in]  center  :   The input scalar or vector size `ndim` or matrix of shape `(ndim, nsam)` of the same type and kind as the input `point`,
    !>                          representing the center of the MVN distribution(s).<br>
    !>                          (**optional**, default = `0`).
    !>
    !>  \return
    !>  `mahalSq`           :   The output of the same type and kind as the input `point` representing
    !>                          the square of the Mahalanobis distance of `point` from the MVN distribution(s) specified by the arguments `center` and `invCov`.<br>
    !>                          On output, `mahalSq` is,
    !>                          <ol>
    !>                              <li>    a scalar if,
    !>                                      <ol>
    !>                                          <li>    the input `point` is a scalar.<br>
    !>                                          <li>    the input `point` is of shape `(1:ndim)` and `invCov` is of shape `(1:ndim, 1:ndim)`.<br>
    !>                                      </ol>
    !>                              <li>    a vector of size `(1:nsam)` if the input `point` is of shape `(1:ndim)` and `invCov` is of shape `(1:ndim, 1:ndim, 1:nsam)`.<br>
    !>                              <li>    a vector of size `(1:npnt)` if the input `point` is of shape `(1:ndim, 1:npnt)` and `invCov` is of shape `(1:ndim, 1:ndim)`.<br>
    !>                              <li>    a matrix of size `(1:nsam, 1:npnt)` if the input `point` is of shape `(1:ndim, 1:npnt)` and `invCov` is of shape `(1:ndim, 1:ndim, 1:nsam)`.<br>
    !>                          </ol>
    !>
    !>  \interface{getDisMahalSq}
    !>  \code{.F90}
    !>
    !>      use pm_distanceMahal, only: getDisMahalSq
    !>
    !>      mahalSq = getDisMahalSq(point, invCov) ! elemental
    !>      mahalSq = getDisMahalSq(point, invCov, mean) ! elemental
    !>
    !>      mahalSq = getDisMahalSq(point(1:ndim), invCov(1:ndim, 1:ndim))
    !>      mahalSq = getDisMahalSq(point(1:ndim), invCov(1:ndim, 1:ndim), center(1:ndim))
    !>
    !>      mahalSq(1:npnt) = getDisMahalSq(point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim))
    !>      mahalSq(1:npnt) = getDisMahalSq(point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim), center(1:ndim))
    !>
    !>      mahalSq(1:nsam) = getDisMahalSq(point(1:ndim), invCov(1:ndim, 1:ndim, 1:nsam))
    !>      mahalSq(1:nsam) = getDisMahalSq(point(1:ndim), invCov(1:ndim, 1:ndim, 1:nsam), center(1:ndim, 1:nsam))
    !>
    !>      mahalSq(1:nsam, 1:npnt) = getDisMahalSq(point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim, 1:nsam))
    !>      mahalSq(1:nsam, 1:npnt) = getDisMahalSq(point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim, 1:nsam), center(1:ndim, 1:nsam))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(point, 1) == size(center, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(point, 1) == [size(invCov, 1), size(invCov, 2)])` must hold for the corresponding input arguments.<br>
    !>  The condition `size(center, rank(center)) == size(invCov, rank(invCov))` must hold for the corresponding input arguments.<br>
    !>  The condition `isMatClass(invCov, posdefmat)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` only when all input
    !>  arguments are scalars or arrays of the same rank as other array-like arguments.<br>
    !>
    !>  \note
    !>  The computation of the Mahalanobis distance for complex arguments
    !>  follows the normal intrinsic Fortran rules for complex arithmetic.
    !>
    !>  \example{getDisMahalSq}
    !>  \include{lineno} example/pm_distanceMahal/getDisMahalSq/main.F90
    !>  \compilef{getDisMahalSq}
    !>  \output{getDisMahalSq}
    !>  \include{lineno} example/pm_distanceMahal/getDisMahalSq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceMahal](@ref test_pm_distanceMahal)<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The runtime checks for the complex input `invCov` must be implemented.<br>
    !>
    !>  \final{getDisMahalSq}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getDisMahalSq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_CK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: point, invCov
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_CK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: point, invCov
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_CK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: point, invCov
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_CK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: point, invCov
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_CK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: point, invCov
        complex(CKG)                            :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_RK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: point, invCov
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_RK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: point, invCov
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_RK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: point, invCov
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_RK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: point, invCov
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getDisMahalSqEleInvDef_D0_RK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvDef_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: point, invCov
        real(RKG)                               :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_CK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                :: point, invCov, center
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_CK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                :: point, invCov, center
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_CK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                :: point, invCov, center
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_CK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                :: point, invCov, center
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_CK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                :: point, invCov, center
        complex(CKG)                            :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_RK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: point, invCov, center
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_RK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: point, invCov, center
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_RK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: point, invCov, center
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_RK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: point, invCov, center
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getDisMahalSqEleInvCen_D0_RK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqEleInvCen_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: point, invCov, center
        real(RKG)                               :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_CK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_CK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_CK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_CK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_CK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:)
        complex(CKG)                            :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_RK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_RK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_RK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_RK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqOneInvDef_D1_RK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:)
        real(RKG)                               :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_CK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_CK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_CK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_CK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_CK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_RK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_RK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_RK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_RK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqOneInvCen_D1_RK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_CK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_CK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_CK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_CK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_CK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_RK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_RK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_RK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_RK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqOneInvDef_D2_RK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvDef_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_CK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_CK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_CK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_CK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_CK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG)                            :: mahalSq(size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_RK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_RK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_RK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_RK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqOneInvCen_D2_RK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqOneInvCen_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)                               :: mahalSq(size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_CK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_CK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_CK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_CK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_CK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_RK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_RK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_RK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_RK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqMixInvDef_D1_RK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_CK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_CK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_CK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_CK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_CK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_RK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_RK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_RK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_RK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqMixInvCen_D1_RK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_CK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_CK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_CK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_CK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_CK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_RK5(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_RK4(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_RK3(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_RK2(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqMixInvDef_D2_RK1(point, invCov) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvDef_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_CK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_CK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_CK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_CK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_CK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG)                            :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_RK5(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_RK4(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_RK3(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_RK2(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisMahalSqMixInvCen_D2_RK1(point, invCov, center) result(mahalSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisMahalSqMixInvCen_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)                               :: mahalSq(size(invCov, 3, IK), size(point, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the square of the Mahalanobis distance of a (set of `npnt`) point(s) from a single (or a set of `nsam` independent)
    !>  sample(s) characterized by a (set of) Multivariate Normal (**MVN**) distribution(s) in `ndim` dimensions.
    !>
    !>  \details
    !>  See [pm_distanceMahal](@ref pm_distanceMahal) for the mathematical definition of the Mahalanobis distance.<br>
    !>
    !>  \param[out] mahalSq :   The output of the same type and kind as the input `point` representing
    !>                          the square of the Mahalanobis distance of `point` from the MVN distribution(s) specified by the arguments `center` and `invCov`.<br>
    !>                          On output, `mahalSq` is,
    !>                          <ol>
    !>                              <li>    a scalar if,
    !>                                      <ol>
    !>                                          <li>    the input `point` is a scalar.<br>
    !>                                          <li>    the input `point` is of shape `(1:ndim)` and `invCov` is of shape `(1:ndim, 1:ndim)`.<br>
    !>                                      </ol>
    !>                              <li>    a vector of size `(1:nsam)` if the input `point` is of shape `(1:ndim)` and `invCov` is of shape `(1:ndim, 1:ndim, 1:nsam)`.<br>
    !>                              <li>    a vector of size `(1:npnt)` if the input `point` is of shape `(1:ndim, 1:npnt)` and `invCov` is of shape `(1:ndim, 1:ndim)`.<br>
    !>                              <li>    a matrix of size `(1:nsam, 1:npnt)` if the input `point` is of shape `(1:ndim, 1:npnt)` and `invCov` is of shape `(1:ndim, 1:ndim, 1:nsam)`.<br>
    !>                          </ol>
    !>  \param[in]  point   :   The input scalar, or vector of size `ndim`, or array of shape `(ndim, npnt)` of either <br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          representing the `ndim`-dimensional coordinates of (`npnt`) point(s)
    !>                          whose distance(s) from the specified MVN distribution(s) should be computed.
    !>  \param[in]  invCov  :   The input scalar or matrix of shape `(ndim, ndim)` or
    !>                          array of shape `(ndim, ndim, nsam)` of the same type and kind as the input `point`,
    !>                          representing the inverse covariance matrix(es) of the MVN distribution(s).<br>
    !>                          (**optional**, default = `1` <b>if and only if</b> `point` is a scalar, otherwise, the Identity matrix of shape `(ndim, ndim)`).
    !>  \param[in]  center  :   The input scalar or vector size `ndim` or matrix of shape `(ndim, nsam)` of the same type and kind as the input `point`,
    !>                          representing the center of the MVN distribution(s).<br>
    !>                          (**optional**, default = `0`).
    !>
    !>
    !>  \interface{setDisMahalSq}
    !>  \code{.F90}
    !>
    !>      use pm_distanceMahal, only: setDisMahalSq
    !>
    !>      call setDisMahalSq(mahalSq, point, invCov) ! elemental
    !>      call setDisMahalSq(mahalSq, point, invCov, mean) ! elemental
    !>
    !>      call setDisMahalSq(mahalSq, point(1:ndim), invCov(1:ndim, 1:ndim))
    !>      call setDisMahalSq(mahalSq, point(1:ndim), invCov(1:ndim, 1:ndim), center(1:ndim))
    !>
    !>      call setDisMahalSq(mahalSq(1:npnt), point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim))
    !>      call setDisMahalSq(mahalSq(1:npnt), point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim), center(1:ndim))
    !>
    !>      call setDisMahalSq(mahalSq(1:nsam), point(1:ndim), invCov(1:ndim, 1:ndim, 1:nsam))
    !>      call setDisMahalSq(mahalSq(1:nsam), point(1:ndim), invCov(1:ndim, 1:ndim, 1:nsam), center(1:ndim, 1:nsam))
    !>
    !>      call setDisMahalSq(mahalSq(1:nsam, 1:npnt), point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim, 1:nsam))
    !>      call setDisMahalSq(mahalSq(1:nsam, 1:npnt), point(1:ndim, 1:npnt), invCov(1:ndim, 1:ndim, 1:nsam), center(1:ndim, 1:nsam))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(point, 1) == size(center, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(point, 1) == [size(invCov, 1), size(invCov, 2)])` must hold for the corresponding input arguments.<br>
    !>  The condition `size(center, rank(center)) == size(invCov, rank(invCov))` must hold for the corresponding input arguments.<br>
    !>  The condition `isMatClass(invCov, posdefmat)` must hold for the corresponding input arguments.<br>
    !>  The size and shape of the input `mahalSq` must be consistent with other input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` only when all input
    !>  arguments are scalars or arrays of the same rank as other array-like arguments.<br>
    !>
    !>  \note
    !>  The computation of the Mahalanobis distance for complex arguments
    !>  follows the normal intrinsic Fortran rules for complex arithmetic.<br>
    !>
    !>  \example{setDisMahalSq}
    !>  \include{lineno} example/pm_distanceMahal/setDisMahalSq/main.F90
    !>  \compilef{setDisMahalSq}
    !>  \output{setDisMahalSq}
    !>  \include{lineno} example/pm_distanceMahal/setDisMahalSq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceMahal](@ref test_pm_distanceMahal)<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The runtime checks for the complex input `invCov` must be implemented.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The performance of the implementation for `complex` input can be improved by using `do_product` on columns of `invCov`
    !>  instead of the current implementation working with rows of `invCov` in `matmul`.<br>
    !>
    !>  \final{setDisMahalSq}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface setDisMahalSq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_CK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                    :: point, invCov
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_CK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                    :: point, invCov
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_CK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                    :: point, invCov
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_CK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                    :: point, invCov
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_CK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                    :: point, invCov
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_RK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: point, invCov
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_RK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: point, invCov
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_RK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: point, invCov
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_RK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: point, invCov
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvDef_D0_RK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvDef_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: point, invCov
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_CK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                    :: point, invCov, center
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK4_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_CK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                    :: point, invCov, center
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK3_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_CK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                    :: point, invCov, center
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK2_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_CK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                    :: point, invCov, center
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK1_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_CK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                    :: point, invCov, center
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_RK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: point, invCov, center
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_RK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: point, invCov, center
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_RK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: point, invCov, center
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_RK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: point, invCov, center
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setDisMahalSqEleInvCen_D0_RK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqEleInvCen_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: point, invCov, center
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_CK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_CK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_CK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_CK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_CK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_RK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_RK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_RK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_RK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D1_RK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_CK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_CK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_CK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_CK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_CK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        complex(CKG), intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_RK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_RK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_RK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_RK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D1_RK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:), center(:)
        real(RKG)   , intent(out)                   :: mahalSq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_CK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_CK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_CK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_CK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_CK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_RK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_RK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_RK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_RK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvDef_D2_RK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvDef_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_CK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_CK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_CK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_CK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_CK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_RK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_RK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_RK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_RK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqOneInvCen_D2_RK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqOneInvCen_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:), center(:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_CK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_CK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_CK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_CK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_CK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_RK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_RK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_RK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_RK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D1_RK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_CK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_CK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_CK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_CK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_CK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_RK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_RK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_RK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_RK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D1_RK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_CK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_CK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_CK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_CK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_CK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_RK5(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_RK4(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_RK3(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_RK2(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvDef_D2_RK1(mahalSq, point, invCov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvDef_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_CK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_CK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_CK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_CK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_CK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        complex(CKG), intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_RK5(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_RK4(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_RK3(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_RK2(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setDisMahalSqMixInvCen_D2_RK1(mahalSq, point, invCov, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setDisMahalSqMixInvCen_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:), invCov(:,:,:), center(:,:)
        real(RKG)   , intent(out)   , contiguous    :: mahalSq(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distanceMahal ! LCOV_EXCL_LINE