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
!>  This module contains procedures and generic interfaces for
!>  computing the additive or multiplicative trace of a given square matrix in arbitrary [packing formats](@ref pm_matrixPack).
!>
!>  \details
!>  In linear algebra, the trace of a square matrix \f$A\f$, denoted \f$tr(A)\f$,
!>  is defined to be the sum of elements on the main diagonal (from the upper left to the lower right) of \f$A\f$.<br>
!>  The trace is only defined for a **square matrix** \f$(n\times n)\f$.<br>
!>  It can be proven that the trace of a matrix is the sum of its (complex) eigenvalues (counted with multiplicities).<br>
!>  It can also be proven that \f$tr(AB) = tr(BA)\f$ for any two matrices \f$A\f$ and \f$B\f$.<br>
!>  This implies that similar matrices have the same trace.<br>
!>  As a consequence one can define the trace of a linear operator mapping a finite-dimensional vector space into itself,
!>  since all matrices describing such an operator with respect to a basis are similar.<br>
!>  The trace is related to the derivative of the determinant (through [Jacobi formula](https://en.wikipedia.org/wiki/Jacobi%27s_formula)).<br>
!>
!>  The trace of an \f$(n\times n)\f$ square matrix \f$A\f$ is defined as,
!>  \f{equation}{
!>      \up{tr} (\mathbf {A} ) = \sum_{i=1}^{n}a_{ii} = a_{11} + a_{22} + \dots + a_{nn} ~,
!>  \f}
!>  where \f$a_{ii}\f$ denotes the entry on the \f$i\f$th row and \f$i\f$th column of \f$A\f$.<br>
!>  The entries of \f$A\f$ can be real numbers or (more generally) complex numbers.<br>
!>
!>  \note
!>  The trace is not defined for non-square matrices.<br>
!>  Expressions like \f$tr(exp(A))\f$, where \f$A\f$ is a square matrix, occur frequently in some fields
!>  (e.g. multivariate statistical theory), that a shorthand notation has become common,
!>  \f{equation}{
!>      \up{tre}(A) = \up{tr}(\exp(A)) ~.
!>  \f}
!>  \f$\up{tre}\f$ is sometimes referred to as the **exponential trace function** and is used in the
!>  [Golden–Thompson inequality](https://en.wikipedia.org/wiki/Golden%E2%80%93Thompson_inequality).<br>
!>
!>  This module additionally also computes the **multiplicative trace function** of a square-matrix which is defined as,
!>  \f{equation}{
!>      \up{trl}(\mathbf{A}) = \prod_{i=1}^{n} a_{ii} = a_{11} \times a_{22} \times \dots \times a_{nn} ~,
!>  \f}
!>  where \f$a_{ii}\f$ denotes the entry on the \f$i\f$th row and \f$i\f$th column of \f$A\f$.<br>
!>  The entries of \f$A\f$ can be real numbers or (more generally) complex numbers.<br>
!>  This multiplicative definition of trace appears in the computation of the determinant of square positive definite matrices.<br>
!>
!>  \see
!>  [pm_matrixDet](@ref pm_matrixDet)<br>
!>
!>  \test[test_pm_matrixIndex](@ref test_pm_matrixTrace)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixTrace

    use pm_kind, only: SK, IK, RKD
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixPack, only: rdpack, rdpack_type
    use pm_matrixPack, only: lfpack, lfpack_type
    use pm_matrixPack, only: rfpack, rfpack_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixPack"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the trace of an input square matrix of type `integer`, `complex`, or `real` of arbitrary kind.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_matrixTrace](@ref pm_matrixTrace) for more details.<br>
    !>
    !>  \param[in]  mat     :   The input non-zero-order matrix of shape `(1 : ndim * (ndim - 1) / 2` or `(1:ndim, 1:ndim)` of,
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL,
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the matrix whose trace is to be computed.<br>
    !>                          It can be of rank `1` (of shape `(1 : ndim * (ndim - 1) / 2`) **if and only if**
    !>                          the input argument `mat` is present and set to [lfpack](@ref pm_matrixPack::lfpack).<br>
    !>  \param[in]  pack    :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [lfpack](@ref pm_matrixPack::lfpack)
    !>                                      signifying the [linear full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rdpack](@ref pm_matrixPack::rdpack)
    !>                                      signifying the [rectangular default packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rfpack](@ref pm_matrixPack::rfpack)
    !>                                      signifying the [rectangular full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                          </ol>
    !>                          (**optional**, default = [rdpack](@ref pm_matrixPack::rdpack))
    !>  \param[in]  subset  :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia)
    !>                                      signifying the storage of the upper-diagonal subset of original matrix in the input `matrix`.<br>
    !>                              <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia)
    !>                                      signifying the storage of the lower-diagonal subset of original matrix in the input `matrix`.<br>
    !>                          </ol>
    !>                          (**optional**. It must be present **if and only if** the input argument is present and set to either
    !>                          [lfpack](@ref pm_matrixPack::lfpack) or [rfpack](@ref pm_matrixPack::rfpack).)
    !>
    !>  \return
    !>  `trace`             :   The output scalar of the same type and kind as the input matrix `mat`,
    !>                          containing the trace of the input square matrix.<br>
    !>
    !>  \interface{getMatTrace}
    !>  \code{.F90}
    !>
    !>      use pm_matrixTrace, only: getMatTrace
    !>      use pm_matrixTrace, only: lfpack, rdpack, rfpack, uppDia, lowDia
    !>
    !>      trace = getMatTrace(mat)
    !>      trace = getMatTrace(mat, pack)          ! pack = rdpack
    !>      trace = getMatTrace(mat, pack, subset)  ! pack = rdpack/lfpack, subset = uppDia/lowDia
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input matrix must correspond to a square matrix or a triangular subset of a non-zero-order square matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatTrace](@ref pm_matrixTrace::getMatTrace)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>  [getMatMulTraceLog](@ref pm_matrixTrace::getMatMulTraceLog)<br>
    !>  [pm_matrixPack](@ref pm_matrixPack)<br>
    !>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
    !>
    !>  \example{getMatTrace}
    !>  \include{lineno} example/pm_matrixTrace/getMatTrace/main.F90
    !>  \compilef{getMatTrace}
    !>  \output{getMatTrace}
    !>  \include{lineno} example/pm_matrixTrace/getMatTrace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixTrace](@ref test_pm_matrixTrace)
    !>
    !>  \final{getMatTrace}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getMatTrace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatTrace_DEF_XXX_IK5(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: trace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_DEF_XXX_IK4(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: trace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_DEF_XXX_IK3(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: trace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_DEF_XXX_IK2(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: trace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_DEF_XXX_IK1(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: trace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_DEF_XXX_CK5(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: trace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_DEF_XXX_CK4(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: trace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_DEF_XXX_CK3(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: trace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_DEF_XXX_CK2(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: trace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_DEF_XXX_CK1(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: trace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_DEF_XXX_RK5(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_DEF_XXX_RK4(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_DEF_XXX_RK3(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_DEF_XXX_RK2(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_DEF_XXX_RK1(mat) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_DEF_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
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

#if IK5_ENABLED
    PURE module function getMatTrace_RDP_XXX_IK5(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: trace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_RDP_XXX_IK4(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: trace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_RDP_XXX_IK3(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: trace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_RDP_XXX_IK2(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: trace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_RDP_XXX_IK1(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: trace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_RDP_XXX_CK5(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: trace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_RDP_XXX_CK4(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: trace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_RDP_XXX_CK3(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: trace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_RDP_XXX_CK2(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: trace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_RDP_XXX_CK1(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: trace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_RDP_XXX_RK5(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: trace
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_RDP_XXX_RK4(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: trace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_RDP_XXX_RK3(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: trace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_RDP_XXX_RK2(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: trace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_RDP_XXX_RK1(mat, pack) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RDP_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: trace
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

#if IK5_ENABLED
    PURE module function getMatTrace_RFP_UXD_IK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_RFP_UXD_IK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_RFP_UXD_IK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_RFP_UXD_IK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_RFP_UXD_IK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_RFP_UXD_CK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_RFP_UXD_CK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_RFP_UXD_CK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_RFP_UXD_CK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_RFP_UXD_CK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_RFP_UXD_RK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_RFP_UXD_RK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_RFP_UXD_RK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_RFP_UXD_RK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_RFP_UXD_RK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatTrace_RFP_XLD_IK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_RFP_XLD_IK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_RFP_XLD_IK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_RFP_XLD_IK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_RFP_XLD_IK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_RFP_XLD_CK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_RFP_XLD_CK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_RFP_XLD_CK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_RFP_XLD_CK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_RFP_XLD_CK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_RFP_XLD_RK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_RFP_XLD_RK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_RFP_XLD_RK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_RFP_XLD_RK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_RFP_XLD_RK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_RFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
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

#if IK5_ENABLED
    PURE module function getMatTrace_LFP_UXD_IK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_LFP_UXD_IK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_LFP_UXD_IK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_LFP_UXD_IK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_LFP_UXD_IK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_LFP_UXD_CK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_LFP_UXD_CK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_LFP_UXD_CK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_LFP_UXD_CK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_LFP_UXD_CK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_LFP_UXD_RK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_LFP_UXD_RK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_LFP_UXD_RK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_LFP_UXD_RK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_LFP_UXD_RK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatTrace_LFP_XLD_IK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatTrace_LFP_XLD_IK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatTrace_LFP_XLD_IK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatTrace_LFP_XLD_IK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatTrace_LFP_XLD_IK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: trace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatTrace_LFP_XLD_CK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatTrace_LFP_XLD_CK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatTrace_LFP_XLD_CK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatTrace_LFP_XLD_CK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatTrace_LFP_XLD_CK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: trace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatTrace_LFP_XLD_RK5(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatTrace_LFP_XLD_RK4(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatTrace_LFP_XLD_RK3(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatTrace_LFP_XLD_RK2(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatTrace_LFP_XLD_RK1(mat, pack, subset) result(trace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatTrace_LFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: trace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

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
    !>  Generate and return the **multiplicative trace** of an input square matrix of type `integer`, `complex`, or `real` of arbitrary kind.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_matrixTrace](@ref pm_matrixTrace) for more details.<br>
    !>
    !>  \param[in]  mat     :   The input non-zero-order matrix of shape `(1 : ndim * (ndim - 1) / 2` or `(1:ndim, 1:ndim)` of,
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL,
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the matrix whose **multiplicative trace** is to be computed.<br>
    !>                          It can be of rank `1` (of shape `(1 : ndim * (ndim - 1) / 2`) **if and only if**
    !>                          the input argument `mat` is present and set to [lfpack](@ref pm_matrixPack::lfpack).<br>
    !>  \param[in]  pack    :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [lfpack](@ref pm_matrixPack::lfpack)
    !>                                      signifying the [linear full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rdpack](@ref pm_matrixPack::rdpack)
    !>                                      signifying the [rectangular default packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rfpack](@ref pm_matrixPack::rfpack)
    !>                                      signifying the [rectangular full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                          </ol>
    !>                          (**optional**, default = [rdpack](@ref pm_matrixPack::rdpack))
    !>  \param[in]  subset  :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia)
    !>                                      signifying the storage of the upper-diagonal subset of original matrix in the input `matrix`.<br>
    !>                              <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia)
    !>                                      signifying the storage of the lower-diagonal subset of original matrix in the input `matrix`.<br>
    !>                          </ol>
    !>                          (**optional**. It must be present **if and only if** the input argument is present and set to either
    !>                          [lfpack](@ref pm_matrixPack::lfpack) or [rfpack](@ref pm_matrixPack::rfpack).)
    !>
    !>  \return
    !>  `mulTrace`          :   The output scalar of the same type and kind as the input matrix `mat`,
    !>                          containing the **multiplicative trace** of the input square matrix.<br>
    !>
    !>  \interface{getMatMulTrace}
    !>  \code{.F90}
    !>
    !>      use pm_matrixTrace, only: getMatMulTrace
    !>      use pm_matrixTrace, only: lfpack, rdpack, rfpack, uppDia, lowDia
    !>
    !>      mulTrace = getMatMulTrace(mat)
    !>      mulTrace = getMatMulTrace(mat, pack)          ! pack = rdpack
    !>      mulTrace = getMatMulTrace(mat, pack, subset)  ! pack = rdpack/lfpack, subset = uppDia/lowDia
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input matrix must correspond to a square matrix or a triangular subset of a non-zero-order square matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatTrace](@ref pm_matrixTrace::getMatTrace)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>  [getMatMulTraceLog](@ref pm_matrixTrace::getMatMulTraceLog)<br>
    !>  [pm_matrixPack](@ref pm_matrixPack)<br>
    !>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
    !>
    !>  \example{getMatMulTrace}
    !>  \include{lineno} example/pm_matrixTrace/getMatMulTrace/main.F90
    !>  \compilef{getMatMulTrace}
    !>  \output{getMatMulTrace}
    !>  \include{lineno} example/pm_matrixTrace/getMatMulTrace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixTrace](@ref test_pm_matrixTrace)
    !>
    !>  \final{getMatMulTrace}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getMatMulTrace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_IK5(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_IK4(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_IK3(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_IK2(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_IK1(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        integer(IKC)                                        :: mulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_CK5(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_CK4(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_CK3(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_CK2(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_CK1(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: mulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_RK5(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_RK4(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_RK3(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_RK2(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_DEF_XXX_RK1(mat) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_DEF_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
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

#if IK5_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_IK5(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_IK4(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_IK3(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_IK2(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: mulTrace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_IK1(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        integer(IKC)                                        :: mulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_CK5(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_CK4(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_CK3(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_CK2(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: mulTrace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_CK1(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: mulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_RK5(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: mulTrace
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_RK4(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_RK3(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_RK2(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: mulTrace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_RDP_XXX_RK1(mat, pack) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RDP_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: mulTrace
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

#if IK5_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_IK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_IK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_IK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_IK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_IK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_CK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_CK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_CK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_CK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_CK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_RK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_RK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_RK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_RK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_RFP_UXD_RK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_IK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_IK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_IK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_IK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_IK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_CK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_CK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_CK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_CK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_CK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_RK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_RK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_RK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_RK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_RFP_XLD_RK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_RFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
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

#if IK5_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_IK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_IK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_IK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_IK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_IK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_CK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_CK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_CK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_CK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_CK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_RK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_RK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_RK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_RK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_LFP_UXD_RK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_IK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_IK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_IK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_IK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_IK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                                        :: mulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_CK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_CK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_CK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_CK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_CK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: mulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_RK5(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_RK4(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_RK3(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_RK2(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTrace_LFP_XLD_RK1(mat, pack, subset) result(mulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTrace_LFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: mulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

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
    !>  Generate and return the natural logarithm of the **multiplicative trace** of an input square matrix of type `integer`, `complex`, or `real` of arbitrary kind.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_matrixTrace](@ref pm_matrixTrace) for more details.<br>
    !>
    !>  \param[in]  mat     :   The input non-zero-order matrix of shape `(1 : ndim * (ndim - 1) / 2` or `(1:ndim, 1:ndim)` of,
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL,
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the matrix whose natural logarithm of the **multiplicative trace** is to be computed.<br>
    !>                          It can be of rank `1` (of shape `(1 : ndim * (ndim - 1) / 2`) **if and only if**
    !>                          the input argument `mat` is present and set to [lfpack](@ref pm_matrixPack::lfpack).<br>
    !>  \param[in]  pack    :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [lfpack](@ref pm_matrixPack::lfpack)
    !>                                      signifying the [linear full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rdpack](@ref pm_matrixPack::rdpack)
    !>                                      signifying the [rectangular default packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                              <li>    the constant [rfpack](@ref pm_matrixPack::rfpack)
    !>                                      signifying the [rectangular full packing](@ref pm_matrixPack) of the input matrix.<br>
    !>                          </ol>
    !>                          (**optional**, default = [rdpack](@ref pm_matrixPack::rdpack))
    !>  \param[in]  subset  :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia)
    !>                                      signifying the storage of the upper-diagonal subset of original matrix in the input `matrix`.<br>
    !>                              <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia)
    !>                                      signifying the storage of the lower-diagonal subset of original matrix in the input `matrix`.<br>
    !>                          </ol>
    !>                          (**optional**. It must be present **if and only if** the input argument is present and set to either
    !>                          [lfpack](@ref pm_matrixPack::lfpack) or [rfpack](@ref pm_matrixPack::rfpack).)
    !>
    !>  \return
    !>  `logMulTrace`          :   The output scalar of the same type and kind as the input matrix `mat`,
    !>                          containing the natural logarithm of the **multiplicative trace** of the input square matrix.<br>
    !>
    !>  \interface{getMatMulTraceLog}
    !>  \code{.F90}
    !>
    !>      use pm_matrixTrace, only: getMatMulTraceLog
    !>      use pm_matrixTrace, only: lfpack, rdpack, rfpack, uppDia, lowDia
    !>
    !>      logMulTrace = getMatMulTraceLog(mat)
    !>      logMulTrace = getMatMulTraceLog(mat, pack)          ! pack = rdpack
    !>      logMulTrace = getMatMulTraceLog(mat, pack, subset)  ! pack = rdpack/lfpack, subset = uppDia/lowDia
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input matrix must correspond to a square matrix or a triangular subset of a non-zero-order square matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatTrace](@ref pm_matrixTrace::getMatTrace)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>  [getMatMulTraceLog](@ref pm_matrixTrace::getMatMulTraceLog)<br>
    !>  [pm_matrixPack](@ref pm_matrixPack)<br>
    !>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
    !>
    !>  \example{getMatMulTraceLog}
    !>  \include{lineno} example/pm_matrixTrace/getMatMulTraceLog/main.F90
    !>  \compilef{getMatMulTraceLog}
    !>  \output{getMatMulTraceLog}
    !>  \include{lineno} example/pm_matrixTrace/getMatMulTraceLog/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixTrace](@ref test_pm_matrixTrace)
    !>
    !>  \final{getMatMulTraceLog}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getMatMulTraceLog

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_IK5(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_IK4(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_IK3(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_IK2(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_IK1(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        real(RKD)                                           :: logMulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_CK5(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_CK4(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_CK3(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_CK2(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_CK1(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        complex(CKC)                                        :: logMulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_RK5(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_RK4(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_RK3(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_RK2(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_DEF_XXX_RK1(mat) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_DEF_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
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

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_IK5(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_IK4(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_IK3(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_IK2(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKD)                                           :: logMulTrace
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_IK1(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKD)                                           :: logMulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_CK5(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_CK4(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_CK3(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_CK2(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: logMulTrace
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_CK1(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        complex(CKC)                                        :: logMulTrace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_RK5(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)                    :: mat(:,:)
        real(RKC)                                           :: logMulTrace
        type(rdpack_type)   , intent(in)                    :: pack
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_RK4(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_RK3(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_RK2(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: logMulTrace
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_RDP_XXX_RK1(mat, pack) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RDP_XXX_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rdpack_type)   , intent(in)                    :: pack
        real(RKC)                                           :: logMulTrace
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

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_IK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_IK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_IK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_IK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_IK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_CK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_CK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_CK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_CK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_CK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_RK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_RK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_RK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_RK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_UXD_RK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_IK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_IK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_IK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_IK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_IK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_CK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_CK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_CK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_CK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_CK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_RK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_RK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_RK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_RK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_RFP_XLD_RK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_RFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:,:)
        type(rfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
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

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_IK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_IK5
#endif
        use pm_kind, only: IKC => IK5
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_IK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_IK4
#endif
        use pm_kind, only: IKC => IK4
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_IK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_IK3
#endif
        use pm_kind, only: IKC => IK3
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_IK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_IK2
#endif
        use pm_kind, only: IKC => IK2
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_IK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_IK1
#endif
        use pm_kind, only: IKC => IK1
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_CK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_CK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_CK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_CK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_CK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_RK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_RK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_RK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_RK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_UXD_RK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_UXD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(uppDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_IK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_IK5
#endif
        use pm_kind, only: IKC => IK5
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_IK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_IK4
#endif
        use pm_kind, only: IKC => IK4
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_IK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_IK3
#endif
        use pm_kind, only: IKC => IK3
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_IK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_IK2
#endif
        use pm_kind, only: IKC => IK2
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_IK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_IK1
#endif
        use pm_kind, only: IKC => IK1
        real(RKD)                                           :: logMulTrace
        integer(IKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_CK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_CK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_CK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_CK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_CK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                                        :: logMulTrace
        complex(CKC)        , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_RK5(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_RK4(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_RK3(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_RK2(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatMulTraceLog_LFP_XLD_RK1(mat, pack, subset) result(logMulTrace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatMulTraceLog_LFP_XLD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                           :: logMulTrace
        real(RKC)           , intent(in)                    :: mat(:)
        type(lfpack_type)   , intent(in)                    :: pack
        type(lowDia_type)   , intent(in)                    :: subset
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixTrace