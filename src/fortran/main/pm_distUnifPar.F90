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
!>  This module contains classes and procedures for setting up and
!>  computing the properties of the MultiVariate Uniform Parallelepiped (MVUP) Distribution.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>MultiVariate Uniform Parallelepiped (MVUP) distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>  </ol>
!>
!>  An \f$\ndim\f$-dimensional parallelepiped \f$P\f$ in vector space \f$\mathbb{R}^{\ndim}\f$
!>  can be defined by a set of arbitrary but independent column vectors \f$v_1, \ldots, v_\ndim\f$ as the set,<br>
!>  \f{equation}{
!>      P = {\sum_{i = 1}^{\ndim} t_i v_i ~~~,~~~ 0 \leq t_i < 1} ~,
!>  \f}
!>  where \f$t_i\f$ are a set of coefficients whose defined range allows full coverage of the parallelepiped.<br>
!>
!>  \image html "pm_distUnifPar@parallelepiped.png" "A 3-dimensional parallelepiped." width=25%
!>
!>  The above parallelepiped can be expressed in the form of a square **representative matrix** of edges of rank \f$\ndim\f$,<br>
!>  \f{equation}{
!>      M_R =
!>      \begin{pmatrix}
!>              v_1 ~,~ \vdots ~,~ v_i ~,~ \vdots ~,~ v_{\ndim}
!>      \end{pmatrix}
!>  \f}
!>  where \f$v_i\f$ is a **column vector** representing the \f$i\f$th edge of the parallelepiped.<br>
!>  The corresponding (positive definite) [**Gramian matrix**](https://en.wikipedia.org/wiki/Gram_matrix) of the parallelepiped is,<br>
!>  \f{equation}{
!>      M_G = M_R^T M_R ~,
!>  \f}
!>  where \f$M_R^T\f$ is the transpose of \f$M_R\f$.<br>
!>  The hyper-volume occupied by the parallelepiped is the given by,<br>
!>  \f{equation}{
!>      \ms{Vol}(P) = |M_R| = \sqrt{|M_G|} ~,
!>  \f}
!>  where \f$|M_G|\f$ represents the determinant of \f$M_G\f$.<br>
!>
!>  The Probability Density Function (**PDF**) of the Uniform Parallelepiped distribution with support \f$P\f$ is given by,<br>
!>  \f{equation}{
!>      \pi(X | P) = \frac{1}{\ms{Vol}(P)} = \frac{1}{|M_R|} = \frac{1}{\sqrt{|M_G|}} ~.
!>  \f}
!>
!>  \note
!>  <ol>
!>      <li>    A parallelepiped with a diagonal representative matrix \f$M_R\f$ is a **hyper-rectangle**.<br>
!>              <ol>
!>                  <li>    A hyper-rectangle is uniquely determined by lengths of its \f$\ndim\f$ edges.<br>
!>                  <li>    The corresponding Gramian matrix \f$M_G\f$ of a hyper-rectangle is also diagonal.<br>
!>              </ol>
!>      <li>    A hyper-rectangle with equal diagonal elements in the representative matrix \f$M_R\f$ is a **hyper-cube**.<br>
!>              <ol>
!>                  <li>    A hyper-cube is uniquely determined by a single scalar number representing the length of its edge that is the same along all axes.<br>
!>                  <li>    The representative matrix \f$M_R\f$ of a hyper-cube is a multiple of the Identity matrix.<br>
!>                  <li>    The corresponding Gramian matrix \f$M_G\f$ a hyper-cube is also a multiple of the Identity matrix.<br>
!>              </ol>
!>  </ol>
!>
!>  \see
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>
!>  \test
!>  [test_pm_distUnifPar](@ref test_pm_distUnifPar)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distUnifPar

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distUnifPar"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type MultiVariate Uniform Parallelepiped (MVUP)
    !>  as defined in the description of [pm_distUnifPar](@ref pm_distUnifPar).
    !>
    !>  \details
    !>  See the documentation of [pm_distUnifPar](@ref pm_distUnifPar) for the definition of the MultiVariate Uniform Parallelepiped (MVUP) distribution.
    !>
    !>  \interface{distUnifPar_type}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifPar, only: distUnifPar_type
    !>      type(distUnifPar_type) :: distUnifPar
    !>
    !>      distUnifPar = distUnifPar_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distUnifPar](@ref test_pm_distUnifPar)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distUnifPar_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distUnifPar_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF)
    !>  of the MultiVariate Uniform Parallelepiped (MVUP) Distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifPar](@ref pm_distUnifPar) for details of the definition of the PDF.
    !>
    !>  \param[in]  logLenEdge  :   The input scalar or `contiguous` array of shape `(1:ndim)`,
    !>                              of the same type and kind as the output `logPDF`,
    !>                              containing the natural logarithm of the lengths of the edges of the `ndim`-dimensional
    !>                              rectangular support of the distribution along each of the `ndim` axes.<br>
    !>                              <ol>
    !>                                  <li>    If `logLenEdge` is a scalar, the support of the distribution is assumed to be a `ndim`-dimensional hyper-cube.<br>
    !>                                          In such a case, the input optional argument `ndim` must be specified.<br>
    !>                                  <li>    If `logLenEdge` is a vector, the support of the distribution is assumed to be a hyper-rectangle.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input `repmat` is missing.)
    !>  \param[in]  ndim        :   The input positive scalar of type `integer` of default kind \IK, containing the number of dimensions of the domain of the distribution.<br>
    !>                              (**optional**. It must be present **if and only if** the input `logLenEdge` is present and is a scalar.)
    !>  \param[in]  repmat      :   The input square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the output `logPDF`,
    !>                              containing the [representative matrix](@ref pm_distUnifPar) \f$M_R\f$ of the parallelepiped.<br>
    !>                              (**optional**. It must be present **if and only if** the input `logLenEdge` is missing.)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the natural logarithm of the PDF of the distribution.
    !>
    !>  \interface{getUnifParLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifPar, only: getUnifParLogPDF
    !>
    !>      logPDF = getUnifParLogPDF(logLenEdge, ndim)        ! Hyper-cube support
    !>      logPDF = getUnifParLogPDF(logLenEdge(1:ndim))      ! Hyper-rectangle support
    !>      logPDF = getUnifParLogPDF(repmat(1:ndim, 1:ndim))  ! Hyper-parallelepiped support
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  The condition `size(repmat, 1) == size(repmat, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition \f$|\ms{repmat}| \neq 0\f$ must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input argument `ndim` is present.<br>
    !>
    !>  \example{getUnifParLogPDF}
    !>  \include{lineno} example/pm_distUnifPar/getUnifParLogPDF/main.F90
    !>  \compilef{getUnifParLogPDF}
    !>  \output{getUnifParLogPDF}
    !>  \include{lineno} example/pm_distUnifPar/getUnifParLogPDF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnifPar](@ref test_pm_distUnifPar)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The implementation of this generic interface with an input representative matrix `repmat` must be improved for better efficiency.<br>
    !>  The current implementation computes the corresponding Gramian matrix which may be slower than directly computing the determinant of
    !>  the representative matrix.<br>
    !>
    !>  \final{getUnifParLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getUnifParLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getUnifCubLogPDF_RK5(logLenEdge, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logLenEdge
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getUnifCubLogPDF_RK4(logLenEdge, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logLenEdge
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getUnifCubLogPDF_RK3(logLenEdge, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logLenEdge
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getUnifCubLogPDF_RK2(logLenEdge, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logLenEdge
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getUnifCubLogPDF_RK1(logLenEdge, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logLenEdge
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifRecLogPDF_RK5(logLenEdge) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: logLenEdge(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifRecLogPDF_RK4(logLenEdge) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: logLenEdge(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifRecLogPDF_RK3(logLenEdge) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: logLenEdge(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifRecLogPDF_RK2(logLenEdge) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: logLenEdge(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifRecLogPDF_RK1(logLenEdge) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: logLenEdge(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifParLogPDF_RK5(repmat) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: repmat(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifParLogPDF_RK4(repmat) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: repmat(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifParLogPDF_RK3(repmat) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: repmat(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifParLogPDF_RK2(repmat) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: repmat(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifParLogPDF_RK1(repmat) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: repmat(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a random vector from the \f$\ndim\f$-dimensional MultiVariate Uniform Parallelepiped (MVUP) Distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifPar](@ref pm_distUnifPar) for details of the definition of the PDF.
    !>
    !>  \param[in]  lb          :   The input scalar or `contiguous` vector of the same type and kind as the output argument `rand`.<br>
    !>                              <ol>
    !>                                  <li>    If `ub` is a scalar, then `lb` must be a scalar representing
    !>                                          the lower bound of the `ndim`-dimensional hyper-cube support of the distribution.
    !>                                  <li>    If `ub` is a vector, then `lb` must be a vector representing
    !>                                          the lower bound of the `ndim`-dimensional hyper-rectangle support of the distribution.
    !>                                  <li>    If `ub` is a matrix, then `lb` must be a vector representing
    !>                                          the origin of the coordinate system with respect to which the `ndim`-dimensional
    !>                                          hyper-parallelepiped support of the distribution is determined.
    !>                              </ol>
    !>                              (**optional**. default = `0`)
    !>  \param[in]  ub          :   The input scalar or `contiguous` vector/matrix of the same type and kind as the output argument `rand`.<br>
    !>                              <ol>
    !>                                  <li>    If `ub` is a scalar, it represents the upper bound of
    !>                                          the `ndim`-dimensional hyper-cube support of the distribution.<br>
    !>                                          In such a case, the input optional argument `ndim` must be specified.<br>
    !>                                  <li>    If `ub` is a vector of shape `(1:ndim)`, it represents the upper bounds of the
    !>                                          `ndim`-dimensional hyper-rectangle support of the distribution along each dimension.<br>
    !>                                  <li>    If `ub` is a square matrix of shape `(1:ndim, 1:ndim)`, it contains the representative
    !>                                          matrix of the `ndim`-dimensional hyper-parallelepiped support of the distribution.<br>
    !>                              </ol>
    !>  \param[in]  ndim        :   The input positive scalar of type `integer` of default kind \IK,
    !>                              containing the number of dimensions of the domain of the distribution.<br>
    !>                              (**optional**. It must be present **if and only if** the input `lb` and `ub` are scalar.)
    !>
    !>  \return
    !>  `rand`                  :   The output vector of shape `(1:ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing a random vector from the support of the distribution.
    !>
    !>  \interface{getUnifParRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifPar, only: getUnifParRand
    !>
    !>      rand(1:ndim) = getUnifParRand(ub(1:ndim))             ! Uniform Hyper-rectangle support
    !>      rand(1:ndim) = getUnifParRand(ub, ndim)               ! Uniform Hyper-cube support
    !>
    !>      rand(1:ndim) = getUnifParRand(lb(1:ndim), ub(1:ndim)) ! Uniform Hyper-rectangle support
    !>      rand(1:ndim) = getUnifParRand(lb, ub, ndim)           ! Uniform Hyper-cube support
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lb /= ub)` must hold for the corresponding input arguments.<br>
    !>  Although there is theoretically no limit on the possible values of `lb` and `ub` with respect to each other,
    !>  the condition `all(lb < ub)` is expected (but not checked) to hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \example{getUnifParRand}
    !>  \include{lineno} example/pm_distUnifPar/getUnifParRand/main.F90
    !>  \compilef{getUnifParRand}
    !>  \output{getUnifParRand}
    !>  \include{lineno} example/pm_distUnifPar/getUnifParRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distUnifPar/getUnifParRand/main.py
    !>  \vis
    !>  \image html example/pm_distUnifPar/getUnifParRand/getUnifParRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifPar](@ref test_pm_distUnifPar)
    !>
    !>  \final{getUnifParRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getUnifParRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifCubRandDU_RK5(ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifCubRandDU_RK4(ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifCubRandDU_RK3(ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifCubRandDU_RK2(ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifCubRandDU_RK1(ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifCubRandLU_RK5(lb, ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: lb, ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifCubRandLU_RK4(lb, ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: lb, ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifCubRandLU_RK3(lb, ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: lb, ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifCubRandLU_RK2(lb, ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: lb, ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifCubRandLU_RK1(lb, ub, ndim) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifCubRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: lb, ub
        real(RKG)                                           :: rand(ndim)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRecRandDU_RK5(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifRecRandDU_RK4(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifRecRandDU_RK3(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifRecRandDU_RK2(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRecRandDU_RK1(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifRecRandLU_RK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifRecRandLU_RK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifRecRandLU_RK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifRecRandLU_RK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifRecRandLU_RK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifRecRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifParRandDU_RK5(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifParRandDU_RK4(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifParRandDU_RK3(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifParRandDU_RK2(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifParRandDU_RK1(ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUnifParRandLU_RK5(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUnifParRandLU_RK4(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUnifParRandLU_RK3(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUnifParRandLU_RK2(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUnifParRandLU_RK1(lb, ub) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifParRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
        real(RKG)                                           :: rand(size(ub, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a random vector from the \f$\ndim\f$-dimensional MultiVariate Uniform Parallelepiped (MVUP) Distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifPar](@ref pm_distUnifPar) for details of the definition of the PDF.
    !>
    !>  \param[inout]   rand    :   The input/output vector of shape `(1:ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL.
    !>                              </ol>
    !>                              On input, it must contain a vector of randomly uniformly-distributed numbers in the range \f$[0, 1)\f$.<br>
    !>                              On output, it will contain a random vector from the support of the target distribution.<br>

    !>  \param[in]  lb          :   The input scalar or `contiguous` vector of the same type and kind as the output argument `rand`.<br>
    !>                              <ol>
    !>                                  <li>    If `ub` is a scalar, then `lb` must be a scalar representing
    !>                                          the lower bound of the `ndim`-dimensional hyper-cube support of the distribution.
    !>                                  <li>    If `ub` is a vector, then `lb` must be a vector representing
    !>                                          the lower bound of the `ndim`-dimensional hyper-rectangle support of the distribution.
    !>                                  <li>    If `ub` is a matrix, then `lb` must be a vector representing
    !>                                          the origin of the coordinate system with respect to which the `ndim`-dimensional
    !>                                          hyper-parallelepiped support of the distribution is determined.
    !>                              </ol>
    !>                              (**optional**. default = `0`)
    !>  \param[in]  ub          :   The input scalar or `contiguous` vector/matrix of the same type and kind as the output argument `rand`.<br>
    !>                              <ol>
    !>                                  <li>    If `ub` is a scalar, it must contain the upper bound of
    !>                                          the `ndim`-dimensional hyper-cube support of the distribution.<br>
    !>                                          In such a case, the input optional argument `ndim` must be specified.<br>
    !>                                  <li>    If `ub` is a vector of shape `(1:ndim)`, it must contain the upper bounds of the
    !>                                          `ndim`-dimensional hyper-rectangle support of the distribution along each dimension.<br>
    !>                                  <li>    If `ub` is a square matrix of shape `(1:ndim, 1:ndim)`, it must contain the representative
    !>                                          matrix of the `ndim`-dimensional hyper-parallelepiped support of the distribution.<br>
    !>                              </ol>
    !>
    !>  \interface{setUnifParRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifPar, only: setUnifParRand
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), ub)                              ! Uniform Hyper-cube support starting at the origin.
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), lb, ub)                          ! Uniform Hyper-cube support with lower bound `lb`.
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), ub(1:ndim))                      ! Uniform Hyper-rectangle support starting at the origin.
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), lb(1:ndim), ub(1:ndim))          ! Uniform Hyper-rectangle support with lower bound `lb`.
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), ub(1:ndim, 1:ndim))              ! Uniform Hyper-parallelepiped support starting at the origin.
    !>
    !>      call random_number(rand(1:ndim))
    !>      call setUnifParRand(rand(1:ndim), lb(1:ndim), ub(1:ndim, 1:ndim))  ! Uniform Hyper-parallelepiped support starting at `lb`.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0. <= rand .and. rand < 1.)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(rand) == size(lb(:)))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(rand) == shape(ub))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(lb /= ub(:))` must hold for the corresponding input arguments.<br>
    !>  Although there is theoretically no limit on the possible values of `lb` and `ub` with respect to each other,
    !>  the condition `all(lb(:) < ub(:))` is expected (but not checked) to hold for the corresponding input arguments.<br>
    !>  The input representative matrix of the parallelepiped support of the distribution `ub(:,:)` must be non-singular.<br>
    !>  In other words, the columns of the matrix must span the \f$\mathbb{R}^{\ndim}\f$ space.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \example{setUnifParRand}
    !>  \include{lineno} example/pm_distUnifPar/setUnifParRand/main.F90
    !>  \compilef{setUnifParRand}
    !>  \output{setUnifParRand}
    !>  \include{lineno} example/pm_distUnifPar/setUnifParRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distUnifPar/setUnifParRand/main.py
    !>  \vis
    !>  \image html example/pm_distUnifPar/setUnifParRand/setUnifParRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifPar](@ref test_pm_distUnifPar)
    !>
    !>  \final{setUnifParRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface setUnifParRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCubRandDU_RK5(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: ub
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCubRandDU_RK4(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: ub
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCubRandDU_RK3(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: ub
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCubRandDU_RK2(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: ub
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCubRandDU_RK1(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifCubRandLU_RK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: lb, ub
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifCubRandLU_RK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: lb, ub
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifCubRandLU_RK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: lb, ub
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifCubRandLU_RK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: lb, ub
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifCubRandLU_RK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifCubRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)                            :: lb, ub
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRecRandDU_RK5(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRecRandDU_RK4(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRecRandDU_RK3(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRecRandDU_RK2(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRecRandDU_RK1(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifRecRandLU_RK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifRecRandLU_RK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifRecRandLU_RK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifRecRandLU_RK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifRecRandLU_RK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifRecRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifParRandDU_RK5(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandDU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifParRandDU_RK4(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandDU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifParRandDU_RK3(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandDU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifParRandDU_RK2(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandDU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifParRandDU_RK1(rand, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandDU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: ub(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setUnifParRandLU_RK5(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setUnifParRandLU_RK4(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setUnifParRandLU_RK3(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setUnifParRandLU_RK2(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setUnifParRandLU_RK1(rand, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUnifParRandLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous            :: rand(:)
        real(RKG)   , intent(in)    , contiguous            :: lb(:), ub(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distUnifPar