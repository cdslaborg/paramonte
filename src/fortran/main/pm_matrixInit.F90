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
!>  This module contains procedures and generic interfaces relevant to generating and initializing matrices of arbitrary shapes `(:, :)`.<br>
!>
!>  \details
!>  The following guidelines describe the various existing methods within the ParaMonte library for initializing matrices.<br>
!>  <ul>
!>      <li>    Use the interface
!>              [getMatInit(..., subset = uppLowDia, ...)](@ref pm_matrixInit::getMatInit) or
!>              [setMatInit(..., subset = uppLowDia, ...)](@ref pm_matrixInit::setMatInit)
!>              to initialize the **upper, diagonal, and lower elements** of a **rectangular subset** of a new or existing matrix.<br>
!>              \htmlonly
!>                  <div style="display: block; margin-left: auto; margin-right: auto; padding: 0 4px 0 4px;">
!>                      <img src="pm_matrixInit@doffneg.png" style="width:30%;">
!>                      <img src="pm_matrixInit@doffzero.png" style="width:30%;">
!>                      <img src="pm_matrixInit@doffpos.png" style="width:30%;">
!>                  </div>
!>              \endhtmlonly
!>              <br>
!>      <li>    Use the interface
!>              [getMatInit(..., subset = uppLow, ...)](@ref pm_matrixInit::getMatInit) or
!>              [setMatInit(..., subset = uppLow, ...)](@ref pm_matrixInit::setMatInit)
!>              to initialize the **lower and upper elements** of a **rectangular subset** of a new or existing matrix.<br>
!>              \htmlonly
!>                  <div style="display: block; margin-left: auto; margin-right: auto; padding: 0 4px 0 4px;">
!>                      <img src="pm_matrixInitUppLow@doffneg.png" style="width:30%;">
!>                      <img src="pm_matrixInitUppLow@doffzero.png" style="width:30%;">
!>                      <img src="pm_matrixInitUppLow@doffpos.png" style="width:30%;">
!>                  </div>
!>              \endhtmlonly
!>              <br>
!>      <li>    Use the interface
!>              [getMatInit(..., subset = lowDia, ...)](@ref pm_matrixInit::getMatInit) or
!>              [setMatInit(..., subset = lowDia, ...)](@ref pm_matrixInit::setMatInit)
!>              to initialize the **diagonal and lower elements** of a **rectangular subset** of a new or existing matrix.<br>
!>              \htmlonly
!>                  <div style="display: block; margin-left: auto; margin-right: auto; padding: 0 4px 0 4px;">
!>                      <img src="pm_matrixInitLowDia@doffpos.png" style="width:30%;">
!>                  </div>
!>              \endhtmlonly
!>              <br>
!>      <li>    Use the interface
!>              [getMatInit(..., subset = uppDia, ...)](@ref pm_matrixInit::getMatInit) or
!>              [setMatInit(..., subset = uppDia, ...)](@ref pm_matrixInit::setMatInit)
!>              to initialize the **diagonal and upper elements** of a **rectangular subset** of a new or existing matrix.<br>
!>              \htmlonly
!>                  <div style="display: block; margin-left: auto; margin-right: auto; padding: 0 4px 0 4px;">
!>                      <img src="pm_matrixInitUppDia@doffneg.png" style="width:30%;">
!>                  </div>
!>              \endhtmlonly
!>              <br>
!>      <li>    Use the interface
!>              [getMatInit(..., subset = vdia, ...)](@ref pm_matrixInit::getMatInit) or
!>              [setMatInit(..., subset = vdia, ...)](@ref pm_matrixInit::setMatInit)
!>              to initialize the **diagonal elements** of a **rectangular subset** of a new or existing matrix.<br>
!>              \htmlonly
!>                  <div style="display: block; margin-left: auto; margin-right: auto; padding: 0 4px 0 4px;">
!>                      <img src="pm_matrixInitDia.png" style="width:30%;">
!>                  </div>
!>              \endhtmlonly
!>              <br>
!>  </ul>
!>
!>  \lapack{3.11}
!>  `SLASET`, `DLASET`, `CLASET`, `ZLASET`.<br>
!>
!>  \see
!>  [pm_arrayInit](@ref pm_arrayInit)<br>
!>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
!>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
!>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getMatInit_vs_setMatInit, The runtime performance of [getMatInit](@ref pm_matrixInit::getMatInit) vs. [setMatInit](@ref pm_matrixInit::setMatInit)}
!>  \include{lineno} benchmark/pm_matrixInit/getMatInit_vs_setMatInit/main.F90
!>  \compilefb{getMatInit_vs_setMatInit}
!>  \postprocb{getMatInit_vs_setMatInit}
!>  \include{lineno} benchmark/pm_matrixInit/getMatInit_vs_setMatInit/main.py
!>  \visb{getMatInit_vs_setMatInit}
!>  \image html benchmark/pm_matrixInit/getMatInit_vs_setMatInit/benchmark.getMatInit_vs_setMatInit.runtime.png width=1000
!>  \image html benchmark/pm_matrixInit/getMatInit_vs_setMatInit/benchmark.getMatInit_vs_setMatInit.runtime.ratio.png width=1000
!>  \moralb{getMatInit_vs_setMatInit}
!>      -#  The procedures under the generic interface [getMatInit](@ref pm_matrixInit::getMatInit) are functions while
!>          the procedures under the generic interface [setMatInit](@ref pm_matrixInit::setMatInit) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs less efficiently than the subroutine interface.<br>
!>          Note that this benchmark does not even include the cost of repeated reallcations, that is, the allocation of matrix happens only once in all tests prior to each benchmark.<br>
!>      -#  Furthermore, the `getMatInit()` implementation with Fortran-intrinsic `reshape()` and implicit-loop appears to be significantly
!>          slower than both the subroutine and function implementations.<br>
!>
!>  \test
!>  [test_pm_matrixInit](@ref test_pm_matrixInit)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixInit

    use pm_kind, only: SK, IK
    use pm_matrixSubset, only: dia, dia_type
    use pm_matrixSubset, only: low, low_type
    use pm_matrixSubset, only: upp, upp_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: uppLow, uppLow_type
    use pm_matrixSubset, only: uppLowDia, uppLowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixInit"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a matrix of shape `(shape(1), shape(2))` with the upper/lower triangle
    !>  and the diagonal elements of the matrix set to the corresponding requested input values.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_matrixInit](@ref pm_matrixInit) for illustrations of possible initialization formats.<br>
    !>
    !>  \param[in]      shape   :   The input `contiguous` vector of size `2` of type `integer` of default kind \IK,
    !>                              containing the shape of the output matrix (as returned by the Fortran intrinsic function `shape()`).<br>
    !>                              By definition, both elements of `shape` must be non-negative.<br>
    !>  \param[in]      subset  :   The input argument that can be either,<br>
    !>                              <ol>
    !>                                  <li>    the constant [upp](@ref pm_matrixSubset::upp), signifying the initialization of the upper-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [low](@ref pm_matrixSubset::low), signifying the initialization of the lower-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [dia](@ref pm_matrixSubset::dia), signifying the initialization of the diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia), signifying the initialization of the lower-triangle and diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia), signifying the initialization of the upper-triangle and diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppLow](@ref pm_matrixSubset::uppLow), signifying the initialization of the upper/lower-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppLowDia](@ref pm_matrixSubset::uppLowDia), signifying the initialization of the upper/lower-triangle and diagonal elements of the requested subset of the matrix.
    !>                              </ol>
    !>                              This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                              Note that setting only the lower/upper triangles of a matrix can be readily converted to the problem of setting the upper/lower triangle and diagonal elements of a matrix.<br>
    !>  \param[in]      vupp    :   The input scalar of the same type and kind as `mat`,<br>
    !>                              containing the value for all upper-triangular elements of the output matrix.<br>
    !>                              (**optional**. It can be present only if the requested input `subset` format is compatible.)
    !>  \param[in]      vlow    :   The input scalar of the same type and kind as `mat`,<br>
    !>                              containing the value for all lower-triangular elements of the output matrix.<br>
    !>                              (**optional**. It can be present only if the requested input `subset` format is compatible.)
    !>  \param[in]      vdia    :   The input scalar or `contiguous` array of rank `1` of the same type and kind as `mat` containing the diagonal elements of the output matrix.<br>
    !>                              (**optional**. It can be present only if the requested input `subset` format is compatible.)
    !>  \param[in]      nrow    :   The input non-negative `integer` of default kind \IK representing the number of rows of the initialized block of the matrix.<br>
    !>                              Setting `nrow = size(mat,1), roff = 0` is equivalent to initializing all rows of the matrix.<br>
    !>                              (**optional**, default = `size(mat,1)`)
    !>  \param[in]      ncol    :   The input non-negative `integer` of default kind \IK representing the number of columns of the initialized block of the matrix.<br>
    !>                              Setting `ncol = size(mat,2), roff = 0` is equivalent to initializing all columns of the matrix.<br>
    !>                              (**optional**, default = `size(mat,2)`)
	!>  \param[in]      ndia	:   The input non-negative integer of default kind \IK representing the number of diagonal elements to write to the initialized block of the matrix.<br>
    !>                              It can be present **only if** the input argument `subset` is set to [dia](@ref pm_matrixSubset::dia) and the input argument `vdia` is a scalar.<br>
    !>                              If `vdia` is a vector of rank `1`, then `ndia` is automatically set to the size of the input vector `vdia`.<br>
    !>                              (**optional**. default = `min(shape(1) - roff, shape(2) - coff)`.)
    !>  \param[in]      roff    :   The input non-negative `integer` of default kind \IK representing the **offset**
    !>                              of the top-left corner of the initialization block from the first row of the matrix.<br>
    !>                              The initialization of the matrix will start at row `1 + roff`.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]      coff    :   The input non-negative `integer` of default kind \IK representing the **offset**
    !>                              of the top-left corner of the initialization block from the first column of the matrix.<br>
    !>                              The initialization of the matrix will start at column `1 + coff`.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]      doff    :   The input `integer` of default kind \IK representing the **offset**
    !>                              of the diagonal of the initialization block from the top-left corner of the initialization block.<br>
    !>                              <ul>
    !>                                  <li>    Setting `doff > 0` implies a diagonal start with column offset `doff` from the top-left corner of the initialization block.<br>
    !>                                  <li>    Setting `doff < 0` implies a diagonal start with row offset `-doff` from the top-left corner of the initialization block.<br>
    !>                                  <li>    Setting `doff = 0` implies the same diagonal start as the top-left corner of the initialization block.<br>
    !>                              </ul>
    !>                              (**optional**, default = `0`. It can be present **only if** the input arguments `vupp` or `vlow` are present.)
    !>
    !>  \return
    !>  `mat`                   :   The output `contiguous` matrix of arbitrary shape `(:, :)` of,
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL with arbitrary `len` type parameter, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ol>
    !>                              On output, the upper/lower triangle and diagonal elements of `mat` (as specified by the input argument `subset`) are set to the corresponding input values.<br>
    !>                              All other matrix elements remain uninitialized on output.<br>
    !>
    !>  \interface{getMatInit}
    !>  \code{.F90}
    !>
    !>      use pm_matrixInit, only: getMatInit, dia, uppDia, lowDia, uppLow, uppLowDia
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vupp, roff = roff, coff = coff, doff = doff) ! subset = upp: Initialize upper-triangular subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vlow, roff = roff, coff = coff, doff = doff) ! subset = low: Initialize lower-triangular subset of matrix.
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vdia   , ndia = ndia, roff = roff, coff = coff) ! Initialize diagonal subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vdia(:)             , roff = roff, coff = coff) ! Initialize diagonal subset of matrix.
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vlow, vdia   , nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vlow, vdia(:), nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vdia   , nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vdia(:), nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vlow, nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-lower subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vlow, nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-lower subset of matrix.
    !>
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vlow, vdia   , nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>      mat(@shape) = getMatInit(shape, subset, vupp, vlow, vdia(:), nrow = nrow, ncol = ncol, roff = roff, coff = coff, doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings associated with [setMatInit](@ref pm_matrixInit::setMatInit) also apply to this generic interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `SLASET`, `DLASET`, `CLASET`, `ZLASET`.<br>
    !>
    !>  \see
    !>  [setMatInit](@ref pm_matrixInit::setMatInit)<br>
    !>  [pm_arrayInit](@ref pm_arrayInit)<br>
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>
    !>  \example{getMatInit}
    !>  \include{lineno} example/pm_matrixInit/getMatInit/main.F90
    !>  \compilef{getMatInit}
    !>  \output{getMatInit}
    !>  \include{lineno} example/pm_matrixInit/getMatInit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixInit](@ref test_pm_matrixInit)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getMatInitXXD_D2XX0_SK5()
    !>                ||| ||||| |||
    !>                ||| ||||| Output type/kind
    !>                ||| ||||The rank of `vdia` argument. `X` implies missing argument.
    !>                ||| |||The rank of `vlow` argument. `X` implies missing argument.
    !>                ||| ||The rank of `vupp` argument. `X` implies missing argument.
    !>                ||| |The rank of the output `mat`.
    !>                ||| `D` stands for dimension.
    !>                ||The `vdia` argument presence: `D` implies presence. `X` implies missing argument.
    !>                |The `vlow` argument presence: `L` implies presence. `X` implies missing argument.
    !>                The `vupp` argument presence: `U` implies presence. `X` implies missing argument.
    !>  \endcode
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.8.0 20221119}
    !>  \desc
    !>  \ifort thows a strange error for compiling this module when [dia](@ref pm_matrixSubset::dia) constant is `use`d within this module.<br>
    !>  Apparently, \ifort confuses this imported constant with the dummy procedure argument
    !>  of the same type in procedures named `getMatInitXXD_D2XX0*` or `getMatInitXXD_D2XX1`.<br>
    !>  \remedy
    !>  There does not appear to exist any way to resolve this intel error currently,
    !>  other than not importing the constant [dia](@ref pm_matrixSubset::dia) in this module.<br>
    !>  Instead, one has to use [dia_type()](@ref pm_matrixSubset::dia_type()) to call the procedures when resolving the generic interface.<br>
    !>  As such, the argument names `upp, low, dia` were renamed to `vupp, vlow, vdia` to emphasize
    !>  the `value` attribute of these arguments and to resolve the Intel name-conflict bug.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  This generic interface should be extended to matrices of different packing formats besides the default.
    !>
    !>  \final{getMatInit}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    interface getMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMatInitXXD_D2XX0_SK5(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitXXD_D2XX0_SK4(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitXXD_D2XX0_SK3(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitXXD_D2XX0_SK2(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitXXD_D2XX0_SK1(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitXXD_D2XX0_IK5(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitXXD_D2XX0_IK4(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitXXD_D2XX0_IK3(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitXXD_D2XX0_IK2(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitXXD_D2XX0_IK1(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitXXD_D2XX0_LK5(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitXXD_D2XX0_LK4(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitXXD_D2XX0_LK3(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitXXD_D2XX0_LK2(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitXXD_D2XX0_LK1(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitXXD_D2XX0_CK5(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitXXD_D2XX0_CK4(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitXXD_D2XX0_CK3(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitXXD_D2XX0_CK2(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitXXD_D2XX0_CK1(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitXXD_D2XX0_RK5(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitXXD_D2XX0_RK4(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitXXD_D2XX0_RK3(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitXXD_D2XX0_RK2(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitXXD_D2XX0_RK1(shape, subset, vdia, ndia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: ndia
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMatInitXXD_D2XX1_SK5(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitXXD_D2XX1_SK4(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitXXD_D2XX1_SK3(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitXXD_D2XX1_SK2(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitXXD_D2XX1_SK1(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitXXD_D2XX1_IK5(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitXXD_D2XX1_IK4(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitXXD_D2XX1_IK3(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitXXD_D2XX1_IK2(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitXXD_D2XX1_IK1(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitXXD_D2XX1_LK5(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitXXD_D2XX1_LK4(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitXXD_D2XX1_LK3(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitXXD_D2XX1_LK2(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitXXD_D2XX1_LK1(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitXXD_D2XX1_CK5(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitXXD_D2XX1_CK4(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitXXD_D2XX1_CK3(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitXXD_D2XX1_CK2(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitXXD_D2XX1_CK1(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitXXD_D2XX1_RK5(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitXXD_D2XX1_RK4(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitXXD_D2XX1_RK3(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitXXD_D2XX1_RK2(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitXXD_D2XX1_RK1(shape, subset, vdia, roff, coff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXXD_D2XX1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: roff, coff
        type(dia_type)              , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
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

#if SK5_ENABLED
    PURE module function getMatInitXLD_D2X00_SK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitXLD_D2X00_SK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitXLD_D2X00_SK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitXLD_D2X00_SK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitXLD_D2X00_SK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitXLD_D2X00_IK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitXLD_D2X00_IK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitXLD_D2X00_IK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitXLD_D2X00_IK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitXLD_D2X00_IK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitXLD_D2X00_LK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitXLD_D2X00_LK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitXLD_D2X00_LK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitXLD_D2X00_LK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitXLD_D2X00_LK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitXLD_D2X00_CK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitXLD_D2X00_CK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitXLD_D2X00_CK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitXLD_D2X00_CK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitXLD_D2X00_CK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitXLD_D2X00_RK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitXLD_D2X00_RK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitXLD_D2X00_RK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitXLD_D2X00_RK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitXLD_D2X00_RK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X00_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMatInitXLD_D2X01_SK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitXLD_D2X01_SK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitXLD_D2X01_SK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitXLD_D2X01_SK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitXLD_D2X01_SK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitXLD_D2X01_IK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitXLD_D2X01_IK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitXLD_D2X01_IK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitXLD_D2X01_IK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitXLD_D2X01_IK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitXLD_D2X01_LK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitXLD_D2X01_LK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitXLD_D2X01_LK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitXLD_D2X01_LK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitXLD_D2X01_LK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitXLD_D2X01_CK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitXLD_D2X01_CK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitXLD_D2X01_CK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitXLD_D2X01_CK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitXLD_D2X01_CK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitXLD_D2X01_RK5(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitXLD_D2X01_RK4(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitXLD_D2X01_RK3(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitXLD_D2X01_RK2(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitXLD_D2X01_RK1(shape, subset, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitXLD_D2X01_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(lowDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
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

#if SK5_ENABLED
    PURE module function getMatInitUXD_D20X0_SK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitUXD_D20X0_SK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitUXD_D20X0_SK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitUXD_D20X0_SK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitUXD_D20X0_SK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitUXD_D20X0_IK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitUXD_D20X0_IK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitUXD_D20X0_IK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitUXD_D20X0_IK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitUXD_D20X0_IK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitUXD_D20X0_LK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitUXD_D20X0_LK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitUXD_D20X0_LK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitUXD_D20X0_LK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitUXD_D20X0_LK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitUXD_D20X0_CK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitUXD_D20X0_CK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitUXD_D20X0_CK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitUXD_D20X0_CK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitUXD_D20X0_CK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitUXD_D20X0_RK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitUXD_D20X0_RK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitUXD_D20X0_RK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitUXD_D20X0_RK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitUXD_D20X0_RK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMatInitUXD_D20X1_SK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitUXD_D20X1_SK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitUXD_D20X1_SK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitUXD_D20X1_SK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitUXD_D20X1_SK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vupp
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitUXD_D20X1_IK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitUXD_D20X1_IK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitUXD_D20X1_IK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitUXD_D20X1_IK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitUXD_D20X1_IK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vupp
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitUXD_D20X1_LK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitUXD_D20X1_LK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitUXD_D20X1_LK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitUXD_D20X1_LK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitUXD_D20X1_LK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vupp
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitUXD_D20X1_CK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitUXD_D20X1_CK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitUXD_D20X1_CK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitUXD_D20X1_CK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitUXD_D20X1_CK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vupp
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitUXD_D20X1_RK5(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitUXD_D20X1_RK4(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitUXD_D20X1_RK3(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitUXD_D20X1_RK2(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitUXD_D20X1_RK1(shape, subset, vupp, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitUXD_D20X1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vupp
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppDia_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
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

#if SK5_ENABLED
    PURE module function getMatInitULX_D200X_SK5(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        character(len(vupp,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitULX_D200X_SK4(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        character(len(vupp,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitULX_D200X_SK3(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        character(len(vupp,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitULX_D200X_SK2(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        character(len(vupp,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitULX_D200X_SK1(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        character(len(vupp,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitULX_D200X_IK5(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitULX_D200X_IK4(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitULX_D200X_IK3(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitULX_D200X_IK2(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitULX_D200X_IK1(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitULX_D200X_LK5(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitULX_D200X_LK4(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitULX_D200X_LK3(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitULX_D200X_LK2(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitULX_D200X_LK1(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitULX_D200X_CK5(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitULX_D200X_CK4(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitULX_D200X_CK3(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitULX_D200X_CK2(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitULX_D200X_CK1(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitULX_D200X_RK5(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitULX_D200X_RK4(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitULX_D200X_RK3(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitULX_D200X_RK2(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitULX_D200X_RK1(shape, subset, vupp, vlow, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULX_D200X_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vupp, vlow
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLow_type)           , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
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

#if SK5_ENABLED
    PURE module function getMatInitULD_D2000_SK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitULD_D2000_SK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitULD_D2000_SK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitULD_D2000_SK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitULD_D2000_SK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitULD_D2000_IK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitULD_D2000_IK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitULD_D2000_IK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitULD_D2000_IK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitULD_D2000_IK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitULD_D2000_LK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitULD_D2000_LK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitULD_D2000_LK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitULD_D2000_LK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitULD_D2000_LK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitULD_D2000_CK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitULD_D2000_CK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitULD_D2000_CK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitULD_D2000_CK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitULD_D2000_CK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitULD_D2000_RK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitULD_D2000_RK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitULD_D2000_RK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitULD_D2000_RK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitULD_D2000_RK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2000_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in)                :: vdia
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMatInitULD_D2001_SK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMatInitULD_D2001_SK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMatInitULD_D2001_SK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMatInitULD_D2001_SK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMatInitULD_D2001_SK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                :: vupp, vlow
        character(*,SKC)            , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        character(len(vdia,IK),SKC)                             :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMatInitULD_D2001_IK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMatInitULD_D2001_IK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMatInitULD_D2001_IK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMatInitULD_D2001_IK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMatInitULD_D2001_IK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)                :: vupp, vlow
        integer(IKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        integer(IKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMatInitULD_D2001_LK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMatInitULD_D2001_LK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMatInitULD_D2001_LK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMatInitULD_D2001_LK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMatInitULD_D2001_LK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)                :: vupp, vlow
        logical(LKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        logical(LKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInitULD_D2001_CK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInitULD_D2001_CK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInitULD_D2001_CK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInitULD_D2001_CK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInitULD_D2001_CK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)                :: vupp, vlow
        complex(CKC)                , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        complex(CKC)                                            :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInitULD_D2001_RK5(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInitULD_D2001_RK4(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInitULD_D2001_RK3(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInitULD_D2001_RK2(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInitULD_D2001_RK1(shape, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff) result(mat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInitULD_D2001_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)                :: vupp, vlow
        real(RKC)                   , intent(in), contiguous    :: vdia(:)
        integer(IK)                 , intent(in)                :: shape(2)
        integer(IK)                 , intent(in), optional      :: nrow, ncol
        integer(IK)                 , intent(in), optional      :: roff, coff, doff
        type(uppLowDia_type)        , intent(in)                :: subset
        real(RKC)                                               :: mat(shape(1), shape(2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the upper/lower triangle and the diagonal elements of the input matrix of arbitrary shape `(:,:)` to the requested input values.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_matrixInit](@ref pm_matrixInit) for illustrations of possible initialization formats.<br>
    !>
    !>  \param[inout]   mat     :   The input/output `contiguous` matrix of arbitrary shape `(:, :)` of,
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL with arbitrary `len` type parameter, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ol>
    !>                              On output, the upper/lower triangle and diagonal elements of `mat` (as specified by the input argument `subset`) are set to the corresponding input values.<br>
    !>                              All other matrix elements remain intact on output.<br>
    !>  \param[in]      subset  :   The input argument that can be either,
    !>                              <ol>
    !>                                  <li>    the constant [upp](@ref pm_matrixSubset::upp), signifying the initialization of the upper-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [low](@ref pm_matrixSubset::low), signifying the initialization of the lower-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [dia](@ref pm_matrixSubset::dia), signifying the initialization of the diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia), signifying the initialization of the lower-triangle and diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia), signifying the initialization of the upper-triangle and diagonal elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppLow](@ref pm_matrixSubset::uppLow), signifying the initialization of the upper/lower-triangle elements of the requested subset of the matrix.
    !>                                  <li>    the constant [uppLowDia](@ref pm_matrixSubset::uppLowDia), signifying the initialization of the upper/lower-triangle and diagonal elements of the requested subset of the matrix.
    !>                              </ol>
    !>                              This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                              Note that setting only the lower/upper triangles of a matrix can be readily converted to the problem of setting the upper/lower triangle and diagonal elements of a matrix.<br>
    !>  \param[in]      vupp    :   The input scalar of the same type and kind as `mat`,<br>
    !>                              containing the value for all upper-triangular elements of the output matrix.<br>
    !>                              (**optional**, It must be present **if and only if** the input argument `subset` specifies an update to the upper-triangle elements.)
    !>  \param[in]      vlow    :   The input scalar of the same type and kind as `mat`,<br>
    !>                              containing the value for all lower-triangular elements of the output matrix.<br>
    !>                              (**optional**, It must be present **if and only if** the input argument `subset` specifies an update to the lower-triangle elements.)
    !>  \param[in]      vdia    :   The input scalar or `contiguous` array of rank `1` of the same type and kind as `mat` containing the diagonal elements of the output matrix.<br>
    !>                              (**optional**, It must be present **if and only if** the input argument `subset` specifies an update to the diagonal elements.)
    !>  \param[in]      nrow    :   The input non-negative `integer` of default kind \IK representing the number of rows of the initialized block of the matrix.<br>
    !>                              Setting `nrow = size(mat,1), roff = 0` is equivalent to initializing all rows of the matrix.<br>
    !>                              (**optional**, default = `size(mat, 1)`. It must be present **if and only if** any of the input arguments `ncol`, `roff`, or `coff` are also present.)
    !>  \param[in]      ncol    :   The input non-negative `integer` of default kind \IK representing the number of columns of the initialized block of the matrix.<br>
    !>                              Setting `ncol = size(mat,2), roff = 0` is equivalent to initializing all columns of the matrix.<br>
    !>                              (**optional**, default = `size(mat, 2)`. It must be present **if and only if** any of the input arguments `nrow`, `roff`, or `coff` are also present.)
	!>  \param[in]      ndia	:   The input non-negative integer of default kind \IK representing the number of diagonal elements to write to the initialized block of the matrix.<br>
    !>                              It can be present **only if** the input argument `subset` is set to [dia](@ref pm_matrixSubset::dia) and the input argument `vdia` is a scalar.<br>
    !>                              If `vdia` is a vector of rank `1`, then `ndia` is automatically set to the size of the input vector `vdia`.<br>
    !>                              (**optional**. default = `minval(shape(mat) - [roff, coff])` if `vdia` is a scalar, otherwise, `size(vdia(:))`.)
    !>  \param[in]      roff    :   The input non-negative `integer` of default kind \IK representing the **offset**
    !>                              of the top-left corner of the initialization block from the first row of the matrix.<br>
    !>                              The initialization of the matrix will start at row `1 + roff`.<br>
    !>                              (**optional**, default = `0`. It must be present **if and only if** any of the input arguments `nrow`, `ncol`, `ndia`, or `coff` are also present.)
    !>  \param[in]      coff    :   The input non-negative `integer` of default kind \IK representing the **offset**
    !>                              of the top-left corner of the initialization block from the first column of the matrix.<br>
    !>                              The initialization of the matrix will start at column `1 + coff`.<br>
    !>                              (**optional**, default = `0`. It must be present **if and only if** any of the input arguments `nrow`, `ncol`, `ndia`, or `roff` are also present.)
    !>  \param[in]      doff    :   The input `integer` of default kind \IK representing the **offset**
    !>                              of the diagonal of the initialization block from the top-left corner of the initialization block.<br>
    !>                              <ul>
    !>                                  <li>    Setting `doff > 0` implies a diagonal start with column offset `doff` from the top-left corner of the initialization block.<br>
    !>                                  <li>    Setting `doff < 0` implies a diagonal start with row offset `-doff` from the top-left corner of the initialization block.<br>
    !>                                  <li>    Setting `doff = 0` implies the same diagonal start as the top-left corner of the initialization block.<br>
    !>                              </ul>
    !>                              (**optional**, default = `0`. It can be present **only if** the input arguments `vupp` or `vlow` are present.)
    !>
    !>  \interface{setMatInit}
    !>  \code{.F90}
    !>
    !>      use pm_matrixInit, only: setMatInit, dia, uppDia, lowDia, uppLow, uppLowDia
    !>
    !>      ! implicit interface
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, doff = doff) ! Initialize upper-triangular subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vlow, doff = doff) ! Initialize lower-triangular subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vdia    , ndia) ! Initialize diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vdia          ) ! Initialize diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vdia(:)       ) ! Initialize diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vlow, vdia   , doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vlow, vdia(:), doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vdia   , doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vdia(:), doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, doff = doff) ! Initialize upper-lower subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, doff = doff) ! Initialize upper-lower subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, vdia   , doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, vdia(:), doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>
    !>      ! explicit interface
    !>
    !>      call setMatInit(mat(:,:), subset, vdia    , ndia  , roff, coff) ! Initialize diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vdia            , roff, coff) ! Initialize diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vdia(:)         , roff, coff) ! Initialize diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, nrow, ncol, roff, coff, doff = doff) ! Initialize upper-triangular subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vlow, nrow, ncol, roff, coff, doff = doff) ! Initialize lower-triangular subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vlow, vdia   , nrow, ncol, roff, coff, doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vlow, vdia(:), nrow, ncol, roff, coff, doff = doff) ! Initialize lower-diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vdia   , nrow, ncol, roff, coff, doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vdia(:), nrow, ncol, roff, coff, doff = doff) ! Initialize upper-diagonal subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, nrow, ncol, roff, coff, doff = doff) ! Initialize upper-lower subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, nrow, ncol, roff, coff, doff = doff) ! Initialize upper-lower subset of matrix.
    !>
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, vdia   , nrow, ncol, roff, coff, doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>      call setMatInit(mat(:,:), subset, vupp, vlow, vdia(:), nrow, ncol, roff, coff, doff = doff) ! Initialize upper-lower-diagonal subset of matrix.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `ndia >= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `nrow >= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `ncol >= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `roff >= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `coff >= 0` must hold for the corresponding input arguments.<br>
    !>  The condition `nrow + roff <= size(mat, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ncol + coff <= size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition `-nrow < doff .and. doff < ncol` must hold for the corresponding input arguments.<br>
    !>  The condition `ndia <= min(size(mat, 1) - roff, size(mat, 2) - coff)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(vdia(:)) == merge(min(nrow + doff, ncol), min(nrow, ncol - doff), doff < 0)` must hold for the corresponding input arguments.<br>
    !>  The condition `len(vupp) <= len(mat)` must hold for the corresponding arguments of type `character` of any kinds.<br>
    !>  The condition `len(vlow) <= len(mat)` must hold for the corresponding arguments of type `character` of any kinds.<br>
    !>  The condition `len(vdia) <= len(mat)` must hold for the corresponding arguments of type `character` of any kinds.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `SLASET`, `DLASET`, `CLASET`, `ZLASET`.<br>
    !>
    !>  \see
    !>  [getMatInit](@ref pm_matrixInit::getMatInit)<br>
    !>  [pm_arrayInit](@ref pm_arrayInit)<br>
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>
    !>  \example{setMatInit}
    !>  \include{lineno} example/pm_matrixInit/setMatInit/main.F90
    !>  \compilef{setMatInit}
    !>  \output{setMatInit}
    !>  \include{lineno} example/pm_matrixInit/setMatInit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixInit](@ref test_pm_matrixInit)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setMatInitXXD_D2XX0_SK5
    !>                ||| ||||| |||
    !>                ||| ||||| Output type/kind
    !>                ||| ||||The rank of `vdia` argument. `X` implies missing argument.
    !>                ||| ||||`F` implies filling the full diagonal scalar `dia` with `ndia` input argument missing.
    !>                ||| |||The rank of `vlow` argument. `X` implies missing argument.
    !>                ||| ||The rank of `vupp` argument. `X` implies missing argument.
    !>                ||| |The rank of the output `mat`.
    !>                ||| `D` stands for dimension.
    !>                ||The `vdia` argument presence: `D` implies presence. `X` implies missing argument.
    !>                |The `vlow` argument presence: `L` implies presence. `X` implies missing argument.
    !>                The `vupp` argument presence: `U` implies presence. `X` implies missing argument.
    !>  \endcode
    !>
    !>  \todo
    !>  \pvhigh
    !>  This generic interface should be extended to matrices of different packing formats besides the default.
    !>
    !>  \final{setMatInit}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! EXP_XLX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_SK5(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_SK4(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_SK3(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_SK2(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_SK1(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_IK5(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_IK4(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_IK3(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_IK2(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_IK1(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_LK5(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_LK4(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_LK3(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_LK2(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_LK1(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_CK5(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_CK4(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_CK3(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_CK2(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_CK1(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_RK5(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_RK4(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_RK3(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_RK2(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLX_D2X0X_RK1(mat, subset, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLX_D2X0X_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_UXX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_SK5(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_SK4(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_SK3(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_SK2(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_SK1(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_IK5(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_IK4(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_IK3(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_IK2(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_IK1(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_LK5(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_LK4(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_LK3(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_LK2(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_LK1(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_CK5(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_CK4(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_CK3(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_CK2(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_CK1(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_RK5(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_RK4(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_RK3(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_RK2(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXX_D20XX_RK1(mat, subset, vupp, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXX_D20XX_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_XXD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_SK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_SK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_SK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_SK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_SK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_IK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_IK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_IK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_IK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_IK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_LK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_LK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_LK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_LK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_LK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_CK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_CK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_CK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_CK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_CK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_RK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_RK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_RK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_RK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XXF_RK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XXF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_SK5(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_SK4(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_SK3(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_SK2(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_SK1(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_IK5(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_IK4(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_IK3(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_IK2(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_IK1(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_LK5(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_LK4(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_LK3(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_LK2(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_LK1(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_CK5(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_CK4(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_CK3(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_CK2(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_CK1(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_RK5(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_RK4(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_RK3(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_RK2(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX0_RK1(mat, subset, vdia, ndia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_SK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_SK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_SK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_SK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_SK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_IK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_IK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_IK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_IK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_IK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_LK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_LK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_LK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_LK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_LK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_CK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_CK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_CK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_CK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_CK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_RK5(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_RK4(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_RK3(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_RK2(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XXD_D2XX1_RK1(mat, subset, vdia, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XXD_D2XX1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        integer(IK)                     , intent(in)                    :: roff, coff
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_XLD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_SK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_SK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_SK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_SK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_SK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_IK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_IK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_IK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_IK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_IK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_LK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_LK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_LK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_LK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_LK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_CK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_CK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_CK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_CK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_CK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_RK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_RK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_RK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_RK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X00_RK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X00_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_SK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_SK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_SK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_SK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_SK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_IK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_IK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_IK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_IK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_IK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_LK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_LK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_LK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_LK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_LK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_CK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_CK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_CK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_CK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_CK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_RK5(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_RK4(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_RK3(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_RK2(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_XLD_D2X01_RK1(mat, subset, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_XLD_D2X01_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_UXD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_SK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_SK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_SK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_SK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_SK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_IK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_IK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_IK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_IK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_IK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_LK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_LK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_LK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_LK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_LK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_CK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_CK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_CK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_CK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_CK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_RK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_RK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_RK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_RK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X0_RK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X0_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_SK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_SK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_SK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_SK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_SK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_IK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_IK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_IK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_IK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_IK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_LK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_LK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_LK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_LK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_LK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_CK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_CK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_CK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_CK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_CK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_RK5(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_RK4(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_RK3(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_RK2(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_UXD_D20X1_RK1(mat, subset, vupp, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_UXD_D20X1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_ULX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_SK5(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_SK4(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_SK3(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_SK2(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_SK1(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_IK5(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_IK4(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_IK3(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_IK2(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_IK1(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_LK5(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_LK4(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_LK3(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_LK2(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_LK1(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_CK5(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_CK4(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_CK3(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_CK2(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_CK1(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_RK5(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_RK4(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_RK3(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_RK2(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULX_D200X_RK1(mat, subset, vupp, vlow, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULX_D200X_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! EXP_ULD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_SK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_SK5
#endif
        use pm_kind, only: SKC => SK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_SK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_SK4
#endif
        use pm_kind, only: SKC => SK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_SK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_SK3
#endif
        use pm_kind, only: SKC => SK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_SK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_SK2
#endif
        use pm_kind, only: SKC => SK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_SK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_SK1
#endif
        use pm_kind, only: SKC => SK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_IK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_IK5
#endif
        use pm_kind, only: IKC => IK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_IK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_IK4
#endif
        use pm_kind, only: IKC => IK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_IK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_IK3
#endif
        use pm_kind, only: IKC => IK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_IK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_IK2
#endif
        use pm_kind, only: IKC => IK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_IK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_IK1
#endif
        use pm_kind, only: IKC => IK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_LK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_LK5
#endif
        use pm_kind, only: LKC => LK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_LK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_LK4
#endif
        use pm_kind, only: LKC => LK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_LK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_LK3
#endif
        use pm_kind, only: LKC => LK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_LK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_LK2
#endif
        use pm_kind, only: LKC => LK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_LK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_LK1
#endif
        use pm_kind, only: LKC => LK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_CK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_CK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_CK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_CK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_CK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_RK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_RK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_RK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_RK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2000_RK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2000_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_SK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_SK5
#endif
        use pm_kind, only: SKC => SK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_SK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_SK4
#endif
        use pm_kind, only: SKC => SK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_SK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_SK3
#endif
        use pm_kind, only: SKC => SK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_SK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_SK2
#endif
        use pm_kind, only: SKC => SK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_SK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_SK1
#endif
        use pm_kind, only: SKC => SK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_IK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_IK5
#endif
        use pm_kind, only: IKC => IK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_IK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_IK4
#endif
        use pm_kind, only: IKC => IK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_IK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_IK3
#endif
        use pm_kind, only: IKC => IK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_IK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_IK2
#endif
        use pm_kind, only: IKC => IK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_IK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_IK1
#endif
        use pm_kind, only: IKC => IK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_LK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_LK5
#endif
        use pm_kind, only: LKC => LK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_LK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_LK4
#endif
        use pm_kind, only: LKC => LK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_LK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_LK3
#endif
        use pm_kind, only: LKC => LK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_LK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_LK2
#endif
        use pm_kind, only: LKC => LK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_LK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_LK1
#endif
        use pm_kind, only: LKC => LK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_CK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_CK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_CK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_CK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_CK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_RK5(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_RK4(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_RK3(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_RK2(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_EXP_ULD_D2001_RK1(mat, subset, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_EXP_ULD_D2001_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IK)                     , intent(in)                    :: nrow, ncol
        integer(IK)                     , intent(in)                    :: roff, coff
        real(RKC)                       , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! IMP_XLX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_SK5(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_SK4(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_SK3(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_SK2(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_SK1(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_IK5(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_IK4(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_IK3(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_IK2(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_IK1(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_LK5(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_LK4(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_LK3(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_LK2(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_LK1(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_CK5(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_CK4(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_CK3(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_CK2(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_CK1(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_RK5(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_RK4(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_RK3(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_RK2(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLX_D2X0X_RK1(mat, subset, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLX_D2X0X_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(low_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_UXX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_SK5(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_SK4(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_SK3(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_SK2(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_SK1(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_IK5(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_IK4(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_IK3(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_IK2(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_IK1(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_LK5(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_LK4(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_LK3(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_LK2(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_LK1(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_CK5(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_CK4(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_CK3(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_CK2(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_CK1(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_RK5(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_RK4(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_RK3(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_RK2(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXX_D20XX_RK1(mat, subset, vupp, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXX_D20XX_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(upp_type)                  , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_XXD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_SK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_SK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_SK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_SK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_SK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_IK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_IK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_IK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_IK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_IK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_LK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_LK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_LK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_LK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_LK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_CK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_CK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_CK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_CK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_CK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_RK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_RK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_RK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_RK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XXF_RK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XXF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_SK5(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_SK4(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_SK3(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_SK2(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_SK1(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_IK5(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_IK4(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_IK3(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_IK2(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_IK1(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_LK5(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_LK4(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_LK3(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_LK2(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_LK1(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_CK5(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_CK4(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_CK3(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_CK2(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_CK1(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_RK5(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_RK4(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_RK3(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_RK2(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX0_RK1(mat, subset, vdia, ndia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX0_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vdia
        integer(IK)                     , intent(in)                    :: ndia
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_SK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_SK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_SK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_SK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_SK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_IK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_IK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_IK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_IK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_IK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_LK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_LK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_LK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_LK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_LK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_CK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_CK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_CK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_CK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_CK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_RK5(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_RK4(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_RK3(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_RK2(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XXD_D2XX1_RK1(mat, subset, vdia)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XXD_D2XX1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        type(dia_type)                  , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_XLD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_SK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_SK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_SK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_SK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_SK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_IK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_IK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_IK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_IK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_IK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_LK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_LK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_LK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_LK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_LK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_CK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_CK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_CK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_CK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_CK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_RK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_RK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_RK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_RK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X00_RK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X00_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_SK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_SK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_SK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_SK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_SK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_IK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_IK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_IK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_IK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_IK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_LK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_LK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_LK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_LK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_LK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_CK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_CK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_CK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_CK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_CK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_RK5(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_RK4(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_RK3(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_RK2(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_XLD_D2X01_RK1(mat, subset, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_XLD_D2X01_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(lowDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_UXD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_SK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_SK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_SK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_SK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_SK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_IK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_IK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_IK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_IK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_IK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_LK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_LK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_LK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_LK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_LK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_CK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_CK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_CK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_CK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_CK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_RK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_RK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_RK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_RK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X0_RK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X0_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_SK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_SK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_SK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_SK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_SK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_IK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_IK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_IK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_IK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_IK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_LK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_LK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_LK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_LK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_LK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_CK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_CK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_CK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_CK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_CK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_RK5(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_RK4(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_RK3(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_RK2(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_UXD_D20X1_RK1(mat, subset, vupp, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_UXD_D20X1_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppDia_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_ULX

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_SK5(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_SK5
#endif
        use pm_kind, only: SKC => SK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_SK4(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_SK4
#endif
        use pm_kind, only: SKC => SK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_SK3(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_SK3
#endif
        use pm_kind, only: SKC => SK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_SK2(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_SK2
#endif
        use pm_kind, only: SKC => SK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_SK1(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_SK1
#endif
        use pm_kind, only: SKC => SK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_IK5(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_IK4(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_IK3(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_IK2(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_IK1(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_LK5(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_LK5
#endif
        use pm_kind, only: LKC => LK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_LK4(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_LK4
#endif
        use pm_kind, only: LKC => LK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_LK3(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_LK3
#endif
        use pm_kind, only: LKC => LK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_LK2(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_LK2
#endif
        use pm_kind, only: LKC => LK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_LK1(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_LK1
#endif
        use pm_kind, only: LKC => LK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_CK5(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_CK4(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_CK3(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_CK2(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_CK1(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_RK5(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_RK4(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_RK3(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_RK2(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULX_D200X_RK1(mat, subset, vupp, vlow, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULX_D200X_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)                     , intent(in)    , optional      :: doff
        type(uppLow_type)               , intent(in)                    :: subset
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! IMP_ULD

    interface setMatInit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_SK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_SK5
#endif
        use pm_kind, only: SKC => SK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_SK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_SK4
#endif
        use pm_kind, only: SKC => SK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_SK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_SK3
#endif
        use pm_kind, only: SKC => SK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_SK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_SK2
#endif
        use pm_kind, only: SKC => SK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_SK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_SK1
#endif
        use pm_kind, only: SKC => SK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
        character(*,SKC)                , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_IK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_IK5
#endif
        use pm_kind, only: IKC => IK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_IK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_IK4
#endif
        use pm_kind, only: IKC => IK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_IK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_IK3
#endif
        use pm_kind, only: IKC => IK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_IK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_IK2
#endif
        use pm_kind, only: IKC => IK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_IK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_IK1
#endif
        use pm_kind, only: IKC => IK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
        integer(IKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_LK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_LK5
#endif
        use pm_kind, only: LKC => LK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_LK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_LK4
#endif
        use pm_kind, only: LKC => LK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_LK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_LK3
#endif
        use pm_kind, only: LKC => LK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_LK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_LK2
#endif
        use pm_kind, only: LKC => LK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_LK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_LK1
#endif
        use pm_kind, only: LKC => LK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
        logical(LKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_CK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_CK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_CK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_CK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_CK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
        complex(CKC)                    , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_RK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_RK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_RK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_RK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2000_RK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2000_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
        real(RKC)                       , intent(in)                    :: vdia
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_SK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_SK5
#endif
        use pm_kind, only: SKC => SK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_SK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_SK4
#endif
        use pm_kind, only: SKC => SK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_SK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_SK3
#endif
        use pm_kind, only: SKC => SK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_SK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_SK2
#endif
        use pm_kind, only: SKC => SK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_SK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_SK1
#endif
        use pm_kind, only: SKC => SK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        character(*,SKC)                , intent(inout) , contiguous    :: mat(:,:)
        character(*,SKC)                , intent(in)    , contiguous    :: vdia(:)
        character(*,SKC)                , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_IK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_IK5
#endif
        use pm_kind, only: IKC => IK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_IK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_IK4
#endif
        use pm_kind, only: IKC => IK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_IK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_IK3
#endif
        use pm_kind, only: IKC => IK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_IK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_IK2
#endif
        use pm_kind, only: IKC => IK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_IK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_IK1
#endif
        use pm_kind, only: IKC => IK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        integer(IKC)                    , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)                    , intent(in)    , contiguous    :: vdia(:)
        integer(IKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_LK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_LK5
#endif
        use pm_kind, only: LKC => LK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_LK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_LK4
#endif
        use pm_kind, only: LKC => LK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_LK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_LK3
#endif
        use pm_kind, only: LKC => LK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_LK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_LK2
#endif
        use pm_kind, only: LKC => LK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_LK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_LK1
#endif
        use pm_kind, only: LKC => LK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        logical(LKC)                    , intent(inout) , contiguous    :: mat(:,:)
        logical(LKC)                    , intent(in)    , contiguous    :: vdia(:)
        logical(LKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_CK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_CK5
#endif
        use pm_kind, only: CKC => CK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_CK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_CK4
#endif
        use pm_kind, only: CKC => CK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_CK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_CK3
#endif
        use pm_kind, only: CKC => CK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_CK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_CK2
#endif
        use pm_kind, only: CKC => CK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_CK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_CK1
#endif
        use pm_kind, only: CKC => CK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        complex(CKC)                    , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)                    , intent(in)    , contiguous    :: vdia(:)
        complex(CKC)                    , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_RK5(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_RK4(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_RK3(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_RK2(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInit_IMP_ULD_D2001_RK1(mat, subset, vupp, vlow, vdia, doff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInit_IMP_ULD_D2001_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(uppLowDia_type)            , intent(in)                    :: subset
        integer(IK)                     , intent(in)    , optional      :: doff
        real(RKC)                       , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)                       , intent(in)    , contiguous    :: vdia(:)
        real(RKC)                       , intent(in)                    :: vupp, vlow
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixInit ! LCOV_EXCL_LINE