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
!>  This module contains procedures and generic interfaces for computing the <b>Discrete Fourier Transform</b>
!>  of a `real` or `complex` sequence using a <b>mixed-radix decimation-in-frequency Fast-Fourier Transform</b>.<br>
!>
!>  \details
!>  The discrete Fourier transform (**DFT**) converts a finite sequence of equally-spaced samples of a function
!>  into a same-length sequence of equally-spaced samples of the discrete-time Fourier transform (**DTFT**),
!>  which is a complex-valued function of frequency.<br>
!>
!>  Since DFT deals with a finite amount of data, it can be implemented
!>  in computers by numerical algorithms or even dedicated hardware.<br>
!>  These implementations usually employ efficient fast Fourier transform (**FFT**) algorithms,
!>  so much so that the terms **FFT** and **DFT** are used interchangeably.<br>
!>
!>  The **DFT** transforms a sequence of \f$N\f$ complex numbers
!>  \f$\{ x_j \} := x_0, x_1, \ldots, x_{N − 1}\f$ into another sequence of complex numbers,
!>  \f$\{z_k\} := z_0, z_1, \ldots, z_{N - 1}\f$ which is defined by,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      x_k
!>      &=& \sum_{j=0}^{N-1} z_j \cdot e^{-\frac {i 2\pi}{N}jk} ~, \\
!>      &=& \sum_{j=0}^{N-1} z_j \cdot \left[\cos\left(\frac{2\pi}{N}jk\right) - i \cdot \sin\left(\frac{2 \pi}{N}jk\right)\right],
!>  \f}
!>
!>  where \f$i\f$ represents the **imaginary unit**.<br>
!>  The naive evaluation of the DFT is a matrix-vector multiplication
!>  that costs \f$\mathcal{O}(n^2)\f$ operations for \f$N\f$ data-points.<br>
!>  Fast Fourier Transform (**FFT**) algorithms use a divide-and-conquer strategy to factorize the matrix into smaller sub-matrices,
!>  corresponding to the integer factors of the length \f$N\f$.<br>
!>  If \f$N\f$ can be factorized into a product of integers \f$f_1 f_2 \ldots f_m\f$,
!>  then the DFT can be computed in \f$\mathcal{O}(n \sum f_i)\f$ operations.<br>
!>  For a radix-2 FFT this gives an operation count of \f$\mathcal{O}(n \log_2 n)\f$.<br>
!>
!>  For every <b>Forward DFT</b> ([setFFTF](@ref pm_fftpack::setFFTF)),
!>  which expresses the data sequence in terms of frequency coefficients),
!>  there is a corresponding <b>Reverse (Backward)</b> DFT.<br>
!>  The **Reverse DFT** ([setFFTR](@ref pm_fftpack::setFFTR)) takes the output of the Forward DFT and
!>  returns the original data sequence given to the Forward DFT, <b>multiplied by the length of the sequence</b>,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      N z_j
!>      &=& \sum_{k=0}^{N-1} x_k \cdot e^{\frac {i 2\pi}{N}jk} ~, \\
!>      &=& \sum_{k=0}^{N-1} x_k \cdot \left[\cos\left(\frac{2\pi}{N}jk\right) + i \cdot \sin\left(\frac{2 \pi}{N}jk\right)\right],
!>  \f}
!>
!>  A third definition, the <b>Inverse DFT</b> ([setFFTI](@ref pm_fftpack::setFFTI)), further normalizes the **Reverse DFT** and returns,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      z_j
!>      &=& \frac{1}{N} \sum_{k=0}^{N-1} x_k \cdot e^{\frac {i 2\pi}{N}jk} ~, \\
!>      &=& \frac{1}{N} \sum_{k=0}^{N-1} x_k \cdot \left[\cos\left(\frac{2\pi}{N}jk\right) + i \cdot \sin\left(\frac{2 \pi}{N}jk\right)\right],
!>  \f}
!>
!>  such that the resulting **normalized Reverse (Backward) DFT** or **Inverse DFT** becomes a true inverse of the **Forward DFT**.<br>
!>  The Reverse DFT becomes relevant when the overall scale of the result is unimportant.<br>
!>  In such cases, it is convenient to use the Reverse (Backward) DFT instead of the Inverse DFT to avoid the redundant divisions in the normalization step.<br>
!>  A 64-bit `real` division is typically on the order of \f$\sim20\f$ CPU cycles on the contemporary hardware.<br>
!>  Avoiding the redundant normalizations can therefore become a significant saving if the DFT is to computed repeatedly for an array of millions of elements.<br>
!>
!>  The Cooley–Tukey algorithm
!>  ==========================
!>
!>  The Cooley–Tukey algorithm, named after J. W. Cooley and John Tukey, is the most common fast Fourier transform (FFT) algorithm.<br>
!>  It re-expresses the discrete Fourier transform (DFT) of an arbitrary composite size \f$N = N_{1} N_{2}\f$ in terms of \f$N_1\f$ smaller DFTs of sizes \f$N_2\f$, recursively,
!>  to reduce the computation time to \f$\mathcal{O}(N log N)\f$ for highly composite \f$N\f$ (smooth numbers).<br>
!>  Because of the algorithm importance, specific variants and implementation styles have become known by their own names, as described below.<br>
!>  Because the Cooley–Tukey algorithm breaks the DFT into smaller DFTs, it can be combined arbitrarily with any other algorithm for the DFT.<br>
!>  For example, The Rader or Bluestein algorithm can be used to handle large prime factors that cannot be decomposed by Cooley–Tukey,
!>  or the prime-factor algorithm can be exploited for greater efficiency in separating out relatively prime factors.<br>
!>
!>  The algorithm, along with its recursive application, was invented by Carl Friedrich Gauss.<br>
!>  Cooley and Tukey independently rediscovered and popularized it 160 years later.<br>
!>
!>  Intuition
!>  ---------
!>
!>  The Cooley–Tukey algorithms recursively re-express a DFT of a composite size \f$N = N_1 N_2\f$ as:<br>
!>  <ol>
!>      <li>    Perform \f$N_1\f$ DFTs of size \f$N_2\f$.
!>      <li>    Multiply by complex roots of unity (often called the twiddle factors).
!>      <li>    Perform \f$N_2\f$ DFTs of size \f$N_1\f$.
!>  </ol>
!>  Typically, either \f$N_1\f$ or \f$N_2\f$ is a small factor (not necessarily prime), called the **radix** (which can differ between stages of the recursion).<br>
!>  If \f$N_1\f$ is the radix, it is called a **decimation in time (DIT) algorithm**, whereas if \f$N_2\f$ is the radix, it is **decimation in frequency (DIF) algorithm**, also called the **Sande–Tukey algorithm**.<br>
!>
!>  \note
!>  The routines of this module are generic **mixed-radix** implementation,
!>  meaning that the can be used to compute the forward and reverse FFT of <b>arbitrary-length</b> `real` and `complex` data sequences.<br>
!>  By contrast, **radix-2 FFT** routines (such as those of the [Numerical Recipes](http://numerical.recipes)), require a potentially-padded sequence (with trailing zeros)
!>  such that the length of the input sequence is always a power of \f$2\f$.<br>
!>  The **mixed-radix** algorithms are generally faster than the `radix-2` algorithms at the cost of requiring an extra memory storage of the same length as the input sequence.<br>
!>  The mixed-radix algorithm uses optimized small length sub-transform FFTs which are combined to create larger FFT of the sequence.<br>
!>  FFTPACK implements efficient sub-transform for factors of 2, 3, 4, and 5.<br>
!>  The computations for other factors fall back to a general length-\f$N\f$ algorithm which uses Singleton method for efficiently computing a DFT.<br>
!>  This module is \f$\mathcal{O}(N^2)\f$, slower than an explicitly dedicated factor implementation. However, it works for arbitrary length sub-transforms.<br>
!>  Arbitrary-length sub-transforms are factorized as much as possible.<br>
!>  For example, a sub-sequence length of \f$143\f$ will be factorized into \f$11\times13\f$.<br>
!>  Large prime factors (e.g., \f$n = 2\times 3\times 3,351,773\f$) are the worst case scenario, because they cannot be further factorized.<br>
!>  The large prime factors in sequence should be avoided as much as possible because their \f$\mathcal{O}(N^2)\f$ computational cost scaling will dominate the run-time.<br>
!>  The mixed-radix initialization procedures under the generic interface [getFactorFFT](@ref pm_fftpack::getFactorFFT) return the list of factors chosen for a given input sequence length \f$N\f$.<br>
!>  The output factors from these procedures can be checked to ensure an efficient factorization and to **estimate the run-time**.<br>
!>  The run-time of the FFT algorithms of this module roughly scale as \f$N\sum f_i\f$, where the \f$f_i\f$ are the factors of the sequence length \f$N\f$.<br>
!>  If specific data lengths appear in the FFT problem that cannot be efficiently factorized using the existing small-prime factor implementations,
!>  the performance can be improved by explicitly implementing the algorithms for the specific factors of interest.<br>
!>
!>  \note
!>  See [this catalog](https://primes.utm.edu/primes/lists/all.txt) for a list of prime numbers.<br>
!>
!>  \attention
!>  There are two possible choices for the sign of the exponential in the forward and inverse FFT equations in the above.<br>
!>  While both conventions are commonly used, the FFTPACK library (as given here) uses the above convention (a negative exponential for the forward transform).<br>
!>  A prominent implementation that uses the opposite convention is provided by the [Numerical Recipes](http://numerical.recipes).<br>
!>
!>  \remark
!>  This reimplementation of the original FFTPACK library achieves the following goals:<br>
!>  <ol>
!>      <li>    Extension of the functionality of the library to arbitrary `complex` and `real` precision kinds (e.g., \RKALL, \CKALL).<br>
!>      <li>    Usage and flexibility improvements to the syntax and interface of the library.<br>
!>      <li>    Performance improvements to the original library.<br>
!>      <li>    Safety improvements to the library, particularly,<br>
!>              <ol>
!>                  <li>    the removal of all assumed-shape dummy arguments,<br>
!>                  <li>    the removal of all implicit conversions between `complex` and `real` actual and dummy arguments,<br>
!>                  <li>    the removal of all implicit conversions between `complex` and `integer` actual and dummy arguments,<br>
!>                  <li>    the addition of abundant runtime array-bound checks to reduce misspecified input arguments.<br>
!>              </ol>
!>      <li>    Removal of all fixed-format FORTRAN77 syntax and their replacement with equivalent free-format syntax.<br>
!>      <li>    Removal of all obsolescent language features such as `goto` statements from the original library.<br>
!>      <li>    Significant reduction in the library size with the help of generics and macros.<br>
!>  </ol>
!>
!>  \note
!>  An acyclic time series may need to be tapered by applying split-cosine-bell tapering to the time series prior to a Fourier Analysis.<br>
!>  See Chapter 5 of the book *Fourier Analysis of Time Series, Peter Bloomfield, Wiley-Interscience, 1976*.<br>
!>
!>  \see
!>  [FFTPACK](http://www.netlib.org/fftpack/index.html), the original mixed-radix FFT library in FORTRAN77 by Paul Swarztrauber of NCAR.<br>
!>  [FFTPACK 5.1](https://github.com/jlokimlin/fftpack5.1), an unofficial mirror of FFTPACK 5.1  by Paul Swarztrauber of NCAR and Richard A. Valent (2011).<br>
!>  [FFTPACK 5.1](https://people.sc.fsu.edu/~jburkardt/f77_src/fftpack5.1/fftpack5.1_reference.html), an unofficial reassembly of FFTPACK 5.1 by the original authors.<br>
!>  [Modernized FFTPACK](https://github.com/fortran-lang/fftpack), a modernized version of the original FFTPACK library in modern Fortran by the Fortran-lang community.<br>
!>  [The GNU FFTPACK in C](https://www.gnu.org/software/gsl/doc/html/fft.html), an incomplete translation of the original FORTRAN77 FFTPACK library in C language.<br>
!>  Paul Swarztrauber, 1982, Vectorizing the Fast Fourier Transforms, Parallel Computations, G. Rodrigue, ed., Academic Press, New York<br>
!>  Paul Swarztrauber, 1984, Fast Fourier Transforms Algorithms for Vector Computers, Parallel Computing, pp.45-63<br>
!>  Clive Temperton, 1983, Self-sorting Mixed-radix FFTs, Journal of Computational Physics, 52, 340-350<br>
!>
!>  \test
!>  [test_pm_fftpack](@ref test_pm_fftpack)
!>
!>  \todo
!>  \pvhigh
!>  The sine and cosine forward and reverse transforms must be added to this module.<br>
!>
!>  \final
!>  This module is a complete modernization and significant extension of the original
!>  [netlib.org](https://netlib.org/fftpack) public-domain FFTPACK library in FORTRAN77 programming language.<br>
!>  If you use any or parts of this module, you should also cite the following reference in addition to the relevant ParaMonte library reference.<br>
!>  <ul>
!>      <li>    P. N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations (G. Rodrigue, ed.), Academic Press, 1982, pp. 51-83.<br>
!>  </ul>
!>  \verbatim
!>                        FFTPACK
!>
!>  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!>
!>                    version 4  april 1985
!>
!>       a package of fortran subprograms for the fast fourier
!>        transform of periodic and other symmetric sequences
!>
!>                           by
!>
!>                    paul n swarztrauber
!>
!>    national center for atmospheric research  boulder,colorado 80307
!>
!>     which is sponsored by the national science foundation
!>
!>  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!>  \endverbatim
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_fftpack

    use pm_kind, only: SK, IK, LK
    use pm_array, only: allocatable, allocatable_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_fftpack"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the factorization vector `factor` of the specified input sequence length and the corresponding vector of trigonometric coefficients `coef`.<br>
    !>
    !>  \details
    !>  The factorization and the trigonometric coefficient vectors are required for computing the Forward or Reverse FFT.<br>
    !>  Note that the computed factoring is **not** necessarily a **complete prime factoring** of the length of the input `data` sequence.<br>
    !>  This is because some of target algorithms corresponding to some higher factors are more efficient than the composition of the lower factors.<br>
    !>  For example, the corresponding algorithms to the composite factors `4` and `6` are faster than combining the corresponding algorithms for `2*2` and `2*3`.<br>
    !>
    !>  \param[in]  data    :   The input `contiguous` vector of either,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, or<br>
    !>                          </ul>
    !>                          containing the data sequence whose Forward or Reverse FFT is to computed.<br>
    !>                          Only the length and the type of input vector is used within the algorithm.<br>
    !>  \param[out] coef    :   The output `contiguous` vector of the same type, kind, and size as the input `data` sequence vector,
    !>                          containing the trigonometric look up table required for Forward or Reverse FFT of the specified `data` sequence.<br>
    !>                          (**optional**. If missing, the trigonometric coefficients will not be computed.)
    !>  \param[in]  attr    :   The input scalar of type [allocatable_type](@ref pm_array::allocatable_type) signifying the `allocatable` status of the input `coef` vector.<br>
    !>                          If present, the input `coef` will be reallocated to `size(data)`.<br>
    !>                          (**optional**. It can be present only if the input argument `coef` has the `allocatable` attribute.)
    !>
    !>  \return
    !>  `factor`            :   The output `allocatable` array of shape `(:)` of type `integer` of default kind \IK,
    !>                          containing the factorization of the length of the `data` sequence whose Forward or Reverse FFT is to be computed.<br>
    !>                          By definition, the condition `product(factor) == size(data)` holds.<br>
    !>
    !>  \interface{getFactorFFT}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: getFactorFFT
    !>      integer(IK), allocatable :: factor(:)
    !>
    !>      factor = getFactorFFT(data(1:lenData))
    !>      factor = getFactorFFT(data(1:lenData), coef(1:lenData))
    !>      factor = getFactorFFT(data(1:lenData), coef(:), attr)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(data) > 1` must hold for the corresponding arguments.<br>
    !>  The condition `size(data) == size(coef)` must hold for the corresponding arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the output argument `coef` is present.<br>
    !>
    !>  \remark
    !>  While the contents of the input `data` sequence is not used within the algorithm, its presence and requirement is to minimize user-mistakes.<br>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>
    !>  \example{getFactorFFT}
    !>  \include{lineno} example/pm_fftpack/getFactorFFT/main.F90
    !>  \compilef{getFactorFFT}
    !>  \output{getFactorFFT}
    !>  \include{lineno} example/pm_fftpack/getFactorFFT/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \todo
    !>  \plow
    !>  Extension to higher order factors may be worthwhile in future.
    !>
    !>  \final{getFactorFFT}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getFactorFFT

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function setFactorFFTDXX_CK5(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function setFactorFFTDXX_CK4(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function setFactorFFTDXX_CK3(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function setFactorFFTDXX_CK2(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function setFactorFFTDXX_CK1(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function setFactorFFTDXX_RK5(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function setFactorFFTDXX_RK4(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function setFactorFFTDXX_RK3(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function setFactorFFTDXX_RK2(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function setFactorFFTDXX_RK1(data) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDXX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function setFactorFFTDCX_CK5(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK4_ENABLED
    impure module function setFactorFFTDCX_CK4(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK3_ENABLED
    impure module function setFactorFFTDCX_CK3(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK2_ENABLED
    impure module function setFactorFFTDCX_CK2(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if CK1_ENABLED
    impure module function setFactorFFTDCX_CK1(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function setFactorFFTDCX_RK5(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK4_ENABLED
    impure module function setFactorFFTDCX_RK4(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK3_ENABLED
    impure module function setFactorFFTDCX_RK3(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK2_ENABLED
    impure module function setFactorFFTDCX_RK2(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

#if RK1_ENABLED
    impure module function setFactorFFTDCX_RK1(data, coef) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCX_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function setFactorFFTDCA_CK5(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if CK4_ENABLED
    impure module function setFactorFFTDCA_CK4(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if CK3_ENABLED
    impure module function setFactorFFTDCA_CK3(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if CK2_ENABLED
    impure module function setFactorFFTDCA_CK2(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if CK1_ENABLED
    impure module function setFactorFFTDCA_CK1(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG), intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function setFactorFFTDCA_RK5(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if RK4_ENABLED
    impure module function setFactorFFTDCA_RK4(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if RK3_ENABLED
    impure module function setFactorFFTDCA_RK3(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if RK2_ENABLED
    impure module function setFactorFFTDCA_RK2(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

#if RK1_ENABLED
    impure module function setFactorFFTDCA_RK1(data, coef, attr) result(factor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFactorFFTDCA_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)   , intent(out)   , allocatable   :: coef(:)
        integer(IK)                 , allocatable   :: factor(:)
        type(allocatable_type)      , intent(in)    :: attr
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Forward Fourier Transform (a.k.a. Fourier Analysis)
    !>  of a periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  See the documentation of [setFFTF](@ref pm_fftpack::setFFTF) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of the same type, kind, and size as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTF}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: getFFTF
    !>
    !>      fft = getFFTF(data(:))
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTF](@ref pm_fftpack::setFFTF).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{getFFTF}
    !>  \include{lineno} example/pm_fftpack/getFFTF/main.F90
    !>  \compilef{getFFTF}
    !>  \output{getFFTF}
    !>  \include{lineno} example/pm_fftpack/getFFTF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{getFFTF}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getFFTF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getFFTF_CK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTF_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTF_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTF_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTF_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getFFTF_RK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTF_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTF_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTF_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTF_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Forward Fourier Transform (or equivalently, the Fourier coefficients) of a periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  For an input `data` sequence of length \f$N\f$ of type `complex`, the output coefficients `fftf` (either stored in `data` or `work`) contain the following,
    !>  \f{equation}{
    !>      \ms{fftf}(k) = \sum_{j = 1}^{N} \ms{data}(j) \exp\left( -i \cdot (k - 1)(j - 1) \frac{2\pi}{N} \right) ~, \\
    !>  \f}
    !>  where \f$i = \sqrt{-1}\f$.<br>
    !>
    !>  For an input `data` sequence of length \f$N\f$ of type `real`, the output coefficients `fftf` (either stored in `data` or `work`) contain the following,
    !>  \f{eqnarray}{
    !>      \ms{fftf}(1)    &=& \sum_{j = 1}^{N} +\ms{data}(j) ~, \\
    !>      \ms{fftf}(2k-2) &=& \sum_{j = 1}^{N} +\ms{data}(j) \cos\left((k - 1)(j - 1) \frac{2\pi}{N}\right) ~~~,~~~ k = \left\lceil\frac{N}{2}\right\rceil ~, \\
    !>      \ms{fftf}(2k-1) &=& \sum_{j = 1}^{N} -\ms{data}(j) \sin\left((k - 1)(j - 1) \frac{2\pi}{N}\right) ~~~,~~~ k = \left\lceil\frac{N}{2}\right\rceil ~, \\
    !>      \ms{fftf}(N)    &=& \sum_{j = 1}^{N} +\ms{data}(j) (-1)^{(j - 1)} ~~~,~~~ \ms{iff} ~ \frac{N}{2} = \left\lceil\frac{N}{2}\right\rceil ~,
    !>  \f}
    !>
    !>  A call to [setFFTF()](@ref pm_fftpack::setFFTF) followed by a call to [setFFTR()](@ref pm_fftpack::setFFTR) will multiply the initial sequence `data` by its length \f$N\f$.<br>
    !>  A call to [setFFTF()](@ref pm_fftpack::setFFTF) followed by a call to [setFFTI()](@ref pm_fftpack::setFFTI) will retrieve the initial sequence `data`.<br>
    !>  See the documentation of [pm_fftpack](@ref pm_fftpack) for more details.<br>
    !>
    !>  \param[in]      factor      :   The input `contiguous` vector of shape `(:)` of type `integer` of default kind \IK,
    !>                                  containing the factorization of the length of the input `data` sequence whose FFT is to be computed.<br>
    !>                                  This input argument along with `coef` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[in]      coef        :   The input `contiguous` vector of shape `(1:size(data))` of the same type and kind as the input argument `data`,
    !>                                  containing the trigonometric look up table required for FFT of the specified `data` sequence.<br>
    !>                                  This input argument along with `factor` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[inout]   data        :   The input/output `contiguous` vector of arbitrary size of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the periodic sequence whose FFT is to be computed.<br>
    !>                                  On output, `data` contains the FFT result if `inwork == .true.`.<br>
    !>                                  Otherwise, the original input data is completely destroyed on return.<br>
    !>  \param[out]     work        :   The output `contiguous` vector of the same type, kind, and size as the input `data`, that is used as a workspace.<br>
    !>                                  On output, `work` contains the FFT result if `inwork == .false.`.<br>
    !>  \param[out]     inwork      :   The output scalar of type `logical` of default kind \LK.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the FFT result is stored in the output `work` argument upon return from the procedure.<br>
    !>                                      <li>    If `.false.`, the FFT result is stored in the output `data` argument upon return from the procedure.<br>
    !>                                  </ol>
    !>
    !>  \interface{setFFTF}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: setFFTF
    !>
    !>      call setFFTF(factor(:), coef(1:size(data)), data(:), work(1:size(data)), inwork)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(data) == size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(data) == size(work)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{setFFTF}
    !>  \include{lineno} example/pm_fftpack/setFFTF/main.F90
    !>  \compilef{setFFTF}
    !>  \output{setFFTF}
    !>  \include{lineno} example/pm_fftpack/setFFTF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{setFFTF}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTF_CK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTF_CK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTF_CK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTF_CK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTF_CK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTF_RK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTF_RK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTF_RK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTF_RK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTF_RK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Reverse (unnormalized) Fourier Transform of a
    !>  periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  See the documentation of [setFFTR](@ref pm_fftpack::setFFTR) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of the same type, kind, and size as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTR}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: getFFTR
    !>
    !>      fft = getFFTR(data(:))
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTR](@ref pm_fftpack::setFFTR).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{getFFTR}
    !>  \include{lineno} example/pm_fftpack/getFFTR/main.F90
    !>  \compilef{getFFTR}
    !>  \output{getFFTR}
    !>  \include{lineno} example/pm_fftpack/getFFTR/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{getFFTR}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getFFTR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getFFTR_CK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTR_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTR_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTR_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTR_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getFFTR_RK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTR_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTR_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTR_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTR_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Reverse (unnormalized) Fourier Transform of a periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  For an input `data` sequence of length \f$N\f$ of type `complex`, the output coefficients `fftr` (either stored in `data` or `work`) contain the following,
    !>  \f{equation}{
    !>      \ms{fftr}(k) = \sum_{j = 1}^{N} \ms{data}(j) \exp\left( +i \cdot (k - 1)(j - 1) \frac{2\pi}{N} \right) ~, \\
    !>  \f}
    !>  where \f$i = \sqrt{-1}\f$.<br>
    !>
    !>  For an input `data` sequence of length \f$N\f$ of type `real`, the output coefficients `fftr` (either stored in `data` or `work`) contain the following,
    !>  \f{eqnarray}{
    !>      \ms{fftr}(k)    &=& \ms{data}(1) \\
    !>                      &+& \sum_{j = 2}^{\left\lceil\frac{N}{2}\right\rceil} 2\left[ \ms{data}(2j - 2) \cos\left((k - 1)(j - 1)\frac{2\pi}{N}\right) - \ms{data}(2j - 1) \sin\left((k - 1)(j - 1)\frac{2\pi}{N}\right) \right] \\
    !>                      &+& (-1)^{(k - 1)} \ms{data}(N) \times 2\left( \frac{n+1}{2} - \left\lceil\frac{N}{2}\right\rceil \right)
    !>      ~,
    !>  \f}
    !>  where the expression \f$2\left( \frac{n+1}{2} - \left\lceil\frac{N}{2}\right\rceil \right)\f$ in the last term ensures
    !>  the last term is non-zero only when \f$N\f$ is even.<br>
    !>
    !>  A call to [setFFTF()](@ref pm_fftpack::setFFTF) followed by a call to [setFFTR()](@ref pm_fftpack::setFFTR) will multiply the initial sequence `data` by its length \f$N\f$.<br>
    !>  A call to [setFFTF()](@ref pm_fftpack::setFFTF) followed by a call to [setFFTI()](@ref pm_fftpack::setFFTI) will retrieve the initial sequence `data`.<br>
    !>  See the documentation of [pm_fftpack](@ref pm_fftpack) for more details.<br>
    !>
    !>  \param[in]      factor      :   The input `contiguous` vector of shape `(:)` of type `integer` of default kind \IK,
    !>                                  containing the factorization of the length of the input `data` sequence whose FFT is to be computed.<br>
    !>                                  This input argument along with `coef` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[in]      coef        :   The input `contiguous` vector of shape `(1:size(data))` of the same type and kind as the input argument `data`,
    !>                                  containing the trigonometric look up table required for FFT of the specified `data` sequence.<br>
    !>                                  This input argument along with `factor` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[inout]   data        :   The input/output `contiguous` vector of arbitrary size of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the periodic sequence whose FFT is to be computed.<br>
    !>                                  On output, `data` contains the FFT result if `inwork == .true.`.<br>
    !>                                  Otherwise, the original input data is completely destroyed on return.<br>
    !>  \param[out]     work        :   The output `contiguous` vector of the same type, kind, and size as the input `data`, that is used as a workspace.<br>
    !>                                  On output, `work` contains the FFT result if `inwork == .false.`.<br>
    !>  \param[out]     inwork      :   The output scalar of type `logical` of default kind \LK.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the FFT result is stored in the output `work` argument upon return from the procedure.<br>
    !>                                      <li>    If `.false.`, the FFT result is stored in the output `data` argument upon return from the procedure.<br>
    !>                                  </ol>
    !>
    !>  \interface{setFFTR}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: setFFTR
    !>
    !>      call setFFTR(factor(:), coef(1:size(data)), data(:), work(1:size(data)), inwork)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(data) == size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(data) == size(work)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{setFFTR}
    !>  \include{lineno} example/pm_fftpack/setFFTR/main.F90
    !>  \compilef{setFFTR}
    !>  \output{setFFTR}
    !>  \include{lineno} example/pm_fftpack/setFFTR/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{setFFTR}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTR_CK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTR_CK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTR_CK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTR_CK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTR_CK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTR_RK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTR_RK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTR_RK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTR_RK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTR_RK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Inverse (normalized) Fourier Transform of a
    !>  periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  See the documentation of [setFFTI](@ref pm_fftpack::setFFTI) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of the same type, kind, and size as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTI}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: getFFTI
    !>
    !>      fft = getFFTI(data(:))
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTR](@ref pm_fftpack::setFFTR).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{getFFTI}
    !>  \include{lineno} example/pm_fftpack/getFFTI/main.F90
    !>  \compilef{getFFTI}
    !>  \output{getFFTI}
    !>  \include{lineno} example/pm_fftpack/getFFTI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{getFFTI}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getFFTI

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getFFTI_CK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTI_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTI_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTI_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTI_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getFFTI_RK5(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTI_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTI_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTI_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTI_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(size(data, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Inverse (normalized) Fourier Transform of a periodic sequence of type `complex` or `real` of arbitrary kind type parameter.
    !>
    !>  \details
    !>  See the documentation of [setFFTR](@ref pm_fftpack::setFFTR) for more details.<br>
    !>
    !>  \param[in]      factor      :   The input `contiguous` vector of shape `(:)` of type `integer` of default kind \IK,
    !>                                  containing the factorization of the length of the input `data` sequence whose FFT is to be computed.<br>
    !>                                  This input argument along with `coef` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[in]      coef        :   The input `contiguous` vector of shape `(1:size(data))` of the same type and kind as the input argument `data`,
    !>                                  containing the trigonometric look up table required for FFT of the specified `data` sequence.<br>
    !>                                  This input argument along with `factor` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[inout]   data        :   The input/output `contiguous` vector of arbitrary size of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the periodic sequence whose FFT is to be computed.<br>
    !>                                  On output, `data` contains the FFT result if `inwork == .true.`.<br>
    !>                                  Otherwise, the original input data is completely destroyed on return.<br>
    !>  \param[out]     work        :   The output `contiguous` vector of the same type, kind, and size as the input `data`, that is used as a workspace.<br>
    !>                                  On output, `work` contains the FFT result if `inwork == .false.`.<br>
    !>  \param[out]     inwork      :   The output scalar of type `logical` of default kind \LK.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the FFT result is stored in the output `work` argument upon return from the procedure.<br>
    !>                                      <li>    If `.false.`, the FFT result is stored in the output `data` argument upon return from the procedure.<br>
    !>                                  </ol>
    !>
    !>  \interface{setFFTI}
    !>  \code{.F90}
    !>
    !>      use pm_fftpack, only: setFFTI
    !>
    !>      call setFFTI(factor(:), coef(1:size(data)), data(:), work(1:size(data)), inwork)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(data) == size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(data) == size(work)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftpack::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftpack::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftpack::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftpack::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftpack::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftpack::setFFTI)<br>
    !>
    !>  \example{setFFTI}
    !>  \include{lineno} example/pm_fftpack/setFFTI/main.F90
    !>  \compilef{setFFTI}
    !>  \output{setFFTI}
    !>  \include{lineno} example/pm_fftpack/setFFTI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftpack](@ref test_pm_fftpack)
    !>
    !>  \final{setFFTI}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTI

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTI_CK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTI_CK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTI_CK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTI_CK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTI_CK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        complex(CKG), intent(in)    , contiguous    :: coef(:)
        complex(CKG), intent(inout) , contiguous    :: data(:)
        complex(CKG), intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTI_RK5(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTI_RK4(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTI_RK3(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTI_RK2(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTI_RK1(factor, coef, data, work, inwork)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)    , contiguous    :: factor(:)
        real(RKG)   , intent(in)    , contiguous    :: coef(:)
        real(RKG)   , intent(inout) , contiguous    :: data(:)
        real(RKG)   , intent(out)   , contiguous    :: work(:)
        logical(LK) , intent(out)                   :: inwork
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_fftpack