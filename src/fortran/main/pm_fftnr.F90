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
!>  of a `real` or `complex` sequence using <b>radix-2 Cooley–Tukey Fast-Fourier Transform</b>.<br>
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
!>  For every <b>Forward DFT</b> ([setFFTF](@ref pm_fftnr::setFFTF)),
!>  which expresses the data sequence in terms of frequency coefficients),
!>  there is a corresponding <b>Reverse (Backward)</b> DFT.<br>
!>  The **Reverse DFT** ([setFFTR](@ref pm_fftnr::setFFTR)) takes the output of the Forward DFT and
!>  returns the original data sequence given to the Forward DFT, <b>multiplied by the length of the sequence</b>,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      N z_j
!>      &=& \sum_{k=0}^{N-1} x_k \cdot e^{\frac {i 2\pi}{N}jk} ~, \\
!>      &=& \sum_{k=0}^{N-1} x_k \cdot \left[\cos\left(\frac{2\pi}{N}jk\right) + i \cdot \sin\left(\frac{2 \pi}{N}jk\right)\right],
!>  \f}
!>
!>  A third definition, the <b>Inverse DFT</b> ([setFFTI](@ref pm_fftnr::setFFTI)), further normalizes the **Reverse DFT** and returns,<br>
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
!>  The radix-2 DIT case
!>  --------------------
!>
!>  A radix-2 **decimation-in-time** (DIT) FFT is the simplest and most common form of the Cooley–Tukey algorithm,
!>  although highly optimized Cooley–Tukey implementations typically use other forms of the algorithm some of which are implemented in [pm_fftpack](@ref pm_fftpack).<br>
!>  Radix-2 DIT divides a DFT of size N into two interleaved DFTs (hence the name **radix-2**) of size \f$N / 2\f$ with each recursive stage.<br>
!>  Consider again the discrete Fourier transform (DFT) defined by the formula:<br>
!>  \f{equation}{
!>    X_{k}=\sum_{n=0}^{N-1}x_{n}e^{-{\frac {2\pi i}{N}}nk},
!>  \f}
!>  where \f$k\f$ is an integer ranging from \f$0\f$ to \f$N − 1\f$.<br>
!>  Radix-2 DIT first computes the DFTs of the even-indexed inputs \f$(x_{2m}=x_{0},x_{2},\ldots ,x_{N-2})\f$ and of the odd-indexed inputs \f$(x_{2m+1}=x_{1},x_{3},\ldots ,x_{N-1})\f$,
!>  and then combines those two results to produce the DFT of the whole sequence.<br>
!>  This idea can then be performed recursively to reduce the overall runtime to \f$\mathcal{O}(N log N)\f$.<br>
!>  This simplified form assumes that \f$N\f$ is a power of two;<br>
!>  since the number of sample points N can usually be chosen freely by the application (e.g. by changing the sample rate or window, **zero-padding**, etcetera), this is often not an important restriction.<br>
!>  The radix-2 DIT algorithm rearranges the DFT of the function \f$x_{n}\f$ into two parts: a sum over the even-numbered indices \f$n = {2m}\f$ and a sum over the odd-numbered indices \f$n = {2m+1}\f$:<br>
!>  \f{equation}{
!>      {\begin{matrix}X_{k}&=&\sum \limits_{m=0}^{N/2-1}x_{2m}e^{-{\frac {2\pi i}{N}}(2m)k}+\sum \limits_{m=0}^{N/2-1}x_{2m+1}e^{-{\frac {2\pi i}{N}}(2m+1)k}\end{matrix}}
!>  \f}
!>  One can factor a common multiplier \f$e^{-{\frac {2\pi i}{N}}k}\f$ out of the second sum, as shown in the equation below.<br>
!>  It is then clear that the two sums are the DFT of the even-indexed part \f$x_{{2m}}\f$ and the DFT of odd-indexed part \f$x_{{2m+1}}\f$ of the function \f$x_{n}\f$.<br>
!>  Denote the DFT of the Even-indexed inputs \f$x_{{2m}}\f$ by \f$E_{k}\f$ and the DFT of the Odd-indexed inputs \f$x_{2m+1}\f$ by \f$O_{k}\f$ and we obtain:<br>
!>  \f{equation}{
!>      \begin{matrix}X_{k}
!>      = \underbrace{\sum \limits_{m=0}^{N/2-1}x_{2m}e^{-{\frac {2\pi i}{N/2}}mk}}_{\mathrm {DFT\;of\;even-indexed\;part\;of\;} x_{n}}{}
!>      + e^{-{\frac {2\pi i}{N}}k} \underbrace{\sum \limits_{m=0}^{N/2-1}x_{2m+1}e^{-{\frac {2\pi i}{N/2}}mk}}_{\mathrm {DFT\;of\;odd-indexed\;part\;of\;} x_{n}}
!>      = E_{k}+e^{-{\frac {2\pi i}{N}}k}O_{k}\qquad {\text{ for }}k=0,\dots ,{\frac {N}{2}}-1.\end{matrix}
!>  \f}
!>  Note that the equalities hold for \f$k = 0,\dots ,N-1\f$ but the crux is that \f$E_{k}\f$ and \f$O_{k}\f$ are calculated in this way for \f$k=0, \dots, {\frac {N}{2}}-1\f$ only.<br>
!>  Thanks to the periodicity of the complex exponential, \f$X_{k+{\frac {N}{2}}}\f$ is also obtained from \f$E_{k}\f$ and \f$O_{k}\f$:<br>
!>  \f{equation}{
!>      \begin{aligned}
!>          X_{k+{\frac {N}{2}}}
!>          &=\sum \limits_{m=0}^{N/2-1}x_{2m}e^{-{\frac {2\pi i}{N/2}}m(k+{\frac {N}{2}})}+e^{-{\frac {2\pi i}{N}}(k+{\frac {N}{2}})}\sum \limits_{m=0}^{N/2-1}x_{2m+1}e^{-{\frac {2\pi i}{N/2}}m(k+{\frac {N}{2}})} \\
!>          &=\sum \limits_{m=0}^{N/2-1}x_{2m}e^{-{\frac {2\pi i}{N/2}}mk}e^{-2\pi mi}+e^{-{\frac {2\pi i}{N}}k}e^{-\pi i}\sum \limits_{m=0}^{N/2-1}x_{2m+1}e^{-{\frac {2\pi i}{N/2}}mk}e^{-2\pi mi} \\
!>          &=\sum \limits_{m=0}^{N/2-1}x_{2m}e^{-{\frac {2\pi i}{N/2}}mk}-e^{-{\frac {2\pi i}{N}}k}\sum \limits_{m=0}^{N/2-1}x_{2m+1}e^{-{\frac {2\pi i}{N/2}}mk}\\
!>          &=E_{k}-e^{-{\frac {2\pi i}{N}}k}O_{k}
!>      \end{aligned}
!>  \f}
!>  We can rewrite \f$X_k\f$ and \f$X_{k+{\frac {N}{2}}}\f$ as:<br>
!>  \f{equation}{
!>      \begin{matrix}
!>          X_{k}
!>          &=& E_{k}+e^{-{\frac {2\pi i}{N}}{k}}O_{k}\\X_{k+{\frac {N}{2}}}
!>          &=& E_{k}-e^{-{\frac {2\pi i}{N}}{k}}O_{k}
!>      \end{matrix}
!>  \f}
!>  This result, expressing the DFT of length \f$N\f$ recursively in terms of two DFTs of size \f$N/2\f$, is the core of the radix-2 DIT fast Fourier transform.<br>
!>  The algorithm gains its speed by re-using the results of intermediate computations to compute multiple DFT outputs.<br>
!>  Note that final outputs are obtained by a \f$+/−\f$ combination of \f$E_{k}\f$ and \f$O_{k}\exp(-2\pi ik/N)\f$, which is simply a size-2 DFT (sometimes called a butterfly in this context);<br>
!>  when this is generalized to larger radices below, the size-2 DFT is replaced by a larger DFT (which itself can be evaluated with an FFT).<br>
!>
!>  \note
!>  The routines of this module provide a **radix-2** implementation of the Cooley–Tukey algorithm,
!>  meaning that the can be used to compute the forward and reverse FFT of `real` and `complex` data sequences whose length is a power of `2`.<br>
!>  This requires a potentially zero-padded sequence (with trailing zeros) such that the length of the input sequence is always a power of \f$2\f$.<br>
!>  One of the most famous such implementation is that of [Numerical Recipes](http://numerical.recipes)), require a potentially-padded sequence (with trailing zeros)
!>  The **mixed-radix** algorithms as implemented in [pm_fftpack](@ref pm_fftpack) are generally faster than the `radix-2` algorithms
!>  at the cost of requiring an extra memory storage of the same length as the input sequence.<br>
!>
!>  \devnote
!>  The naming of this module is in honor of the Numerical Recipes book in Fortran
!>  that contains one of the most widely used implementations of radix-2 algorithm.<br>
!>
!>  \see
!>  [pm_fftpack](@ref pm_fftpack).<br>
!>  Numerical Recipes in Fortran, Press et al. 1992.<br>
!>
!>  \test
!>  [test_pm_fftnr](@ref test_pm_fftnr)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_fftnr

    use pm_kind, only: SK, IK, LK, RKB

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_fftnr"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Forward Fourier Transform (a.k.a. Fourier Analysis)
    !>  of a periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of size [getExpNext(size(data), 2)](@ref pm_mathExp::getExpNext)
    !>                              of the same type and kind as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTF}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: getFFTF
    !>      type_of(data) :: fft(1 : getExpNext(size(data)))
    !>
    !>      fft(:) = getFFTF(data(:))
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTF](@ref pm_fftnr::setFFTF).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{getFFTF}
    !>  \include{lineno} example/pm_fftnr/getFFTF/main.F90
    !>  \compilef{getFFTF}
    !>  \output{getFFTF}
    !>  \include{lineno} example/pm_fftnr/getFFTF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
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
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTF_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTF_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTF_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTF_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
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
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTF_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTF_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTF_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTF_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Inverse (normalized by `2 / size(data)`) Fourier Transform
    !>  of a periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of size [getExpNext(size(data), 2)](@ref pm_mathExp::getExpNext)
    !>                              of the same type and kind as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTI}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: getFFTI
    !>      type_of(data) :: fft(1 : getExpNext(size(data)))
    !>
    !>      fft(:) = getFFTI(data(:))
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTR](@ref pm_fftnr::setFFTR).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{getFFTI}
    !>  \include{lineno} example/pm_fftnr/getFFTI/main.F90
    !>  \compilef{getFFTI}
    !>  \output{getFFTI}
    !>  \include{lineno} example/pm_fftnr/getFFTI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
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
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTI_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTI_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTI_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTI_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
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
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTI_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTI_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTI_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTI_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Reverse (unnormalized) Fourier Transform of a
    !>  periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[in]      data    :   The input `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>
    !>  \return
    !>  `fft`                   :   The output vector of size [getExpNext(size(data), 2)](@ref pm_mathExp::getExpNext)
    !>                              of the same type and kind as the input `data`, containing the FFT result.<br>
    !>
    !>  \interface{getFFTR}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: getFFTR
    !>      type_of(data) :: fft(1 : getExpNext(size(data)))
    !>
    !>      fft(:) = getFFTR(data(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(data)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \note
    !>  This functional generic interface is simply a more flexible but slower wrapper
    !>  around the subroutine generic interface [setFFTR](@ref pm_fftnr::setFFTR).<br>
    !>  As such, this functional interface can be significantly slower than the corresponding subroutine interface.<br>
    !>
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{getFFTR}
    !>  \include{lineno} example/pm_fftnr/getFFTR/main.F90
    !>  \compilef{getFFTR}
    !>  \output{getFFTR}
    !>  \include{lineno} example/pm_fftnr/getFFTR/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
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
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK4_ENABLED
    impure module function getFFTR_CK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK3_ENABLED
    impure module function getFFTR_CK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK2_ENABLED
    impure module function getFFTR_CK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if CK1_ENABLED
    impure module function getFFTR_CK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: data(:)
        complex(CKG)                                :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
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
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK4_ENABLED
    impure module function getFFTR_RK4(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK3_ENABLED
    impure module function getFFTR_RK3(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK2_ENABLED
    impure module function getFFTR_RK2(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

#if RK1_ENABLED
    impure module function getFFTR_RK1(data) result(fft)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFFTR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: data(:)
        real(RKG)                                   :: fft(2**ceiling(log(real(size(data, 1, IK), RKB)) / log(2._RKB)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Forward Fourier Transform of a periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[inout]   data    :   The input/output `contiguous` vector of arbitrary size of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the periodic sequence whose FFT is to be computed.<br>
    !>                              On output, `data` contains the FFT result.<br>
    !>
    !>  \interface{setFFTF}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: setFFTF
    !>
    !>      call setFFTF(data(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(data) .and. isIntPow(size(data))` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{setFFTF}
    !>  \include{lineno} example/pm_fftnr/setFFTF/main.F90
    !>  \compilef{setFFTF}
    !>  \output{setFFTF}
    !>  \include{lineno} example/pm_fftnr/setFFTF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
    !>
    !>  \final{setFFTF}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTF_CK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTF_CK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTF_CK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTF_CK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTF_CK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTF_RK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTF_RK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTF_RK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTF_RK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTF_RK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Inverse (normalized by `2 / size(data)`) Fourier Transform of a periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[inout]   data        :   The input/output `contiguous` vector of arbitrary size of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the periodic sequence whose FFT is to be computed.<br>
    !>                                  On output, `data` contains the FFT result.<br>
    !>
    !>  \interface{setFFTI}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: setFFTI
    !>
    !>      call setFFTI(data(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(data) .and. isIntPow(size(data))` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{setFFTI}
    !>  \include{lineno} example/pm_fftnr/setFFTI/main.F90
    !>  \compilef{setFFTI}
    !>  \output{setFFTI}
    !>  \include{lineno} example/pm_fftnr/setFFTI/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
    !>
    !>  \final{setFFTI}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTI

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTI_CK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTI_CK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTI_CK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTI_CK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTI_CK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTI_RK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTI_RK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTI_RK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTI_RK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTI_RK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Reverse (unnormalized) Fourier Transform of a periodic sequence of type `complex` or `real` of arbitrary kind parameter.
    !>
    !>  \details
    !>  See the documentation of [pm_fftnr](@ref pm_fftnr) for more details.<br>
    !>
    !>  \param[inout]   data        :   The input/output `contiguous` vector of arbitrary size of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the periodic sequence whose FFT is to be computed.<br>
    !>                                  On output, `data` contains the FFT result.<br>
    !>
    !>  \interface{setFFTR}
    !>  \code{.F90}
    !>
    !>      use pm_fftnr, only: setFFTR
    !>
    !>      call setFFTR(data(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(data) .and. isIntPow(size(data))` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFFTF](@ref pm_fftnr::getFFTF)<br>
    !>  [getFFTR](@ref pm_fftnr::getFFTR)<br>
    !>  [getFFTI](@ref pm_fftnr::getFFTI)<br>
    !>  [setFFTF](@ref pm_fftnr::setFFTF)<br>
    !>  [setFFTR](@ref pm_fftnr::setFFTR)<br>
    !>  [setFFTI](@ref pm_fftnr::setFFTI)<br>
    !>
    !>  \example{setFFTR}
    !>  \include{lineno} example/pm_fftnr/setFFTR/main.F90
    !>  \compilef{setFFTR}
    !>  \output{setFFTR}
    !>  \include{lineno} example/pm_fftnr/setFFTR/main.out.F90
    !>
    !>  \test
    !>  [test_pm_fftnr](@ref test_pm_fftnr)
    !>
    !>  \final{setFFTR}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setFFTR

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setFFTR_CK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setFFTR_CK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setFFTR_CK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setFFTR_CK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setFFTR_CK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setFFTR_RK5(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setFFTR_RK4(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setFFTR_RK3(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setFFTR_RK2(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setFFTR_RK1(data)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setFFTR_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(inout) , contiguous    :: data(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_fftnr