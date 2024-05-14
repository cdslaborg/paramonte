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
!>  This module contains classes and procedures for computing properties related to the cross correlation of random samples.
!>
!>  \details
!>
!>  Cross Correlation
!>  =================
!>
!>  Cross-correlation is a measure of similarity of two series as a function of the displacement of one relative to the other.<br>
!>  This is also known as a **sliding dot product** or **sliding inner-product**.<br>
!>  It is commonly used for searching a long signal for a shorter, known feature.<br>
!>  It has applications in pattern recognition, single particle analysis, electron tomography, averaging, cryptanalysis, and neurophysiology.<br>
!>  The cross-correlation is similar in nature to the *convolution of two functions*.<br>
!>  In an **autocorrelation**, which is the cross-correlation of a signal with itself, there will always be a peak at a lag of zero, and its size will be the *signal energy*.<br>
!>
!>  The term cross-correlation also refers to the correlations between the entries of two random vectors \f$\mathbf{X}\f$ and \f$\mathbf {Y}\f$,
!>  while the correlations of a random vector \f$\mathbf{X}\f$ are the correlations between the entries of \f$\mathbf{X}\f$ itself,
!>  those forming the correlation matrix of \f$\mathbf{X}\f$.<br>
!>  If each of \f$\mathbf{X}\f$ and \f$\mathbf{Y}\f$ is a scalar random variable which is realized repeatedly in a time series,
!>  then the correlations of the various temporal instances of \f$\mathbf{X}\f$ are known as **autocorrelations** of \f$\mathbf{X}\f$,
!>  and the cross-correlations of \f$\mathbf{X}\f$ with \f$\mathbf{Y}\f$ across time are temporal cross-correlations.<br>
!>  In probability and statistics, the definition of correlation always includes a standardizing factor in such a way that correlations have values between −1 and +1.
!>
!>  Interpretation
!>  --------------
!>
!>  If \f$X\f$ and \f$Y\f$ are two independent random variables with probability density functions \f$f\f$ and \f$g\f$, respectively,
!>  then the probability density of the difference \f$Y − X\f$ is formally given by the **cross-correlation** \f$f\star g\f$.<br>
!>  Equivalently, the convolution \f$f * g\f$ (equivalent to the cross-correlation of \f$f(t)\f$ and \f$g(-t)\f$ gives the probability density function of the sum \f$X + Y\f$.<br>
!>
!>  Definition
!>  ----------
!>
!>  For continuous functions \f$f\f$ and \f$g\f$, the cross-correlation is defined as:<br>
!>  \f{equation}{
!>      (f\star g)(\tau) = \int_{-\infty }^{\infty }{\overline {f(t)}}g(t+\tau )\,dt
!>  \f}
!>  which is equivalent to,
!>  \f{equation}{
!>      (f\star g)(\tau ) = \int_{-\infty }^{\infty }{\overline {f(t-\tau )}}g(t)\, dt
!>  \f}
!>  where \f$\overline{f(t)}\f$ denotes the complex conjugate of \f$f(t)\f$, and \f$\tau\f$ is called **displacement** or **lag**.<br>
!>  For highly-correlated \f$f\f$ and \f$g\f$ which have a maximum cross-correlation at a particular \f$\tau\f$,
!>  a feature in \f$f\f$ at \f$t\f$ also occurs later in \f$g\f$ at \f$t + \tau\f$, hence \f$g\f$ **is said to lag** \f$f\f$ by \f$\tau\f$.<br>
!>
!>  Properties
!>  ----------
!>
!>  <ol>
!>      <li>    The cross-correlation of functions \f$f(t)\f$ and \f$g(t)\f$ is equivalent to the convolution (denoted by \f$*\f$) of \f$\overline{f(-t)}\f$ and \f$g(t)\f$.<br>
!>              That is:<br>
!>              \f{equation}{
!>                  [f(t)\star g(t)](t) = [{\overline{f(-t)}} * g(t)](t) ~.
!>              \f}
!>      <li>    \f$[f(t)\star g(t)](t) = [{\overline{g(t)}}\star{\overline {f(t)}}](-t)\f$ ~.
!>      <li>    If \f$f\f$ is a Hermitian function, then \f$f\star g = f*g\f$.
!>      <li>    If both \f$f\f$ and \f$g\f$ are Hermitian, then \f$f\star g = g\star f\f$.
!>      <li>    \f$\left(f\star g\right)\star \left(f\star g\right) = \left(f\star f\right)\star \left(g\star g\right)\f$.
!>      <li>    Analogous to the convolution theorem, the cross-correlation satisfies,
!>              \f{equation}{
!>                  \mathcal{F} \left\{f\star g\right\} = {\overline{{\mathcal {F}}\left\{f\right\}}} \cdot {\mathcal{F}}\left\{g\right\} ~,
!>              \f}
!>              where \f${\mathcal{F}}\f$ denotes the Fourier transform, and an \f${\overline {f}}\f$ again indicates the complex conjugate of \f$f\f$,
!>              since \f$\mathcal{F}\left\{{\overline {f(-t)}}\right\} = {\overline{{\mathcal {F}}\left\{f(t)\right\}}}\f$.<br>
!>              Coupled with fast Fourier transform algorithms, this property is often exploited for the efficient numerical computation of cross-correlations.<br>
!>      <li>    The cross-correlation is related to the spectral density (see Wiener–Khinchin theorem).<br>
!>      <li>    The cross-correlation of a convolution of \f$f\f$ and \f$h\f$ with a function \f$g\f$ is the convolution of the cross-correlation of \f$g\f$ and \f$f\f$ with the kernel \f$h\f$:<br>
!>              \f{equation}{
!>                  g\star \left(f * h\right) = \left(g\star f\right) * h.
!>              \f}
!>  </ol>
!>
!>  Computation
!>  ===========
!>
!>  The fastest methods for computing for Cross-correlation of large sequences rely on the fast-fourier-transform and the correlation/convolution theorems.<br>
!>  <ol>
!>      <li>    Convert the input signals \f$f\f$ and \f$g\f$ signals to signals of type `complex`.<br>
!>      <li>    Left-pad the input signals \f$f\f$ and \f$g\f$ with zeros such that the sequences become of minimum length `size(f) + size(g) - 1`.<br>
!>              Such padding is to ensure the periodic nature of the FFT method does not affect the computed correlations by overlapping the sequences circularly.<br>
!>      <li>    [Compute the FFT](@ref pm_fftpack::getFFTF) of both signals and multiply the results of the first FFT with the conjugate of the second.<br>
!>      <li>    [Compute the inverse FFT](@ref pm_fftpack::getFFTI) of the resulting elemental multiplication.<br>
!>  </ol>
!>  The above steps can be summarized as `getFFTI(getFFTF(getPaddedl(f, size(f) + size(g) - 1, (0., 0.))) + conjg(getFFTF(getPaddedl(f, size(f) + size(g) - 1, (0., 0.)))))`
!>  assuming both `f` and `g` are already of type `complex`.<br>
!>  The resulting slice `ccf(1 : size(f))` contains the CCF for lags `[(lag, lag = 0, size(f) - 1)]` and
!>  the slice `ccf(size(f) + 1 : size(f) + size(g) - 1)` contains the CCF for lags `[(lag, lag = -size(g), -1)]`.<br>
!>
!>  The convolution theorem
!>  -----------------------
!>
!>  The convolution theorem states that under suitable conditions the Fourier transform of a
!>  convolution of two functions (or signals) is the pointwise product of their Fourier transforms.<br>
!>  More generally, convolution in one domain (e.g., time domain) equals point-wise multiplication in the other domain (e.g., frequency domain).<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>  [Cross Correlation](https://en.wikipedia.org/wiki/Cross-correlation)<br>
!>
!>  \test
!>  [test_pm_sampleCCF](@ref test_pm_sampleCCF)
!>
!>  \finmain{pm_sampleCCF}
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 01:45 AM, August 21, 2018, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleCCF

    use pm_kind, only: SK, IK, LK
    use pm_container, only: css_type
    use pm_array, only: nothing, nothing_type
    use pm_sampleNorm, only: zscore, zscore_type
    use pm_sampleScale, only: stdscale, stdscale_type
    use pm_sampleShift, only: meanshift, meanshift_type
    use pm_matrixSubset, only: uppDia_type, uppDia
    use pm_matrixSubset, only: lowDia_type, lowDia
    use pm_matrixSubset, only: upp_type, upp
    use pm_matrixSubset, only: low_type, low
    use pm_fftpack, only: getFactorFFT

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleCCF"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the auto-correlation function (ACF) \f$(f\star f)(\tau)\f$
    !>  of the discrete signal \f$f\f$ lagging itself for a range of lags.<br>
    !>
    !>  \details
    !>  This generic interface a flexible wrapper for the lower-level potentially faster generic interface [setACF](@ref pm_sampleCCF::setACF).<br>
    !>  See the documentation of the parent module [pm_sampleCCF](@ref pm_sampleCCF) for algorithmic details and auto-correlation definition.<br>
    !>
    !>  \note
    !>  By default, the input sequence `f` is right-padded with zeros such that the resulting new sequence
    !>  has the length `max(lag(size(lag)), size(f) - 1) - min(lag(1), 1 - size(f)) + 1` before computing the auto-correlation.<br>
    !>
    !>  \param[in]  f           :   The input `contiguous` vector of arbitrary size (of minimum `2`) of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the first sequence in the auto-correlation computation.<br>
    !>  \param[in]  lag         :   The input vector of type `integer` of default kind \IK of arbitrary size,
    !>                              containing the set of lags in **ascending order** for which the correlation must be computed.<br>
    !>                              (**optional**. default = `getRange(0, size(f) - 1)`)
    !>  \param[in]  norm        :   The input scalar constant that can be:<br>
    !>                              <ol>
    !>                                  <li>    the constant [nothing](@ref pm_array::nothing) or an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                          implying that the input `f` must be used as is without any normalization.<br>
    !>                                  <li>    the constant [meanshift](@ref pm_sampleShift::meanshift) or an object of type [meanshift_type](@ref pm_sampleShift::meanshift_type),
    !>                                          implying that the input `f` must be mean-shifted before computing the auto-correlation.<br>
    !>                                  <li>    the constant [stdscale](@ref pm_sampleScale::stdscale) or an object of type [stdscale_type](@ref pm_sampleScale::stdscale_type),
    !>                                          implying that the input `f` must be scaled by a factor of `1 / sum(abs(f)**2)`.<br>
    !>                                          This option is particularly useful when the input sequence is already mean-shifted but not properly scaled.<br>
    !>                                  <li>    the constant [zscore](@ref pm_sampleNorm::zscore) or an object of type [zscore_type](@ref pm_sampleNorm::zscore_type),
    !>                                          implying that the input `f` must be mean-shifted and further scaled by a factor of `1 / sum(abs(f)**2)`.<br>
    !>                              </ol>
    !>                              (**optional**, default = [zscore](@ref pm_sampleNorm::zscore).)
    !>
    !>  \return
    !>  `acf`                   :   The output `allocatable` vector of the same type and kind as the input `f` of size `size(lag)`
    !>                              containing the auto-correlation of the input `f` for the specified input lags.<br>
    !>                              Note that `acf` is always real-valued at lag `0` by definition, even when the input `f` is of type `complex`.<br>
    !>                              For all non-zero lags, the imaginary component of the auto-correlation function is an odd function.<br>
    !>
    !>  \interface{getACF}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCCF, only: getACF
    !>
    !>      acf(1:size(lag)) = getACF(f(:), lag = lag, norm = norm)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `isAscending(lag)` must hold for the corresponding input arguments.<br>
    !>  This generic interface is merely a flexible wrapper around the generic `subroutine` interface [setACF](@ref pm_sampleCCF::setACF).<br>
    !>  As such, all conditions and warnings associated with [setACF](@ref pm_sampleCCF::setACF) equally apply to this generic interface.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getACF](@ref pm_sampleCCF::getACF)<br>
    !>  [setACF](@ref pm_sampleCCF::setACF)<br>
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getACF}
    !>  \include{lineno} example/pm_sampleCCF/getACF/main.F90
    !>  \compilef{getACF}
    !>  \output{getACF}
    !>  \include{lineno} example/pm_sampleCCF/getACF/main.out.F90
    !>  \postproc{getACF}
    !>  \include{lineno} example/pm_sampleCCF/getACF/main.py
    !>  \vis{getACF}
    !>  \image html pm_sampleCCF/getACF/getACF.crd.sin.RK.png width=700
    !>  \image html pm_sampleCCF/getACF/getACF.acf.sin.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleCCF](@ref test_pm_sampleCCF)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getACF_D1_CK5()
    !>         ||| || |||
    !>         ||| || The type and kind parameters of the input sequence.
    !>         ||| The dimension of the input sequence `f`.
    !>         ACF: Cross-Correlation Function.
    !>  \endcode
    !>
    !>  \finmain{getACF}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! D1

    interface getACF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getACF_D1_CK5(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:)
        complex(CKC)            , allocatable                           :: acf(:)
    end function
#endif

#if CK4_ENABLED
    module function getACF_D1_CK4(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:)
        complex(CKC)            , allocatable                           :: acf(:)
    end function
#endif

#if CK3_ENABLED
    module function getACF_D1_CK3(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:)
        complex(CKC)            , allocatable                           :: acf(:)
    end function
#endif

#if CK2_ENABLED
    module function getACF_D1_CK2(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:)
        complex(CKC)            , allocatable                           :: acf(:)
    end function
#endif

#if CK1_ENABLED
    module function getACF_D1_CK1(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:)
        complex(CKC)            , allocatable                           :: acf(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getACF_D1_RK5(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:)
        real(RKC)               , allocatable                           :: acf(:)
    end function
#endif

#if RK4_ENABLED
    module function getACF_D1_RK4(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:)
        real(RKC)               , allocatable                           :: acf(:)
    end function
#endif

#if RK3_ENABLED
    module function getACF_D1_RK3(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:)
        real(RKC)               , allocatable                           :: acf(:)
    end function
#endif

#if RK2_ENABLED
    module function getACF_D1_RK2(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:)
        real(RKC)               , allocatable                           :: acf(:)
    end function
#endif

#if RK1_ENABLED
    module function getACF_D1_RK1(f, lag, norm) result(acf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getACF_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:)
        real(RKC)               , allocatable                           :: acf(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getACF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the auto-correlation function (ACF) \f$(f\star f)(\tau)\f$ of the discrete signal \f$f\f$
    !>  lagging itself for a range of lags that spans the sequence length.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_sampleCCF](@ref pm_sampleCCF) for algorithmic details and sample correlation matrix definition.<br>
    !>  Note that `acf` is always real-valued at lag `0` by definition, even when the input `f` is of type `complex`.<br>
    !>  For all non-zero lags, the imaginary component of the auto-correlation function is an odd function.<br>
    !>
    !>  \warning
    !>  This low-level generic interface neither right-pads nor normalizes the input sequence in any form.<br>
    !>  The input sequence is used as is for computing the ACF.<br>
    !>  In general, the following steps should be taken before and after computing the ACF using this routine,
    !>  <ol>
    !>      <li>    If the original (unpadded) sequence is not already mean-shifted (with zero mean),
    !>              it is highly recommended to mean-shift the sequence first before computing the ACF.<br>
    !>              This can be readily done by the routines of [pm_sampleMean](@ref pm_sampleMean) and [pm_sampleShift](@ref pm_sampleShift).<br>
    !>      <li>    While the input mean-shifted sequence can be passed to this generic interface as is,
    !>              it is highly recommended to right-pad the mean-shifted sequence by zeros to the minimum length of `size(f) * 2 - 1`.<br>
    !>              This means that the length of the input sequence should always be an odd number, although not necessarily, if you know what you are doing.<br>
    !>      <li>    The output `acf` must from the generic interface must be multiplied by `1 / f(1)` to properly normalize the output `acf` to the range `[-1, +1]`.<br>
    !>              While easy, this can also be readily done via the generic interfaces of [pm_sampleScale](@ref pm_sampleScale).<br>
    !>  </ol>
    !>
    !>  \note
    !>  When the input sequence is right-padded such that `size(acf) = size(f) * 2 - 1`, then,
    !>  <ol>
    !>      <li>    the resulting slice `acf(1 : size(f))` by this generic interface contains the auto-correlation corresponding to lags `[(i, i = 0, size(f))]`.
    !>      <li>    the resulting slice `acf(size(f) + 1 : size(acf))` contains the auto-correlation corresponding to lags `[(i, i = -size(f) + 1, -1)]`.
    !>  </ol>
    !>
    !>  \param[in]      factor  :   The input `contiguous` vector of shape `(:)` of type `integer` of default kind \IK,
    !>                              containing the factorization of the length of the input sequence `f` whose FFT is to be computed.<br>
    !>                              This input argument along with the input argument `coef` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[in]      coef    :   The input `contiguous` vector of shape `(1:size(data))` of the same type and kind as the input argument `f`,
    !>                              containing the trigonometric look up table required for FFT of the specified sequence.<br>
    !>                              This input argument along with `factor` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[inout]   f       :   The input `contiguous` vector of arbitrary size (of minimum `2`) of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the first sequence in the auto-correlation computation.<br>
    !>                              On output, the contents of `f` are destroyed.<br>
    !>                              If the output `inf = .true.`, then `f` contains the resulting unnormalized ACF of the input sequence.<br>
    !>  \param[out]     work    :   The output `contiguous` vector of the same type, kind, and size as the input `f`, used as a workspace.<br>
    !>                              If the condition `inf = .false.` holds on output, then `work` contains the resulting unnormalized ACF of the input sequence `f`.<br>
    !>  \param[out]     inf     :   The output scalar of type `logical` of default kind \LK.<br>
    !>                              <ol>
    !>                                  <li>    If `.true.`, the resulting ACF is stored in the output argument `f` upon return from the procedure.<br>
    !>                                  <li>    If `.false.`, the resulting ACF is stored in the output argument `work` upon return from the procedure.<br>
    !>                              </ol>
    !>
    !>  \interface{getACF}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCCF, only: setACF
    !>
    !>      call setACF(factor(:), coef(1:nseq), f(1:nseq), work(1:nseq), inf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `2 < size(f)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(f) == size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(f) == size(work)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getACF](@ref pm_sampleCCF::getACF)<br>
    !>  [setACF](@ref pm_sampleCCF::setACF)<br>
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{setACF}
    !>  \include{lineno} example/pm_sampleCCF/setACF/main.F90
    !>  \compilef{setACF}
    !>  \output{setACF}
    !>  \include{lineno} example/pm_sampleCCF/setACF/main.out.F90
    !>  \postproc{setACF}
    !>  \include{lineno} example/pm_sampleCCF/setACF/main.py
    !>  \vis{setACF}
    !>  \image html pm_sampleCCF/setACF/setACF.crd.sin.RK.png width=700
    !>  \image html pm_sampleCCF/setACF/setACF.acf.sin.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleCCF](@ref test_pm_sampleCCF)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setACF_FP_D1_CK5()
    !>         ||| || || |||
    !>         ||| || || The type and kind parameters of the input sequence.
    !>         ||| || The dimension of the input sequence `f`.
    !>         ||| The method used: FP => fftpack.
    !>         ACF: Cross-Correlation Function.
    !>  \endcode
    !>
    !>  \finmain{setACF}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin<br>

    ! FP - D1

    interface setACF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setACF_FP_D1_CK5(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setACF_FP_D1_CK4(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setACF_FP_D1_CK3(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setACF_FP_D1_CK2(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setACF_FP_D1_CK1(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setACF_FP_D1_RK5(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setACF_FP_D1_RK4(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setACF_FP_D1_RK3(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setACF_FP_D1_RK2(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setACF_FP_D1_RK1(factor, coef, f, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setACF_FP_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setACF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cross-correlation function (CCF) \f$(f\star g)(\tau)\f$
    !>  of the discrete signal \f$g\f$ lagging the signal the discrete signal \f$f\f$ for a range of specified lags.<br>
    !>
    !>  \details
    !>  This generic interface a flexible wrapper for the lower-level potentially faster generic interface [setCCF](@ref pm_sampleCCF::setCCF).<br>
    !>  See the documentation of the parent module [pm_sampleCCF](@ref pm_sampleCCF) for algorithmic details and cross-correlation definition.<br>
    !>
    !>  \note
    !>  The input sequences `f` and `g` are right-padded with zeros to a minimum length of `max(lag(size(lag)), size(f) - 1) - min(lag(1), 1 - size(g)) + 1`
    !>  before computing the cross-correlation.<br>
    !>
    !>  \param[in]  f           :   The input `contiguous` vector of arbitrary size (of minimum `2`) of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the first sequence in the cross-correlation computation.<br>
    !>  \param[in]  g           :   The input `contiguous` vector of arbitrary size (of minimum `2`) of the same type and kind as the input argument `f`,
    !>                              containing the second sequence in the cross-correlation computation.<br>
    !>  \param[in]  lag         :   The input vector of type `integer` of default kind \IK of arbitrary size,
    !>                              containing the set of lags in **ascending order** for which the correlation must be computed.<br>
    !>                              (**optional**. default = `getRange(1 - size(g), size(f) - 1)`)
    !>  \param[in]  norm        :   The input scalar constant that can be:<br>
    !>                              <ol>
    !>                                  <li>    the constant [nothing](@ref pm_array::nothing) or an object of type [nothing_type](@ref pm_array::nothing_type),
    !>                                          implying that the input `f` and `g` must be used as is without any normalization.<br>
    !>                                  <li>    the constant [meanshift](@ref pm_sampleShift::meanshift) or an object of type [meanshift_type](@ref pm_sampleShift::meanshift_type),
    !>                                          implying that the input `f` and `g` must be mean-shifted before computing the cross-correlation.<br>
    !>                                  <li>    the constant [stdscale](@ref pm_sampleScale::stdscale) or an object of type [stdscale_type](@ref pm_sampleScale::stdscale_type),
    !>                                          implying that the input `f` and `g` must be scaled by a factor of `1 / sqrt(sum(abs(f)**2)) / sqrt(sum(abs(g)**2))`.<br>
    !>                                          This option is particularly useful when the input sequences are already mean-shifted but not properly scaled.<br>
    !>                                  <li>    the constant [zscore](@ref pm_sampleNorm::zscore) or an object of type [zscore_type](@ref pm_sampleNorm::zscore_type),
    !>                                          implying that the input `f` and `g` must be mean-shifted and further scaled by a
    !>                                          factor of `1 / sqrt(sum(abs(f)**2)) / sqrt(sum(abs(g)**2))`.<br>
    !>                              </ol>
    !>                              (**optional**, default = [zscore](@ref pm_sampleNorm::zscore).)
    !>
    !>  \return
    !>  `ccf`                   :   The output `allocatable` vector of the same type and kind as the input `f` of size `max(lag(size(lag)), size(f) - 1) - min(lag(1), 1 - size(g)) + 1`,
    !>                              containing the cross-correlation of the input sequences `f` and `g`.
    !>
    !>  \interface{getCCF}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCCF, only: getCCF
    !>
    !>      ccf(:) = getCCF(f(:), g(:), lag = lag, norm = norm)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `isAscending(lag)` must hold for the corresponding input arguments.<br>
    !>  This generic interface is merely a flexible wrapper around the generic `subroutine` interface [setCCF](@ref pm_sampleCCF::setCCF).<br>
    !>  As such, all conditions and warnings associated with [setCCF](@ref pm_sampleCCF::setCCF) equally apply to this generic interface.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getACF](@ref pm_sampleCCF::getACF)<br>
    !>  [setACF](@ref pm_sampleCCF::setACF)<br>
    !>  [getCCF](@ref pm_sampleCCF::getCCF)<br>
    !>  [setCCF](@ref pm_sampleCCF::setCCF)<br>
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getCCF}
    !>  \include{lineno} example/pm_sampleCCF/getCCF/main.F90
    !>  \compilef{getCCF}
    !>  \output{getCCF}
    !>  \include{lineno} example/pm_sampleCCF/getCCF/main.out.F90
    !>  \postproc{getCCF}
    !>  \include{lineno} example/pm_sampleCCF/getCCF/main.py
    !>  \vis{getCCF}
    !>  \image html pm_sampleCCF/getCCF/getCCF.crd.sin.cos.RK.png width=700
    !>  \image html pm_sampleCCF/getCCF/getCCF.ccf.sin.cos.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleCCF](@ref test_pm_sampleCCF)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setCCF_FG_CK5()
    !>         ||| || |||
    !>         ||| || The type and kind parameters of the input sequences.
    !>         ||| The input sequences `f` and `g`.
    !>         CCF: Cross-Correlation Function.
    !>  \endcode
    !>
    !>  \finmain{getCCF}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! FG

    interface getCCF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getCCF_FG_CK5(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_CK5
#endif
        use pm_kind, only: CKC => CK5
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:), g(:)
        complex(CKC)            , allocatable                           :: ccf(:)
    end function
#endif

#if CK4_ENABLED
    module function getCCF_FG_CK4(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_CK4
#endif
        use pm_kind, only: CKC => CK4
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:), g(:)
        complex(CKC)            , allocatable                           :: ccf(:)
    end function
#endif

#if CK3_ENABLED
    module function getCCF_FG_CK3(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_CK3
#endif
        use pm_kind, only: CKC => CK3
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:), g(:)
        complex(CKC)            , allocatable                           :: ccf(:)
    end function
#endif

#if CK2_ENABLED
    module function getCCF_FG_CK2(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_CK2
#endif
        use pm_kind, only: CKC => CK2
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:), g(:)
        complex(CKC)            , allocatable                           :: ccf(:)
    end function
#endif

#if CK1_ENABLED
    module function getCCF_FG_CK1(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_CK1
#endif
        use pm_kind, only: CKC => CK1
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        complex(CKC)            , intent(in)    , contiguous            :: f(:), g(:)
        complex(CKC)            , allocatable                           :: ccf(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getCCF_FG_RK5(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_RK5
#endif
        use pm_kind, only: RKC => RK5
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:), g(:)
        real(RKC)               , allocatable                           :: ccf(:)
    end function
#endif

#if RK4_ENABLED
    module function getCCF_FG_RK4(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_RK4
#endif
        use pm_kind, only: RKC => RK4
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:), g(:)
        real(RKC)               , allocatable                           :: ccf(:)
    end function
#endif

#if RK3_ENABLED
    module function getCCF_FG_RK3(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_RK3
#endif
        use pm_kind, only: RKC => RK3
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:), g(:)
        real(RKC)               , allocatable                           :: ccf(:)
    end function
#endif

#if RK2_ENABLED
    module function getCCF_FG_RK2(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_RK2
#endif
        use pm_kind, only: RKC => RK2
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:), g(:)
        real(RKC)               , allocatable                           :: ccf(:)
    end function
#endif

#if RK1_ENABLED
    module function getCCF_FG_RK1(f, g, lag, norm) result(ccf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCCF_FG_RK1
#endif
        use pm_kind, only: RKC => RK1
        class(*)                , intent(in)                , optional  :: norm
        integer(IK)             , intent(in)    , contiguous, optional  :: lag(:)
        real(RKC)               , intent(in)    , contiguous            :: f(:), g(:)
        real(RKC)               , allocatable                           :: ccf(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the cross-correlation function (CCF) \f$(f\star g)(\tau)\f$
    !>  of the discrete signal \f$g\f$ lagging the signal the discrete signal \f$f\f$ for a range of
    !>  lags that spans the maximum of the lengths of the two sequences.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_sampleCCF](@ref pm_sampleCCF) for algorithmic details and sample correlation matrix definition.<br>
    !>
    !>  \warning
    !>  This low-level generic interface neither right-pads nor normalizes the input sequences in any form.<br>
    !>  The input sequences are used as is for computing the CCF.<br>
    !>  In general, the following steps should be taken before and after computing the CCF using this routine,
    !>  <ol>
    !>      <li>    If the original (unpadded) sequences are not already mean-shifted (with zero mean),
    !>              it is highly recommended to mean-shift the sequences first before computing the CCF.<br>
    !>              This can be readily done by the routines of [pm_sampleMean](@ref pm_sampleMean) and [pm_sampleShift](@ref pm_sampleShift).<br>
    !>      <li>    While the input mean-shifted sequences can be passed to this generic interface as is,
    !>              it is highly recommended to right-pad the mean-shifted sequences by zeros to the minimum length of `size(f) + size(g) - 1`.<br>
    !>      <li>    The output CCF must from the generic interface must be multiplied by `1 / size(ccf) / sqrt(sum(abs(f)**2) + sum(abs(g)**2))`
    !>              to properly normalize the output CCF to the range `[-1, +1]`.<br>
    !>              If the input sequences have unit-variances, the normalization factor reduces to `1 / size(ccf)`.<br>
    !>              While easy, this rescaling can also be readily done via the generic interfaces of [pm_sampleScale](@ref pm_sampleScale).<br>
    !>  </ol>
    !>
    !>  \note
    !>  When the input sequences are right-padded to the size `size(f) + size(g) - 1`, then,
    !>  <ol>
    !>      <li>    the resulting slice `ccf(1 : size(g))` by this generic interface contains the cross-correlation corresponding to lags `[(i, i = 0, size(g))]`.
    !>      <li>    the resulting slice `ccf(size(g) + 1 : size(ccf))` contains the cross-correlation corresponding to lags `[(i, i = -size(f) + 1, -1)]`.
    !>  </ol>
    !>
    !>  \param[in]      factor  :   The input `contiguous` vector of shape `(:)` of type `integer` of default kind \IK,
    !>                              containing the factorization of the length of the input sequence `f` whose FFT is to be computed.<br>
    !>                              This input argument along with the input argument `coef` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[in]      coef    :   The input `contiguous` vector of shape `(1:size(data))` of the same type and kind as the input argument `f`,
    !>                              containing the trigonometric look up table required for FFT of the specified sequences.<br>
    !>                              This input argument along with `factor` is the direct output of [getFactorFFT](@ref pm_fftpack::getFactorFFT).<br>
    !>  \param[inout]   f       :   The input `contiguous` vector of arbitrary size (of minimum `2`) of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the first sequence in the cross-correlation computation.<br>
    !>                              On output, the contents of `f` are destroyed.<br>
    !>                              If the output `inf = .true.`, then `f` contains the resulting unnormalized CCF of the input sequences.<br>
    !>  \param[inout]   g       :   The input `contiguous` vector of the same type and kind and size as the input argument `f`,
    !>                              containing the second sequence in the cross-correlation computation.<br>
    !>                              On output, the contents of `g` are destroyed.<br>
    !>                              If the condition `inf = .false.` holds on output, then `g` contains the resulting unnormalized CCF of the input sequences.<br>
    !>  \param[out]     work    :   The output `contiguous` vector of the same type, kind, and size as the input `f`, used as a workspace.<br>
    !>  \param[out]     inf     :   The output scalar of type `logical` of default kind \LK.<br>
    !>                              <ol>
    !>                                  <li>    If `.true.`, the resulting `ccf` is stored in the output argument `f` upon return from the procedure.<br>
    !>                                  <li>    If `.false.`, the resulting `ccf` is stored in the output argument `g` upon return from the procedure.<br>
    !>                              </ol>
    !>
    !>  \interface{getCCF}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCCF, only: setCCF
    !>
    !>      call setCCF(factor(:), coef(1:nseq), f(1:nseq), g(1:nseq), work(1:nseq), inf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `2 < size(f)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(f) == size(g)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(f) == size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(f) == size(work)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getACF](@ref pm_sampleCCF::getACF)<br>
    !>  [setACF](@ref pm_sampleCCF::setACF)<br>
    !>  [getCCF](@ref pm_sampleCCF::getCCF)<br>
    !>  [setCCF](@ref pm_sampleCCF::setCCF)<br>
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{setCCF}
    !>  \include{lineno} example/pm_sampleCCF/setCCF/main.F90
    !>  \compilef{setCCF}
    !>  \output{setCCF}
    !>  \include{lineno} example/pm_sampleCCF/setCCF/main.out.F90
    !>  \postproc{setCCF}
    !>  \include{lineno} example/pm_sampleCCF/setCCF/main.py
    !>  \vis{setCCF}
    !>  \image html pm_sampleCCF/setCCF/setCCF.crd.sin.cos.RK.png width=700
    !>  \image html pm_sampleCCF/setCCF/setCCF.ccf.sin.cos.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleCCF](@ref test_pm_sampleCCF)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setCCF_FP_FG_CK5()
    !>         ||| || || |||
    !>         ||| || || The type and kind parameters of the input sequences.
    !>         ||| || The input sequences `f` and `g`.
    !>         ||| The method used: FP => fftpack.
    !>         CCF: Cross-Correlation Function.
    !>  \endcode
    !>
    !>  \finmain{setCCF}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin<br>

    ! FP - FG

    interface setCCF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCCF_FP_FG_CK5(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_CK5
#endif
        use pm_kind, only: CKC => CK5
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:), g(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCCF_FP_FG_CK4(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_CK4
#endif
        use pm_kind, only: CKC => CK4
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:), g(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCCF_FP_FG_CK3(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_CK3
#endif
        use pm_kind, only: CKC => CK3
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:), g(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCCF_FP_FG_CK2(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_CK2
#endif
        use pm_kind, only: CKC => CK2
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:), g(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCCF_FP_FG_CK1(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_CK1
#endif
        use pm_kind, only: CKC => CK1
        logical(LK)             , intent(out)                   :: inf
        complex(CKC)            , intent(out)   , contiguous    :: work(:)
        complex(CKC)            , intent(inout) , contiguous    :: f(:), g(:)
        complex(CKC)            , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCCF_FP_FG_RK5(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_RK5
#endif
        use pm_kind, only: RKC => RK5
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:), g(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCCF_FP_FG_RK4(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_RK4
#endif
        use pm_kind, only: RKC => RK4
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:), g(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCCF_FP_FG_RK3(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_RK3
#endif
        use pm_kind, only: RKC => RK3
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:), g(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCCF_FP_FG_RK2(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_RK2
#endif
        use pm_kind, only: RKC => RK2
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:), g(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCCF_FP_FG_RK1(factor, coef, f, g, work, inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCCF_FP_FG_RK1
#endif
        use pm_kind, only: RKC => RK1
        logical(LK)             , intent(out)                   :: inf
        real(RKC)               , intent(out)   , contiguous    :: work(:)
        real(RKC)               , intent(inout) , contiguous    :: f(:), g(:)
        real(RKC)               , intent(in)    , contiguous    :: coef(:)
        integer(IK)             , intent(in)    , contiguous    :: factor(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleCCF ! LCOV_EXCL_LINE