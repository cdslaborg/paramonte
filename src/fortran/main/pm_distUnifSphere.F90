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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Uniform Spherical distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of
!>  the uniform distribution of points on the surface of a an arbitrary `n`-dimensional sphere:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>  </ol>
!>
!>  n-sphere
!>  --------
!>
!>  An \f$n\f$-sphere or a hypersphere is a topological space that is homeomorphic to a standard \f$n\f$-sphere,
!>  which is the set of points in \f$(n + 1)\f$-dimensional Euclidean space that are situated at a constant distance \f$r\f$ from a fixed point, called the **center**.<br>
!>  It is the generalization of an ordinary sphere in the ordinary three-dimensional space.<br>
!>  The *radius* of a sphere is the constant distance of its points to the center.<br>
!>  When the sphere has unit radius, it is usual to call it the unit \f$n\f$-sphere or simply the \f$n\f$-sphere for brevity.<br>
!>  In terms of the standard norm, the \f$n\f$-sphere is defined as,<br>
!>  \f{equation}{
!>      S^{n} = \left\{x\in\mathbb{R}^{n+1}:\left\|x\right\|=1\right\} ~,
!>  \f}
!>  and an n-sphere of radius \f$r\f$ can be defined as
!>  \f{equation}{
!>      S^{n}(r)=\left\{x\in \mathbb {R} ^{n+1}:\left\|x\right\|=r\right\} ~.
!>  \f}
!>  The dimension of \f$n\f$-sphere is \f$n\f$, and must not be confused with the dimension \f$(n + 1)\f$ of the Euclidean space in which it is naturally embedded.<br>
!>  An \f$n\f$-sphere is the surface or boundary of an \f$(n + 1)\f$-dimensional ball.<br>
!>  In particular:<br>
!>  <ol>
!>      <li>    The pair of points at the ends of a (one-dimensional) line segment is a 0-sphere.
!>      <li>    A circle, which is the one-dimensional circumference of a (two-dimensional) ball, is a 1-sphere.
!>      <li>    The two-dimensional surface of a three-dimensional ball is a 2-sphere, commonly called a sphere.
!>      <li>    The three-dimensional boundary of a (four-dimensional) 4-ball is a 3-sphere.
!>      <li>    The \f$(n – 1)\f$-dimensional boundary of a (\f$n\f$-dimensional) \f$\ndim\f$-ball is an \f$(n – 1)\f$-sphere.
!>  </ol>
!>
!>  Probability Density Function (PDF)
!>  ----------------------------------
!>
!>  The PDF of an \f$n\f$-sphere can be readily computed as the inverse of the surface area of the sphere, which has the closed form,
!>  \f{equation}{
!>      S_{n-1} = {\frac {n\pi^{\frac{n}{2}}}{\Gamma\left({\frac{n}{2}} + 1\right)}} r^{n-1} = {\frac{2\pi^{\frac{n}{2}}}{\Gamma\left({\frac{n}{2}}\right)}} r^{n-1} = n\frac{V_{n}}{r}~.
!>  \f}
!>  for the \f$(n−1)\f$-dimensional surface of the sphere \f$S_{n−1}\f$ of radius \f$r\f$, where \f$V_{n + 1}\f$ is the volume of the corresponding \f$n\f$-ball.<br>
!>  Unlike the volume and the surface area of an \f$n\f$-sphere, the computation of the surface area of an arbitrary hyper-ellipsoid is very much involved.<br>
!>  See, for example, the discussion [here](https://analyticphysics.com/Higher%20Dimensions/Ellipsoids%20in%20Higher%20Dimensions.htm) for more information.<br>
!>
!>  Random Number Generation
!>  ------------------------
!>
!>  The **RNG** generic interfaces within this module generate uniformly-distributed random vectors on the surface an \f$n\f$-sphere
!>  using a generalization of the proposed approach of Marsaglia (1972) for choosing a point from the surface of a sphere.<br>
!>  <ol>
!>      <li>    Define \f$\ndim = n + 1\f$ as the dimension of the space within which the \f$n\f$-sphere is embedded.
!>      <li>    Generate a normalized (unit) \f$\ndim\f$-dimensional [Multivariate Normal random vector](@ref pm_distMultiNorm),
!>              \f{equation}{
!>                  \unit{r} = \frac{1}{\sqrt{\sum_{i=1}^\ndim g_i^2}} \sum_{j=1}^{\ndim} g_j \unit{r}_j ~,
!>              \f}
!>              where \f$\unit{r}_j\f$ are the unit vectors representing the orthogonal basis the of \f$\ndim\f$-space.<br>
!>              This unit vector \f$\unit{r}\f$ has uniform random orientation and distribution on the surface of an \f$n\f$-sphere
!>              of unit radius centered at the origin of the coordinate system.<br>
!>      <li>    **Optionally**, this random uniform surface point can be further transformed to the surface of an \f$\ndim\f$-ellipsoid via an affine transformation.<br>
!>              However, the transformation of an \f$\ndim\f$-ball to an \f$\ndim\f$-ellipsoid does not preserve the uniformity of the random points on the ellipsoid *surface*.<br>
!>              In other words, the resulting random points on the surface of the \f$\ndim\f$-ellipsoid are not uniformly distributed.<br>
!>              Instead, the random points will be clustered toward the sharp corners of the \f$\ndim\f$-ellipsoid, i.e., along the directions with the largest eigenvalues.<br>
!>              To transform the uniform \f$n\f$-sphere to the a non-uniform distribution on the surface of an \f$\ndim\f$-ellipsoid:<br>
!>              <ol>
!>                  <li>    First compute the [Cholesky decomposition](@ref pm_matrixChol) of the representative Gramian matrix \f$\gramian\f$ of the ellipsoid,
!>                          \f{equation}{
!>                              \gramian = \mat{L}\mat{L}^H ~,
!>                          \f}
!>                          where \f$\mat{L}\f$ is the left triangular matrix resulting from the Cholesky factorization, and \f$\mat{L}^T\f$ is its Hermitian transpose.<br>
!>                  <li>    Then the vector,
!>                          \f{equation}{
!>                              \bs{r}_\ell = \mat{L} ~ \bs{r}_\sphere + \mu_\ell ~,
!>                          \f}
!>                          is **non-uniformly** distributed on the surface of ellipsoid \f$\ell\f$ centered at the original \f$n\f$-sphere.<br>
!>                          The actual distribution of the resulting points on the ellipsoidal surface is a complicated
!>                          function of the curvature of ellipsoidal surface at different points.<br>
!>                          Unlike the case of [uniformly-distributed points within a hyper-ellipsoid](@ref pm_distUnifEll),
!>                          there is straightforward solution for generating uniform points on the surface an arbitrary dimensional hyper-ellipsoid.<br>
!>                          The exact equations and numerical integrations for computing the surface curvature of a hyper-ellipsoid are,
!>                          for example, discussed [here](https://analyticphysics.com/Higher%20Dimensions/Ellipsoids%20in%20Higher%20Dimensions.htm).<br>
!>              </ol>
!>  </ol>
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [pm_distUnifSphere](@ref pm_distUnifSphere)<br>
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  Marsaglia G., et al. (1972), Choosing a point from the surface of a sphere. The Annals of Mathematical Statistics 43(2):645-646.<br>
!>
!>  \test
!>  [test_pm_distUnifSphere](@ref test_pm_distUnifSphere)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distUnifSphere

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distUnifSphere"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type MultiVariate Uniform Spherical
    !>  as defined in the description of [pm_distUnifSphere](@ref pm_distUnifSphere).
    !>
    !>  \details
    !>  See the documentation of [pm_distUnifSphere](@ref pm_distUnifSphere) for the definition of the MultiVariate Uniform Spherical distribution.
    !>
    !>  \interface{distUnifSphere_type}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifSphere, only: distUnifSphere_type
    !>      type(distUnifSphere_type) :: distUnifSphere
    !>
    !>      distUnifSphere = distUnifSphere_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distUnifSphere](@ref test_pm_distUnifSphere)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distUnifSphere_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distUnifSphere_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF)
    !>  of the uniform distribution on an \f$n\f$-sphere embedded in an \f$\ndim\f$-dimensional space.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifSphere](@ref pm_distUnifSphere) for details of the definition of the PDF.<br>
    !>
    !>  \note
    !>  Note that the procedures of this generic interface do not require an
    !>  input `X` value representing a location within the domain of the density function.<br>
    !>  This is fine and intentional by design because the PDF is uniform across the entire support of the PDF.<br>
    !>
    !>  \param[in]  logRadius   :   The input scalar of the same type and kind as the output `logPDF`,
    !>                              containing the natural logarithm of the radius of the \f$n\f$-sphere embedded in the \f$\ndim\f$-dimensional space.<br>
    !>                              A value of `1` corresponds to unit \f$n\f$-sphere.<br>
    !>  \param[in]  ndim        :   The input positive scalar of type `integer` of default kind \IK,
    !>                              containing the dimension of space within which the \f$n\f$-sphere is **embedded**.<br>
    !>                              In other words, \f$\ndim = n + 1\f$.<br>
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the natural logarithm of the PDF of the uniform distribution on the surface of \f$n\f$-sphere.<br>
    !>
    !>  \interface{getUnifSphereLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifSphere, only: getUnifSphereLogPDF
    !>
    !>      logPDF = getUnifSphereLogPDF(logRadius, ndim
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \example{getUnifSphereLogPDF}
    !>  \include{lineno} example/pm_distUnifSphere/getUnifSphereLogPDF/main.F90
    !>  \compilef{getUnifSphereLogPDF}
    !>  \output{getUnifSphereLogPDF}
    !>  \include{lineno} example/pm_distUnifSphere/getUnifSphereLogPDF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnifSphere](@ref test_pm_distUnifSphere)
    !>
    !>  \final{getUnifSphereLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getUnifSphereLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getUnifSphereLogPDF_D0_RK5(logRadius, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifSphereLogPDF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logRadius
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getUnifSphereLogPDF_D0_RK4(logRadius, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifSphereLogPDF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logRadius
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getUnifSphereLogPDF_D0_RK3(logRadius, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifSphereLogPDF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logRadius
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getUnifSphereLogPDF_D0_RK2(logRadius, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifSphereLogPDF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logRadius
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getUnifSphereLogPDF_D0_RK1(logRadius, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifSphereLogPDF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logRadius
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a (collection) of random vector(s) of size `ndim` uniformly distributed on the surface of an \f$n\f$-sphere,
    !>  optionally with the specified input `mean(1:ndim)` and optionally affine-transformed to a non-uniform distribution on
    !>  the surface of an \f$(n+1)\f$-ellipsoid represented by the Cholesky Factorization of its Gramian matrix.
    !>
    !>  \details
    !>  The procedures of this generic interface are merely wrappers around the subroutine interface [setUnifRand](@ref pm_distUnifSphere::setUnifSphereRand).<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type).)
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the center of the sphere.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`. It must be present if the input argument `chol` is missing.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains the [Cholesky Factorization](@ref pm_matrixChol)
    !>                              of the Gramian matrix of the corresponding hyper-ellipsoid on which the output random vectors must be distributed proportional to the ellipsoidal surface curvature.<br>
    !>                              (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` is also present.)
    !>  \param[in]      subset  :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the input `chol` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the input `chol` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `chol` is present.)
    !>  \param[in]      nsam    :   The input scalar `integer` of default kind \IK containing the number of random vectors to generate.<br>
    !>                              (**optional**. If present, the output `rand` is of rank `2`, otherwise is of rank `1`.)
    !>
    !>  \return
    !>  `rand`                  :   The output vector of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the random output vector(s):<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `nsam` is missing, then `rand` shall be of shape `(1:ndim)`.<br>
    !>                                  <li>    If the input argument `nsam` is present, then `rand` shall be of shape `(1:ndim, 1:nsam)`.<br>
    !>                              </ol>
    !>
    !>  \interface{getUnifSphereRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifSphere, only: getUnifSphereRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      rand(1:ndim) = getUnifSphereRand(mean(1:ndim))
    !>      rand(1:ndim) = getUnifSphereRand(chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getUnifSphereRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      rand(1:ndim) = getUnifSphereRand(rng, mean(1:ndim))
    !>      rand(1:ndim) = getUnifSphereRand(rng, chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getUnifSphereRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(rng, mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(rng, chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifSphereRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifSphereRand](@ref pm_distUnifSphere::getUnifSphereRand)<br>
    !>  [setUnifSphereRand](@ref pm_distUnifSphere::setUnifSphereRand)<br>
    !>  [getUnifEllRand](@ref pm_distUnifEll::getUnifEllRand)<br>
    !>  [setUnifEllRand](@ref pm_distUnifEll::setUnifEllRand)<br>
    !>  [getUnifParRand](@ref pm_distUnifPar::getUnifParRand)<br>
    !>  [setUnifParRand](@ref pm_distUnifPar::setUnifParRand)<br>
    !>
    !>  \example{getUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/getUnifSphereRand/main.F90
    !>  \compilef{getUnifSphereRand}
    !>  \output{getUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/getUnifSphereRand/main.out.F90
    !>  \postproc{getUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/getUnifSphereRand/main.py
    !>  \vis{getUnifSphereRand}
    !>  \image html pm_distUnifSphere/getUnifSphereRand/getUnifSphereRandMean.RK.png width=700
    !>  \image html pm_distUnifSphere/getUnifSphereRand/getUnifSphereRandChol.RK.png width=700
    !>  \image html pm_distUnifSphere/getUnifSphereRand/getUnifSphereRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifSphere](@ref test_pm_distUnifSphere)
    !>
    !>  \final{getUnifSphereRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D1_RK5(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D1_RK4(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D1_RK3(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D1_RK2(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D1_RK1(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGF

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGX

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGD

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D2_RK5(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D2_RK4(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D2_RK3(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D2_RK2(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_DC_XXX_D2_RK1(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_DM_AC_UXD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_DM_AC_XLD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_AC_UXD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGD_AM_AC_XLD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGD_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGF

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGF_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGF_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGX

    interface getUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getUSR_RNGX_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUSR_RNGX_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a (collection) of random vector(s) of size `ndim` uniformly distributed on the surface of an \f$n\f$-sphere,
    !>  optionally with the specified input `mean(1:ndim)` and optionally affine-transformed to a non-uniform distribution on
    !>  the surface of an \f$(n+1)\f$-ellipsoid represented by the Cholesky Factorization of its Gramian matrix.
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type).)
    !>  \param[out]     rand    :   The output `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:nsam)` of<br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the random output vector.<br>
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the center of the sphere.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`. It must be present if the input argument `chol` is missing.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains the [Cholesky Factorization](@ref pm_matrixChol)
    !>                              of the Gramian matrix of the corresponding hyper-ellipsoid on which the output random vectors must be distributed proportional to the ellipsoidal surface curvature.<br>
    !>                              (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` is also present.)
    !>  \param[in]      subset  :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the input `chol` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the input `chol` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `chol` is present.)
    !>
    !>  \interface{setUnifSphereRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifSphere, only: setUnifSphereRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      call setUnifSphereRand(rand(1:ndim), mean(1:ndim))
    !>      call setUnifSphereRand(rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifSphereRand(rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      call setUnifSphereRand(rng, rand(1:ndim), mean(1:ndim))
    !>      call setUnifSphereRand(rng, rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifSphereRand(rng, rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      call setUnifSphereRand(rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setUnifSphereRand(rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifSphereRand(rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      call setUnifSphereRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setUnifSphereRand(rng, rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifSphereRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mean, 1) == size(rand, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(chol) == size(rand, 1))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface are `pure` when the input argument `rng` is set to
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) and the compile-time macro `CHECK_ENABLED` is set to `0` or is undefined.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getUnifSphereRand](@ref pm_distUnifSphere::getUnifSphereRand)<br>
    !>  [setUnifSphereRand](@ref pm_distUnifSphere::setUnifSphereRand)<br>
    !>  [getUnifEllRand](@ref pm_distUnifEll::getUnifEllRand)<br>
    !>  [setUnifEllRand](@ref pm_distUnifEll::setUnifEllRand)<br>
    !>  [getUnifParRand](@ref pm_distUnifPar::getUnifParRand)<br>
    !>  [setUnifParRand](@ref pm_distUnifPar::setUnifParRand)<br>
    !>
    !>  \example{setUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/setUnifSphereRand/main.F90
    !>  \compilef{setUnifSphereRand}
    !>  \output{setUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/setUnifSphereRand/main.out.F90
    !>  \postproc{setUnifSphereRand}
    !>  \include{lineno} example/pm_distUnifSphere/setUnifSphereRand/main.py
    !>  \vis{setUnifSphereRand}
    !>  \image html pm_distUnifSphere/setUnifSphereRand/setUnifSphereRand.RK.png width=700
    !>  \image html pm_distUnifSphere/setUnifSphereRand/setUnifSphereRandMean.RK.png width=700
    !>  \image html pm_distUnifSphere/setUnifSphereRand/setUnifSphereRandChol.RK.png width=700
    !>  \image html pm_distUnifSphere/setUnifSphereRand/setUnifSphereRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifSphere](@ref test_pm_distUnifSphere)
    !>
    !>  \final{setUnifSphereRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D1_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D1_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D1_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D1_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D1_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D1_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D1_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D1_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D1_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D1_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGF

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGX

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGD

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D2_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D2_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D2_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D2_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_DC_XXX_D2_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D2_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D2_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D2_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D2_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_DC_XXX_D2_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_UXD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_DM_AC_XLD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_UXD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGD_AM_AC_XLD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGD_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGF

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGF_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGF_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGX

    interface setUnifSphereRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setUSR_RNGX_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setUSR_RNGX_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distUnifSphere