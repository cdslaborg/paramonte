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
!>  This module contains procedures and routines for the creating test datasets for clustering algorithms.
!>
!>  \details
!>  The [mmvue_type](@ref pm_clusTest::mmvue_type) generates and stores information about a randomly generated
!>  [multiple multivariate uniform ellipsoidal distribution](@ref pm_distUnifElls) in default library precision \RK.<br>
!>  The choice of default library precision \RK is currently enforced by the lack of widespread support for Parameterized Derived Types in Fortran compilers (particularly gfortran).<br>
!>
!>  \see
!>  [pm_clusKmeans](@ref pm_clusKmeans)<br>
!>  [pm_clusPartition](@ref pm_clusPartition)<br>
!>
!>  \test
!>  [test_pm_clusTest](@ref test_pm_clusTest)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 03, 2017, 2:16 PM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_clusTest

    use pm_err, only: err_type
    use pm_kind, only: SK, IK, LK, RKG => RK
    use pm_distUnif, only: rngf, rngf_type, xoshiro256ssw_type, rngx_type => xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_clusTest"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating objects containing the range of specifications of an [MMVUE](@ref pm_distUnifElls)
    !>  distribution as necessary by the derived type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>
    !>
    !>  \details
    !>  See also the documentation of the type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>
    !>
    !>  \interface{range_type}
    !>  \code{.F90}
    !>
    !>      use pm_clusTest, only: range_type
    !>      type(range_type) :: range
    !>
    !>      range = range_type(ndim = ndim(1:2), nell = nell(1:2), nsam = nsam(1:2), mean = mean(1:2), std = std(1:2))
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [mmvue_type](@ref pm_clusTest::mmvue_type)<br>
    !>  [mmvue_typer](@ref pm_clusTest::mmvue_typer)<br>
    !>
    !>  \test
    !>  [test_pm_clusTest](@ref test_pm_clusTest)
    !>
    !>  \final{range_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type :: range_type
        integer(IK) :: ndim(2) = [2, 2]                     !<  \public The vector of size `2` of type `integer` of default kind \IK, containing the positive-integer lower and upper limits of the randomly chosen number of dimensions of the domain of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK) :: nell(2) = [1, 20]                    !<  \public The vector of size `2` of type `integer` of default kind \IK, containing the positive-integer lower and upper limits of the randomly chosen number of ellipsoids in the domain of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK) :: nsam(2) = [5000, 5000]               !<  \public The vector of size `2` of type `integer` of default kind \IK, containing the positive-integer lower and upper limits of the randomly chosen number of random points sampled from the support of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        real(RKG)   :: mean(2) = [-1._RKG, 1._RKG]          !<  \public The vector of size `2` of type `real` of default kind \RK, containing the lower and upper limits of the randomly chosen means of the ellipsoids of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        real(RKG)   :: std(2) = [0.001_RKG, 1._RKG]         !<  \public The vector of size `2` of type `real` of default kind \RK, containing the positive-real lower and upper limits of the randomly chosen means of the ellipsoids of the [MMVUE distribution](@ref pm_distUnifElls).<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating objects containing the specifications of a realization of an
    !>  [MMUVE distribution](@ref pm_distUnifElls) and a collection of randomly sampled points from it.<br>
    !>
    !>  \details
    !>  See also the documentation of the type constructor [mmvue_typer](@ref pm_clusTest::mmvue_typer).<br>
    !>
    !>  \interface{mmvue_type}
    !>  \code{.F90}
    !>
    !>      use pm_clusTest, only: mmvue_type
    !>      type(mmvue_type) :: mmvue
    !>
    !>      mmvue = mmvue_type(rng, range = range_type(), nsim = nsim)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [mmvue_typer](@ref pm_clusTest::mmvue_typer)<br>
    !>
    !>  \example{mmvue_type}
    !>  \include{lineno} example/pm_clusTest/mmvue_type/main.F90
    !>  \compilef{mmvue_type}
    !>  \output{mmvue_type}
    !>  \include{lineno} example/pm_clusTest/mmvue_type/main.out.F90
    !>  \postproc{mmvue_type}
    !>  \include{lineno} example/pm_clusTest/mmvue_type/main.py
    !>  \vis{mmvue_type}
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.0.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.1.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.2.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.3.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.4.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.5.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.6.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.7.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.8.svg width=700
    !>  <br><br>
    !>  \image html pm_clusTest/mmvue_type/mmvue_type.9.svg width=700
    !>
    !>  \test
    !>  [test_pm_clusTest](@ref test_pm_clusTest)
    !>
    !>  \todo
    !>  \phigh
    !>  This type must be generalized to accept `real` components of arbitrary kind, once GNU compiler support for PDTs is strong.<br>
    !>
    !>  \final{mmvue_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    type :: mmvue_type
        type(range_type)                :: range                    !<  \public The scalar object of type [range_type](@ref pm_clusTest::range_type) containing the specifications of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK)                     :: ndim                     !<  \public The number of dimensions of the domain of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK)                     :: nell                     !<  \public The number of ellipsoids in the domain of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK)                     :: nsam                     !<  \public The number of random points sampled from the domain of the [MMVUE distribution](@ref pm_distUnifElls).<br>
        integer(IK)                     :: nsim = 10000             !<  \public The scalar of type `integer` of default kind \IK, representing the number of simulated random vectors from
                                                                    !!          the ellipsoids of the distribution to estimate the effective volumes of the ellipsoids (i.e., the inverse of the PDF of the distribution).<br>
        real(RKG)                       :: logVolUnitBall           !<  \public The natural logarithm of the volume of a unit ball in `ndim` dimensions.<br>
        real(RKG)                       :: logDensity               !<  \public The scalar of type `real` containing the natural logarithm of the inverse of the *effective* volume
                                                                    !!          of a single sampled point from the distribution domain, normalized by the volume of a unit ball.<br>
                                                                    !!          In other words, `logDensity = log(mmvue_type%%nsam) - mmvue_type%%logSumVolNormedEff`.<br>
        real(RKG)                       :: logSumVolNormedEff       !<  \public The natural logarithm of the sum of the effective volumes of the distribution ellipsoids, normalized by `logVolUnitBall`.<br>
                                                                    !!          This measure includes shared volumes between ellipsoids only once from one individual ellipsoid.<br>
                                                                    !!          Because this measure does not have closed expression, it must be approximated using Monte Carlo approach.<br>
                                                                    !!          By definition, this sum is always smaller than or equal to the sum of the individual normalized ellipsoid volumes.<br>
                                                                    !!          The equality happens only when no two ellipsoids overlap.<br>
        real(RKG)       , allocatable   :: std(:,:)                 !<  \public The matrix of shape `(1:ndim, 1:nell)` containing the standard deviations of the ellipsoids in the distribution.<br>
        real(RKG)       , allocatable   :: sample(:,:)              !<  \public The matrix of shape `(1:ndim, 1:nsam)` containing the set of uniformly-sampled random points from the distribution domain.<br>
        real(RKG)       , allocatable   :: mean(:,:)                !<  \public The matrix of shape `(1:ndim, 1:nell)` containing the set of `nell` ellipsoid centers of the distribution domain.<br>
        real(RKG)       , allocatable   :: invGram(:,:,:)           !<  \public The matrix of shape `(1:ndim, 1:ndim, 1:nell)` containing the set of `nell` ellipsoid inverse-Gramian matrices of the distribution domain.<br>
        real(RKG)       , allocatable   :: choLowGramUpp(:,:,:)     !<  \public The matrix of shape `(1:ndim, 1:ndim + 1, 1:nell)` containing the set of `nell` ellipsoid Cholesky-Gramian half (triangular) matrices of the distribution domain.<br>
        real(RKG)       , allocatable   :: cumPropVolNormed(:)      !<  \public The vector of shape `(0:nell)` containing the cumulative proportion of volumes of individual ellipsoids in the distribution, used for generating the random `sample`.<br>
        real(RKG)       , allocatable   :: logVolNormed(:)          !<  \public The vector of shape `(1:nell)` containing the natural logarithm of the set of `nell` volumes of the Gramian matrices of the distribution domain, normalized by the volume of a unit ball.<br>
        real(RKG)       , allocatable   :: mahalSq(:,:)             !<  \public The matrix of shape `(1:nell, 1:nsam)` containing the square of the Mahalanobis distances the randomly sampled points from each ellipsoid center in the distribution.<br>
        real(RKG)       , allocatable   :: invmul(:)                !<  \public The vector of shape `(1:nsam)` containing the inverse of sample multiplicity, that is, the ellipsoidal membership count of each randomly-sampled point from the distribution.<br>
                                                                    !!          For example, a multiplicity of `3` means that the `isam` sampled point falls within exactly three ellipsoids in the distribution.<br>
        integer(IK)     , allocatable   :: size(:)                  !<  \public The vector of shape `(1:nell)` containing the set of `nell` counts of randomly-sampled points from within each ellipsoid in the distribution.<br>
        integer(IK)     , allocatable   :: membership(:)            !<  \public The vector of shape `(1:nsam)` containing the ellipsoidal membership of each randomly-sampled point from the distribution.<br>
        type(err_type)                  :: err
    contains
        procedure, pass :: write => mmvue_type_write
    end type

    !>  \cond excluded
    interface mmvue_type
        module procedure :: mmvue_typer_rngf, mmvue_typer_rngx
    end interface
    !>  \endcond

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an object of type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>
    !>
    !>  \brief
    !>  This procedure is the constructor of the derived type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>.
    !>  See the documentation of [mmvue_type](@ref pm_clusTest::mmvue_type) for example usage.<br>
    !>
    !>  \param[inout]   rng         :   The input/output scalar that can be an object of,
    !>                                  <ol>
    !>                                      <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                              implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                      <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                              implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                                  </ol>
    !>  \param[in]      range       :   The input scalar object of type [range_type](@ref pm_clusTest::range_type) containing the optionally adjustable settings of the [MMVUE distribution](@ref pm_distUnifElls).<br>
    !>                                  containing the lower and upper limits for the random number of dimensions of the domain of the target [MMUVE distribution](@ref pm_distUnifElls).<br>
    !>                                  (**optional**, default = [range_type()](@ref pm_clusTest::range_type).)
    !>  \param[in]      nsim        :   The input scalar of type `integer` of default kind \IK,
    !>                                  representing the number of simulated random vectors from the [MMUVE distribution](@ref pm_distUnifElls) for computing the effective volume of the distribution ellipsoids.<br>
    !>                                  The larger the value of `nsim`, the more accurate the volume of the potentially-overlapping multi-ellipsoidal support of the distribution will be computed.<br>
    !>                                  (**optional**. The default value is set by the corresponding `nsim` component of [mmvue_type](@ref pm_clusTest::mmvue_type).)
    !>
    !>  \return
    !>  `self`                      :   The output scalar object of type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>
    !>
    !>  \interface{mmvue_typer}
    !>  \code{.F90}
    !>
    !>      use pm_clusTest, only: mmvue_type
    !>      type(mmvue_type) :: mmvue
    !>
    !>      mmvue = mmvue_type(rng, range = range_type(), nsim = nsim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < range%ndim(1) .and. range%ndim(1) <= range%ndim(2)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < range%nell(1) .and. range%nell(1) <= range%nell(2)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < range%nsam(1) .and. range%nsam(1) <= range%nsam(2)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < range%std(1) .and. range%std(1) <= range%std(2)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < nsim` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \see
    !>  [mmvue_type](@ref pm_clusTest::mmvue_type)<br>
    !>
    !>  \test
    !>  [test_pm_clusTest](@ref test_pm_clusTest)
    !>
    !>  \final{mmvue_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface mmvue_typer

    impure module function mmvue_typer_rngf(rng, range, nsim) result(self)
#if INTEL_COMPILER_ENABLED && DLL_ENABLED && (WINDOWS_ENABLED || DARWIN_ENABLED)
        !DEC$ ATTRIBUTES DLLEXPORT :: mmvue_typer_rngf
#endif
        type(range_type), intent(in), optional :: range
        integer(IK), intent(in), optional :: nsim
        type(rngf_type), intent(in) :: rng
        type(mmvue_type) :: self
    end function

    impure module function mmvue_typer_rngx(rng, range, nsim) result(self)
#if INTEL_COMPILER_ENABLED && DLL_ENABLED && (WINDOWS_ENABLED || DARWIN_ENABLED)
        !DEC$ ATTRIBUTES DLLEXPORT :: mmvue_typer_rngx
#endif
        type(range_type), intent(in), optional :: range
        integer(IK), intent(in), optional :: nsim
        type(rngx_type), intent(inout) :: rng
        type(mmvue_type) :: self
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Write the specifications of an object of type [mmvue_type](@ref pm_clusTest::mmvue_type) to an
    !>  external `sequential` formatted filed with the specified input file unit `funit`.<br>
    !>
    !>  \brief
    !>  This procedure is a dynamic method of the derived type [mmvue_type](@ref pm_clusTest::mmvue_type).<br>.
    !>  See the documentation of [mmvue_type](@ref pm_clusTest::mmvue_type) for example usage.<br>
    !>
    !>  \param[in]  self    :   The input scalar object of type [mmvue_type](@ref pm_clusTest::mmvue_type) that is **implicitly** passed to the procedure.<br>
    !>  \param[in]  funit   :   The input scalar `integer` of default kind \IK, containing the unit number associated with the external formatted `sequential`
    !>                          file to which the specifications of the object of type [mmvue_type](@ref pm_clusTest::mmvue_type) must be appended.<br>
    !>
    !>  \interface{mmvue_type_write}
    !>  \code{.F90}
    !>
    !>      use pm_clusTest, only: mmvue_type
    !>      type(mmvue_type) :: mmvue
    !>
    !>      call mmvue%write(funit)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [mmvue_type](@ref pm_clusTest::mmvue_type)<br>
    !>  [mmvue_type_write](@ref pm_clusTest::mmvue_type_write)<br>
    !>
    !>  \test
    !>  [test_pm_clusTest](@ref test_pm_clusTest)
    !>
    !>  \final{mmvue_type_write}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface mmvue_type_write

    impure module subroutine mmvue_type_write(self, funit)
#if INTEL_COMPILER_ENABLED && DLL_ENABLED && (WINDOWS_ENABLED || DARWIN_ENABLED)
        !DEC$ ATTRIBUTES DLLEXPORT :: mmvue_type_write
#endif
        class(mmvue_type), intent(in) :: self
        integer(IK), intent(in) :: funit
    end subroutine

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_clusTest
