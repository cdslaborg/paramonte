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
!>  This module contains the derived types for generating `allocatable` containers of scalar,
!>  vector, matrix, or cube of integer, real, complex, logical, and string values of arbitrary kinds.<br>
!>
!>  \details
!>  This module offer two broad categories of container derived types:<br>
!>  <ol>
!>      <li>    containers that end with \f$\ms{*_type}\f$ are **derived type** containers of *intrinsic types* of **default kind**: \SK, \IK, \LK, \CK, \RK.<br>
!>      <li>    containers that end with \f$\ms{*_pdt}\f$ are **parameterized derived** type containers of *intrinsic types* of **arbitrary kind**.<br>
!>  </ol>
!>
!>  Each of the above categories contains the following subcategories based on their component (`val`) rank:<br>
!>  <ol>
!>      <li>    containers that end with \f$\ms{*_type}\f$ are **derived type** containers of *intrinsic types* of **default kind**: \SK, \IK, \LK, \CK, \RK.<br>
!>      <li>    containers that end with \f$\ms{*_pdt}\f$ are **parameterized derived** type containers of *intrinsic types* of **arbitrary kind type parameter**.<br>
!>  </ol>
!>
!>  Naming convention
!>  -----------------
!>
!>  The containers of this module use a concise letter-based naming convention:<br>
!>  <ol>
!>      <li>    The **first letter** in all container names \f$\ms{c**_type}\f$ or \f$\ms{c**_pdt}\f$ stands for <b>c</b>ontainer.<br>
!>      <li>    The **second letter** in all container names \f$\ms{*x*_type}\f$ or \f$\ms{*x*_pdt}\f$ stands for the **rank** of the component of the container:<br>
!>              <ol>
!>                  <li>    The letter \f$\ms{s}\f$ as in \f$\ms{*s*_type}\f$ or \f$\ms{*s*_pdt}\f$ stands for a <b>s</b>calar-component container (`container%%val`).<br>
!>                  <li>    The letter \f$\ms{v}\f$ as in \f$\ms{*v*_type}\f$ or \f$\ms{*v*_pdt}\f$ stands for a <b>v</b>ector-component container (`container%%val(:)`).<br>
!>                  <li>    The letter \f$\ms{m}\f$ as in \f$\ms{*m*_type}\f$ or \f$\ms{*m*_pdt}\f$ stands for a <b>m</b>atrix-component container (`container%%val(:,:)`).<br>
!>                  <li>    The letter \f$\ms{c}\f$ as in \f$\ms{*c*_type}\f$ or \f$\ms{*c*_pdt}\f$ stands for a <b>c</b>ube-component container (`container%%val(:,:,:)`).<br>
!>              </ol>
!>      <li>    The **third letter** in all container names \f$\ms{**x_type}\f$ or \f$\ms{**x_pdt}\f$ stands for the **type** of the component of the container:<br>
!>              <ol>
!>                  <li>    The letter \f$\ms{s}\f$ as in \f$\ms{**s_type}\f$ or \f$\ms{**s_pdt}\f$ stands for a <b>s</b>tring-component (character) container.<br>
!>                  <li>    The letter \f$\ms{i}\f$ as in \f$\ms{**i_type}\f$ or \f$\ms{**i_pdt}\f$ stands for a <b>i</b>nteger-component container.<br>
!>                  <li>    The letter \f$\ms{l}\f$ as in \f$\ms{**l_type}\f$ or \f$\ms{**l_pdt}\f$ stands for a <b>l</b>ogical-component container.<br>
!>                  <li>    The letter \f$\ms{c}\f$ as in \f$\ms{**c_type}\f$ or \f$\ms{**c_pdt}\f$ stands for a <b>c</b>omplex-component container.<br>
!>                  <li>    The letter \f$\ms{r}\f$ as in \f$\ms{**r_type}\f$ or \f$\ms{**r_pdt}\f$ stands for a <b>r</b>eal-component container.<br>
!>              </ol>
!>      <li>    The **final suffix** \f$\ms{***_type}\f$ and \f$\ms{***_pdt}\f$ stand for the *derived type* and *Parameterized Derived Type* (PDT) containers.<br>
!>  </ol>
!>
!>  **For example:**
!>  <ol>
!>      <li>    Read [css_type](@ref pm_container::css_type) as <b>c</b>ontainer of <b>s</b>calar of <b>string</b> (`character`) derived type (of default kind \SK).
!>              This container is widely known as **varying string** or simply **string** in the Fortran community.<br>
!>      <li>    Read [css_pdt](@ref pm_container::css_pdt) as <b>c</b>ontainer of <b>s</b>calar of <b>string</b> (`character`) derived type (of kind \SKALL).
!>              This container is widely known as **varying string** or simply **string** in the Fortran community.<br>
!>      <li>    Read [csi_type](@ref pm_container::csi_type) as <b>c</b>ontainer of <b>s</b>calar of <b>i</b>nteger derived type (of default kind \IK).
!>      <li>    Read [cvr_type](@ref pm_container::cvr_type) as <b>c</b>ontainer of <b>v</b>ector of <b>r</b>eal derived type (of default kind \RK).
!>  </ol>
!>
!>  What is the use of containers?
!>  ------------------------------
!>
!>  Containers are essential for creating [jagged arrays](https://en.wikipedia.org/wiki/Jagged_array).<br>
!>  While string containers are the most popular container type in almost all languages,
!>  other container types and ranks offered by this module find heavy use in Machine Learning and Data Science or even simple computational tasks.<br>
!>
!>  \warning
!>  The `elemental` versions of the generic interfaces of this procedure are known to yield incorrect result if compiled by gfortran version 10-12.<br>
!>  The Intel Fortran compiler yields the correct `elemental` results.<br>
!>
!>  \devnote
!>  Currently, only the scalar containers (matching names \f$\ms{cs*_pdt}\f$ and \f$\ms{cs*_type}\f$) have custom `elemental` constructors.<br>
!>  The was no perceived usage and consensus on how to design custom constructors for containers
!>  of higher-rank components (e.g., vector, matrix, and cube containers).<br>
!>
!>  \see
!>  [pm_kind](@ref pm_kind)<br>
!>
!>  \todo
!>  \pmed
!>  Currently, only the following containers are exemplified:<br>
!>  <ol>
!>      <li>    [css_pdt](@ref pm_container::css_pdt)<br>
!>      <li>    [css_type](@ref pm_container::css_type)<br>
!>  </ol>
!>  Examples for the following generic interfaces must be added:<br>
!>  <ol>
!>      <li>    [csi_pdt](@ref pm_container::csi_pdt)<br>
!>      <li>    [csl_pdt](@ref pm_container::csl_pdt)<br>
!>      <li>    [csc_pdt](@ref pm_container::csc_pdt)<br>
!>      <li>    [csr_pdt](@ref pm_container::csr_pdt)<br>
!>      <li>    [csi_type](@ref pm_container::csi_type)<br>
!>      <li>    [csl_type](@ref pm_container::csl_type)<br>
!>      <li>    [csc_type](@ref pm_container::csc_type)<br>
!>      <li>    [csr_type](@ref pm_container::csr_type)<br>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_container

    use pm_kind, only: SK, IK, LK, CK, RK

    implicit none

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [css_type](@ref pm_container::css_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>s</b>tring objects.<br>
    !>
    !>  \details
    !>  The [css_type](@ref pm_container::css_type) container is widely known as **varying string** or **string container** within the Fortran language.<br>
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `character` of default kind \SK.<br>
    !>  \param[in]  trimmed :   The input scalar or array of the same shape as the input `val`, of type `logical` of default kind \LK.<br>
    !>                          If `.true.`, the corresponding elements of the input `val` will **not** be trimmed.<br>
    !>                          Setting this argument to `.true.` will lead to better runtime performance.<br>
    !>                          (**optional**, default = `.false.`, i.e., the trailing blanks in all input `val` will be trimmed).
    !>
    !>  \return
    !>  `container`         :   The output object of type [css_type](@ref pm_container::css_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{css_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_container, only: css_type
    !>      type(css_type), allocatable  :: containers(:)
    !>      type(css_type) :: container
    !>
    !>      container = css_type(val, trimmed = trimmed)
    !>
    !>      ! example
    !>
    !>      container = css_type(SK_"aaa")
    !>      containers = css_type([character(3, SK) :: "a", "aa", "aaa"])
    !>      containers = css_type([character(3, SK) :: "a", "aa", "aaa"], trimmed = .false._LK)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [css_pdt](@ref pm_container::css_pdt)<br>
    !>  [cvs_pdt](@ref pm_container::cvs_pdt)<br>
    !>  [cms_pdt](@ref pm_container::cms_pdt)<br>
    !>  [ccs_pdt](@ref pm_container::ccs_pdt)<br>
    !>  [css_type](@ref pm_container::css_type)<br>
    !>  [cvs_type](@ref pm_container::cvs_type)<br>
    !>  [cms_type](@ref pm_container::cms_type)<br>
    !>  [ccs_type](@ref pm_container::ccs_type)<br>
    !>
    !>  \example{css_type}
    !>  \include{lineno} example/pm_container/css_type/main.F90
    !>  \compilef{css_type}
    !>  \output{css_type}
    !>  \include{lineno} example/pm_container/css_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{css_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: css_type
        character(:, SK), allocatable :: val
    end type

    !>  \cond excluded
    interface css_type
    pure elemental module function css_typer_D0(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: css_typer_D0
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csi_type](@ref pm_container::csi_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>i</b>nteger objects.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `integer` of default kind \LK.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csi_type](@ref pm_container::csi_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csi_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_container, only: csi_type
    !>      type(csi_type) :: container
    !>
    !>      container(..) = csi_type(val(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [csi_pdt](@ref pm_container::csi_pdt)<br>
    !>  [cvi_pdt](@ref pm_container::cvi_pdt)<br>
    !>  [cmi_pdt](@ref pm_container::cmi_pdt)<br>
    !>  [cci_pdt](@ref pm_container::cci_pdt)<br>
    !>  [csi_type](@ref pm_container::csi_type)<br>
    !>  [cvi_type](@ref pm_container::cvi_type)<br>
    !>  [cmi_type](@ref pm_container::cmi_type)<br>
    !>  [cci_type](@ref pm_container::cci_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csi_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csi_type
        integer(IK), allocatable :: val
    end type

    !>  \cond excluded
    interface csi_type
    pure elemental module function csi_typer_D0(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: csi_typer_D0
#endif
        use pm_kind, only: IKG => IK
        integer(IKG)        , intent(in)                :: val
        type(csi_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csl_type](@ref pm_container::csl_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>l</b>ogical objects.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `logical` of default kind \LK.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csl_type](@ref pm_container::csl_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csl_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_container, only: csl_type
    !>      type(csl_type) :: container
    !>
    !>      container(..) = csl_type(val(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [csl_pdt](@ref pm_container::csl_pdt)<br>
    !>  [cvl_pdt](@ref pm_container::cvl_pdt)<br>
    !>  [cml_pdt](@ref pm_container::cml_pdt)<br>
    !>  [ccl_pdt](@ref pm_container::ccl_pdt)<br>
    !>  [csl_type](@ref pm_container::csl_type)<br>
    !>  [cvl_type](@ref pm_container::cvl_type)<br>
    !>  [cml_type](@ref pm_container::cml_type)<br>
    !>  [ccl_type](@ref pm_container::ccl_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csl_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csl_type
        logical(LK), allocatable :: val
    end type

    !>  \cond excluded
    interface csl_type
    pure elemental module function csl_typer_D0(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: csl_typer_D0
#endif
        use pm_kind, only: LKG => LK
        logical(LKG)        , intent(in)                :: val
        type(csl_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csc_type](@ref pm_container::csc_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>l</b>ogical objects.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `complex` of default kind \CK.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csc_type](@ref pm_container::csc_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csc_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: CK
    !>      use pm_container, only: csc_type
    !>      type(csc_type) :: container
    !>
    !>      container(..) = csc_type(val(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [csc_pdt](@ref pm_container::csc_pdt)<br>
    !>  [cvc_pdt](@ref pm_container::cvc_pdt)<br>
    !>  [cmc_pdt](@ref pm_container::cmc_pdt)<br>
    !>  [ccc_pdt](@ref pm_container::ccc_pdt)<br>
    !>  [csc_type](@ref pm_container::csc_type)<br>
    !>  [cvc_type](@ref pm_container::cvc_type)<br>
    !>  [cmc_type](@ref pm_container::cmc_type)<br>
    !>  [ccc_type](@ref pm_container::ccc_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csc_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csc_type
        complex(CK), allocatable :: val
    end type

    !>  \cond excluded
    interface csc_type
    pure elemental module function csc_typer_D0(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: csc_typer_D0
#endif
        use pm_kind, only: CKG => CK
        complex(CKG)        , intent(in)                :: val
        type(csc_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csr_type](@ref pm_container::csr_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>l</b>ogical objects.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `real` of default kind \RK.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csr_type](@ref pm_container::csr_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csr_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: RK
    !>      use pm_container, only: csr_type
    !>      type(csr_type) :: container
    !>
    !>      container(..) = csr_type(val(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [csr_pdt](@ref pm_container::csr_pdt)<br>
    !>  [cvr_pdt](@ref pm_container::cvr_pdt)<br>
    !>  [cmr_pdt](@ref pm_container::cmr_pdt)<br>
    !>  [ccr_pdt](@ref pm_container::ccr_pdt)<br>
    !>  [csr_type](@ref pm_container::csr_type)<br>
    !>  [cvr_type](@ref pm_container::cvr_type)<br>
    !>  [cmr_type](@ref pm_container::cmr_type)<br>
    !>  [ccr_type](@ref pm_container::ccr_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csr_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csr_type
        real(RK), allocatable :: val
    end type

    !>  \cond excluded
    interface csr_type
    pure elemental module function csr_typer_D0(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: csr_typer_D0
#endif
        use pm_kind, only: RKG => RK
        real(RKG)           , intent(in)                :: val
        type(csr_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csp_type](@ref pm_container::csp_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar unlimited <b>p</b>olymorphic objects.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of any type of any kind.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csp_type](@ref pm_container::csp_type) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csp_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: RK
    !>      use pm_container, only: csp_type
    !>      type(csp_type) :: container
    !>
    !>      container(..) = csp_type(val(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    !>  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    !>  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    !>
    !>  \see
    !>  [csp_type](@ref pm_container::csp_type)<br>
    !>  [cvp_type](@ref pm_container::cvp_type)<br>
    !>  [cmp_type](@ref pm_container::cmp_type)<br>
    !>  [ccp_type](@ref pm_container::ccp_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csp_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csp_type
        class(*), allocatable :: val
    end type

    !>  \cond excluded
    interface csp_type
    pure elemental module function csp_typer_D0(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: csp_typer_D0
#endif
        use pm_kind, only: RKG => RK
        class(*)            , intent(in)                :: val
        type(csp_type)                                  :: container
    end function
    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type `character` of default kind \SK of arbitrary length type parameter.<br>
    !>
    !>  \copydetails css_type
    type :: cvs_type
        character(:, SK)    , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type `integer` of default kind \IK.<br>
    !>
    !>  \copydetails css_type
    type :: cvi_type
        integer(IK)         , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type `logical` of default kind \LK.<br>
    !>
    !>  \copydetails css_type
    type :: cvl_type
        logical(LK)         , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type `complex` of default kind \CK.<br>
    !>
    !>  \copydetails css_type
    type :: cvc_type
        complex(CK)         , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type `real` of default kind \RK.<br>
    !>
    !>  \copydetails css_type
    type :: cvr_type
        real(RK)            , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a vector component of type unlimited polymorphic.<br>
    !>
    !>  \final{cvp_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    type :: cvp_type
        class(*)            , allocatable   :: val(:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type `character` of default kind \SK of arbitrary length type parameter.<br>
    !>
    !>  \copydetails css_type
    type :: cms_type
        character(:, SK)    , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type `integer` of default kind \IK.<br>
    !>
    !>  \copydetails css_type
    type :: cmi_type
        integer(IK)         , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type `logical` of default kind \LK.<br>
    !>
    !>  \copydetails css_type
    type :: cml_type
        logical(LK)         , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type `complex` of default kind \CK.<br>
    !>
    !>  \copydetails css_type
    type :: cmc_type
        complex(CK)         , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type `real` of default kind \RK.<br>
    !>
    !>  \copydetails css_type
    type :: cmr_type
        real(RK)            , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a matrix component of type unlimited polymorphic.<br>
    !>
    !>  \final{cmp_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    type :: cmp_type
        class(*)            , allocatable   :: val(:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type `character` of default kind \SK of arbitrary length type parameter.<br>
    !>
    !>  \copydetails css_type
    type :: ccs_type
        character(:, SK)    , allocatable   :: val(:,:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type `integer` of default kind \IK.<br>
    !>
    !>  \copydetails css_type
    type :: cci_type
        integer(IK)         , allocatable   :: val(:,:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type `logical` of default kind \LK.<br>
    !>
    !>  \copydetails css_type
    type :: ccl_type
        logical(LK)         , allocatable   :: val(:,:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type `complex` of default kind \CK.<br>
    !>
    !>  \copydetails css_type
    type :: ccc_type
        complex(CK)         , allocatable   :: val(:,:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type `real` of default kind \RK.<br>
    !>
    !>  \copydetails css_type
    type :: ccr_type
        real(RK)            , allocatable   :: val(:,:,:)
    end type

    !>  \brief
    !>  This is the **derived type** for generating a container of a cube component of type unlimited polymorphic.<br>
    !>
    !>  \final{ccp_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    type :: ccp_type
        class(*)            , allocatable   :: val(:,:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

    !>  \brief
    !>  This is the [css_pdt](@ref pm_container::css_pdt) parameterized type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>s</b>tring objects of kind \SKALL.<br>
    !>
    !>  \details
    !>  The [css_pdt](@ref pm_container::css_pdt) container is widely known as **varying string** or **string container** within the Fortran language.<br>
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>  \param[in]  trimmed :   The input scalar or array of the same shape as the input `val`, of type `logical` of default kind \LK.<br>
    !>                          If `.true.`, the corresponding elements of the input `val` will **not** be trimmed.<br>
    !>                          Setting this argument to `.true.` will lead to better runtime performance.<br>
    !>                          (**optional**, default = `.false.`, i.e., all input `val` will be left-adjusted and trimmed).
    !>
    !>  \return
    !>  `container`         :   The output object of type [css_pdt(kind(val)](@ref pm_container::css_pdt) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{css_pdt}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: css_pdt
    !>      type(css_pdt(kind("a"))), allocatable  :: containers(:)
    !>      type(css_pdt(kind("a"))) :: container
    !>
    !>      container = css_pdt(val, trimmed = trimmed)
    !>
    !>      ! example
    !>
    !>      container = css_pdt("aaa")
    !>      containers = css_pdt([character(3) :: "a", "aa", "aaa"])
    !>      containers = css_pdt([character(3) :: "a", "aa", "aaa"], trimmed = .false._LK)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [css_pdt](@ref pm_container::css_pdt)<br>
    !>  [cvs_pdt](@ref pm_container::cvs_pdt)<br>
    !>  [cms_pdt](@ref pm_container::cms_pdt)<br>
    !>  [ccs_pdt](@ref pm_container::ccs_pdt)<br>
    !>  [css_type](@ref pm_container::css_type)<br>
    !>  [cvs_type](@ref pm_container::cvs_type)<br>
    !>  [cms_type](@ref pm_container::cms_type)<br>
    !>  [ccs_type](@ref pm_container::ccs_type)<br>
    !>
    !>  \example{css_pdt}
    !>  \include{lineno} example/pm_container/css_pdt/main.F90
    !>  \compilef{css_pdt}
    !>  \output{css_pdt}
    !>  \include{lineno} example/pm_container/css_pdt/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{css_pdt}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: css_pdt(kind)
        integer             , kind          :: kind = SK
        character(:, kind)  , allocatable   :: val
    end type

    !>  \cond excluded
    interface css_pdt

#if SK5_ENABLED
    pure elemental module function constructCon_D0_PSSK5(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_pdt(SKG))                              :: container
    end function
#endif

#if SK4_ENABLED
    pure elemental module function constructCon_D0_PSSK4(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_pdt(SKG))                              :: container
    end function
#endif

#if SK3_ENABLED
    pure elemental module function constructCon_D0_PSSK3(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_pdt(SKG))                              :: container
    end function
#endif

#if SK2_ENABLED
    pure elemental module function constructCon_D0_PSSK2(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_pdt(SKG))                              :: container
    end function
#endif

#if SK1_ENABLED
    pure elemental module function constructCon_D0_PSSK1(val, trimmed) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: val
        logical(LK)         , intent(in)    , optional  :: trimmed
        type(css_pdt(SKG))                              :: container
    end function
#endif

    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a vector component of type `character` of arbitrary kind of `allocatable` length.
    type :: cvs_pdt(kind)
        integer             , kind          :: kind = SK
        character(:, kind)  , allocatable   :: val(:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a matrix component of type `character` of arbitrary kind of `allocatable` length.
    type :: cms_pdt(kind)
        integer             , kind          :: kind = SK
        character(:, kind)  , allocatable   :: val(:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a cube component of type `character` of arbitrary kind of `allocatable` length.
    type :: ccs_pdt(kind)
        integer             , kind          :: kind = SK
        character(:, kind)  , allocatable   :: val(:,:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csi_pdt](@ref pm_container::csi_pdt) parameterized type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>i</b>nteger objects of kind \IKALL.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `integer` of kind \IKALL.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csi_pdt(kind(val))](@ref pm_container::csi_pdt) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csi_pdt}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: csi_pdt
    !>      type(csi_pdt(kind(val))) :: container
    !>
    !>      container = csi_pdt(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [csi_pdt](@ref pm_container::csi_pdt)<br>
    !>  [cvi_pdt](@ref pm_container::cvi_pdt)<br>
    !>  [cmi_pdt](@ref pm_container::cmi_pdt)<br>
    !>  [cci_pdt](@ref pm_container::cci_pdt)<br>
    !>  [csi_type](@ref pm_container::csi_type)<br>
    !>  [cvi_type](@ref pm_container::cvi_type)<br>
    !>  [cmi_type](@ref pm_container::cmi_type)<br>
    !>  [cci_type](@ref pm_container::cci_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csi_pdt}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csi_pdt(kind)
        integer             , kind          :: kind = IK
        integer(kind)       , allocatable   :: val
    end type

    !>  \cond excluded
    interface csi_pdt

#if IK5_ENABLED
    pure elemental module function constructCon_D0_PSIK5(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)        , intent(in)                :: val
        type(csi_pdt(IKG))                              :: container
    end function
#endif

#if IK4_ENABLED
    pure elemental module function constructCon_D0_PSIK4(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)        , intent(in)                :: val
        type(csi_pdt(IKG))                              :: container
    end function
#endif

#if IK3_ENABLED
    pure elemental module function constructCon_D0_PSIK3(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)        , intent(in)                :: val
        type(csi_pdt(IKG))                              :: container
    end function
#endif

#if IK2_ENABLED
    pure elemental module function constructCon_D0_PSIK2(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)        , intent(in)                :: val
        type(csi_pdt(IKG))                              :: container
    end function
#endif

#if IK1_ENABLED
    pure elemental module function constructCon_D0_PSIK1(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)        , intent(in)                :: val
        type(csi_pdt(IKG))                              :: container
    end function
#endif

    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a vector component of type `integer` of arbitrary kind.
    type :: cvi_pdt(kind)
        integer             , kind          :: kind = IK
        integer(kind)       , allocatable   :: val(:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a matrix component of type `integer` of arbitrary kind.
    type :: cmi_pdt(kind)
        integer             , kind          :: kind = IK
        integer(kind)       , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a cube component of type `integer` of arbitrary kind.
    type :: cci_pdt(kind)
        integer             , kind          :: kind = IK
        integer(kind)       , allocatable   :: val(:,:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csl_pdt](@ref pm_container::csl_pdt) parameterized type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>l</b>ogical objects of kind \LKALL.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `integer` of kind \LKALL.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csl_pdt(kind(val))](@ref pm_container::csl_pdt) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csl_pdt}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: csl_pdt
    !>      type(csl_pdt(kind(val))) :: container
    !>
    !>      container = csl_pdt(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [csl_pdt](@ref pm_container::csl_pdt)<br>
    !>  [cvl_pdt](@ref pm_container::cvl_pdt)<br>
    !>  [cml_pdt](@ref pm_container::cml_pdt)<br>
    !>  [ccl_pdt](@ref pm_container::ccl_pdt)<br>
    !>  [csl_type](@ref pm_container::csl_type)<br>
    !>  [cvl_type](@ref pm_container::cvl_type)<br>
    !>  [cml_type](@ref pm_container::cml_type)<br>
    !>  [ccl_type](@ref pm_container::ccl_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csl_pdt}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csl_pdt(kind)
        integer             , kind          :: kind = LK
        integer(kind)       , allocatable   :: val
    end type

    !>  \cond excluded
    interface csl_pdt

#if LK5_ENABLED
    pure elemental module function constructCon_D0_PSLK5(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        integer(LKG)        , intent(in)                :: val
        type(csl_pdt(LKG))                              :: container
    end function
#endif

#if LK4_ENABLED
    pure elemental module function constructCon_D0_PSLK4(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        integer(LKG)        , intent(in)                :: val
        type(csl_pdt(LKG))                              :: container
    end function
#endif

#if LK3_ENABLED
    pure elemental module function constructCon_D0_PSLK3(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        integer(LKG)        , intent(in)                :: val
        type(csl_pdt(LKG))                              :: container
    end function
#endif

#if LK2_ENABLED
    pure elemental module function constructCon_D0_PSLK2(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        integer(LKG)        , intent(in)                :: val
        type(csl_pdt(LKG))                              :: container
    end function
#endif

#if LK1_ENABLED
    pure elemental module function constructCon_D0_PSLK1(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        integer(LKG)        , intent(in)                :: val
        type(csl_pdt(LKG))                              :: container
    end function
#endif

    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of an `allocatable` vector component of type `logical` of arbitrary kind.
    type :: cvl_pdt(kind)
        integer             , kind          :: kind = LK
        logical(kind)       , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a matrix component of type `logical` of arbitrary kind.
    type :: cml_pdt(kind)
        integer             , kind          :: kind = LK
        logical(kind)       , allocatable   :: val(:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a cube component of type `logical` of arbitrary kind.
    type :: ccl_pdt(kind)
        integer             , kind          :: kind = LK
        logical(kind)       , allocatable   :: val(:,:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csc_pdt](@ref pm_container::csc_pdt) parameterized type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>c</b>omplex objects of kind \CKALL.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `complex` of kind \CKALL.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csc_pdt(kind(val))](@ref pm_container::csc_pdt) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csc_pdt}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: csc_pdt
    !>      type(csc_pdt(kind(val))) :: container
    !>
    !>      container = csc_pdt(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [csc_pdt](@ref pm_container::csc_pdt)<br>
    !>  [cvc_pdt](@ref pm_container::cvc_pdt)<br>
    !>  [cmc_pdt](@ref pm_container::cmc_pdt)<br>
    !>  [ccc_pdt](@ref pm_container::ccc_pdt)<br>
    !>  [csc_type](@ref pm_container::csc_type)<br>
    !>  [cvc_type](@ref pm_container::cvc_type)<br>
    !>  [cmc_type](@ref pm_container::cmc_type)<br>
    !>  [ccc_type](@ref pm_container::ccc_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csc_pdt}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csc_pdt(kind)
        integer             , kind          :: kind = CK
        complex(kind)       , allocatable   :: val
    end type

    !>  \cond excluded
    interface csc_pdt

#if CK5_ENABLED
    pure elemental module function constructCon_D0_PSCK5(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)                :: val
        type(csc_pdt(CKG))                              :: container
    end function
#endif

#if CK4_ENABLED
    pure elemental module function constructCon_D0_PSCK4(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)                :: val
        type(csc_pdt(CKG))                              :: container
    end function
#endif

#if CK3_ENABLED
    pure elemental module function constructCon_D0_PSCK3(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)                :: val
        type(csc_pdt(CKG))                              :: container
    end function
#endif

#if CK2_ENABLED
    pure elemental module function constructCon_D0_PSCK2(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)                :: val
        type(csc_pdt(CKG))                              :: container
    end function
#endif

#if CK1_ENABLED
    pure elemental module function constructCon_D0_PSCK1(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)                :: val
        type(csc_pdt(CKG))                              :: container
    end function
#endif

    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of an `allocatable` vector component of type `complex` of arbitrary kind.
    type :: cvc_pdt(kind)
        integer             , kind          :: kind = CK
        complex(kind)       , allocatable   :: val(:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a matrix component of type `complex` of arbitrary kind.
    type :: cmc_pdt(kind)
        integer             , kind          :: kind = CK
        complex(kind)       , allocatable   :: val(:,:)
    end type

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a cube component of type `complex` of arbitrary kind.
    type :: ccc_pdt(kind)
        integer             , kind          :: kind = CK
        complex(kind)       , allocatable   :: val(:,:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [csr_pdt](@ref pm_container::csr_pdt) parameterized type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>r</b>eal objects of kind \RKALL.<br>
    !>
    !>  \details
    !>  The constructor of this derived type takes the following arguments.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of type `complex` of kind \RKALL.<br>
    !>
    !>  \return
    !>  `container`         :   The output object of type [csr_pdt(kind(val))](@ref pm_container::csr_pdt) of the same kind, rank, and shape
    !>                          as the input `val`, each element of which contains the corresponding element of the input `val`.
    !>
    !>  \interface{csr_pdt}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: csr_pdt
    !>      type(csr_pdt(kind(val))) :: container
    !>
    !>      container = csr_pdt(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [csr_pdt](@ref pm_container::csr_pdt)<br>
    !>  [cvr_pdt](@ref pm_container::cvr_pdt)<br>
    !>  [cmr_pdt](@ref pm_container::cmr_pdt)<br>
    !>  [ccr_pdt](@ref pm_container::ccr_pdt)<br>
    !>  [csr_type](@ref pm_container::csr_type)<br>
    !>  [cvr_type](@ref pm_container::cvr_type)<br>
    !>  [cmr_type](@ref pm_container::cmr_type)<br>
    !>  [ccr_type](@ref pm_container::ccr_type)<br>
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \final{csc_pdt}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: csr_pdt(kind)
        integer             , kind          :: kind = RK
        real(kind)          , allocatable   :: val
    end type

    !>  \cond excluded
    interface csr_pdt

#if RK5_ENABLED
    pure elemental module function constructCon_D0_PSRK5(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                :: val
        type(csr_pdt(RKG))                              :: container
    end function
#endif

#if RK4_ENABLED
    pure elemental module function constructCon_D0_PSRK4(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                :: val
        type(csr_pdt(RKG))                              :: container
    end function
#endif

#if RK3_ENABLED
    pure elemental module function constructCon_D0_PSRK3(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                :: val
        type(csr_pdt(RKG))                              :: container
    end function
#endif

#if RK2_ENABLED
    pure elemental module function constructCon_D0_PSRK2(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                :: val
        type(csr_pdt(RKG))                              :: container
    end function
#endif

#if RK1_ENABLED
    pure elemental module function constructCon_D0_PSRK1(val) result(container)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructCon_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                :: val
        type(csr_pdt(RKG))                              :: container
    end function
#endif

    end interface
    !>  \endcond excluded

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of an `allocatable` vector component of type `real` of arbitrary kind.
    type :: cvr_pdt(kind)
        integer             , kind          :: kind = RK
        real(kind)          , allocatable   :: val(:)
    end type

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a matrix component of type `real` of arbitrary kind.
    type :: cmr_pdt(kind)
        integer             , kind          :: kind = RK
        real(kind)          , allocatable   :: val(:,:)
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the **parameterized derived type** for generating a container of a cube component of type `real` of arbitrary kind.
    type :: ccr_pdt(kind)
        integer             , kind          :: kind = RK
        real(kind)          , allocatable   :: val(:,:,:)
    end type

#endif
!PDT_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_isless
    !>  Generate and return the result of comparing the values of two input containers of scalar values using the `<` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same shape as the input array-like arguments of
    !>                          type `logical` of default kind \LK containing the result of the comparison of the
    !>                          values of the two input containers of scalar values via the `<` operator.
    !>
    !>  \interface{isless}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: operator(<)
    !>      logical(LK) :: itis
    !>
    !>      itis = con1 <= con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{isless}
    !>  \include{lineno} example/pm_container/isless/main.F90
    !>  \compilef{isless}
    !>  \output{isless}
    !>  \include{lineno} example/pm_container/isless/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{isless}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(<)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function isless_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isless_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isless_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isless_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isless_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isless_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isless_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isless_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isless_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isless_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isless_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isless_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isless_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isless_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isless_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isless_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isless_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isless_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isless_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isless_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isless_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isless_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isless_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isless_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isless_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isless_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isless_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isless_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isless_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isless_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isless_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_ismore
    !>  Generate and return the result of comparing the values of two input containers of scalar values using the `>` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same shape as the input array-like arguments of
    !>                          type `logical` of default kind \LK containing the result of the comparison of the
    !>                          values of the two input containers of scalar values via the `>` operator.
    !>
    !>  \interface{ismore}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: operator(>)
    !>      logical(LK) :: itis
    !>
    !>      itis = con1 > con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{ismore}
    !>  \include{lineno} example/pm_container/ismore/main.F90
    !>  \compilef{ismore}
    !>  \output{ismore}
    !>  \include{lineno} example/pm_container/ismore/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with \gfortran{10.3}.
    !>  \code{.F90}
    !>
    !>      program elementalBugGfortran10
    !>
    !>          implicit none
    !>
    !>          type :: strc
    !>              character(:, SK), allocatable :: value
    !>          end type strc
    !>
    !>          write(*,*) "ismore([strc('A'), strc('B'), strc('C')], strc('B'))"
    !>          write(*,*)  ismore([strc('A'), strc('B'), strc('C')], strc('B'))
    !>
    !>          write(*,*) "['A', 'B', 'C'] > 'B'"
    !>          write(*,*)  ['A', 'B', 'C'] > 'B'
    !>
    !>      contains
    !>
    !>          pure elemental function ismore(con1, con2) result(itis)
    !>              type(strc), intent(in) :: con1
    !>              type(strc), intent(in) :: con2
    !>              logical(LK) :: itis
    !>              itis = con1%val > con2%val
    !>          end function
    !>
    !>      end program elementalBugGfortran10
    !>
    !>  \endcode
    !>  The above code yields a correct answer when compiled with \ifort{2021.5}.
    !>  \code
    !>      ismore([strc('A'), strc('B'), strc('C')], strc('B'))
    !>      F F T
    !>      ['A', 'B', 'C'] > 'B'
    !>      F F T
    !>  \endcode
    !>  However, it yields an incorrect answer with \gfortran{10.3-11}.
    !>  \code
    !>      ismore([strc('A'), strc('B'), strc('C')], strc('B'))
    !>      F T T
    !>      ['A', 'B', 'C'] > 'B'
    !>      F F T
    !>  \endcode
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{ismore}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(>)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function ismore_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismore_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismore_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismore_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismore_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function ismore_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function ismore_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function ismore_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function ismore_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function ismore_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function ismore_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function ismore_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function ismore_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function ismore_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function ismore_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function ismore_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function ismore_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function ismore_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function ismore_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function ismore_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function ismore_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function ismore_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function ismore_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function ismore_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function ismore_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function ismore_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function ismore_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function ismore_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function ismore_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function ismore_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismore_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_isleq
    !>  Generate and return the result of comparing the values of two input containers of scalar values using the `<=` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same shape as the input array-like arguments of
    !>                          type `logical` of default kind \LK containing the result of the comparison of the
    !>                          values of the two input containers of scalar values via the `<=` operator.
    !>
    !>  \interface{isleq}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: operator(<=)
    !>      logical(LK) :: itis
    !>
    !>      itis = con1 <= con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{isleq}
    !>  \include{lineno} example/pm_container/isleq/main.F90
    !>  \compilef{isleq}
    !>  \output{isleq}
    !>  \include{lineno} example/pm_container/isleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{isleq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(<=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function isleq_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isleq_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isleq_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isleq_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isleq_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isleq_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isleq_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isleq_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isleq_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isleq_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isleq_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isleq_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isleq_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isleq_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isleq_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isleq_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isleq_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isleq_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isleq_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isleq_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isleq_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isleq_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isleq_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isleq_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isleq_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isleq_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isleq_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isleq_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isleq_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isleq_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isleq_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_ismeq
    !>  Generate and return the result of comparing the values of two input containers of scalar values using the `>=` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same shape as the input array-like arguments of
    !>                          type `logical` of default kind \LK containing the result of the comparison of the
    !>                          values of the two input containers of scalar values via the `>=` operator.
    !>
    !>  \interface{ismeq}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: operator(>=)
    !>      logical(LK) :: itis
    !>
    !>      itis = con1 >= con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{ismeq}
    !>  \include{lineno} example/pm_container/ismeq/main.F90
    !>  \compilef{ismeq}
    !>  \output{ismeq}
    !>  \include{lineno} example/pm_container/ismeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{ismeq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(>=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function ismeq_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismeq_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismeq_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismeq_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function ismeq_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function ismeq_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function ismeq_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function ismeq_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function ismeq_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function ismeq_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function ismeq_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function ismeq_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function ismeq_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function ismeq_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function ismeq_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function ismeq_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function ismeq_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function ismeq_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function ismeq_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function ismeq_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function ismeq_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function ismeq_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function ismeq_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function ismeq_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function ismeq_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function ismeq_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function ismeq_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function ismeq_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function ismeq_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function ismeq_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ismeq_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_isneq
    !>  Generate and return the result of comparing the values of two input containers of scalar values using the `/=` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same shape as the input array-like arguments of
    !>                          type `logical` of default kind \LK containing the result of the comparison of the
    !>                          values of the two input containers of scalar values via the `/=` operator.
    !>
    !>  \interface{isneq}
    !>  \code{.F90}
    !>
    !>      use pm_container, only: operator(/=)
    !>      logical(LK) :: itis
    !>
    !>      itis = con1 /= con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{isneq}
    !>  \include{lineno} example/pm_container/isneq/main.F90
    !>  \compilef{isneq}
    !>  \output{isneq}
    !>  \include{lineno} example/pm_container/isneq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{isneq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(/=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function isneq_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isneq_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isneq_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isneq_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function isneq_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isneq_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isneq_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isneq_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isneq_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isneq_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isneq_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isneq_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isneq_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isneq_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isneq_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function isneq_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function isneq_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function isneq_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function isneq_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function isneq_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isneq_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isneq_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isneq_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isneq_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isneq_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isneq_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isneq_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isneq_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isneq_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isneq_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isneq_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_iseq
    !>  Generate and return the result of comparing the values of the two input containers of scalar of arbitrary type and kind using the `==` operator.
    !>
    !>  \param[in]  con1    :   The input scalar or array of the same rank and shape as the input array-like `con2` of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>  \param[in]  con2    :   The input scalar, or array of the same rank as the input array-like `con1`, of the same type and kind as `con1`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar or array of the same rank and shape as the input array-like arguments, of type `logical` of default kind \LK,
    !>                          containing the result of the comparison of the values of the two input containers via the `==` operator.
    !>
    !>  \interface{iseq}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_container, only: operator(==), css_pdt
    !>      type(css_pdt) :: con1, con2
    !>      logical(LK) :: equal
    !>
    !>      equal = con1 == con2
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The ranks and sizes of all array-like input arguments must match.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{iseq}
    !>  \include{lineno} example/pm_container/iseq/main.F90
    !>  \compilef{iseq}
    !>  \output{iseq}
    !>  \include{lineno} example/pm_container/iseq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{iseq}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface operator(==)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function iseq_D0_D0_BSSK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_BSSK
#endif
        type(css_type)          , intent(in)                :: con1
        type(css_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function iseq_D0_D0_BSIK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_BSIK
#endif
        type(csi_type)          , intent(in)                :: con1
        type(csi_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function iseq_D0_D0_BSLK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_BSLK
#endif
        type(csl_type)          , intent(in)                :: con1
        type(csl_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function iseq_D0_D0_BSCK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_BSCK
#endif
        type(csc_type)          , intent(in)                :: con1
        type(csc_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

    pure elemental module function iseq_D0_D0_BSRK(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_BSRK
#endif
        type(csr_type)          , intent(in)                :: con1
        type(csr_type)          , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function iseq_D0_D0_PSSK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK4_ENABLED
    pure elemental module function iseq_D0_D0_PSSK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK3_ENABLED
    pure elemental module function iseq_D0_D0_PSSK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK2_ENABLED
    pure elemental module function iseq_D0_D0_PSSK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if SK1_ENABLED
    pure elemental module function iseq_D0_D0_PSSK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(in)                :: con1
        type(css_pdt(SKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function iseq_D0_D0_PSIK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK4_ENABLED
    pure elemental module function iseq_D0_D0_PSIK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK3_ENABLED
    pure elemental module function iseq_D0_D0_PSIK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK2_ENABLED
    pure elemental module function iseq_D0_D0_PSIK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if IK1_ENABLED
    pure elemental module function iseq_D0_D0_PSIK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(in)                :: con1
        type(csi_pdt(IKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function iseq_D0_D0_PSLK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK4_ENABLED
    pure elemental module function iseq_D0_D0_PSLK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK3_ENABLED
    pure elemental module function iseq_D0_D0_PSLK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK2_ENABLED
    pure elemental module function iseq_D0_D0_PSLK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if LK1_ENABLED
    pure elemental module function iseq_D0_D0_PSLK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(in)                :: con1
        type(csl_pdt(LKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function iseq_D0_D0_PSCK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK4_ENABLED
    pure elemental module function iseq_D0_D0_PSCK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK3_ENABLED
    pure elemental module function iseq_D0_D0_PSCK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK2_ENABLED
    pure elemental module function iseq_D0_D0_PSCK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if CK1_ENABLED
    pure elemental module function iseq_D0_D0_PSCK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(in)                :: con1
        type(csc_pdt(CKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function iseq_D0_D0_PSRK5(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK4_ENABLED
    pure elemental module function iseq_D0_D0_PSRK4(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK3_ENABLED
    pure elemental module function iseq_D0_D0_PSRK3(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK2_ENABLED
    pure elemental module function iseq_D0_D0_PSRK2(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

#if RK1_ENABLED
    pure elemental module function iseq_D0_D0_PSRK1(con1, con2) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(in)                :: con1
        type(csr_pdt(RKG))      , intent(in)                :: con2
        logical(LK)                                         :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_container_assign
    !>  Assign the contents of the input `source` to the output `destin`.
    !>
    !>  \param[out] destin  :   The output scalar or array of the arbitrary rank and shape that can be of,
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type)
    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt)
    !>                              <li>    type [csi_pdt](@ref pm_container::csi_pdt)
    !>                              <li>    type [csl_pdt](@ref pm_container::csl_pdt)
    !>                              <li>    type [csc_pdt](@ref pm_container::csc_pdt)
    !>                              <li>    type [csr_pdt](@ref pm_container::csr_pdt)
    !>                          </ol>
    !>                          whose contents will be the taken from the contents of the input `source` argument upon return.
    !>  \param[in]  source  :   The input scalar, or array of the same rank as the input array-like `destin`,
    !>                          of the same type and kind as `destin`, whose contents will be assigned to the output `destin`.
    !>
    !>  \interface{assignment}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_container, only: assignment(=)
    !>
    !>      destin = source
    !>      destin(..) = source
    !>      destin(..) = source(..)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `rank(source) == 0 .or. rank(destin) == rank(source)` must hold (i.e., only `elemental` assignments are possible).<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  This generic primarily exists to circumvent the gfortran bug as of version 11 in intrinsic assignment of arrays of type [css_pdt](@ref pm_container::css_pdt).<br>
    !>
    !>  \see
    !>  [operator(==)](@ref pm_container_iseq)<br>
    !>  [operator(<)](@ref pm_container_isless)<br>
    !>  [operator(>)](@ref pm_container_ismore)<br>
    !>  [operator(>=)](@ref pm_container_ismeq)<br>
    !>  [operator(<=)](@ref pm_container_isleq)<br>
    !>  [operator(/=)](@ref pm_container_isneq)<br>
    !>  [assignment(=)](@ref pm_container_assign)<br>
    !>
    !>  \example{assignment}
    !>  \include{lineno} example/pm_container/assign/main.F90
    !>  \compilef{assignment}
    !>  \output{assignment}
    !>  \include{lineno} example/pm_container/assign/main.out.F90
    !>
    !>  \test
    !>  [test_pm_container](@ref test_pm_container)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The elemental implementations of the procedures under this generic interface yield incorrect results with gfortran 10.3.
    !>  \remedy
    !>  Currently unknown.
    !>
    !>  \todo
    !>  \pvlow The functionality of this generic interface can be extended to input arrays of higher rank.
    !>
    !>  \final{assignment}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface assignment(=)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module subroutine assign_D0_D0_BSSK(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_BSSK
#endif
        type(css_type)          , intent(out)               :: destin
        type(css_type)          , intent(in)                :: source
    end subroutine

    pure elemental module subroutine assign_D0_D0_BSIK(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_BSIK
#endif
        type(csi_type)          , intent(out)               :: destin
        type(csi_type)          , intent(in)                :: source
    end subroutine

    pure elemental module subroutine assign_D0_D0_BSLK(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_BSLK
#endif
        type(csl_type)          , intent(out)               :: destin
        type(csl_type)          , intent(in)                :: source
    end subroutine

    pure elemental module subroutine assign_D0_D0_BSCK(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_BSCK
#endif
        type(csc_type)          , intent(out)               :: destin
        type(csc_type)          , intent(in)                :: source
    end subroutine

    pure elemental module subroutine assign_D0_D0_BSRK(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_BSRK
#endif
        type(csr_type)          , intent(out)               :: destin
        type(csr_type)          , intent(in)                :: source
    end subroutine

#if PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine assign_D0_D0_PSSK5(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))      , intent(out)               :: destin
        type(css_pdt(SKG))      , intent(in)                :: source
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine assign_D0_D0_PSSK4(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))      , intent(out)               :: destin
        type(css_pdt(SKG))      , intent(in)                :: source
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine assign_D0_D0_PSSK3(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))      , intent(out)               :: destin
        type(css_pdt(SKG))      , intent(in)                :: source
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine assign_D0_D0_PSSK2(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))      , intent(out)               :: destin
        type(css_pdt(SKG))      , intent(in)                :: source
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine assign_D0_D0_PSSK1(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))      , intent(out)               :: destin
        type(css_pdt(SKG))      , intent(in)                :: source
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module subroutine assign_D0_D0_PSIK5(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSIK5
#endif
        use pm_kind, only: IKG => IK5
        type(csi_pdt(IKG))      , intent(out)               :: destin
        type(csi_pdt(IKG))      , intent(in)                :: source
    end subroutine
#endif

#if IK4_ENABLED
    pure elemental module subroutine assign_D0_D0_PSIK4(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSIK4
#endif
        use pm_kind, only: IKG => IK4
        type(csi_pdt(IKG))      , intent(out)               :: destin
        type(csi_pdt(IKG))      , intent(in)                :: source
    end subroutine
#endif

#if IK3_ENABLED
    pure elemental module subroutine assign_D0_D0_PSIK3(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSIK3
#endif
        use pm_kind, only: IKG => IK3
        type(csi_pdt(IKG))      , intent(out)               :: destin
        type(csi_pdt(IKG))      , intent(in)                :: source
    end subroutine
#endif

#if IK2_ENABLED
    pure elemental module subroutine assign_D0_D0_PSIK2(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSIK2
#endif
        use pm_kind, only: IKG => IK2
        type(csi_pdt(IKG))      , intent(out)               :: destin
        type(csi_pdt(IKG))      , intent(in)                :: source
    end subroutine
#endif

#if IK1_ENABLED
    pure elemental module subroutine assign_D0_D0_PSIK1(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSIK1
#endif
        use pm_kind, only: IKG => IK1
        type(csi_pdt(IKG))      , intent(out)               :: destin
        type(csi_pdt(IKG))      , intent(in)                :: source
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module subroutine assign_D0_D0_PSLK5(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSLK5
#endif
        use pm_kind, only: LKG => LK5
        type(csl_pdt(LKG))      , intent(out)               :: destin
        type(csl_pdt(LKG))      , intent(in)                :: source
    end subroutine
#endif

#if LK4_ENABLED
    pure elemental module subroutine assign_D0_D0_PSLK4(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSLK4
#endif
        use pm_kind, only: LKG => LK4
        type(csl_pdt(LKG))      , intent(out)               :: destin
        type(csl_pdt(LKG))      , intent(in)                :: source
    end subroutine
#endif

#if LK3_ENABLED
    pure elemental module subroutine assign_D0_D0_PSLK3(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSLK3
#endif
        use pm_kind, only: LKG => LK3
        type(csl_pdt(LKG))      , intent(out)               :: destin
        type(csl_pdt(LKG))      , intent(in)                :: source
    end subroutine
#endif

#if LK2_ENABLED
    pure elemental module subroutine assign_D0_D0_PSLK2(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSLK2
#endif
        use pm_kind, only: LKG => LK2
        type(csl_pdt(LKG))      , intent(out)               :: destin
        type(csl_pdt(LKG))      , intent(in)                :: source
    end subroutine
#endif

#if LK1_ENABLED
    pure elemental module subroutine assign_D0_D0_PSLK1(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSLK1
#endif
        use pm_kind, only: LKG => LK1
        type(csl_pdt(LKG))      , intent(out)               :: destin
        type(csl_pdt(LKG))      , intent(in)                :: source
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module subroutine assign_D0_D0_PSCK5(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSCK5
#endif
        use pm_kind, only: CKG => CK5
        type(csc_pdt(CKG))      , intent(out)               :: destin
        type(csc_pdt(CKG))      , intent(in)                :: source
    end subroutine
#endif

#if CK4_ENABLED
    pure elemental module subroutine assign_D0_D0_PSCK4(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSCK4
#endif
        use pm_kind, only: CKG => CK4
        type(csc_pdt(CKG))      , intent(out)               :: destin
        type(csc_pdt(CKG))      , intent(in)                :: source
    end subroutine
#endif

#if CK3_ENABLED
    pure elemental module subroutine assign_D0_D0_PSCK3(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSCK3
#endif
        use pm_kind, only: CKG => CK3
        type(csc_pdt(CKG))      , intent(out)               :: destin
        type(csc_pdt(CKG))      , intent(in)                :: source
    end subroutine
#endif

#if CK2_ENABLED
    pure elemental module subroutine assign_D0_D0_PSCK2(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSCK2
#endif
        use pm_kind, only: CKG => CK2
        type(csc_pdt(CKG))      , intent(out)               :: destin
        type(csc_pdt(CKG))      , intent(in)                :: source
    end subroutine
#endif

#if CK1_ENABLED
    pure elemental module subroutine assign_D0_D0_PSCK1(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSCK1
#endif
        use pm_kind, only: CKG => CK1
        type(csc_pdt(CKG))      , intent(out)               :: destin
        type(csc_pdt(CKG))      , intent(in)                :: source
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine assign_D0_D0_PSRK5(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSRK5
#endif
        use pm_kind, only: RKG => RK5
        type(csr_pdt(RKG))      , intent(out)               :: destin
        type(csr_pdt(RKG))      , intent(in)                :: source
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine assign_D0_D0_PSRK4(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSRK4
#endif
        use pm_kind, only: RKG => RK4
        type(csr_pdt(RKG))      , intent(out)               :: destin
        type(csr_pdt(RKG))      , intent(in)                :: source
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine assign_D0_D0_PSRK3(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSRK3
#endif
        use pm_kind, only: RKG => RK3
        type(csr_pdt(RKG))      , intent(out)               :: destin
        type(csr_pdt(RKG))      , intent(in)                :: source
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine assign_D0_D0_PSRK2(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSRK2
#endif
        use pm_kind, only: RKG => RK2
        type(csr_pdt(RKG))      , intent(out)               :: destin
        type(csr_pdt(RKG))      , intent(in)                :: source
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine assign_D0_D0_PSRK1(destin, source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: assign_D0_D0_PSRK1
#endif
        use pm_kind, only: RKG => RK1
        type(csr_pdt(RKG))      , intent(out)               :: destin
        type(csr_pdt(RKG))      , intent(in)                :: source
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
!PDT_ENABLED

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! >  \brief
    ! >  This is the [csis_type](@ref pm_container::csis_type) type for generating instances of <b>c</b>ontainer of <b>s</b>calar of <b>i</b>nteger and <b>s</b>tring.<br>
    ! >
    ! >  \details
    ! >  The [csis_type](@ref pm_container::csis_type) container is frequently used within the ParaMonte library for encapsulate an `integer` number and its conversion to string.<br>
    ! >  The constructor of this derived type takes the following arguments.<br>
    ! >
    ! >  \param[in]  int     :   The input scalar or array of arbitrary rank of type `integer` of default kind \IK.<br>
    ! >  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of default kind \SK.<br>
    ! >
    ! >  \return
    ! >  `container`         :   The output object of type [csis_type](@ref pm_container::csis_type) of the same kind, rank, and shape
    ! >                          as the input array-like arguments `val`, each element of which contains the corresponding element of the input `val`.
    ! >
    ! >  \interface{csis_type}
    ! >  \code{.F90}
    ! >
    ! >      use pm_kind, only: SK
    ! >      use pm_container, only: csis_type
    ! >      type(csis_type), allocatable  :: containers(:)
    ! >      type(csis_type) :: container
    ! >
    ! >      container = csis_type(val, trimmed = trimmed)
    ! >
    ! >      ! example
    ! >
    ! >      container = csis_type(SK_"aaa")
    ! >      containers = csis_type([character(3, SK) :: "a", "aa", "aaa"])
    ! >      containers = csis_type([character(3, SK) :: "a", "aa", "aaa"], trimmed = .false._LK)
    ! >
    ! >  \endcode
    ! >
    ! >  \pure
    ! >
    ! >  \elemental
    ! >
    ! >  \remark
    ! >  The existence of this **derived type** is primarily motivated by the persistent bugs in the
    ! >  implementations of `allocatable` parameterized derived types in gfortran compilers < 12.<br>
    ! >  As such, this definition might be redundant once gfortran PDT bugs are resolved.<br>
    ! >
    ! >  \see
    ! >  [css_pdt](@ref pm_container::css_pdt)<br>
    ! >  [cvs_pdt](@ref pm_container::cvs_pdt)<br>
    ! >  [cms_pdt](@ref pm_container::cms_pdt)<br>
    ! >  [ccs_pdt](@ref pm_container::ccs_pdt)<br>
    ! >  [css_type](@ref pm_container::css_type)<br>
    ! >  [cvs_type](@ref pm_container::cvs_type)<br>
    ! >  [cms_type](@ref pm_container::cms_type)<br>
    ! >  [ccs_type](@ref pm_container::ccs_type)<br>
    ! >
    ! >  \example{csis_type}
    ! >  \include{lineno} example/pm_container/csis_type/main.F90
    ! >  \compile
    ! >  \output
    ! >  \include{lineno} example/pm_container/csis_type/main.out.F90
    ! >
    ! >  \test
    ! >  [test_pm_container](@ref test_pm_container)
    ! >
    ! >  \final{csis_type}
    ! >
    ! >  \author
    ! >  \FatemehBagheri, Tuesday April 30, 2019, 12:58 PM, SEIR, UTA
    ! >  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    !type :: csis_type
    !    integer(IK) :: int
    !    character(:, SK), allocatable :: str
    !end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Generate and return the extraction of the `val` component of the input scalar container.
!    !>
!    !>  \details
!    !>  In simple words, this generic interface extracts the value of the input container(s).<br>
!    !>  When the input container(s) is of type [css_type](@ref pm_container::css_type),
!    !>  the output values will have a length type parameter that
!    !>  is the maximum of the input `val` components.<br>
!    !>
!    !>  \param[in]  con     :   The input scalar or array of arbitrary rank `1` of,
!    !>                          <ol>
!    !>                              <li>    type [css_type](@ref pm_container::css_type)
!    !>                              <li>    type [csi_type](@ref pm_container::csi_type)
!    !>                              <li>    type [csl_type](@ref pm_container::csl_type)
!    !>                              <li>    type [csc_type](@ref pm_container::csc_type)
!    !>                              <li>    type [csr_type](@ref pm_container::csr_type)
!    !>                          </ol>
!    !>
!    !>  \return
!    !>  `val`               :   The output scalar or array of the same rank and shape as the input array-like argument `con`,
!    !>                          of the same type, kind, and value as the `val` component of the input container(s).
!    !>
!    !>  \interface{getVal}
!    !>  \code{.F90}
!    !>
!    !>      use pm_container, only: getVal
!    !>
!    !>      val(@shape(con)) = getVal(con(..))
!    !>
!    !>  \endcode
!    !>
!    !>  \note
!    !>  If the `val` component of the input argument `con` is unallocated,
!    !>  the output `val` will remain uninitialized on output.
!    !>
!    !>  \pure
!    !>
!    !>  \elemental
!    !>
!    !>  \see
!    !>  [getCharSeq](@ref pm_str::getCharSeq)<br>
!    !>  [getCharVec](@ref pm_str::getCharVec)<br>
!    !>
!    !>  \example{getVal}
!    !>  \include{lineno} example/pm_container/getVal/main.F90
!    !>  \compilef{getVal}
!    !>  \output{getVal}
!    !>  \include{lineno} example/pm_container/getVal/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_container](@ref test_pm_container)
!    !>
!    !>  \todo
!    !>  \pvlow The functionality of this generic interface can be extended to input containers of higher rank component.<br>
!    !>
!    !>  \final{getVal}
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
!    interface getVal
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    pure elemental module function getVal_D0_BSSK(con) result(val)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getVal_D0_BSSK
!#endif
!        use pm_kind, only: SKG => SK
!        type(css_type)          , intent(in)                :: con
!        character(len(con%val, IK),SKG)                     :: val
!    end function
!
!    pure elemental module function getVal_D0_BSIK(con) result(val)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getVal_D0_BSIK
!#endif
!        use pm_kind, only: IKG => IK
!        type(csi_type)          , intent(in)                :: con
!        integer(IKG)                                        :: val
!    end function
!
!    pure elemental module function getVal_D0_BSLK(con) result(val)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getVal_D0_BSLK
!#endif
!        use pm_kind, only: LKG => LK
!        type(csl_type)          , intent(in)                :: con
!        logical(LKG)                                        :: val
!    end function
!
!    pure elemental module function getVal_D0_BSCK(con) result(val)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getVal_D0_BSCK
!#endif
!        use pm_kind, only: CKG => CK
!        type(csc_type)          , intent(in)                :: con
!        complex(CKG)                                        :: val
!    end function
!
!    pure elemental module function getVal_D0_BSRK(con) result(val)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getVal_D0_BSRK
!#endif
!        use pm_kind, only: RKG => RK
!        type(csr_type)          , intent(in)                :: con
!        real(RKG)                                           :: val
!    end function
!
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_container ! LCOV_EXCL_LINE