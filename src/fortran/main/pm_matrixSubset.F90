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
!>  This module contains abstract and concrete derived types that are required for compile-time
!>  resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.<br>
!>  Such procedures frequently need to work on the upper/lower-diagonal triangular blocks of some of their input matrix arguments.<br>
!>
!>  \details
!>  Within this module, a **lower**-triangular storage for a subset of an arbitrary matrix has the form,
!>  \f{equation}{
!>      L = {\begin{bmatrix}\ell _{1,1}&&&&\\\ell _{2,1}&\ell _{2,2}&&&\\\ell _{3,1}&\ell _{3,2}&\ddots &&\\\vdots &\vdots &\ddots &\ddots &\\\ell _{n,1}&\ell _{n,2}&\ldots &\ell _{n,n-1}&\ell _{n,n}\end{bmatrix}} ~,
!>  \f}
!>  while an **upper**-triangular storage for a subset of an arbitrary matrix has the form,
!>  \f{equation}{
!>      \ms{U} = {\begin{bmatrix}u_{1,1}&u_{1,2}&u_{1,3}&\ldots &u_{1,n}\\&u_{2,2}&u_{2,3}&\ldots &u_{2,n}\\&&\ddots &\ddots &\vdots \\&&&\ddots &u_{n-1,n}\\0&&&&u_{n,n}\end{bmatrix}} ~.
!>  \f}
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixSubset

    use pm_kind, only: SK
    use pm_array, only: nothing_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixSubset"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of storage (upper-diagonal triangular, lower-diagonal triangular, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{subset_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: subset_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  upper-triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [upp](@ref pm_matrixSubset::upp)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{upp_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: upp_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [upp_type](@ref pm_matrixSubset::upp_type) that is exclusively used to request
    !>  upper-triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{upp}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(upp_type), parameter :: upp = upp_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: upp
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  lower-triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [low](@ref pm_matrixSubset::low)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{low_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: low_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [low_type](@ref pm_matrixSubset::low_type) that is exclusively used to request
    !>  lower-triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{low}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(low_type), parameter :: low = low_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: low
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  unit (or Identity or diagonal) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [dia](@ref pm_matrixSubset::dia)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{dia_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: dia_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [dia_type](@ref pm_matrixSubset::dia_type) that is exclusively used to request
    !>  unit (or Identity or diagonal) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{dia}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(dia_type), parameter :: dia = dia_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: dia
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  upper-lower triangular (<b>excluding diagonal</b>) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [uppLow](@ref pm_matrixSubset::uppLow)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppLow_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: uppLow_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [uppLow_type](@ref pm_matrixSubset::uppLow_type) that is exclusively used to request
    !>  upper-lower triangular (<b>excluding diagonal</b>) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppLow}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(uppLow_type), parameter :: uppLow = uppLow_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: uppLow
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  full diagonal and upper-lower triangular (<b>excluding diagonal</b>) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [uppLowDia](@ref pm_matrixSubset::uppLowDia)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppLowDia_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: uppLowDia_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type) that is exclusively used to request
    !>  full diagonal and upper-lower triangular (<b>excluding diagonal</b>) storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppLowDia}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(uppLowDia_type), parameter :: uppLowDia = uppLowDia_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: uppLowDia
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  upper-diagonal triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [uppDia](@ref pm_matrixSubset::uppDia)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppDia_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: uppDia_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type) that is exclusively used to request
    !>  upper-diagonal triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{uppDia}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(uppDia_type), parameter :: uppDia = uppDia_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: uppDia
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  lower-diagonal triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lowDia](@ref pm_matrixSubset::lowDia)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{lowDia_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(subset_type) :: lowDia_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type) that is exclusively used to request
    !>  lower-diagonal triangular storage format of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [upp](@ref pm_matrixSubset::upp)<br>
    !>  [low](@ref pm_matrixSubset::low)<br>
    !>  [dia](@ref pm_matrixSubset::dia)<br>
    !>  [uppLow](@ref pm_matrixSubset::uppLow)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppLowDia](@ref pm_matrixSubset::uppLowDia)<br>
    !>  [dia_type](@ref pm_matrixSubset::dia_type)<br>
    !>  [uppLow_type](@ref pm_matrixSubset::uppLow_type)<br>
    !>  [uppDia_type](@ref pm_matrixSubset::uppDia_type)<br>
    !>  [lowDia_type](@ref pm_matrixSubset::lowDia_type)<br>
    !>  [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type)<br>
    !>  [subset_type](@ref pm_matrixSubset::subset_type)<br>
    !>
    !>  \final{lowDia}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(lowDia_type), parameter :: lowDia = lowDia_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lowDia
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the objects representing the complementary
    !>  subset of the input matrix subset `sub` with respect to the main diagonal of the matrix.
    !>
    !>  \param[in]  sub  :   The input scalar (or array of arbitrary rank and shape) of,
    !>                          <ol>
    !>                              <li>    type [upp_type](@ref pm_matrixSubset::upp_type),
    !>                              <li>    type [low_type](@ref pm_matrixSubset::low_type),
    !>                              <li>    type [dia_type](@ref pm_matrixSubset::dia_type),
    !>                              <li>    type [uppDia_type](@ref pm_matrixSubset::uppDia_type),
    !>                              <li>    type [lowDia_type](@ref pm_matrixSubset::lowDia_type),
    !>                              <li>    type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                              <li>    type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                          </ol>
    !>
    !>  \return
    !>  `subComp`           :   The output object of the same rank and shape as the input subset `sub` whose
    !>                          type represents the matrix subset corresponding to the complementary subset of
    !>                          the input `sub` with respect to the main diagonal of the matrix.<br>
    !>                          <ol>
    !>                              <li>    If `sub` is of type [upp_type](@ref pm_matrixSubset::upp_type), then the output `subComp` is of type [lowDia_type](@ref pm_matrixSubset::lowDia_type).
    !>                              <li>    If `sub` is of type [low_type](@ref pm_matrixSubset::low_type), then the output `subComp` is of type [uppDia_type](@ref pm_matrixSubset::uppDia_type).
    !>                              <li>    If `sub` is of type [dia_type](@ref pm_matrixSubset::dia_type), then the output `subComp` is of type [uppLow_type](@ref pm_matrixSubset::uppLow_type).
    !>                              <li>    If `sub` is of type [uppDia_type](@ref pm_matrixSubset::uppDia_type), then the output `subComp` is of type [low_type](@ref pm_matrixSubset::low_type).
    !>                              <li>    If `sub` is of type [lowDia_type](@ref pm_matrixSubset::lowDia_type), then the output `subComp` is of type [upp_type](@ref pm_matrixSubset::upp_type).
    !>                              <li>    If `sub` is of type [uppLow_type](@ref pm_matrixSubset::uppLow_type), then the output `subComp` is of type [dia_type](@ref pm_matrixSubset::dia_type).
    !>                              <li>    If `sub` is of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type), then the output `subComp` is of type [nothing_type](@ref pm_array::nothing_type).
    !>                          </ol>
    !>
    !>  \interface{getSubComp}
    !>  \code{.F90}
    !>
    !>      use pm_matrixSubset, only: getSubComp
    !>
    !>      subComp(..) =  getSubComp(sub(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getSubComp](@ref pm_matrixSubset::getSubComp)<br>
    !>  [getSubSymm](@ref pm_matrixSubset::getSubSymm)<br>
    !>  [getSubUnion](@ref pm_matrixSubset::getSubUnion)<br>
    !>
    !>  \final{getSubComp}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface getSubComp

    pure elemental module function getSubCompUXX(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompUXX
#endif
        type(upp_type)      , intent(in)    :: sub
        type(lowDia_type)                   :: subComp
    end function

    pure elemental module function getSubCompXLX(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompXLX
#endif
        type(low_type)      , intent(in)    :: sub
        type(uppDia_type)                   :: subComp
    end function

    pure elemental module function getSubCompXXD(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompXXD
#endif
        type(dia_type)      , intent(in)    :: sub
        type(uppLow_type)                   :: subComp
    end function

    pure elemental module function getSubCompUXD(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompUXD
#endif
        type(uppDia_type)   , intent(in)    :: sub
        type(low_type)                      :: subComp
    end function

    pure elemental module function getSubCompXLD(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompXLD
#endif
        type(lowDia_type)   , intent(in)    :: sub
        type(upp_type)                      :: subComp
    end function

    pure elemental module function getSubCompULX(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompULX
#endif
        type(uppLow_type)   , intent(in)    :: sub
        type(dia_type)                      :: subComp
    end function

    pure elemental module function getSubCompULD(sub) result(subComp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubCompULD
#endif
        type(uppLowDia_type), intent(in)    :: sub
        type(nothing_type)                  :: subComp
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the objects representing the symmetric mirror
    !>  subset of the input matrix subset `sub` with respect to the main diagonal of the matrix.
    !>
    !>  \param[in]  sub     :   The input scalar (or array of arbitrary rank and shape) of,
    !>                          <ol>
    !>                              <li>    type [upp_type](@ref pm_matrixSubset::upp_type),
    !>                              <li>    type [low_type](@ref pm_matrixSubset::low_type),
    !>                              <li>    type [dia_type](@ref pm_matrixSubset::dia_type),
    !>                              <li>    type [uppDia_type](@ref pm_matrixSubset::uppDia_type),
    !>                              <li>    type [lowDia_type](@ref pm_matrixSubset::lowDia_type),
    !>                              <li>    type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                              <li>    type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                          </ol>
    !>
    !>  \return
    !>  `subSymm`           :   The output object of the same rank and shape as the input subset `sub` whose
    !>                          type represents the matrix subset corresponding to the symmetric mirror subset of
    !>                          the input subset `sub` with respect to the main diagonal of the matrix.<br>
    !>                          <ol>
    !>                              <li>    If `sub` is of type [upp_type](@ref pm_matrixSubset::upp_type), then the output `subSymm` is of type [low_type](@ref pm_matrixSubset::low_type).
    !>                              <li>    If `sub` is of type [low_type](@ref pm_matrixSubset::low_type), then the output `subSymm` is of type [upp_type](@ref pm_matrixSubset::upp_type).
    !>                              <li>    If `sub` is of type [dia_type](@ref pm_matrixSubset::dia_type), then the output `subSymm` is of type [dia_type](@ref pm_matrixSubset::dia_type).
    !>                              <li>    If `sub` is of type [uppDia_type](@ref pm_matrixSubset::uppDia_type), then the output `subSymm` is of type [lowDia_type](@ref pm_matrixSubset::lowDia_type).
    !>                              <li>    If `sub` is of type [lowDia_type](@ref pm_matrixSubset::lowDia_type), then the output `subSymm` is of type [uppDia_type](@ref pm_matrixSubset::uppDia_type).
    !>                              <li>    If `sub` is of type [uppLow_type](@ref pm_matrixSubset::uppLow_type), then the output `subSymm` is of type [uppLow_type](@ref pm_matrixSubset::uppLow_type).
    !>                              <li>    If `sub` is of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type), then the output `subSymm` is of type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type).
    !>                          </ol>
    !>
    !>  \interface{getSubSymm}
    !>  \code{.F90}
    !>
    !>      use pm_matrixSubset, only: getSubSymm
    !>
    !>      subSymm(..) =  getSubSymm(sub(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getSubComp](@ref pm_matrixSubset::getSubComp)<br>
    !>  [getSubSymm](@ref pm_matrixSubset::getSubSymm)<br>
    !>  [getSubUnion](@ref pm_matrixSubset::getSubUnion)<br>
    !>
    !>  \final{getSubSymm}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface getSubSymm

    pure elemental module function getSubSymmUXX(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmUXX
#endif
        type(upp_type)      , intent(in)    :: sub
        type(low_type)                      :: subSymm
    end function

    pure elemental module function getSubSymmXLX(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmXLX
#endif
        type(low_type)      , intent(in)    :: sub
        type(upp_type)                      :: subSymm
    end function

    pure elemental module function getSubSymmXXD(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmXXD
#endif
        type(dia_type)      , intent(in)    :: sub
        type(dia_type)                      :: subSymm
    end function

    pure elemental module function getSubSymmUXD(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmUXD
#endif
        type(uppDia_type)   , intent(in)    :: sub
        type(lowDia_type)                   :: subSymm
    end function

    pure elemental module function getSubSymmXLD(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmXLD
#endif
        type(lowDia_type)   , intent(in)    :: sub
        type(uppDia_type)                   :: subSymm
    end function

    pure elemental module function getSubSymmULX(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmULX
#endif
        type(uppLow_type)   , intent(in)    :: sub
        type(uppLow_type)                   :: subSymm
    end function

    pure elemental module function getSubSymmULD(sub) result(subSymm)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubSymmULD
#endif
        type(uppLowDia_type), intent(in)    :: sub
        type(uppLowDia_type)                :: subSymm
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the object representing the union of two input matrix subsets.
    !>
    !>  \param[in]  sub1    :   The input scalar (or array of the same rank as the rank of other array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type [upp_type](@ref pm_matrixSubset::upp_type),
    !>                              <li>    type [low_type](@ref pm_matrixSubset::low_type),
    !>                              <li>    type [dia_type](@ref pm_matrixSubset::dia_type),
    !>                              <li>    type [uppDia_type](@ref pm_matrixSubset::uppDia_type),
    !>                              <li>    type [lowDia_type](@ref pm_matrixSubset::lowDia_type),
    !>                              <li>    type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                              <li>    type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                          </ol>
    !>  \param[in]  sub2    :   The input scalar (or array of the same rank as the rank of other array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type [upp_type](@ref pm_matrixSubset::upp_type),
    !>                              <li>    type [low_type](@ref pm_matrixSubset::low_type),
    !>                              <li>    type [dia_type](@ref pm_matrixSubset::dia_type),
    !>                              <li>    type [uppDia_type](@ref pm_matrixSubset::uppDia_type),
    !>                              <li>    type [lowDia_type](@ref pm_matrixSubset::lowDia_type),
    !>                              <li>    type [uppLow_type](@ref pm_matrixSubset::uppLow_type),
    !>                              <li>    type [uppLowDia_type](@ref pm_matrixSubset::uppLowDia_type),
    !>                          </ol>
    !>
    !>  \return
    !>  `subUnion`          :   The output scalar (or array of the same rank as the rank of other array-like arguments) of
    !>                          type resulting from the union of the two input matrix subsets.<br>
    !>                          For example,<br>
    !>                          <ol>
    !>                              <li>    the union of [upp](@ref pm_matrixSubset::upp) and [low](@ref pm_matrixSubset::low) yields [uppLow](@ref pm_matrixSubset::uppLow).
    !>                              <li>    the union of [upp](@ref pm_matrixSubset::upp) and [dia](@ref pm_matrixSubset::dia) yields [uppDia](@ref pm_matrixSubset::uppDia).
    !>                              <li>    the union of [uppDia](@ref pm_matrixSubset::uppDia) and [dia](@ref pm_matrixSubset::dia) yields [uppDia](@ref pm_matrixSubset::uppDia).
    !>                              <li>    the union of [uppDia](@ref pm_matrixSubset::uppDia) and [nothing](@ref pm_array::nothing) yields [uppDia](@ref pm_matrixSubset::uppDia).
    !>                          </ol>
    !>
    !>  \interface{getSubUnion}
    !>  \code{.F90}
    !>
    !>      use pm_matrixSubset, only: getSubUnion
    !>
    !>      subUnion(..) =  getSubUnion(sub1(..), sub2(..))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getSubComp](@ref pm_matrixSubset::getSubComp)<br>
    !>  [getSubSymm](@ref pm_matrixSubset::getSubSymm)<br>
    !>  [getSubUnion](@ref pm_matrixSubset::getSubUnion)<br>
    !>
    !>  \final{getSubUnion}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! XXX

    interface getSubUnion

    pure elemental module function getSubUnion_XXX_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_XXX
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(nothing_type)                  :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_UXX
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(upp_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_XLX
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(low_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_XXD
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(dia_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_UXD
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_XLD
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_ULX
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXX_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXX_ULD
#endif
        type(nothing_type)  , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! UXX

    interface getSubUnion

    pure elemental module function getSubUnion_UXX_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_XXX
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(upp_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_UXX
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(upp_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_XLX
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_XXD
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_UXD
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_XLD
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_ULX
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXX_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXX_ULD
#endif
        type(upp_type)      , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! XLX

    interface getSubUnion

    pure elemental module function getSubUnion_XLX_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_XXX
#endif
        type(low_type)      , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(low_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_UXX
#endif
        type(low_type)      , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_XLX
#endif
        type(low_type)      , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(low_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_XXD
#endif
        type(low_type)      , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_UXD
#endif
        type(low_type)      , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_XLD
#endif
        type(low_type)      , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_ULX
#endif
        type(low_type)      , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLX_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLX_ULD
#endif
        type(low_type)      , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! XXD

    interface getSubUnion

    pure elemental module function getSubUnion_XXD_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_XXX
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(dia_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_UXX
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_XLX
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_XXD
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(dia_type)                      :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_UXD
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_XLD
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_ULX
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_XXD_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XXD_ULD
#endif
        type(dia_type)      , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! UXD

    interface getSubUnion

    pure elemental module function getSubUnion_UXD_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_XXX
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_UXX
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_XLX
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_XXD
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_UXD
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_XLD
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_ULX
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_UXD_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_UXD_ULD
#endif
        type(uppDia_type)   , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! XLD

    interface getSubUnion

    pure elemental module function getSubUnion_XLD_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_XXX
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_UXX
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_XLX
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_XXD
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_UXD
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_XLD
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(lowDia_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_ULX
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_XLD_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_XLD_ULD
#endif
        type(lowDia_type)   , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! ULX

    interface getSubUnion

    pure elemental module function getSubUnion_ULX_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_XXX
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_UXX
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_XLX
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_XXD
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_UXD
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_XLD
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_ULX
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLow_type)                   :: subUnion
    end function

    pure elemental module function getSubUnion_ULX_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULX_ULD
#endif
        type(uppLow_type)   , intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

    ! ULD

    interface getSubUnion

    pure elemental module function getSubUnion_ULD_XXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_XXX
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(nothing_type)  , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_UXX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_UXX
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(upp_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_XLX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_XLX
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(low_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_XXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_XXD
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(dia_type)      , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_UXD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_UXD
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(uppDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_XLD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_XLD
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(lowDia_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_ULX(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_ULX
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(uppLow_type)   , intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    pure elemental module function getSubUnion_ULD_ULD(sub1, sub2) result(subUnion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubUnion_ULD_ULD
#endif
        type(uppLowDia_type), intent(in)    :: sub1
        type(uppLowDia_type), intent(in)    :: sub2
        type(uppLowDia_type)                :: subUnion
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixSubset ! LCOV_EXCL_LINE