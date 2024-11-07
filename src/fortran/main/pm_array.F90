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
!>  resolution of procedures within the generic interfaces of the ParaMonte library for various array operations.<br>
!>
!>  \details
!>  Such procedures frequently need to work in a specific direction of some of their input array arguments.<br>
!>  While the English words, **inverse**, **reverse**, **backward** are sometimes used synonymously,
!>  there are subtle differences between them that make different in practice.<br>
!>  <ol>
!>      <li>    The words **forward/backward** are adjectives that frequently refer to a **direction** along which an action is performed.<br>
!>              **Example**: Walking backward...<br>
!>      <li>    The words **inverse/reverse** are adjective/nouns/verbs that frequently refer to an **action** that is performed.<br>
!>              <ol>
!>                  <li>    The word **reverse** frequently refers to performing an **action** in the **opposite order**.<br>
!>                          **Example**: Walking in reverse...<br>
!>                  <li>    The word **inverse** frequently refers to performing an **action** of **inversion** on something (e.g., turning something upside down).<br>
!>                          **Example**: The inverse of `5` is `1/5`.<br>
!>              </ol>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_array

    use pm_kind, only: SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_array"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of operation (reverse, inverse, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{action_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: action_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request **no action** on a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [nothing](@ref pm_array::nothing)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{nothing_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(action_type) :: nothing_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [nothing_type](@ref pm_array::nothing_type) that is exclusively used
    !>  to request **no action** on a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{nothing}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(nothing_type), parameter :: nothing = nothing_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: nothing
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [reverse](@ref pm_array::reverse)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{reverse_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(action_type) :: reverse_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [reverse_type](@ref pm_array::reverse_type) that is exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{reverse}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(reverse_type), parameter :: reverse = reverse_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: reverse
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [inverse](@ref pm_array::inverse)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{inverse_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(action_type) :: inverse_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [inverse_type](@ref pm_array::inverse_type) that is exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{inverse}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(inverse_type), parameter :: inverse = inverse_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: inverse
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of operation (reverse, inverse, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{direction_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: direction_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [forward](@ref pm_array::forward)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{forward_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(direction_type) :: forward_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [forward_type](@ref pm_array::forward_type) that is exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{forward}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(forward_type), parameter :: forward = forward_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: forward
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [backward](@ref pm_array::backward)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{backward_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(direction_type) :: backward_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [backward_type](@ref pm_array::backward_type) that is exclusively used
    !>  to request reversal of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{backward}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(backward_type), parameter :: backward = backward_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: backward
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of array sides (left, right, leftRight, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{side_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: side_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request left side of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [left](@ref pm_array::left)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{left_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(side_type) :: left_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [left_type](@ref pm_array::left_type) that is exclusively used
    !>  to request the left side of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{left}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(left_type), parameter :: left = left_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: left
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request the right side of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [right](@ref pm_array::right)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{right_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(side_type) :: right_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [right_type](@ref pm_array::right_type) that is exclusively used
    !>  to request the right side of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{right}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(right_type), parameter :: right = right_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: right
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request left-and-right sides of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [leftRight](@ref pm_array::leftRight)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{leftRight_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(side_type) :: leftRight_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [leftRight_type](@ref pm_array::leftRight_type) that is exclusively used
    !>  to request the left-and-right sides of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [left](@ref pm_array::left)<br>
    !>  [right](@ref pm_array::right)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [reverse](@ref pm_array::reverse)<br>
    !>  [inverse](@ref pm_array::inverse)<br>
    !>  [forward](@ref pm_array::forward)<br>
    !>  [backward](@ref pm_array::backward)<br>
    !>  [leftRight](@ref pm_array::leftRight)<br>
    !>  [left_type](@ref pm_array::left_type)<br>
    !>  [right_type](@ref pm_array::right_type)<br>
    !>  [nothing_type](@ref pm_array::nothing_type)<br>
    !>  [reverse_type](@ref pm_array::reverse_type)<br>
    !>  [inverse_type](@ref pm_array::inverse_type)<br>
    !>  [forward_type](@ref pm_array::forward_type)<br>
    !>  [backward_type](@ref pm_array::backward_type)<br>
    !>  [leftRight_type](@ref pm_array::leftRight_type)<br>
    !>  [direction_type](@ref pm_array::direction_type)<br>
    !>  [action_type](@ref pm_array::action_type)<br>
    !>  [side_type](@ref pm_array::side_type)<br>
    !>
    !>  \final{leftRight}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(leftRight_type), parameter :: leftRight = leftRight_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: leftRight
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different sequence border patterns (e.g., adjacent (contiguous), discrete, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [adjacent](@ref pm_array::adjacent)<br>
    !>  [discrete](@ref pm_array::discrete)<br>
    !>  [adjacent_type](@ref pm_array::adjacent_type)<br>
    !>  [discrete_type](@ref pm_array::discrete_type)<br>
    !>  [border_type](@ref pm_array::border_type)<br>
    !>
    !>  \final{border_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: border_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the adjacent sequence border within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [adjacent](@ref pm_array::adjacent)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [adjacent](@ref pm_array::adjacent)<br>
    !>  [discrete](@ref pm_array::discrete)<br>
    !>  [adjacent_type](@ref pm_array::adjacent_type)<br>
    !>  [discrete_type](@ref pm_array::discrete_type)<br>
    !>  [border_type](@ref pm_array::border_type)<br>
    !>
    !>  \final{adjacent_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(border_type) :: adjacent_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [adjacent_type](@ref pm_array::adjacent_type) that is exclusively used
    !>  to signify the adjacent sequence border within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [adjacent](@ref pm_array::adjacent)<br>
    !>  [discrete](@ref pm_array::discrete)<br>
    !>  [adjacent_type](@ref pm_array::adjacent_type)<br>
    !>  [discrete_type](@ref pm_array::discrete_type)<br>
    !>  [border_type](@ref pm_array::border_type)<br>
    !>
    !>  \final{adjacent}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(adjacent_type), parameter :: adjacent = adjacent_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: adjacent
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the discrete sequence border within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [discrete](@ref pm_array::discrete)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [adjacent](@ref pm_array::adjacent)<br>
    !>  [discrete](@ref pm_array::discrete)<br>
    !>  [adjacent_type](@ref pm_array::adjacent_type)<br>
    !>  [discrete_type](@ref pm_array::discrete_type)<br>
    !>  [border_type](@ref pm_array::border_type)<br>
    !>
    !>  \final{discrete_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(border_type) :: discrete_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [discrete_type](@ref pm_array::discrete_type) that is exclusively used
    !>  to signify the discrete sequence border within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [adjacent](@ref pm_array::adjacent)<br>
    !>  [discrete](@ref pm_array::discrete)<br>
    !>  [adjacent_type](@ref pm_array::adjacent_type)<br>
    !>  [discrete_type](@ref pm_array::discrete_type)<br>
    !>  [border_type](@ref pm_array::border_type)<br>
    !>
    !>  \final{discrete}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(discrete_type), parameter :: discrete = discrete_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: discrete
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different sequence order patterns (e.g., sorted, ascending, descending, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{order_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: order_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the sorted sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [sorted](@ref pm_array::sorted)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{sorted_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(order_type) :: sorted_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [sorted_type](@ref pm_array::sorted_type) that is exclusively used
    !>  to signify the sorted sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{sorted}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(sorted_type), parameter :: sorted = sorted_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: sorted
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the monotonic sequence order (e.g., **strictly** ascending, descending, equal, ...) within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [monotonic](@ref pm_array::monotonic)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{monotonic_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(order_type) :: monotonic_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [monotonic_type](@ref pm_array::monotonic_type) that is exclusively used
    !>  to signify the monotonic sequence order (e.g., **strictly** ascending, descending, equal, ...) within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{monotonic}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(monotonic_type), parameter :: monotonic = monotonic_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: monotonic
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the ascending sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [ascending](@ref pm_array::ascending)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{ascending_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(sorted_type) :: ascending_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [ascending_type](@ref pm_array::ascending_type) that is exclusively used
    !>  to signify the ascending sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [monotonic](@ref pm_array::monotonic)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [monotonic_type](@ref pm_array::monotonic_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{ascending}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(ascending_type), parameter :: ascending = ascending_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ascending
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the descending sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [descending](@ref pm_array::descending)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{descending_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(sorted_type) :: descending_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [descending_type](@ref pm_array::descending_type) that is exclusively used
    !>  to signify the descending sequence order within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [sorted](@ref pm_array::sorted)<br>
    !>  [ascending](@ref pm_array::ascending)<br>
    !>  [descending](@ref pm_array::descending)<br>
    !>  [sorted_type](@ref pm_array::sorted_type)<br>
    !>  [ascending_type](@ref pm_array::ascending_type)<br>
    !>  [descending_type](@ref pm_array::descending_type)<br>
    !>  [order_type](@ref pm_array::order_type)<br>
    !>
    !>  \final{descending}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(descending_type), parameter :: descending = descending_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: descending
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify that an array of arbitrary rank has the `allocatable` attribute.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [allocatable](@ref pm_array::allocatable)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [allocatable](@ref pm_array::allocatable)<br>
    !>  [allocatable_type](@ref pm_array::allocatable_type)<br>
    !>
    !>  \final{allocatable_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: allocatable_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [allocatable_type](@ref pm_array::allocatable_type) that is exclusively used
    !>  to signify that an array of arbitrary rank has the `allocatable` attribute.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.
    !>
    !>  \see
    !>  [allocatable](@ref pm_array::allocatable)<br>
    !>  [allocatable_type](@ref pm_array::allocatable_type)<br>
    !>
    !>  \final{allocatable}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(allocatable_type), parameter :: allocatable = allocatable_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: allocatable
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_array ! LCOV_EXCL_LINE