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
!>  Such procedures frequently have multiple implementations determined by the specific control flow used in the algorithm.<br>
!>
!>  \details
!>  In computer science, control flow (or flow of control) is the order in which individual statements,
!>  instructions or function calls of an imperative program are executed or evaluated.<br>
!>  There are several types of programming control flows including,
!>  <ol>
!>      <li>    **Sequence control flow** refers to performing evaluations or statement executions one after another, all serially and in the specified order.
!>      <li>    **Selection control flow** refers to performing evaluations or statement executions conditionally using *if*, *unless*, *switch*, *case*, or other similar programming constructs.
!>      <li>    **Iteration control flow** (or **looping**) refers to performing evaluations or statement executions repeatedly using *for*, *while*, *repeat*, *until*, or other similar programming constructs.
!>      <li>    **Procedural Abstraction control flow** (or subroutine call) refers to performing evaluations or statement executions by calling a sequence of procedures.
!>      <li>    **Recursion control flow** refers to performing evaluations or statement executions repeatedly using procedure self calls until a condition is met.
!>      <li>    **Nondeterminacy control flow** refers to performing evaluations or statement executions randomly from a set of predefined possible flows in the code.
!>      <li>    **Concurrency control flow** refers to performing evaluations or statement executions concurrently with any dependence of the individual tasks on each other.
!>  </ol>
!>
!>  This module contains an incomplete list of derived types that help resolve the various implementations of the same algorithms using different control flows at compile time.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_control

    use pm_kind, only: SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_control"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of iterations (e.g., looping, recursion, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{control_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: control_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request **sequence control flow** within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [sequence](@ref pm_control::sequence)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{sequence_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(control_type) :: sequence_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [sequence_type](@ref pm_control::sequence_type) that is exclusively used
    !>  to request **sequence control flow** within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{sequence}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(sequence_type), parameter :: sequence = sequence_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: sequence
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request **recursive procedure** interfaces within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [recursion](@ref pm_control::recursion)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{recursion_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(control_type) :: recursion_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [recursion_type](@ref pm_control::recursion_type) that is exclusively used
    !>  to request **recursive procedure** interfaces within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{recursion}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(recursion_type), parameter :: recursion = recursion_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: recursion
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request **looping procedure** interfaces within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [iteration](@ref pm_control::iteration)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{iteration_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(control_type) :: iteration_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [iteration_type](@ref pm_control::iteration_type) that is exclusively used
    !>  to request **looping procedure** interfaces within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{iteration}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(iteration_type), parameter :: iteration = iteration_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: iteration
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request **selection control flow** within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [selection](@ref pm_control::selection)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{selection_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(control_type) :: selection_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [selection_type](@ref pm_control::selection_type) that is exclusively used
    !>  to request **selection control flow** within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sequence](@ref pm_control::sequence)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [selection](@ref pm_control::selection)<br>
    !>  [sequence_type](@ref pm_control::sequence_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>  [selection_type](@ref pm_control::selection_type)<br>
    !>  [control_type](@ref pm_control::control_type)<br>
    !>
    !>  \final{selection}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(selection_type), parameter :: selection = selection_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: selection
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_control ! LCOV_EXCL_LINE