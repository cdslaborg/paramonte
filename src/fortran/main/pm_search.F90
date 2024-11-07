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
!>  resolution of procedures within the generic interfaces of the ParaMonte library for various search operations.<br>
!>  Such procedures frequently need to work with a specific searching method within the algorithm implementation.<br>
!>
!>  \details
!>  In computer science, a search algorithm is an algorithm designed to solve a search problem.<br>
!>  Search algorithms work to retrieve information stored within particular data structure,
!>  or calculated in the search space of a problem domain, with either discrete or continuous values.<br>
!>
!>  Search algorithms can be classified based on their mechanism of searching into three types of algorithms:<br>
!>  Search algorithms are commonly evaluated by their computational complexity, or maximum theoretical run time.<br>
!>  Comparison search algorithms improve on linear searching by successively eliminating records based on comparisons 
!>  of the keys until the target record is found, and can be applied on data structures with a defined order.<br>
!>  Digital search algorithms work based on the properties of digits in data structures by using numerical keys.<br>
!>  Finally, hashing directly maps keys to records based on a hash function.<br>
!>  <ol>
!>      <li>    **Linear search algorithms check every record for the one associated with a target key in a linear fashion.<br>
!>      <li>    **Binary**, or **half-interval** search algorithms repeatedly target the center of the search structure and divide the search space in half.<br>
!>              <ol>
!>                  <li>    Binary search algorithms have a maximum complexity of \f$\mathcal{O}(\log n)\f$, or logarithmic time.<br>
!>                  <li>    In simple terms, the maximum number of operations needed to find the search target is a logarithmic function of the size of the search space.<br>
!>              </ol>
!>      <li>    **Hashing**, search algorithms directly map keys to records based on a [hash function](https://en.wikipedia.org/wiki/Hash_function).<br>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_search

    use pm_kind, only: SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_search"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different search algorithms (linear, binary, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [linear](@ref pm_search::linear)<br>
    !>  [binary](@ref pm_search::binary)<br>
    !>  [linear_type](@ref pm_search::linear_type)<br>
    !>  [binary_type](@ref pm_search::binary_type)<br>
    !>  [search_type](@ref pm_search::search_type)<br>
    !>
    !>  \final{search_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: search_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request linear search algorithm within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [linear](@ref pm_search::linear)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [linear](@ref pm_search::linear)<br>
    !>  [binary](@ref pm_search::binary)<br>
    !>  [linear_type](@ref pm_search::linear_type)<br>
    !>  [binary_type](@ref pm_search::binary_type)<br>
    !>  [search_type](@ref pm_search::search_type)<br>
    !>
    !>  \final{linear_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(search_type) :: linear_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [linear_type](@ref pm_search::linear_type) that is exclusively used
    !>  to request linear search algorithm within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [linear](@ref pm_search::linear)<br>
    !>  [binary](@ref pm_search::binary)<br>
    !>  [linear_type](@ref pm_search::linear_type)<br>
    !>  [binary_type](@ref pm_search::binary_type)<br>
    !>  [search_type](@ref pm_search::search_type)<br>
    !>
    !>  \final{linear}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(linear_type), parameter :: linear = linear_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: linear
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request binary search algorithm within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [binary](@ref pm_search::binary)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [linear](@ref pm_search::linear)<br>
    !>  [binary](@ref pm_search::binary)<br>
    !>  [linear_type](@ref pm_search::linear_type)<br>
    !>  [binary_type](@ref pm_search::binary_type)<br>
    !>  [search_type](@ref pm_search::search_type)<br>
    !>
    !>  \final{binary_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(search_type) :: binary_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [binary_type](@ref pm_search::binary_type) that is exclusively used
    !>  to request binary search algorithm within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [linear](@ref pm_search::linear)<br>
    !>  [binary](@ref pm_search::binary)<br>
    !>  [linear_type](@ref pm_search::linear_type)<br>
    !>  [binary_type](@ref pm_search::binary_type)<br>
    !>  [search_type](@ref pm_search::search_type)<br>
    !>
    !>  \final{binary}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(binary_type), parameter :: binary = binary_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: binary
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_search ! LCOV_EXCL_LINE