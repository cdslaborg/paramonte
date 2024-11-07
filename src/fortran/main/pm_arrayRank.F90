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
!>  This module contains procedures and generic interfaces for obtaining **various rankings of elements** of arrays of various types.
!>
!>  \details
!>  Depending on the applications, the rank of the elements of an array can be defined in different ways:
!>
!>  <ol>
!>      <li>    <b>Ordinal ranking (`1234`) ranking:</b> [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal) or [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
!>              This kind of ranking of values is widely known as ordinal (`1234`) ranking.<br>
!>              In ordinal ranking, all items receive distinct ordinal numbers, including items that compare equal.<br>
!>              The assignment of distinct ordinal numbers to items that compare equal can be done at random, or arbitrarily,
!>              but it is generally preferable to use a system that is arbitrary but consistent,
!>              as this gives stable results if the ranking is done multiple times.<br>
!>              In computer data processing, ordinal ranking is also referred to as <b>row numbering</b>.
!>              That is, if `A < B == C < D`, then the sequence `ABCD` has the <b>ordinal ranking</b> `1234`.<br>
!>
!>      <li>    <b>Standard competition (`1224`) ranking:</b> [getRankStandard](@ref pm_arrayRank::getRankStandard) or [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
!>              This kind of ranking of values is widely known as Standard Competition (`1224`) ranking.<br>
!>              In Standard Competition ranking, items that compare equal receive the same ranking number,
!>              and then a gap is left in the ranking numbers. The number of ranking numbers that are left out
!>              in this gap is one less than the number of items that compared equal.<br>
!>              Equivalently, the ranking number of each item is `1` plus the number of items ranked above it.<br>
!>              This ranking strategy is frequently adopted for competitions, as it means that if two (or more) competitors
!>              tie for a position in the ranking, and the position of all those ranked below them is unaffected
!>              (i.e., a competitor only comes second if exactly one person scores better than them,
!>              third if exactly two people score better than them, fourth if exactly three people score better than them, etc.).<br>
!>              Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
!>              *B* gets ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `4` (*fourth*).<br>
!>              That is, if `A < B == C < D`, then the sequence `ABCD` has the Standard Competition ranking `1224`.<br>
!>
!>      <li>    <b>Modified competition (`1334`) ranking:</b> [getRankModified](@ref pm_arrayRank::getRankModified) or [setRankModified](@ref pm_arrayRank::setRankModified)<br>
!>              This kind of ranking of values is widely known as Modified Competition (`1334`) ranking.<br>
!>              Sometimes, competition ranking is done by leaving the gaps in the ranking numbers before the sets of equal-ranking items
!>              (rather than after them as in Standard Competition ranking).<br>
!>              The number of ranking numbers that are left out in this gap
!>              remains one less than the number of items that compared equal.<br>
!>              Equivalently, the ranking number of each item is equal to the number of items ranked equal to it or above it.<br>
!>              This ranking ensures that a competitor only comes second if they score higher than all but one of their opponents,
!>              third if they score higher than all but two of their opponents, etc.<br>
!>              Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked head of *D*, then *A* gets ranking
!>              number `1` (*first*), *B* gets ranking number `3` (*joint third*), *C* also gets ranking number `3` (*joint third*)
!>              and *D* gets ranking number `4` (*fourth*). In this case, nobody would get ranking number `2` (*second*) (left as a gap).<br>
!>              That is, if `A < B == C < D`, then the sequence `ABCD` has the Modified Competition ranking `1334`.<br>
!>
!>      <li>    <b>Dense (`1223`) ranking:</b> [getRankDense](@ref pm_arrayRank::getRankDense) or [setRankDense](@ref pm_arrayRank::setRankDense)<br>
!>              This kind of ranking of values is widely known as dense (`1223`) ranking.<br>
!>              In Dense ranking, items that compare equally receive the same ranking number, and the next items receive the immediately following ranking number.<br>
!>              Equivalently, the ranking number of each item is `1` plus the number of items ranked above it that are distinct with respect to the ranking order.<br>
!>              Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number 1 (*first*), *B* gets
!>              ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `3` (*Third*).<br>
!>              That is, if `A < B == C < D`, then the sequence `ABCD` has the Dense ranking `1223`.<br>
!>              Dense ranking effective factorizes the array into classes of unique values.<br>
!>              Therefore, the Dense rank of each element of the array is simply its class <b>level</b>.
!>
!>      <li>    <b>Fractional (`1 2.5 2.5 4`) ranking:</b> [getRankFractional](@ref pm_arrayRank::getRankFractional) or [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
!>              This kind of ranking of values is widely known as fractional (`1 2.5 2.5 4`) ranking.<br>
!>              In Fractional ranking, items that compare equal receive the same ranking number, which is the mean of what they would have under ordinal rankings;<br>
!>              Equivalently, the ranking number of 1 plus the number of items ranked above it plus half the number of items equal to it.<br>
!>              This strategy has the property that the sum of the ranking numbers is the same as under ordinal ranking.<br>
!>              For this reason, it is used in computing Borda counts and ranking statistics (e.g., Spearman Correlation).<br>
!>              Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
!>              *B* and *C* each get ranking number `2.5` (average of *joint second/third*) and *D* gets ranking number `4` (*fourth*).<br>
!>              That is, if `A < B == C < D`, then the sequence `ABCD` has the Fractional ranking `1223`.<br>
!>              <b>Example:</b><br>
!>              Suppose the data set is `1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0, 5.0`.<br>
!>              The ordinal ranks are `1, 2, 3, 4, 5, 6, 7, 8, 9`.<br>
!>              For `v = 1.0`, the Fractional rank is the average of the ordinal ranks: `(1 + 2) / 2 = 1.5`.<br>
!>              In a similar manner, for `v = 5.0`, the Fractional rank is `(7 + 8 + 9) / 3 = 8.0`.<br>
!>              Thus the Fractional ranks are: `1.5, 1.5, 3.0, 4.5, 4.5, 6.0, 8.0, 8.0, 8.0`<br>
!>  </ol>
!>
!>  \warning
!>  The support for ranking of string containers is disabled when the library is built with
!>  the GNU Fortran compiler because of the lack of support for Parameterized Derived Types (PDTs) in \gfortran.
!>
!>  \note
!>  Obtaining the **ordinal ranking** of an array is very similar to obtaining the sorted indices of the array.<br>
!>  For more information see this [article](https://en.wikipedia.org/wiki/Ranking).<br>
!>
!>  \todo
!>  \plow
!>  The relevant benchmarks comparing the functional and subroutine interfaces should be added here.<br>
!>
!>  \todo
!>  \pvhigh
!>  Support for ranking of arrays of PDTs must be enabled again as soon as \gfortran supports PDTs.<br>
!>
!>  \todo
!>  \phigh
!>  An optional argument `sorted` must be added to all interfaces within this module to
!>  allow fast computation of the rankings of previously sorted arrays without redundant resorting.<br>
!>
!>  \test
!>  [test_pm_arrayRank](@ref test_pm_arrayRank)<br>
!>
!>  \todo
!>  \pmed
!>  The generic interfaces of this module must be extended to support rankings of matrices along a specified dimension.<br>
!>  Such a pattern occurs, for example, in computing the [Spearman rank correlation matrix](@ref pm_sampleCor).<br>
!>
!>  \todo
!>  \pmed
!>  The current implementation of the generic ranking interfaces of this module are separated from each other.<br>
!>  However, it may be preferrable to merge all generic interfaces into single interface bindings `getRank()` and `setRank()`.<br>
!>  Consequently, an extra argument of class [rank_type](@ref pm_arrayRank::rank_type) must be added to all procedure interfaces to make them distinguishable.<br>
!>  For now, the procedures for different ranking methods were are under separate generic interface names, because of the complexity in the merging of the
!>  fractional ranking procedures (which output `real` ranks of default kind \RK) with the rest (which output `integer` ranks of default kind \IK).<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!   \copydetails pm_arrayRank

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayRank

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayRank"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of array ranking (dense, fractional, modified, ordinal, standard, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{rank_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: rank_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request dense ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [dense](@ref pm_arrayRank::dense)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{dense_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rank_type) :: dense_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [dense_type](@ref pm_arrayRank::dense_type) that is exclusively used
    !>  to request the dense ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{dense}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(dense_type), parameter :: dense = dense_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request the ordinal ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [ordinal](@ref pm_arrayRank::ordinal)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{ordinal_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rank_type) :: ordinal_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [ordinal_type](@ref pm_arrayRank::ordinal_type) that is exclusively used
    !>  to request the ordinal ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{ordinal}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(ordinal_type), parameter :: ordinal = ordinal_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request the modified ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [modified](@ref pm_arrayRank::modified)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{modified_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rank_type) :: modified_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [modified_type](@ref pm_arrayRank::modified_type) that is exclusively used
    !>  to request the modified ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{modified}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(modified_type), parameter :: modified = modified_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request the standard ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [standard](@ref pm_arrayRank::standard)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{standard_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rank_type) :: standard_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [standard_type](@ref pm_arrayRank::standard_type) that is exclusively used
    !>  to request the standard ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{standard}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(standard_type), parameter :: standard = standard_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request the fractional ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [fractional](@ref pm_arrayRank::fractional)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{fractional_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(rank_type) :: fractional_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [fractional_type](@ref pm_arrayRank::fractional_type) that is exclusively used
    !>  to request the fractional ranking of a given array within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [dense](@ref pm_arrayRank::dense)<br>
    !>  [ordinal](@ref pm_arrayRank::ordinal)<br>
    !>  [modified](@ref pm_arrayRank::modified)<br>
    !>  [standard](@ref pm_arrayRank::standard)<br>
    !>  [fractional](@ref pm_arrayRank::fractional)<br>
    !>  [dense_type](@ref pm_arrayRank::dense_type)<br>
    !>  [ordinal_type](@ref pm_arrayRank::ordinal_type)<br>
    !>  [modified_type](@ref pm_arrayRank::modified_type)<br>
    !>  [standard_type](@ref pm_arrayRank::standard_type)<br>
    !>  [fractional_type](@ref pm_arrayRank::fractional_type)<br>
    !>  [rank_type](@ref pm_arrayRank::rank_type)<br>
    !>
    !>  \final{fractional}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(fractional_type), parameter :: fractional = fractional_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Dense rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>dense (`1223`) ranking</b>.<br>
    !>  In Dense ranking, items that compare equally receive the same ranking number, and the next items receive the immediately following ranking number.<br>
    !>  Equivalently, the ranking number of each item is `1` plus the number of items ranked above it that are distinct with respect to the ranking order.<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number 1 (*first*), *B* gets
    !>  ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `3` (*Third*).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Dense ranking** `1223`.<br>
    !>  Dense ranking effective factorizes the array into classes of unique values.<br>
    !>  Therefore, the Dense rank of each element of the array is simply its class **level**.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `rank(1:size(array)`    :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` matches that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Dense rank of the `i`th element of `array`.**
    !>
    !>  \interface{getRankDense}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: getRankDense
    !>
    !>      rank(1:size(array)) = getRankDense(array)
    !>      rank(1:size(array)) = getRankDense(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{getRankDense}
    !>  \include{lineno} example/pm_arrayRank/getRankDense/main.F90
    !>  \compilef{getRankDense}
    !>  \output{getRankDense}
    !>  \include{lineno} example/pm_arrayRank/getRankDense/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{getRankDense}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRankDense

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankDenseDefCom_D0_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseDefCom_D0_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseDefCom_D0_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseDefCom_D0_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseDefCom_D0_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankDenseDefCom_D1_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseDefCom_D1_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseDefCom_D1_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseDefCom_D1_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseDefCom_D1_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankDenseDefCom_D1_IK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankDenseDefCom_D1_IK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankDenseDefCom_D1_IK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankDenseDefCom_D1_IK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankDenseDefCom_D1_IK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankDenseDefCom_D1_LK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankDenseDefCom_D1_LK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankDenseDefCom_D1_LK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankDenseDefCom_D1_LK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankDenseDefCom_D1_LK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankDenseDefCom_D1_CK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankDenseDefCom_D1_CK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankDenseDefCom_D1_CK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankDenseDefCom_D1_CK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankDenseDefCom_D1_CK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankDenseDefCom_D1_RK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankDenseDefCom_D1_RK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankDenseDefCom_D1_RK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankDenseDefCom_D1_RK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankDenseDefCom_D1_RK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankDenseDefCom_D1_PSSK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseDefCom_D1_PSSK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseDefCom_D1_PSSK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseDefCom_D1_PSSK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseDefCom_D1_PSSK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankDenseDefCom_D1_BSSK(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankDenseCusCom_D0_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseCusCom_D0_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseCusCom_D0_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseCusCom_D0_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseCusCom_D0_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankDenseCusCom_D1_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseCusCom_D1_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseCusCom_D1_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseCusCom_D1_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseCusCom_D1_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankDenseCusCom_D1_IK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankDenseCusCom_D1_IK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankDenseCusCom_D1_IK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankDenseCusCom_D1_IK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankDenseCusCom_D1_IK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankDenseCusCom_D1_LK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankDenseCusCom_D1_LK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankDenseCusCom_D1_LK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankDenseCusCom_D1_LK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankDenseCusCom_D1_LK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankDenseCusCom_D1_CK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankDenseCusCom_D1_CK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankDenseCusCom_D1_CK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankDenseCusCom_D1_CK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankDenseCusCom_D1_CK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankDenseCusCom_D1_RK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankDenseCusCom_D1_RK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankDenseCusCom_D1_RK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankDenseCusCom_D1_RK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankDenseCusCom_D1_RK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankDenseCusCom_D1_PSSK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankDenseCusCom_D1_PSSK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankDenseCusCom_D1_PSSK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankDenseCusCom_D1_PSSK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankDenseCusCom_D1_PSSK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
        procedure(logical(LK))                                      :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankDenseCusCom_D1_BSSK(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankDenseCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **Dense rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>dense (`1223`) ranking</b>.<br>
    !>  In Dense ranking, items that compare equally receive the same ranking number, and the next items receive the immediately following ranking number.<br>
    !>  Equivalently, the ranking number of each item is `1` plus the number of items ranked above it that are distinct with respect to the ranking order.<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number 1 (*first*), *B* gets
    !>  ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `3` (*Third*).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Dense ranking** `1223`.<br>
    !>  Dense ranking effective factorizes the array into classes of unique values. Therefore, the Dense rank of each element of the array
    !>  is simply its class **level**.
    !>
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Dense rank of the `i`th element of `array`.**
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setRankDense}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: setRankDense
    !>
    !>      call setRankDense(rank, array)
    !>      call setRankDense(rank, array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Dense ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{setRankDense}
    !>  \include{lineno} example/pm_arrayRank/setRankDense/main.F90
    !>  \compilef{setRankDense}
    !>  \output{setRankDense}
    !>  \include{lineno} example/pm_arrayRank/setRankDense/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{setRankDense}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRankDense

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D0_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D0_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D0_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D0_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D0_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_IK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_IK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_IK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_IK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_IK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_LK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_LK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_LK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_LK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_LK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_CK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_CK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_CK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_CK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_CK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_RK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_RK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_RK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_RK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_RK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_PSSK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_PSSK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_PSSK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_PSSK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankDenseDefCom_D1_PSSK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRankDenseDefCom_D1_BSSK(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankDenseCusCom_D0_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankDenseCusCom_D0_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankDenseCusCom_D0_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankDenseCusCom_D0_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankDenseCusCom_D0_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankDenseCusCom_D1_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankDenseCusCom_D1_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankDenseCusCom_D1_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankDenseCusCom_D1_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankDenseCusCom_D1_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRankDenseCusCom_D1_IK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRankDenseCusCom_D1_IK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRankDenseCusCom_D1_IK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRankDenseCusCom_D1_IK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRankDenseCusCom_D1_IK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRankDenseCusCom_D1_LK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRankDenseCusCom_D1_LK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRankDenseCusCom_D1_LK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRankDenseCusCom_D1_LK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRankDenseCusCom_D1_LK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRankDenseCusCom_D1_CK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRankDenseCusCom_D1_CK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRankDenseCusCom_D1_CK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRankDenseCusCom_D1_CK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRankDenseCusCom_D1_CK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRankDenseCusCom_D1_RK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRankDenseCusCom_D1_RK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRankDenseCusCom_D1_RK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRankDenseCusCom_D1_RK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRankDenseCusCom_D1_RK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setRankDenseCusCom_D1_PSSK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankDenseCusCom_D1_PSSK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankDenseCusCom_D1_PSSK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankDenseCusCom_D1_PSSK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankDenseCusCom_D1_PSSK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setRankDenseCusCom_D1_BSSK(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankDenseCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Fractional rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>fractional (`1 2.5 2.5 4`) ranking</b>.<br>
    !>  In Fractional ranking, items that compare equal receive the same ranking number, which is the mean of what they would have under ordinal
    !>  rankings; equivalently, the ranking number of 1 plus the number of items ranked above it plus half the number of items equal to it.<br>
    !>  This strategy has the property that the sum of the ranking numbers is the same as under ordinal ranking.<br>
    !>  For this reason, it is used in computing Borda counts and ranking statistics (e.g., Spearman Correlation).<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
    !>  *B* and *C* each get ranking number `2.5` (average of *joint second/third*) and *D* gets ranking number `4` (*fourth*).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Fractional ranking** `1223`.<br>
    !>  **Example:**<br>
    !>  Suppose the data set is `1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0, 5.0`.<br>
    !>  The ordinal ranks are `1, 2, 3, 4, 5, 6, 7, 8, 9`.<br>
    !>  For `v = 1.0`, the Fractional rank is the average of the ordinal ranks: `(1 + 2) / 2 = 1.5`.<br>
    !>  In a similar manner, for `v = 5.0`, the Fractional rank is `(7 + 8 + 9) / 3 = 8.0`.<br>
    !>  Thus the Fractional ranks are: `1.5, 1.5, 3.0, 4.5, 4.5, 6.0, 8.0, 8.0, 8.0`
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `rank(1:size(array)`    :   The output `contiguous` array of rank `1` of type `real` of default kind \RK containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Fractional rank of the `i`th element of `array`.**
    !>
    !>  \interface{getRankFractional}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: getRankFractional
    !>
    !>      rank(1:size(array)) = getRankFractional(array)
    !>      rank(1:size(array)) = getRankFractional(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{getRankFractional}
    !>  \include{lineno} example/pm_arrayRank/getRankFractional/main.F90
    !>  \compilef{getRankFractional}
    !>  \output{getRankFractional}
    !>  \include{lineno} example/pm_arrayRank/getRankFractional/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{getRankFractional}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRankFractional

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankFractionalDefCom_D0_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalDefCom_D0_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalDefCom_D0_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalDefCom_D0_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalDefCom_D0_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankFractionalDefCom_D1_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalDefCom_D1_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalDefCom_D1_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalDefCom_D1_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalDefCom_D1_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankFractionalDefCom_D1_IK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankFractionalDefCom_D1_IK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankFractionalDefCom_D1_IK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankFractionalDefCom_D1_IK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankFractionalDefCom_D1_IK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankFractionalDefCom_D1_LK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => RK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankFractionalDefCom_D1_LK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => RK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankFractionalDefCom_D1_LK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => RK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankFractionalDefCom_D1_LK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => RK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankFractionalDefCom_D1_LK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => RK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankFractionalDefCom_D1_CK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => RK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankFractionalDefCom_D1_CK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => RK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankFractionalDefCom_D1_CK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => RK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankFractionalDefCom_D1_CK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => RK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankFractionalDefCom_D1_CK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => RK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankFractionalDefCom_D1_RK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => RK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankFractionalDefCom_D1_RK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => RK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankFractionalDefCom_D1_RK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => RK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankFractionalDefCom_D1_RK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => RK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankFractionalDefCom_D1_RK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => RK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankFractionalDefCom_D1_PSSK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalDefCom_D1_PSSK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalDefCom_D1_PSSK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalDefCom_D1_PSSK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalDefCom_D1_PSSK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankFractionalDefCom_D1_BSSK(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankFractionalCusCom_D0_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalCusCom_D0_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalCusCom_D0_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalCusCom_D0_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalCusCom_D0_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankFractionalCusCom_D1_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalCusCom_D1_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalCusCom_D1_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalCusCom_D1_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalCusCom_D1_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankFractionalCusCom_D1_IK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankFractionalCusCom_D1_IK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankFractionalCusCom_D1_IK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankFractionalCusCom_D1_IK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankFractionalCusCom_D1_IK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankFractionalCusCom_D1_LK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => RK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankFractionalCusCom_D1_LK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => RK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankFractionalCusCom_D1_LK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => RK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankFractionalCusCom_D1_LK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => RK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankFractionalCusCom_D1_LK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => RK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankFractionalCusCom_D1_CK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => RK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankFractionalCusCom_D1_CK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => RK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankFractionalCusCom_D1_CK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => RK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankFractionalCusCom_D1_CK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => RK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankFractionalCusCom_D1_CK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => RK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankFractionalCusCom_D1_RK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => RK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankFractionalCusCom_D1_RK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => RK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankFractionalCusCom_D1_RK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => RK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankFractionalCusCom_D1_RK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => RK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankFractionalCusCom_D1_RK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => RK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankFractionalCusCom_D1_PSSK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankFractionalCusCom_D1_PSSK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankFractionalCusCom_D1_PSSK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankFractionalCusCom_D1_PSSK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankFractionalCusCom_D1_PSSK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                                                   :: rank(size(array, kind = IK))
        procedure(logical(LK))                                      :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankFractionalCusCom_D1_BSSK(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankFractionalCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        real(TKR)                                                   :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **Fractional rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>fractional (`1 2.5 2.5 4`) ranking</b>.<br>
    !>  In Fractional ranking, items that compare equal receive the same ranking number, which is the mean of what they would have under ordinal
    !>  rankings; equivalently, the ranking number of 1 plus the number of items ranked above it plus half the number of items equal to it.<br>
    !>  This strategy has the property that the sum of the ranking numbers is the same as under ordinal ranking.<br>
    !>  For this reason, it is used in computing Borda counts and ranking statistics (e.g., Spearman Correlation).<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
    !>  *B* and *C* each get ranking number `2.5` (average of *joint second/third*) and *D* gets ranking number `4` (*fourth*).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Fractional ranking** `1223`.<br>
    !>  **Example:**<br>
    !>  Suppose the data set is `1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0, 5.0`.<br>
    !>  The ordinal ranks are `1, 2, 3, 4, 5, 6, 7, 8, 9`.<br>
    !>  For `v = 1.0`, the Fractional rank is the average of the ordinal ranks: `(1 + 2) / 2 = 1.5`.<br>
    !>  In a similar manner, for `v = 5.0`, the Fractional rank is `(7 + 8 + 9) / 3 = 8.0`.<br>
    !>  Thus the Fractional ranks are: `1.5, 1.5, 3.0, 4.5, 4.5, 6.0, 8.0, 8.0, 8.0`
    !>
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `real` of default kind \RK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Fractional rank of the `i`th element of `array`.**
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired. In such cases,
    !>                              user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setRankFractional}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: setRankFractional
    !>
    !>      call setRankFractional(rank, array)
    !>      call setRankFractional(rank, array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Fractional ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{setRankFractional}
    !>  \include{lineno} example/pm_arrayRank/setRankFractional/main.F90
    !>  \compilef{setRankFractional}
    !>  \output{setRankFractional}
    !>  \include{lineno} example/pm_arrayRank/setRankFractional/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{setRankFractional}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRankFractional

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D0_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D0_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D0_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D0_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D0_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_IK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_IK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_IK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_IK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_IK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_LK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => RK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_LK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => RK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_LK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => RK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_LK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => RK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_LK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => RK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_CK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => RK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_CK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => RK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_CK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => RK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_CK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => RK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_CK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => RK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_RK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => RK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_RK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => RK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_RK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => RK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_RK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => RK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_RK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => RK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_PSSK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_PSSK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_PSSK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_PSSK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankFractionalDefCom_D1_PSSK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRankFractionalDefCom_D1_BSSK(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankFractionalCusCom_D0_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankFractionalCusCom_D0_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankFractionalCusCom_D0_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankFractionalCusCom_D0_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankFractionalCusCom_D0_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_IK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_IK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_IK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_IK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_IK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_LK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => RK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_LK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => RK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_LK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => RK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_LK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => RK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_LK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => RK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_CK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => RK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_CK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => RK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_CK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => RK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_CK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => RK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_CK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => RK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_RK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => RK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_RK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => RK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_RK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => RK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_RK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => RK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_RK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => RK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setRankFractionalCusCom_D1_PSSK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankFractionalCusCom_D1_PSSK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankFractionalCusCom_D1_PSSK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankFractionalCusCom_D1_PSSK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankFractionalCusCom_D1_PSSK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setRankFractionalCusCom_D1_BSSK(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankFractionalCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        real(TKR)                   , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Modified rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>Modified Competition (`1334`) ranking</b>.<br>
    !>  Sometimes, competition ranking is done by leaving the gaps in the ranking numbers before the sets
    !>  of equal-ranking items (rather than after them as in Standard Competition ranking).<br>
    !>  The number of ranking numbers that are left out in this gap remains one less than the number of items that compared equal.<br>
    !>  Equivalently, the ranking number of each item is equal to the number of items ranked equal to it or above it.<br>
    !>  This ranking ensures that a competitor only comes second if they score higher than all but one of their opponents,
    !>  third if they score higher than all but two of their opponents, etc.<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked head of *D*, then *A* gets ranking
    !>  number `1` (*first*), *B* gets ranking number `3` (*joint third*), *C* also gets ranking number `3` (*joint third*)
    !>  and *D* gets ranking number `4` (*fourth*).<br>
    !>  In this case, nobody would get ranking number `2` (*second*) (left as a gap).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Modified Competition ranking** `1334`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Modified rank of the `i`th element of `array`.**
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `rank(1:size(array)`    :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` matches that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Modified rank of the `i`th element of `array`.**
    !>
    !>  \interface{getRankModified}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: getRankModified
    !>
    !>      rank(1:size(array)) = getRankModified(array)
    !>      rank(1:size(array)) = getRankModified(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Modified ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{getRankModified}
    !>  \include{lineno} example/pm_arrayRank/getRankModified/main.F90
    !>  \compilef{getRankModified}
    !>  \output{getRankModified}
    !>  \include{lineno} example/pm_arrayRank/getRankModified/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{getRankModified}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRankModified

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankModifiedDefCom_D0_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedDefCom_D0_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedDefCom_D0_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedDefCom_D0_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedDefCom_D0_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankModifiedDefCom_D1_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedDefCom_D1_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedDefCom_D1_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedDefCom_D1_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedDefCom_D1_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankModifiedDefCom_D1_IK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankModifiedDefCom_D1_IK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankModifiedDefCom_D1_IK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankModifiedDefCom_D1_IK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankModifiedDefCom_D1_IK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankModifiedDefCom_D1_LK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankModifiedDefCom_D1_LK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankModifiedDefCom_D1_LK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankModifiedDefCom_D1_LK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankModifiedDefCom_D1_LK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankModifiedDefCom_D1_CK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankModifiedDefCom_D1_CK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankModifiedDefCom_D1_CK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankModifiedDefCom_D1_CK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankModifiedDefCom_D1_CK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankModifiedDefCom_D1_RK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankModifiedDefCom_D1_RK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankModifiedDefCom_D1_RK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankModifiedDefCom_D1_RK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankModifiedDefCom_D1_RK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankModifiedDefCom_D1_PSSK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedDefCom_D1_PSSK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedDefCom_D1_PSSK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedDefCom_D1_PSSK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedDefCom_D1_PSSK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankModifiedDefCom_D1_BSSK(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankModifiedCusCom_D0_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedCusCom_D0_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedCusCom_D0_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedCusCom_D0_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedCusCom_D0_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankModifiedCusCom_D1_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedCusCom_D1_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedCusCom_D1_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedCusCom_D1_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedCusCom_D1_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankModifiedCusCom_D1_IK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankModifiedCusCom_D1_IK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankModifiedCusCom_D1_IK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankModifiedCusCom_D1_IK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankModifiedCusCom_D1_IK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankModifiedCusCom_D1_LK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankModifiedCusCom_D1_LK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankModifiedCusCom_D1_LK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankModifiedCusCom_D1_LK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankModifiedCusCom_D1_LK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankModifiedCusCom_D1_CK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankModifiedCusCom_D1_CK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankModifiedCusCom_D1_CK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankModifiedCusCom_D1_CK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankModifiedCusCom_D1_CK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankModifiedCusCom_D1_RK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankModifiedCusCom_D1_RK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankModifiedCusCom_D1_RK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankModifiedCusCom_D1_RK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankModifiedCusCom_D1_RK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankModifiedCusCom_D1_PSSK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankModifiedCusCom_D1_PSSK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankModifiedCusCom_D1_PSSK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankModifiedCusCom_D1_PSSK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankModifiedCusCom_D1_PSSK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
        procedure(logical(LK))                                      :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankModifiedCusCom_D1_BSSK(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankModifiedCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **Modified rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>Modified Competition (`1334`) ranking</b>.<br>
    !>  Sometimes, competition ranking is done by leaving the gaps in the ranking numbers before the sets of equal-ranking items
    !>  (rather than after them as in Standard Competition ranking).<br>
    !>  The number of ranking numbers that are left out in this gap remains one less than the number of items that compared equal.<br>
    !>  Equivalently, the ranking number of each item is equal to the number of items ranked equal to it or above it.<br>
    !>  This ranking ensures that a competitor only comes second if they score
    !>  higher than all but one of their opponents, third if they score higher than all but two of their opponents, etc.<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked head of *D*, then *A* gets ranking
    !>  number `1` (*first*), *B* gets ranking number `3` (*joint third*), *C* also gets ranking number `3` (*joint third*)
    !>  and *D* gets ranking number `4` (*fourth*). In this case, nobody would get ranking number `2` (*second*) (left as a gap).<br>
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Modified Competition ranking** `1334`.<br>
    !>
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Modified rank of the `i`th element of `array`.**
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setRankModified}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: setRankModified
    !>
    !>      call setRankModified(rank, array)
    !>      call setRankModified(rank, array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Modified ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{setRankModified}
    !>  \include{lineno} example/pm_arrayRank/setRankModified/main.F90
    !>  \compilef{setRankModified}
    !>  \output{setRankModified}
    !>  \include{lineno} example/pm_arrayRank/setRankModified/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{setRankModified}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRankModified

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D0_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D0_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D0_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D0_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D0_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_IK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_IK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_IK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_IK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_IK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_LK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_LK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_LK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_LK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_LK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_CK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_CK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_CK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_CK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_CK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_RK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_RK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_RK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_RK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_RK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_PSSK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_PSSK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_PSSK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_PSSK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankModifiedDefCom_D1_PSSK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRankModifiedDefCom_D1_BSSK(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankModifiedCusCom_D0_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankModifiedCusCom_D0_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankModifiedCusCom_D0_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankModifiedCusCom_D0_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankModifiedCusCom_D0_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_IK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_IK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_IK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_IK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_IK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_LK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_LK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_LK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_LK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_LK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_CK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_CK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_CK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_CK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_CK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_RK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_RK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_RK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_RK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_RK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setRankModifiedCusCom_D1_PSSK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankModifiedCusCom_D1_PSSK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankModifiedCusCom_D1_PSSK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankModifiedCusCom_D1_PSSK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankModifiedCusCom_D1_PSSK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setRankModifiedCusCom_D1_BSSK(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankModifiedCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **ordinal rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>ordinal (`1234`) ranking</b>.<br>
    !>  In ordinal ranking, all items receive distinct ordinal numbers, including items that compare equal.<br>
    !>  The assignment of distinct ordinal numbers to items that compare equal can be done at random, or arbitrarily,
    !>  but it is generally preferable to use a system that is arbitrary but consistent,
    !>  as this gives stable results if the ranking is done multiple times.<br>
    !>  In computer data processing, ordinal ranking is also referred to as **row numbering**.
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **ordinal ranking** `1234`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `rank(1:size(array)`    :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` matches that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the ordinal rank of the `i`th element of `array`.**
    !>
    !>  \interface{getRankOrdinal}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: getRankOrdinal
    !>
    !>      rank(1:size(array)) = getRankOrdinal(array)
    !>      rank(1:size(array)) = getRankOrdinal(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{getRankOrdinal}
    !>  \include{lineno} example/pm_arrayRank/getRankOrdinal/main.F90
    !>  \compilef{getRankOrdinal}
    !>  \output{getRankOrdinal}
    !>  \include{lineno} example/pm_arrayRank/getRankOrdinal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{getRankOrdinal}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRankOrdinal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankOrdinalDefCom_D0_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalDefCom_D0_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalDefCom_D0_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalDefCom_D0_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalDefCom_D0_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankOrdinalDefCom_D1_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalDefCom_D1_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalDefCom_D1_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalDefCom_D1_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalDefCom_D1_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankOrdinalDefCom_D1_IK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankOrdinalDefCom_D1_IK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankOrdinalDefCom_D1_IK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankOrdinalDefCom_D1_IK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankOrdinalDefCom_D1_IK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankOrdinalDefCom_D1_LK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankOrdinalDefCom_D1_LK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankOrdinalDefCom_D1_LK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankOrdinalDefCom_D1_LK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankOrdinalDefCom_D1_LK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankOrdinalDefCom_D1_CK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankOrdinalDefCom_D1_CK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankOrdinalDefCom_D1_CK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankOrdinalDefCom_D1_CK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankOrdinalDefCom_D1_CK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankOrdinalDefCom_D1_RK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankOrdinalDefCom_D1_RK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankOrdinalDefCom_D1_RK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankOrdinalDefCom_D1_RK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankOrdinalDefCom_D1_RK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankOrdinalDefCom_D1_PSSK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalDefCom_D1_PSSK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalDefCom_D1_PSSK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalDefCom_D1_PSSK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalDefCom_D1_PSSK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankOrdinalDefCom_D1_BSSK(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankOrdinalCusCom_D0_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalCusCom_D0_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalCusCom_D0_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalCusCom_D0_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalCusCom_D0_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankOrdinalCusCom_D1_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalCusCom_D1_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalCusCom_D1_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalCusCom_D1_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalCusCom_D1_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankOrdinalCusCom_D1_IK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankOrdinalCusCom_D1_IK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankOrdinalCusCom_D1_IK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankOrdinalCusCom_D1_IK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankOrdinalCusCom_D1_IK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankOrdinalCusCom_D1_LK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankOrdinalCusCom_D1_LK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankOrdinalCusCom_D1_LK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankOrdinalCusCom_D1_LK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankOrdinalCusCom_D1_LK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankOrdinalCusCom_D1_CK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankOrdinalCusCom_D1_CK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankOrdinalCusCom_D1_CK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankOrdinalCusCom_D1_CK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankOrdinalCusCom_D1_CK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankOrdinalCusCom_D1_RK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankOrdinalCusCom_D1_RK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankOrdinalCusCom_D1_RK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankOrdinalCusCom_D1_RK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankOrdinalCusCom_D1_RK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankOrdinalCusCom_D1_PSSK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankOrdinalCusCom_D1_PSSK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankOrdinalCusCom_D1_PSSK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankOrdinalCusCom_D1_PSSK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankOrdinalCusCom_D1_PSSK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
        procedure(logical(LK))                                      :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankOrdinalCusCom_D1_BSSK(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankOrdinalCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **ordinal rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>ordinal (`1234`) ranking</b>.<br>
    !>  In ordinal ranking, all items receive distinct ordinal numbers, including items that compare equal.<br>
    !>  The assignment of distinct ordinal numbers to items that compare equal can be done at random, or arbitrarily,
    !>  but it is generally preferable to use a system that is arbitrary but consistent,
    !>  as this gives stable results if the ranking is done multiple times.<br>
    !>  In computer data processing, ordinal ranking is also referred to as **row numbering**.
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **ordinal ranking** `1234`.<br>
    !>
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the ordinal rank of the `i`th element of `array`.**
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setRankOrdinal}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: setRankOrdinal
    !>
    !>      call setRankOrdinal(rank, array)
    !>      call setRankOrdinal(rank, array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{setRankOrdinal}
    !>  \include{lineno} example/pm_arrayRank/setRankOrdinal/main.F90
    !>  \compilef{setRankOrdinal}
    !>  \output{setRankOrdinal}
    !>  \include{lineno} example/pm_arrayRank/setRankOrdinal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{setRankOrdinal}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRankOrdinal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D0_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D0_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D0_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D0_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D0_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_IK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_IK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_IK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_IK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_IK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_LK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_LK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_LK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_LK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_LK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_CK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_CK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_CK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_CK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_CK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_RK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_RK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_RK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_RK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_RK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_PSSK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_PSSK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_PSSK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_PSSK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankOrdinalDefCom_D1_PSSK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRankOrdinalDefCom_D1_BSSK(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankOrdinalCusCom_D0_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankOrdinalCusCom_D0_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankOrdinalCusCom_D0_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankOrdinalCusCom_D0_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankOrdinalCusCom_D0_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_IK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_IK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_IK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_IK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_IK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_LK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_LK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_LK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_LK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_LK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_CK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_CK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_CK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_CK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_CK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_RK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_RK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_RK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_RK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_RK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setRankOrdinalCusCom_D1_PSSK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankOrdinalCusCom_D1_PSSK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankOrdinalCusCom_D1_PSSK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankOrdinalCusCom_D1_PSSK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankOrdinalCusCom_D1_PSSK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setRankOrdinalCusCom_D1_BSSK(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankOrdinalCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Standard rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>Standard Competition (`1224`) ranking</b>.<br>
    !>  In Standard Competition ranking, items that compare equal receive the same ranking number,
    !>  and then a gap is left in the ranking numbers. The number of ranking numbers that are left out
    !>  in this gap is one less than the number of items that compared equal.<br>
    !>  Equivalently, the ranking number of each item is `1` plus the number of items ranked above it.<br>
    !>  This ranking strategy is frequently adopted for competitions, as it means that if two (or more) competitors
    !>  tie for a position in the ranking, and the position of all those ranked below them is unaffected
    !>  (i.e., a competitor only comes second if exactly one person scores better than them,
    !>  third if exactly two people score better than them, fourth if exactly three people score better than them, etc.).<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
    !>  *B* gets ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `4` (*fourth*).
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Standard Competition ranking** `1224`.<br>
    !>
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Standard rank of the `i`th element of `array`.**
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `rank(1:size(array)`    :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` matches that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Standard rank of the `i`th element of `array`.**
    !>
    !>  \interface{getRankStandard}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: getRankStandard
    !>
    !>      rank(1:size(array)) = getRankStandard(array)
    !>      rank(1:size(array)) = getRankStandard(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{getRankStandard}
    !>  \include{lineno} example/pm_arrayRank/getRankStandard/main.F90
    !>  \compilef{getRankStandard}
    !>  \output{getRankStandard}
    !>  \include{lineno} example/pm_arrayRank/getRankStandard/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{getRankStandard}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getRankStandard

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankStandardDefCom_D0_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardDefCom_D0_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardDefCom_D0_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardDefCom_D0_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardDefCom_D0_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankStandardDefCom_D1_SK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardDefCom_D1_SK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardDefCom_D1_SK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardDefCom_D1_SK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardDefCom_D1_SK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankStandardDefCom_D1_IK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankStandardDefCom_D1_IK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankStandardDefCom_D1_IK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankStandardDefCom_D1_IK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankStandardDefCom_D1_IK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankStandardDefCom_D1_LK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankStandardDefCom_D1_LK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankStandardDefCom_D1_LK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankStandardDefCom_D1_LK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankStandardDefCom_D1_LK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankStandardDefCom_D1_CK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankStandardDefCom_D1_CK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankStandardDefCom_D1_CK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankStandardDefCom_D1_CK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankStandardDefCom_D1_CK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankStandardDefCom_D1_RK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankStandardDefCom_D1_RK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankStandardDefCom_D1_RK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankStandardDefCom_D1_RK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankStandardDefCom_D1_RK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankStandardDefCom_D1_PSSK5(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardDefCom_D1_PSSK4(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardDefCom_D1_PSSK3(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardDefCom_D1_PSSK2(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardDefCom_D1_PSSK1(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankStandardDefCom_D1_BSSK(array) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankStandardCusCom_D0_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardCusCom_D0_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardCusCom_D0_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardCusCom_D0_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardCusCom_D0_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(len(array, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getRankStandardCusCom_D1_SK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardCusCom_D1_SK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardCusCom_D1_SK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardCusCom_D1_SK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardCusCom_D1_SK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getRankStandardCusCom_D1_IK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK4_ENABLED
    module function getRankStandardCusCom_D1_IK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK3_ENABLED
    module function getRankStandardCusCom_D1_IK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK2_ENABLED
    module function getRankStandardCusCom_D1_IK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if IK1_ENABLED
    module function getRankStandardCusCom_D1_IK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getRankStandardCusCom_D1_LK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK4_ENABLED
    module function getRankStandardCusCom_D1_LK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK3_ENABLED
    module function getRankStandardCusCom_D1_LK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK2_ENABLED
    module function getRankStandardCusCom_D1_LK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if LK1_ENABLED
    module function getRankStandardCusCom_D1_LK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getRankStandardCusCom_D1_CK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK4_ENABLED
    module function getRankStandardCusCom_D1_CK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK3_ENABLED
    module function getRankStandardCusCom_D1_CK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK2_ENABLED
    module function getRankStandardCusCom_D1_CK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if CK1_ENABLED
    module function getRankStandardCusCom_D1_CK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getRankStandardCusCom_D1_RK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK4_ENABLED
    module function getRankStandardCusCom_D1_RK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK3_ENABLED
    module function getRankStandardCusCom_D1_RK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK2_ENABLED
    module function getRankStandardCusCom_D1_RK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if RK1_ENABLED
    module function getRankStandardCusCom_D1_RK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getRankStandardCusCom_D1_PSSK5(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK4_ENABLED
    module function getRankStandardCusCom_D1_PSSK4(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK3_ENABLED
    module function getRankStandardCusCom_D1_PSSK3(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK2_ENABLED
    module function getRankStandardCusCom_D1_PSSK2(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function
#endif

#if SK1_ENABLED
    module function getRankStandardCusCom_D1_PSSK1(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                                                :: rank(size(array, kind = IK))
        procedure(logical(LK))                                      :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getRankStandardCusCom_D1_BSSK(array, isSorted) result(rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRankStandardCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(TKR)                                                :: rank(size(array, kind = IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the **Standard rank** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the order specified by the input procedure `isSorted()` using the Quicksort algorithm such that `array(rank)`
    !>  will be in ascending order (or in the requested order as specified by `isSorted()`.
    !>
    !>  \details
    !>  This kind of ranking of values is widely known as <b>Standard Competition (`1224`) ranking</b>.<br>
    !>  In Standard Competition ranking, items that compare equal receive the same ranking number,
    !>  and then a gap is left in the ranking numbers. The number of ranking numbers that are left out
    !>  in this gap is one less than the number of items that compared equal.<br>
    !>  Equivalently, the ranking number of each item is `1` plus the number of items ranked above it.<br>
    !>  This ranking strategy is frequently adopted for competitions, as it means that if two (or more) competitors
    !>  tie for a position in the ranking, and the position of all those ranked below them is unaffected
    !>  (i.e., a competitor only comes second if exactly one person scores better than them,
    !>  third if exactly two people score better than them, fourth if exactly three people score better than them, etc.).<br>
    !>  Thus if *A* ranks ahead of *B* and *C* (which compare equal) which are both ranked ahead of *D*, then *A* gets ranking number `1` (*first*),
    !>  *B* gets ranking number `2` (*joint second*), *C* also gets ranking number `2` (*joint second*) and *D* gets ranking number `4` (*fourth*).
    !>  That is, if `A < B == C < D`, then the sequence `ABCD` has the **Standard Competition ranking** `1224`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized container of string of kind \SKALL) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) (container of string of default kind \SK) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements rankings will be computed and returned.
    !>  \param[out] rank        :   The output `contiguous` array of rank `1` of type `integer` of default kind \IK
    !>                              containing the ranks of the corresponding elements of `array`.<br>
    !>                              The size of `rank` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                              **Read `rank(i)` as the Standard rank of the `i`th element of `array`.**
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                              and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined sorting condition
    !>                              within `isSorted()`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a Fortran string (i.e., a scalar `character`),
    !>                              then both input arguments to `isSorted()` are single `character(1,SKG)` where `SKG` is the kind of `array`.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the rank of the input `array` is `1`,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      TYPE(KIND)  , intent(in)    :: a, b
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK)    , intent(in)    :: a, b
    !>                                      integer(IK)         , intent(in)    :: a, b
    !>                                      logical(LK)         , intent(in)    :: a, b
    !>                                      complex(CK)         , intent(in)    :: a, b
    !>                                      real(RK)            , intent(in)    :: a, b
    !>                                      type(css_type)      , intent(in)    :: a, b
    !>                                      type(css_pdt(SK))   , intent(in)    :: a, b
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `isSorted()` when the input `array` is a scalar string,
    !>                              \code{.F90}
    !>                                  function isSorted(a,b) result (sorted)
    !>                                      character(1,SKG), intent(in)    :: a, b
    !>                                      logical(LK)                     :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `SKG` represents the kind of the input string argument `array`.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setRankStandard}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRank, only: setRankStandard
    !>
    !>      call setRankStandard(rank, array)
    !>      call setRankStandard(rank, array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that the definition of `isSorted()`, if present, must be such that `isSorted() .and. .not. isSorted()`
    !>  is equivalent to an equality check for two elements of the input `array`. This equality check is used to
    !>  identify ties within the Standard ranking of the input `array`.
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted` is present.
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{setRankStandard}
    !>  \include{lineno} example/pm_arrayRank/setRankStandard/main.F90
    !>  \compilef{setRankStandard}
    !>  \output{setRankStandard}
    !>  \include{lineno} example/pm_arrayRank/setRankStandard/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRank](@ref test_pm_arrayRank)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with Intel ifort 2021.5 that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \todo
    !>  \plow The current bypass for the PDT name aliasing bug can be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{setRankStandard}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setRankStandard

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D0_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D0_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D0_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D0_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D0_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_SK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_SK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_SK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_SK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_SK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_IK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_IK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_IK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_IK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_IK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_LK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_LK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_LK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_LK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_LK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_CK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_CK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_CK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_CK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_CK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_RK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_RK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_RK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_RK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_RK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_PSSK5(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_PSSK4(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_PSSK3(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_PSSK2(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRankStandardDefCom_D1_PSSK1(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRankStandardDefCom_D1_BSSK(rank, array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardDefCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankStandardCusCom_D0_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D0_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankStandardCusCom_D0_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D0_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankStandardCusCom_D0_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D0_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankStandardCusCom_D0_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D0_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankStandardCusCom_D0_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D0_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setRankStandardCusCom_D1_SK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_SK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankStandardCusCom_D1_SK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_SK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankStandardCusCom_D1_SK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_SK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankStandardCusCom_D1_SK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_SK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankStandardCusCom_D1_SK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_SK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setRankStandardCusCom_D1_IK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_IK5
#endif
        use pm_kind, only: TKR => IK, IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setRankStandardCusCom_D1_IK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_IK4
#endif
        use pm_kind, only: TKR => IK, IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setRankStandardCusCom_D1_IK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_IK3
#endif
        use pm_kind, only: TKR => IK, IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setRankStandardCusCom_D1_IK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_IK2
#endif
        use pm_kind, only: TKR => IK, IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setRankStandardCusCom_D1_IK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_IK1
#endif
        use pm_kind, only: TKR => IK, IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setRankStandardCusCom_D1_LK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_LK5
#endif
        use pm_kind, only: TKR => IK, LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setRankStandardCusCom_D1_LK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_LK4
#endif
        use pm_kind, only: TKR => IK, LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setRankStandardCusCom_D1_LK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_LK3
#endif
        use pm_kind, only: TKR => IK, LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setRankStandardCusCom_D1_LK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_LK2
#endif
        use pm_kind, only: TKR => IK, LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setRankStandardCusCom_D1_LK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_LK1
#endif
        use pm_kind, only: TKR => IK, LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setRankStandardCusCom_D1_CK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_CK5
#endif
        use pm_kind, only: TKR => IK, CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setRankStandardCusCom_D1_CK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_CK4
#endif
        use pm_kind, only: TKR => IK, CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setRankStandardCusCom_D1_CK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_CK3
#endif
        use pm_kind, only: TKR => IK, CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setRankStandardCusCom_D1_CK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_CK2
#endif
        use pm_kind, only: TKR => IK, CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setRankStandardCusCom_D1_CK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_CK1
#endif
        use pm_kind, only: TKR => IK, CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setRankStandardCusCom_D1_RK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_RK5
#endif
        use pm_kind, only: TKR => IK, RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setRankStandardCusCom_D1_RK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_RK4
#endif
        use pm_kind, only: TKR => IK, RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setRankStandardCusCom_D1_RK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_RK3
#endif
        use pm_kind, only: TKR => IK, RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setRankStandardCusCom_D1_RK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_RK2
#endif
        use pm_kind, only: TKR => IK, RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setRankStandardCusCom_D1_RK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_RK1
#endif
        use pm_kind, only: TKR => IK, RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setRankStandardCusCom_D1_PSSK5(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_PSSK5
#endif
        use pm_kind, only: TKR => IK, SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setRankStandardCusCom_D1_PSSK4(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_PSSK4
#endif
        use pm_kind, only: TKR => IK, SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setRankStandardCusCom_D1_PSSK3(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_PSSK3
#endif
        use pm_kind, only: TKR => IK, SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setRankStandardCusCom_D1_PSSK2(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_PSSK2
#endif
        use pm_kind, only: TKR => IK, SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setRankStandardCusCom_D1_PSSK1(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_PSSK1
#endif
        use pm_kind, only: TKR => IK, SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setRankStandardCusCom_D1_BSSK(rank, array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRankStandardCusCom_D1_BSSK
#endif
        use pm_kind, only: TKR => IK, SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(TKR)                , intent(out)   , contiguous    :: rank(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayRank