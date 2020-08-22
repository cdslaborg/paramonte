!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

#if defined PARADRAM
    character(*), parameter :: MODULE_NAME = "@ParaDRAMProposal_mod"
#elif defined PARADISE
    character(*), parameter :: MODULE_NAME = "@ParaDISEProposal_mod"
#endif

    type(Err_type), save :: ProposalErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, abstract :: ProposalAbstract_type
    contains
        procedure(getNew_proc)                  , nopass    , deferred  :: getNew
        procedure(getLogProb_proc)              , nopass    , deferred  :: getLogProb
        procedure(doAdaptation_proc)            , nopass    , deferred  :: doAdaptation
       !procedure(readRestartFileAscii_proc)    , nopass    , deferred  :: readRestartFileAscii
       !procedure(writeRestartFileAscii_proc)   , nopass    , deferred  :: writeRestartFileAscii
#if defined CAF_ENABLED || defined MPI_ENABLED
        procedure(bcastAdaptation_proc)         , nopass    , deferred  :: bcastAdaptation
#endif
    end type ProposalAbstract_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined CAF_ENABLED || defined MPI_ENABLED
    abstract interface
        subroutine bcastAdaptation_proc()
        end subroutine bcastAdaptation_proc
    end interface
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    function getNew_proc( nd            &
                        , counterDRS    &
                        , StateOld      &
                        ) result (StateNew)
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
       !class(ProposalAbstract_type), intent(inout) :: Proposal
        integer(IK), intent(in) :: nd
        integer(IK), intent(in) :: counterDRS
        real(RK)   , intent(in) :: StateOld(nd)
        real(RK)                :: StateNew(nd)
    end function getNew_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    function getLogProb_proc( nd                &
                            , counterDRS        &
                            , StateOld          &
                            , StateNew          &
                            ) result (logProb)
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
        integer(IK), intent(in) :: nd
        integer(IK), intent(in) :: counterDRS
        real(RK)   , intent(in) :: StateOld(nd)
        real(RK)   , intent(in) :: StateNew(nd)
        real(RK)                :: logProb
    end function getLogProb_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    abstract interface
    subroutine doAdaptation_proc( nd                        &
                                , chainSize                 &
                                , Chain                     &
                                , ChainWeight               &
                                , isFreshRun                &
                                , samplerUpdateIsGreedy     &
                                , meanAccRateSinceStart     &
                                , samplerUpdateSucceeded    &
                                , adaptationMeasure         &
                                )
        use Constants_mod, only: IK, RK
        import :: ProposalAbstract_type
       !class(ProposalAbstract_type), intent(inout) :: Proposal
        integer(IK), intent(in)     :: nd
        integer(IK), intent(in)     :: chainSize
        real(RK)   , intent(in)     :: Chain(nd,chainSize)
        integer(IK), intent(in)     :: ChainWeight(chainSize)
        logical    , intent(in)     :: isFreshRun
        logical    , intent(in)     :: samplerUpdateIsGreedy
        real(RK)   , intent(inout)  :: meanAccRateSinceStart
        logical    , intent(out)    :: samplerUpdateSucceeded
        real(RK)   , intent(out)    :: adaptationMeasure
    end subroutine doAdaptation_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !abstract interface
    !subroutine readRestartFileAscii_proc()
    !end subroutine readRestartFileAscii_proc
    !end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !abstract interface
    !subroutine writeRestartFileAscii_proc()
    !end subroutine writeRestartFileAscii_proc
    !end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


