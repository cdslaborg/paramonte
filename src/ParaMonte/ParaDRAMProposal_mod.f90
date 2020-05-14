!***********************************************************************************************************************************
!***********************************************************************************************************************************
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
!   we ask you to acknowledge the ParaMonte library's usage
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module ParaDRAMProposal_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaDRAMProposal_mod"

    type(Err_type), save :: ProposalErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    type, abstract :: Proposal_type
    contains
        procedure(getNew_proc)          , nopass    , deferred  :: getNew
        procedure(doAdaptation_proc)    , nopass    , deferred  :: doAdaptation
        procedure(readRestartFile_proc) , nopass    , deferred  :: readRestartFile
        procedure(writeRestartFile_proc), nopass    , deferred  :: writeRestartFile
#if defined CAF_ENABLED || defined MPI_ENABLED
        procedure(getAdaptation_proc)   , nopass    , deferred  :: getAdaptation
#endif
    end type Proposal_type

    !*******************************************************************************************************************************

#if defined CAF_ENABLED || defined MPI_ENABLED
    abstract interface
        subroutine getAdaptation_proc()
        end subroutine getAdaptation_proc
    end interface
#endif

    abstract interface
        function getNew_proc( nd            &
                            , counterDRS    &
                            , StateOld      &
                            ) result (StateNew)
            use Constants_mod, only: IK, RK
            import :: Proposal_type
           !class(Proposal_type), intent(inout) :: Proposal
            integer(IK), intent(in)             :: nd
            integer(IK), intent(in)             :: counterDRS
            real(RK)   , intent(in)             :: StateOld(nd)
            real(RK)                            :: StateNew(nd)
        end function getNew_proc
    end interface

    abstract interface
        subroutine doAdaptation_proc( nd                        &
                                    , chainSize                 &
                                    , Chain                     &
                                    , ChainWeight               &
                                    , samplerUpdateIsGreedy     &
                                    , meanAccRateSinceStart     &
                                    , samplerUpdateSucceeded    &
                                    , adaptationMeasure         &
                                    )
            use Constants_mod, only: IK, RK
            import :: Proposal_type
           !class(Proposal_type), intent(inout) :: Proposal
            integer(IK), intent(in)             :: nd
            integer(IK), intent(in)             :: chainSize
            real(RK)   , intent(in)             :: Chain(nd,chainSize)
            integer(IK), intent(in)             :: ChainWeight(chainSize)
            logical    , intent(in)             :: samplerUpdateIsGreedy
            real(RK)   , intent(in)             :: meanAccRateSinceStart
            logical    , intent(out)            :: samplerUpdateSucceeded
            real(RK)   , intent(out)            :: adaptationMeasure
        end subroutine doAdaptation_proc
    end interface

    abstract interface
        subroutine readRestartFile_proc()
        end subroutine readRestartFile_proc
    end interface

    abstract interface
        subroutine writeRestartFile_proc()
        end subroutine writeRestartFile_proc
    end interface

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module ParaDRAMProposal_mod