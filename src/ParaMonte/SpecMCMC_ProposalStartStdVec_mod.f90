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
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module SpecMCMC_ProposalStartStdVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ProposalStartStdVec_mod"

    real(RK), allocatable           :: proposalStartStdVec(:) ! namelist input

    type                            :: ProposalStartStdVec_type
        real(RK), allocatable       :: val(:)
        real(RK), allocatable       :: def(:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalStartCorMat, checkForSanity, nullifyNameListVar
    end type ProposalStartStdVec_type

    interface ProposalStartStdVec_type
        module procedure            :: constructProposalStartStdVec
    end interface ProposalStartStdVec_type

    private :: constructProposalStartStdVec, setProposalStartCorMat, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructProposalStartStdVec(nd,methodName) result(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalStartStdVec
#endif
        use Constants_mod, only: IK, NULL_RK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)         :: nd
        character(*), intent(in)        :: methodName
        type(ProposalStartStdVec_type)  :: self
        integer(IK)                     :: i
        allocate( self%def(nd) )
        self%def = 0._RK
        do i = 1,nd
            self%def(i) = 1._RK
        end do
        self%null   = NULL_RK
        self%desc   = &
        "proposalStartStdVec is a real-valued positive vector of length ndim, where ndim is the dimension of the sampling space. &
        &It serves as the best-guess starting Standard Deviation of each of the components of the proposal distribution. &
        &If the initial covariance matrix (ProposalStartCovMat) is missing as an input variable to " // &
        methodName // ", then proposalStartStdVec (along with the input variable ProposalStartCorMat) will be used to construct &
        &the initial covariance matrix of the proposal distribution of the MCMC sampler. &
        &However, if ProposalStartCovMat is present as an input argument to " // methodName // ", then the input proposalStartStdVec &
        &along with the input ProposalStartCorMat will be completely ignored and the input value for ProposalStartCovMat &
        &will be used to construct the initial covariance matrix of the proposal distribution of " // methodName // ". &
        &The default value of proposalStartStdVec is a vector of unit values (i.e., ones) of length ndim."
    end function constructProposalStartStdVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(self,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProposalStartStdVec_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        if (allocated(proposalStartStdVec)) deallocate(proposalStartStdVec)
        allocate(proposalStartStdVec(nd))
        proposalStartStdVec = self%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setProposalStartCorMat(self,proposalStartStdVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalStartCorMat
#endif
        use Constants_mod, only: RK
        implicit none
        class(ProposalStartStdVec_type), intent(inout)  :: self
        real(RK), intent(in)                            :: proposalStartStdVec(:)
        self%val = proposalStartStdVec
        where (self%val==self%null)
                self%val = self%def
        end where
    end subroutine setProposalStartCorMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(self,Err,methodName,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(ProposalStartStdVec_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
        integer(IK)                                 :: i
        do i = 1,nd
            if (self%val(i)<=0._RK) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                            &The input requested value (" // num2str(self%val(i)) // ") for the component " // &
                            num2str(i) // " of the variable proposalStartStdVec for the proposal distribution of " // &
                            methodName // " must be a positive real number.\n\n"
            end if
        end do
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_ProposalStartStdVec_mod