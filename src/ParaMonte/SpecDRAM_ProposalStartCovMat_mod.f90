!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module SpecDRAM_ProposalStartCovMat_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_ProposalStartCovMat_mod"

    real(RK), allocatable           :: ProposalStartCovMat(:,:) ! namelist input

    type                            :: ProposalStartCovMat_type
        logical                     :: isPresent
        real(RK), allocatable       :: Val(:,:)
        real(RK), allocatable       :: Def(:,:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalStartCorMat, checkForSanity, nullifyNameListVar
    end type ProposalStartCovMat_type

    interface ProposalStartCovMat_type
        module procedure            :: constructProposalStartCovMat
    end interface ProposalStartCovMat_type

    private :: constructProposalStartCovMat, setProposalStartCorMat, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructProposalStartCovMat(nd,methodName) result(ProposalStartCovMatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalStartCovMat
#endif
        use Constants_mod, only: IK, NULL_RK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)         :: nd
        character(*), intent(in)        :: methodName
        type(ProposalStartCovMat_type)  :: ProposalStartCovMatObj
        integer(IK)                     :: i
        ProposalStartCovMatObj%isPresent = .false.
        allocate( ProposalStartCovMatObj%Def(nd,nd) )
        ProposalStartCovMatObj%Def = 0._RK
        do i = 1,nd
            ProposalStartCovMatObj%Def(i,i) = 1._RK
        end do
        ProposalStartCovMatObj%null   = NULL_RK
        ProposalStartCovMatObj%desc   = &
        "ProposalStartCovMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the &
        &sampling space. It serves as the best-guess starting covariance matrix of the proposal distribution. &
        &To bring the sampling efficiency of " // methodName // " to within the desired requested range, the covariance matrix will &
        &be adaptively updated throughout the simulation, according to the user's requested schedule. If ProposalStartCovMat &
        &is not provided by the user, its value will be automatically computed from the input variables ProposalStartCorMat and &
        &ProposalStartStdVec. The default value of ProposalStartCovMat is an ndim-by-ndim Identity matrix."
    end function constructProposalStartCovMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(ProposalStartCovMatObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProposalStartCovMat_type), intent(in) :: ProposalStartCovMatObj
        integer(IK), intent(in)                     :: nd
        if (allocated(ProposalStartCovMat)) deallocate(ProposalStartCovMat)
        allocate(ProposalStartCovMat(nd,nd))
        ProposalStartCovMat = ProposalStartCovMatObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setProposalStartCorMat(ProposalStartCovMatObj,ProposalStartCovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalStartCorMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        class(ProposalStartCovMat_type), intent(inout)  :: ProposalStartCovMatObj
        real(RK), intent(in)                            :: ProposalStartCovMat(:,:)
        integer(IK)                                     :: i,j,nd
        nd = size(ProposalStartCovMat(:,1))
        ProposalStartCovMatObj%Val = ProposalStartCovMat
        ProposalStartCovMatObj%isPresent = .false.
        do i = 1, nd
            do j = 1, nd
                if (ProposalStartCovMatObj%Val(j,i)==ProposalStartCovMatObj%null) then
                    ProposalStartCovMatObj%Val(j,i) = ProposalStartCovMatObj%Def(j,i)
                else
                    ProposalStartCovMatObj%isPresent = .true.
                end if
            end do
        end do
    end subroutine setProposalStartCorMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(ProposalStartCovMatObj,Err,methodName,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Matrix_mod, only: isPosDef
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(ProposalStartCovMat_type), intent(in) :: ProposalStartCovMatObj
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
        if (.not.isPosDef(nd,ProposalStartCovMatObj%Val)) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested ProposalStartCovMat for the proposal of " // methodName // &
                        " is not a positive-definite matrix.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecDRAM_ProposalStartCovMat_mod