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

module SpecMCMC_ProposalStartCorMat_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ProposalStartCorMat_mod"

    real(RK), allocatable           :: ProposalStartCorMat(:,:) ! namelist input

    type                            :: ProposalStartCorMat_type
        real(RK), allocatable       :: Val(:,:)
        real(RK), allocatable       :: Def(:,:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalStartCorMat, checkForSanity, nullifyNameListVar
    end type ProposalStartCorMat_type

    interface ProposalStartCorMat_type
        module procedure            :: constructProposalStartCorMat
    end interface ProposalStartCorMat_type

    private :: constructProposalStartCorMat, setProposalStartCorMat, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructProposalStartCorMat(nd,methodName) result(ProposalStartCorMatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalStartCorMat
#endif
        use Constants_mod, only: IK, NULL_RK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)         :: nd
        character(*), intent(in)        :: methodName
        type(ProposalStartCorMat_type)  :: ProposalStartCorMatObj
        integer(IK)                     :: i
        allocate( ProposalStartCorMatObj%Def(nd,nd) )
        ProposalStartCorMatObj%Def = 0._RK
        do i = 1,nd
            ProposalStartCorMatObj%Def(i,i) = 1._RK
        end do
        ProposalStartCorMatObj%null   = NULL_RK
        ProposalStartCorMatObj%desc   = &
        "ProposalStartCorMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the &
        &sampling space. It serves as the best-guess starting correlation matrix of the proposal distribution used by " &
        // methodName // ". &
        &It is used (along with the input vector ProposalStartStdVec) to construct the covariance matrix of the proposal &
        &distribution when the input covariance matrix is missing in the input list of variables. &
        &If the covariance matrix is given as input to "//methodName//", any input values for ProposalStartCorMat, &
        &as well as ProposalStartStdVec, will be automatically ignored by "//methodName//". As input to " // methodName // &
        ", the variable ProposalStartCorMat along with ProposalStartStdVec is especially useful in situations where &
        &obtaining the best-guess covariance matrix is not trivial. The default value of ProposalStartCorMat is an ndim-by-ndim Identity matrix."
    end function constructProposalStartCorMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(ProposalStartCorMatObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProposalStartCorMat_type), intent(in) :: ProposalStartCorMatObj
        integer(IK), intent(in)                     :: nd
        if (allocated(ProposalStartCorMat)) deallocate(ProposalStartCorMat)
        allocate(ProposalStartCorMat(nd,nd))
        ProposalStartCorMat = ProposalStartCorMatObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setProposalStartCorMat(ProposalStartCorMatObj,ProposalStartCorMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalStartCorMat
#endif
        use Constants_mod, only: RK
        implicit none
        class(ProposalStartCorMat_type), intent(inout)  :: ProposalStartCorMatObj
        real(RK), intent(in)                            :: ProposalStartCorMat(:,:)
        ProposalStartCorMatObj%Val = ProposalStartCorMat
        where (ProposalStartCorMatObj%Val==ProposalStartCorMatObj%null)
                ProposalStartCorMatObj%Val = ProposalStartCorMatObj%Def
        end where
    end subroutine setProposalStartCorMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(ProposalStartCorMatObj,Err,methodName,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Matrix_mod, only: isPosDef
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(ProposalStartCorMat_type), intent(in) :: ProposalStartCorMatObj
        integer(IK), intent(in)             :: nd
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if (.not.isPosDef(nd,ProposalStartCorMatObj%Val)) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested ProposalStartCorMat for the proposal of " // methodName // &
                        " is not a positive-definite matrix.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_ProposalStartCorMat_mod