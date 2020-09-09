!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SpecBase_ParallelizationModel_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_ParallelizationModel_mod"

    integer(IK), parameter          :: MAX_LEN_PARALLELIZATION_MODEL = 63

    character(MAX_LEN_PARALLELIZATION_MODEL) :: parallelizationModel

    type                            :: ParallelizationModel_type
        logical                     :: isSinglChain
        logical                     :: isMultiChain
        logical                     :: isForkJoin
        character(10)               :: multiChain
        character(11)               :: singlChain
        character(:), allocatable   :: def
        character(:), allocatable   :: val
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setParallelizationModel, checkForSanity, nullifyNameListVar
    end type ParallelizationModel_type

    interface ParallelizationModel_type
        module procedure            :: constructParallelizationModel
    end interface ParallelizationModel_type

    private :: constructParallelizationModel, setParallelizationModel, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructParallelizationModel(methodName) result(ParallelizationModelObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructParallelizationModel
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: NULL_SK, IK, PMSM
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)        :: methodName
        type(ParallelizationModel_type) :: ParallelizationModelObj

        ParallelizationModelObj%isSinglChain = .false.
        ParallelizationModelObj%isMultiChain = .false.
        ParallelizationModelObj%singlChain = "singleChain"
        ParallelizationModelObj%multiChain = "multiChain"
        ParallelizationModelObj%def = ParallelizationModelObj%singlChain
        ParallelizationModelObj%null = repeat(NULL_SK, MAX_LEN_PARALLELIZATION_MODEL)

        ParallelizationModelObj%desc = &
        "parallelizationModel is a string variable that represents the parallelization method to be used in "// methodName //". &
        &The string value must be enclosed by either single or double quotation marks when provided as input. "
        if (methodName==PMSM%ParaDRAM .or. methodName==PMSM%ParaDISE) then
            ParallelizationModelObj%desc = ParallelizationModelObj%desc // &
            "Two options are currently supported:\n\n&
            &    parallelizationModel = '" // ParallelizationModelObj%multiChain // "'\n\n&
            &            This method uses the Prefect Parallelism scheme in which multiple MCMC chains are generated &
                        &independently of each other. In this case, multiple output MCMC chain files will also be generated.\n\n&
            &    parallelizationModel = '" // ParallelizationModelObj%singlChain // "'\n\n&
            &            This method uses the fork-style parallelization scheme. &
                        &A single MCMC chain file will be generated in this case. At each MCMC step multiple proposal steps &
                        &will be checked in parallel until one proposal is accepted.\n\n&
            &Note that in serial mode, there is no parallelism. Therefore, this option does not affect non-parallel simulations &
            &and its value is ignored. The serial mode is equivalent to either of the parallelism methods with only one simulation &
            &image (processor, core, or thread). &
            &The default value is parallelizationModel = '" // ParallelizationModelObj%def // "'. &
            &Note that the input values are case-insensitive and white-space characters are ignored."
        else
            block
                use Err_mod, only: Err_type, abort
                type(Err_type) :: Err
                Err%occurred = .true.
                Err%msg = MODULE_NAME//": Catastrophic internal error occurred. The simulation method name is not recognized."
                call abort(Err)
            end block
        end if
    end function constructParallelizationModel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(ParallelizationModelObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ParallelizationModel_type), intent(in) :: ParallelizationModelObj
        parallelizationModel = ParallelizationModelObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setParallelizationModel(ParallelizationModelObj,parallelizationModel)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setParallelizationModel
#endif
        use String_mod, only: getLowerCase, replaceStr
        implicit none
        class(ParallelizationModel_type), intent(inout) :: ParallelizationModelObj
        character(*), intent(in)                        :: parallelizationModel
        ParallelizationModelObj%val = trim(adjustl(replaceStr(parallelizationModel," ", "")))
        if (ParallelizationModelObj%val==trim(adjustl(ParallelizationModelObj%null))) ParallelizationModelObj%val = trim(adjustl(ParallelizationModelObj%def))
        if (getLowerCase(ParallelizationModelObj%val)==getLowerCase(ParallelizationModelObj%singlChain)) ParallelizationModelObj%isSinglChain = .true.
        if (getLowerCase(ParallelizationModelObj%val)==getLowerCase(ParallelizationModelObj%multiChain)) ParallelizationModelObj%isMultiChain = .true.
    end subroutine setParallelizationModel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(ParallelizationModel,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(ParallelizationModel_type), intent(in)    :: ParallelizationModel
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        if ( .not.(ParallelizationModel%isSinglChain .or. ParallelizationModel%isMultiChain) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested parallelization method (" // ParallelizationModel%val // &
                        ") represented by variable parallelizationModel cannot be anything other than &
                        &'singleChain' or 'multiChain'. If you don't know an appropriate value &
                        &for ParallelizationModel, drop it from the input list. " // methodName // &
                        " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_ParallelizationModel_mod