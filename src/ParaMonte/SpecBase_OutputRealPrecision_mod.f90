module SpecBase_OutputRealPrecision_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_OutputRealPrecision_mod"

    integer(IK), parameter          :: OUTPUT_REAL_PRECISION = 8

    integer(IK)                     :: outputRealPrecision ! namelist input

    type                            :: OutputRealPrecision_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: str
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setOutputRealPrecision, checkForSanity, nullifyNameListVar
    end type OutputRealPrecision_type

    interface OutputRealPrecision_type
        module procedure                :: constructOutputRealPrecision
    end interface OutputRealPrecision_type

    private :: constructOutputRealPrecision, setOutputRealPrecision

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructOutputRealPrecision(methodName) result(OutputRealPrecisionObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructOutputRealPrecision
#endif
        use Constants_mod, only: NULL_IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)       :: methodName
        type(OutputRealPrecision_type) :: OutputRealPrecisionObj
        OutputRealPrecisionObj%def = OUTPUT_REAL_PRECISION
        OutputRealPrecisionObj%null = NULL_IK
        OutputRealPrecisionObj%desc = &
        "The variable outputRealPrecision is a 32-bit integer number that determines the precision - that is, the number of &
        &significant digits - of the real numbers in the output files of " // methodName // ". Any positive integer is acceptable &
        &as the input value of outputRealPrecision. However, any digits of the output real numbers beyond the accuracy of 64-bit &
        &real numbers (approximately 16 digits of significance) will be meaningless and random. Set this variable to 16 (or larger) &
        &if full reproducibility of the simulation is needed in the future. But keep in mind that larger precisions will result in &
        &larger-size output files. This variable is ignored for binary output (if any occurs during the simulation). &
        &The default value is " // num2str(OutputRealPrecisionObj%def) // "."
    end function constructOutputRealPrecision

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(OutputRealPrecisionObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(OutputRealPrecision_type), intent(in)  :: OutputRealPrecisionObj
        outputRealPrecision = OutputRealPrecisionObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setOutputRealPrecision(OutputRealPrecisionObj,outputRealPrecision)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setOutputRealPrecision
#endif
        use String_mod, only: num2str
        use Constants_mod, only: IK
        implicit none
        class(OutputRealPrecision_type), intent(inout)  :: OutputRealPrecisionObj
        integer(IK), intent(in)                         :: outputRealPrecision
        OutputRealPrecisionObj%val = outputRealPrecision
        if (OutputRealPrecisionObj%val==OutputRealPrecisionObj%null) then
            OutputRealPrecisionObj%val = OutputRealPrecisionObj%def
        end if
        OutputRealPrecisionObj%str = num2str(OutputRealPrecisionObj%val)
    end subroutine setOutputRealPrecision

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(OutputRealPrecisionObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(OutputRealPrecision_type), intent(in) :: OutputRealPrecisionObj
        character(*), intent(in)               :: methodName
        type(Err_type), intent(inout)          :: Err
        character(*), parameter                :: PROCEDURE_NAME = "@checkForSanity()"
        if ( OutputRealPrecisionObj%val<1_IK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable outputRealPrecision must be a positive integer < 16. If you are not sure &
                        &about the appropriate value for this variable, simply drop it from the input. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_OutputRealPrecision_mod