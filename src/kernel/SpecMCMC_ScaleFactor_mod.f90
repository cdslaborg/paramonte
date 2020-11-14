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

module SpecMCMC_ScaleFactor_mod

    use Constants_mod, only: RK, IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ScaleFactor_mod"
    integer(IK), parameter          :: MAX_LEN_STRING_SCALE_FACTOR = 127

    character(:), allocatable       :: scaleFactor

    type                            :: ScaleFactor_type
        real(RK)                    :: val, defVal
        character(:), allocatable   :: str, defStr, null, desc
    contains
        procedure, pass             :: set => setScaleFactor, checkForSanity, nullifyNameListVar
    end type ScaleFactor_type

    interface ScaleFactor_type
        module procedure            :: constructScaleFactor
    end interface ScaleFactor_type

    private :: constructScaleFactor, setScaleFactor, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructScaleFactor(nd,methodName) result(ScaleFactorObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructScaleFactor
#endif
        use Constants_mod, only: RK, IK, NULL_SK
        use String_mod, only: num2str
        use Decoration_mod, only: TAB
        implicit none
        integer(IK), intent(in)     :: nd
        character(*), intent(in)    :: methodName
        type(ScaleFactor_type)      :: ScaleFactorObj

        ScaleFactorObj%defStr = "gelman"
        ScaleFactorObj%defVal = 2.38_RK/sqrt(real(nd,kind=RK))  ! Gelman, Roberts, Gilks (1996): Efficient Metropolis Jumping Rules
        ScaleFactorObj%null = repeat(NULL_SK, MAX_LEN_STRING_SCALE_FACTOR)
        ScaleFactorObj%desc = &
        "scaleFactor is a real-valued positive number (which must be given as string), by the square of which the &
        &covariance matrix of the proposal distribution of the MCMC sampler is scaled. In other words, &
        &the proposal distribution will be scaled in every direction by the value of scaleFactor. &
        &It can also be given in units of the string keyword 'gelman' (which is case-INsensitive) after the paper:\n\n" &
        // TAB // "Gelman, Roberts, and Gilks (1996): 'Efficient Metropolis Jumping Rules'.\n\n&
        &The paper finds that the optimal scaling factor for a Multivariate Gaussian proposal distribution for the &
        &Metropolis-Hastings Markov Chain Monte Carlo sampling of a target Multivariate Normal Distribution &
        &of dimension ndim is given by:\n\n&
        &    scaleFactor = 2.38/sqrt(ndim)  ,  in the limit of ndim -> Infinity.\n\n&
        &Multiples of the gelman scale factors are also acceptable as input and can be specified like the following examples:\n\n&
        &    scaleFactor = '1'\n\n&
        &            multiplies the ndim-dimensional proposal covariance matrix by 1, essentially no change occurs to &
                    &the covariance matrix.\n\n" // &
        '    scaleFactor = "1"\n\n' // &
        "            same as the previous example. The double-quotation marks act the same way as single-quotation marks.\n\n&
        &    scaleFactor = '2.5'\n\n&
        &            multiplies the ndim-dimensional proposal covariance matrix by 2.5.\n\n&
        &    scaleFactor = '2.5*Gelman'\n\n&
        &            multiplies the ndim-dimensional proposal covariance matrix by 2.5 * 2.38/sqrt(ndim).\n\n" // &
        '    scaleFactor = "2.5 * gelman"\n\n' // &
        "            same as the previous example, but with double-quotation marks. space characters are ignored.\n\n" // &
        '    scaleFactor = "2.5 * gelman*gelman*2"\n\n' // &
        "            equivalent to gelmanFactor-squared multiplied by 5.\n\n&
        &Note, however, that the result of Gelman et al. paper applies only to multivariate normal proposal distributions, in &
        &the limit of infinite dimensions. Therefore, care must be taken when using Gelman's scaling factor with non-Gaussian &
        &proposals and target objective functions. Note that only the product symbol (*) can be parsed &
        &in the string value of scaleFactor. The presence of other mathematical symbols or multiple appearances of the product &
        &symbol will lead to a simulation crash. Also, note that the prescription of an acceptance range specified by the input &
        &variable 'targetAcceptanceRate' will lead to dynamic modification of the initial input value of scaleFactor throughout sampling &
        &for adaptiveUpdateCount times. &
        &The default scaleFactor string-value is 'gelman' (for all proposals), which is subsequently converted to 2.38/sqrt(ndim)."
    end function constructScaleFactor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(ScaleFactorObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ScaleFactor_type), intent(in) :: ScaleFactorObj
        scaleFactor = ScaleFactorObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setScaleFactor(ScaleFactorObj,scaleFactor)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setScaleFactor
#endif
        use Constants_mod, only: RK
        implicit none
        class(ScaleFactor_type), intent(inout)  :: ScaleFactorObj
        character(*), intent(in)                :: scaleFactor
        ScaleFactorObj%str = trim(adjustl(scaleFactor))
        if (ScaleFactorObj%str==ScaleFactorObj%null) then
            ScaleFactorObj%str = ScaleFactorObj%defStr
        end if
    end subroutine setScaleFactor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: This subroutine also assigns the value of ScaleFactor. It MUST be executed by all images.
    subroutine checkForSanity(ScaleFactorObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str, String_type
        use Err_mod, only: Err_type
        implicit none
        class(ScaleFactor_type), intent(inout)  :: ScaleFactorObj
        character(*), intent(in)                :: methodName
        type(Err_type), intent(inout)           :: Err
        character(*), parameter                 :: PROCEDURE_NAME = "@checkForSanity()"
        type(String_type)                       :: String
        integer(IK)                             :: i

        ! First convert the scaleFactor string to real value:

        String%value = String%replaceStr(ScaleFactorObj%str," ","") ! remove the white spaces
        if (len_trim(adjustl(String%value))==0) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input string value (" // ScaleFactorObj%str // ") for the variable scaleFactor &
                        &is empty. Make sure the input string follows the syntax rules of " &
                        // methodName // " for this variable. Otherwise drop it from the input list. " &
                        // methodName // " will automatically assign an appropriate value to it.\n\n"
            return
        end if

        ! Now split the string by "*" to real coefficient and character (gelman) parts for further evaluations

        String%Parts = String%splitStr( string = String%value, delimiter = "*", nPart = String%nPart )
        ScaleFactorObj%val = 1._RK
        do i = 1, String%nPart
            if ( String%getLowerCase( String%Parts(i)%record ) == "gelman" ) then
                ScaleFactorObj%val = ScaleFactorObj%val * ScaleFactorObj%defVal
            else
                ScaleFactorObj%val = ScaleFactorObj%val * String%str2real64( str=String%Parts(i)%record, iostat=Err%stat )
                if ( Err%stat/=0 ) then
                    Err%occurred = .true.
                    Err%msg = Err%msg // &
                    MODULE_NAME // PROCEDURE_NAME // ": Error occurred while reading real number.\n&
                    &The input string value for the variable scaleFactor (" // ScaleFactorObj%str // ") does not appear to follow &
                    &the standard syntax rules of "// methodName // " for this variable. '" // String%Parts(i)%record // &
                    "' cannot be parsed into any meaningful token. Please correct the input value, or drop it from the input list, &
                    &in which case, " // methodName // " will automatically assign an appropriate value to it.\n\n"
                    return
                end if
            end if
        end do

        ! Now check if the real value is positive

        if (ScaleFactorObj%val<=0) then
            Err%occurred = .true.
            Err%msg = Err%msg // &
            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
            &The input string value (" // ScaleFactorObj%str // ") translates to a negative real value: " // &
            num2str(ScaleFactorObj%val) // ". &
            &Make sure the input string follows the syntax rules of " // methodName // " for this variable. &
            &Otherwise drop it from the input list. " // methodName // &
            " will automatically assign an appropriate value to it.\n\n"
            return
        end if

    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_ScaleFactor_mod