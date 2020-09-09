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

#if defined CFI_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The procedural C interface to ParaDRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED)
function &
#else
subroutine &
#endif
    runParaDRAM ( ndim          &
                , getLogFunc4C  &
                , InputFileVec  &
                , lenInputFile  &
                ) bind(C, name="runParaDRAM")
#if defined DLL_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
#endif
    use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_f_procpointer !, c_null_char
    use ParaMonteLogFunc_mod, only: getLogFunc_proc
    use Constants_mod, only: IK, RK
    use ParaDRAM_mod, only: ParaDRAM_type

    implicit none

    integer(IK), intent(in), value                          :: ndim
    integer(IK), intent(in), value                          :: lenInputFile
    type(c_funptr), intent(in), value                       :: getLogFunc4C
    character(len=1,kind=c_char), dimension(*), intent(in)  :: InputFileVec
    procedure(getLogFunc_proc), pointer                     :: getLogFunc
    character(:), allocatable                               :: inputFileStr
    type(ParaDRAM_type)                                     :: self
    integer                                                 :: i
#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED)
    integer(IK)                                             :: runParaDRAM
    runParaDRAM = 0_IK
#endif

    ! reconstruct the input file

    !i = 1
    !do
    !    if (InputFileVec(i)==c_null_char) exit
    !    i = i + 1_IK
    !end do

    if (lenInputFile==0_IK) then
        inputFileStr = ""
    else
        allocate(character(lenInputFile) :: inputFileStr)
        do i = 1, lenInputFile
            inputFileStr(i:i) = InputFileVec(i)
        end do
    end if

    ! associate the input C procedure pointer to a Fortran procedure pointer

    call c_f_procpointer(cptr=getLogFunc4C, fptr=getLogFunc)

    ! call runParaDRAM

    if (ndim>0_IK) then
        call self%runSampler( ndim = ndim               &
                            , getLogFunc = getLogFunc   &
                            , inputFile = inputFileStr  &
                            )
    end if
    nullify(getLogFunc)

#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED)
    if (self%Err%occurred) runParaDRAM = -1_IK
end function runParaDRAM
#else
end subroutine runParaDRAM
#endif

#else

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The procedural Fortran interface to ParaDRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine runParaDRAM  ( ndim          &
                        , getLogFunc    &
                        , inputFile     &
                        ) !bind(C, name="runParaDRAM")
#if defined DLL_ENABLED
    !DEC$ ATTRIBUTES DLLEXPORT :: runParaDRAM
#endif
    use ParaMonteLogFunc_mod, only: getLogFunc_proc
    use Constants_mod, only: IK, RK
    use ParaDRAM_mod, only: ParaDRAM_type

    implicit none

    integer(IK) , intent(in)            :: ndim
    procedure(getLogFunc_proc)          :: getLogFunc
    character(*), intent(in), optional  :: inputFile

    type(ParaDRAM_type)                 :: self

    ! call runParaDRAM

    call self%runSampler( ndim = ndim               &
                        , getLogFunc = getLogFunc   &
                        , inputFile = inputFile     &
                        )

end subroutine runParaDRAM

#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
