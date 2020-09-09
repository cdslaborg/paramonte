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

module SpecBase_ChainFileFormat_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_ChainFileFormat_mod"

    integer(IK), parameter          :: MAX_LEN_CHAIN_FILE_FORMAT = 63_IK

    character(MAX_LEN_CHAIN_FILE_FORMAT) :: chainFileFormat

    type                            :: ChainFileFormat_type
        logical                     :: isCompact
        logical                     :: isVerbose
        logical                     :: isBinary
        character(7)                :: compact
        character(7)                :: verbose
        character(6)                :: binary
        character(:), allocatable   :: def
        character(:), allocatable   :: val
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setChainFileFormat, checkForSanity, nullifyNameListVar
    end type ChainFileFormat_type

    interface ChainFileFormat_type
        module procedure            :: constructChainFileFormat
    end interface ChainFileFormat_type

    private :: constructChainFileFormat, setChainFileFormat, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructChainFileFormat(methodName) result(ChainFileFormatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructChainFileFormat
#endif
        use Constants_mod, only: NULL_SK, FILE_EXT, FILE_TYPE
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: methodName
        type(ChainFileFormat_type)  :: ChainFileFormatObj

        ChainFileFormatObj%isCompact = .false.
        ChainFileFormatObj%isVerbose = .false.
        ChainFileFormatObj%isBinary  = .false.
        ChainFileFormatObj%compact = "compact"
        ChainFileFormatObj%verbose = "verbose"
        ChainFileFormatObj%binary = FILE_TYPE%binary
        ChainFileFormatObj%def = ChainFileFormatObj%compact

        ChainFileFormatObj%null = repeat(NULL_SK, MAX_LEN_CHAIN_FILE_FORMAT)

        ChainFileFormatObj%desc = &
        "chainFileFormat is a string variable that represents the format of the output chain file(s) of "// methodName // &
        " simulation. The string value must be enclosed by either single or double quotation marks when provided as input. &
        &Three values are possible:\n\n&
        &    chainFileFormat = 'compact'\n\n&
        &            This is the ASCII (text) file format which is human-readable but does not preserve the full accuracy of the &
                    &output values. It is also a significantly slower mode of chain file generation, compared to the binary file format (see below). &
                    &If the compact format is specified, each of the repeating MCMC states will be condensed into a single entry (row) in &
                    &the output MCMC chain file. Each entry will be then assigned a sample-weight that is equal to the number of repetitions of &
                    &that state in the MCMC chain. Thus, each row in the output chain file will represent a unique sample from the objective function. &
                    &This will lead to a significantly smaller ASCII chain file and faster output size compared to the verbose chain file format (see below).\n\n&
        &    chainFileFormat = 'verbose'\n\n&
        &            This is the ASCII (text) file format which is human-readable but does not preserve the full accuracy of &
                    &the output values. It is also a significantly slower mode of chain file generation, &
                    &compared to both compact and binary chain file formats (see above and below). &
                    &If the verbose format is specified, all MCMC states will have equal sample-weights of 1 in the output chain file. &
                    &The verbose format can lead to much larger chain file sizes than the compact and binary file formats. &
                    &This is especially true if the target objective function has a very high-dimensional state space.\n\n&
        &    chainFileFormat = '" // ChainFileFormatObj%binary // "'\n\n&
        &            This is the binary file format which is not human-readable, but preserves the exact values in the output &
                    &MCMC chain file. It is also often the fastest mode of chain file generation. If the binary file format is chosen, the chain &
                    &will be automatically output in the compact format (but as binary) to ensure the production of the smallest-possible output chain file. &
                    &Binary chain files will have the " // FILE_EXT%binary // " file extensions. Use the binary format if you need full accuracy representation &
                    &of the output values while having the smallest-size output chain file in the shortest time possible.\n\n&
        &The default value is chainFileFormat = '" // ChainFileFormatObj%def // "' as it provides a reasonable trade-off between &
        &speed and output file size while generating human-readable chain file contents. Note that the input values are case-insensitive."
    end function constructChainFileFormat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(ChainFileFormatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ChainFileFormat_type), intent(in) :: ChainFileFormatObj
        chainFileFormat = ChainFileFormatObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setChainFileFormat(ChainFileFormatObj,chainFileFormat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setChainFileFormat
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(ChainFileFormat_type), intent(inout)    :: ChainFileFormatObj
        character(*), intent(in)                      :: chainFileFormat
        ChainFileFormatObj%val = trim(adjustl(chainFileFormat))
        if ( ChainFileFormatObj%val==trim(adjustl(ChainFileFormatObj%null)) ) then
            ChainFileFormatObj%val = trim(adjustl(ChainFileFormatObj%def))
        end if
        if (getLowerCase(ChainFileFormatObj%val)==getLowerCase(ChainFileFormatObj%compact)) ChainFileFormatObj%iscompact = .true.
        if (getLowerCase(ChainFileFormatObj%val)==getLowerCase(ChainFileFormatObj%verbose)) ChainFileFormatObj%isverbose = .true.
        if (getLowerCase(ChainFileFormatObj%val)==getLowerCase(ChainFileFormatObj%binary)) ChainFileFormatObj%isBinary = .true.
    end subroutine setChainFileFormat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(ChainFileFormat,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(ChainFileFormat_type), intent(in)   :: ChainFileFormat
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
        if ( .not.(ChainFileFormat%isCompact .or. ChainFileFormat%isVerbose .or. ChainFileFormat%isBinary) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested chain file format ('" // ChainFileFormat%val // &
                        "') represented by the variable chainFileFormat cannot be anything other than '" // &
                        ChainFileFormat%compact // "' or '" // ChainFileFormat%verbose // "' or '" // ChainFileFormat%binary // "'. &
                        &If you don't know an appropriate value for chainFileFormat, drop it from the input list. " // methodName // &
                        " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_ChainFileFormat_mod