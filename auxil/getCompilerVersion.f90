!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   MIT License
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   Permission is hereby granted, free of charge, to any person obtaining a 
!   copy of this software and associated documentation files (the "Software"), 
!   to deal in the Software without restriction, including without limitation 
!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!   and/or sell copies of the Software, and to permit persons to whom the 
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be 
!   included in all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!   ACKNOWLEDGMENT
!
!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!   As per the ParaMonte library license agreement terms, if you use any parts of 
!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!   work (education/research/industry/development/...) by citing the ParaMonte 
!   library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program getCompilerVersion

    use iso_fortran_env, only: compiler_version, IK => int32
    implicit none
    integer(IK) , parameter     :: INTEL_VERSION(3) = [18_IK,0_IK,0_IK]
    integer(IK) , parameter     :: GNU_VERSION(3) = [7_IK,3_IK,0_IK]
    integer(IK) , allocatable   :: MinVersion(:)
    integer(IK)                 :: fileunit, startindex, endindex
    logical                     :: isIntel = .false.
    logical                     :: isGNU = .false.
    character(:), allocatable   :: string, isParaMonteCompatibleCompiler

    type :: CharVec_type
        character(:)    , allocatable   :: record
    end type CharVec_type

    type(CharVec_type)  , allocatable   :: VersionParts(:)

    isParaMonteCompatibleCompiler = "false"
    string = trim(adjustl(compiler_version()))

    if ( index(string,"GCC") > 0 ) then ! it's Gfortran
        startindex = index(string,"version") + 8
        endindex = len(string)
        isGNU = .true.
    elseif ( index(string,"Intel") > 0 ) then ! it's Intel
        startindex = index(string,"Version") + 8
        endindex = index(string,"Build") - 2
        isIntel = .true.
    end if

    string = trim(adjustl(string(startindex:endindex)))
    open(newunit=fileunit,file="getCompilerVersion.tmp")
    write(fileunit,"(A)") string
    close(fileunit)

    ! check for ParaMonteCompatibleCompiler

    !if ( lge(string(1:3),"7.3") ) isParaMonteCompatibleCompiler = "true"
    !if ( lge(string(1:6),"18.0.0") ) isParaMonteCompatibleCompiler = "true"

    VersionParts = splitStr(string, ".")
    if (isIntel) then
        MinVersion = INTEL_VERSION
    elseif (isGNU) then
        MinVersion = GNU_VERSION
    end if
    if  ( str2int(VersionParts(1)%record) > MinVersion(1) .or. &
        ( str2int(VersionParts(1)%record)== MinVersion(1) .and. str2int(VersionParts(2)%record) >= MinVersion(2) ) ) &
        isParaMonteCompatibleCompiler = "true"

    open(newunit=fileunit,file="isParaMonteCompatibleCompiler.tmp")
    write(fileunit,"(A)") isParaMonteCompatibleCompiler
    close(fileunit)

contains

    function splitStr(string,delimiter,nPart) result(Parts)

        implicit none
        character(len=*)    , intent(in)            :: string, delimiter
        integer(IK)         , intent(out), optional :: nPart
        character(len=:)    , allocatable           :: dummyStr
        type(CharVec_type)  , allocatable           :: Parts(:)
        integer(IK)                                 :: maxNumSplit
        integer(IK)                                 :: stringLen, delimLen, splitCounter, currentPos

        dummyStr  = string
        delimLen  = len(delimiter)
        stringLen = len(dummyStr)

        if (delimLen==0) then
            allocate(Parts(1))
            Parts(1)%record = string
            return
        end if

        maxNumSplit = 1 + stringLen / delimLen
        allocate(Parts(maxNumSplit))
        splitCounter = 1
        loopParseString: do
            if (stringLen<delimLen) then
                Parts(splitCounter)%record = dummyStr
                exit loopParseString
            elseif (stringLen==delimLen) then
                if (dummyStr==delimiter) then
                    Parts(splitCounter)%record = ""
                else
                    Parts(splitCounter)%record = dummyStr
                end if
                exit loopParseString
            elseif (dummyStr(1:delimLen)==delimiter) then
                dummyStr = dummyStr(delimLen+1:stringLen)
                stringLen = len(dummyStr)
                cycle loopParseString
            else
                currentPos = 2
                loopSearchString: do
                    if (dummyStr(currentPos:currentPos+delimLen-1)==delimiter) then
                        Parts(splitCounter)%record = dummyStr(1:currentPos-1)
                        if (currentPos+delimLen>stringLen) then
                            exit loopParseString
                        else
                            splitCounter = splitCounter + 1
                            dummyStr = dummyStr(currentPos+delimLen:stringLen)
                            stringLen = len(dummyStr)
                            cycle loopParseString
                        end if
                    else
                        currentPos = currentPos + 1
                        if (stringLen<currentPos+delimLen-1) then
                            Parts(splitCounter)%record = dummyStr
                            exit loopParseString
                        end if
                        cycle loopSearchString
                    end if
                end do loopSearchString
            end if
        end do loopParseString
        Parts = Parts(1:splitCounter)
        if (present(nPart)) nPart = splitCounter

    end function splitStr

    function str2int(str,iostat)
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        integer(IK)                     :: str2int
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2int
        else
            read(str,*) str2int
        endif
    end function str2int

end program getCompilerVersion
