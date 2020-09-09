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

#if defined MEXPRINT_ENABLED
#include "fintrf.h"
#endif

submodule (Decoration_mod) Routines_mod

    implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function constructDecoration(tabStr,symbol,text,List) result(Decoration)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDecoration
#endif
        use JaggedArray_mod, only: CharVec_type
        implicit none
        character(*), intent(in), optional          :: tabStr, symbol, text
        type(CharVec_type), intent(in), optional    :: List
        type(Decoration_type) :: Decoration
        if (present(tabStr)) then
            Decoration%tab = tabStr
        else
            Decoration%tab = TAB
        end if
        if (present(symbol)) then
            Decoration%symbol = symbol
        else
            Decoration%symbol = STAR
        end if
        if (present(text)) Decoration%text = text
        if (present(List)) Decoration%List = List
    end function constructDecoration 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine writeDecoratedText(text,symbol,width,thicknessHorz,thicknessVert,marginTop,marginBot,outputUnit,newLine)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: writeDecoratedText
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)            :: text
        character(*), intent(in), optional  :: symbol,newLine
        integer(IK) , intent(in), optional  :: width,thicknessHorz,thicknessVert,marginTop,marginBot
        integer(IK) , intent(in), optional  :: outputUnit
        integer(IK)                         :: thicknessVertDefault
        if (present(thicknessVert)) then
            thicknessVertDefault = thicknessVert
        else
            thicknessVertDefault = DECORATION_THICKNESS_VERT
        end if
        if (present(newLine)) then
            call writeDecoratedList( getListOfLines(text,newLine) &
                                   , symbol, width, thicknessHorz, thicknessVert, marginTop, marginBot, outputUnit )
        else
            call write(outputUnit,marginTop,0,thicknessVertDefault, drawLine(symbol,width) )
            call write(outputUnit,0,0,1, sandwich(text,symbol,width,thicknessHorz) )
            call write(outputUnit,0,marginBot,thicknessVertDefault, drawLine(symbol,width) )
        end if
    end subroutine writeDecoratedText

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine writeDecoratedList(List,symbol,width,thicknessHorz,thicknessVert,marginTop,marginBot,outputUnit)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: writeDecoratedList
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK
        implicit none
        type(CharVec_type), allocatable , intent(in)    :: List(:)
        character(*)        , intent(in), optional      :: symbol
        integer(IK)         , intent(in), optional      :: width,thicknessHorz,thicknessVert,marginTop,marginBot
        integer(IK)         , intent(in), optional      :: outputUnit
        integer(IK)                                     :: i
        integer(IK)                                     :: thicknessVertDefault
        if (present(thicknessVert)) then
            thicknessVertDefault = thicknessVert
        else
            thicknessVertDefault = DECORATION_THICKNESS_VERT
        end if
        call write(outputUnit,marginTop,0,thicknessVertDefault, drawLine(symbol,width) )
        do i = 1,size(List)
            call write(outputUnit,0,0,1, sandwich(List(i)%record,symbol,width,thicknessHorz) )
        end do
        call write(outputUnit,0,marginBot,thicknessVertDefault, drawLine(symbol,width) )
    end subroutine writeDecoratedList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pure module function drawLine(symbol,width) result(line)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: drawLine
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in), optional  :: symbol
        integer(IK), intent(in) , optional  :: width

        character(:), allocatable           :: line
        integer(IK)                         :: decorationWidth, symbolLen, symbolIndex, i
        character(:), allocatable           :: decorationSymbol

        if (present(symbol)) then
            if (len(symbol)<1) then
                decorationSymbol = " "
            else
                decorationSymbol = symbol
            end if
        else
            decorationSymbol = STAR
        end if
        symbolLen = len(decorationSymbol)

        if (present(width)) then
            decorationWidth = width
        else
            decorationWidth = DECORATION_WIDTH
        end if

        symbolIndex = 1
        allocate(character(decorationWidth) :: line)
        do i=1,decorationWidth
            line(i:i) = decorationSymbol(symbolIndex:symbolIndex)
            symbolIndex = symbolIndex + 1
            if (symbolIndex>symbolLen) symbolIndex = 1
        end do

    end function drawLine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module function sandwich(text,symbol,width,thicknessHorz) result(sandwichedText)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sandwich
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in), optional  :: text, symbol
        integer(IK), intent(in) , optional  :: width,thicknessHorz
        character(:), allocatable           :: sandwichedText
        integer(IK)                         :: decorationWidth, decorationThicknessHorz
        character(:), allocatable           :: decorationText, decorationSymbol
        integer(IK)                         :: i,decorationTextLen, symbolLen, symbolIndex, leftLimit, rightLimit
        integer(IK)                         :: sandwichedTextStart,decorationTextStart,decorationTextLenCounter

        if (present(symbol)) then
            decorationSymbol = symbol
        else
            decorationSymbol = STAR
        end if

        if (present(width)) then
            decorationWidth = width
        else
            decorationWidth = DECORATION_WIDTH
        end if

        if (present(thicknessHorz)) then
            decorationThicknessHorz = thicknessHorz
        else
            decorationThicknessHorz = DECORATION_THICKNESS_HORZ
        end if

        if (present(text)) then
            decorationText = trim(adjustl(text))
        else
            decorationText = ""
        end if

        if (decorationWidth<1) then
            sandwichedText = ""
            return
        end if

        allocate( character(decorationWidth) :: sandwichedText )
        decorationTextLen = len(decorationText)

        symbolLen = len(symbol)
        symbolIndex = 1
        leftLimit   = decorationThicknessHorz + 1
        rightLimit  = decorationWidth - decorationThicknessHorz + 1

        if (decorationTextLen<1) then
            do i = 1, decorationWidth
                if (i<leftLimit) then
                    sandwichedText(i:i) = decorationSymbol(symbolIndex:symbolIndex)
                    symbolIndex = symbolIndex + 1
                    if (symbolIndex>symbolLen) symbolIndex = 1
                elseif (i<rightLimit) then
                    sandwichedText(i:i) = " "
                    symbolIndex = symbolIndex + 1
                    if (symbolIndex>symbolLen) symbolIndex = 1
                else
                    sandwichedText(i:i) = decorationSymbol(symbolIndex:symbolIndex)
                    symbolIndex = symbolIndex + 1
                    if (symbolIndex>symbolLen) symbolIndex = 1
                end if
            end do
            return
        end if

        sandwichedTextStart = max( leftLimit  , ( decorationWidth - decorationTextLen ) / 2 + 1 )
        decorationTextStart = max( 1 , leftLimit - ( decorationWidth - decorationTextLen ) / 2 )
        decorationTextLenCounter = decorationTextStart

        do i=1,decorationWidth
            if (i<leftLimit) then
                sandwichedText(i:i) = decorationSymbol(symbolIndex:symbolIndex)
                symbolIndex = symbolIndex + 1
                if (symbolIndex>symbolLen) symbolIndex = 1
            elseif (i<rightLimit) then
                if (i<sandwichedTextStart) then
                    sandwichedText(i:i) = " "
                else if ( decorationTextLenCounter<=decorationTextLen ) then
                    sandwichedText(i:i) = decorationText(decorationTextLenCounter:decorationTextLenCounter)
                    decorationTextLenCounter = decorationTextLenCounter + 1
                else
                    sandwichedText(i:i) = " "
                end if
                symbolIndex = symbolIndex + 1
                if (symbolIndex>symbolLen) symbolIndex = 1
            else
                sandwichedText(i:i) = decorationSymbol(symbolIndex:symbolIndex)
                symbolIndex = symbolIndex + 1
                if (symbolIndex>symbolLen) symbolIndex = 1
            end if
        end do

        !! initialize empty container
        !do i=1,decorationWidth
        !    sandwichedText(i:i) = " "
        !end do

        !! add margin
        !do i=1,decorationThicknessHorz
        !    sandwichedText(i:i) = decorationSymbol
        !    sandwichedText(decorationWidth-i+1:decorationWidth-i+1) = decorationSymbol
        !end do

        !! add decorationText in between
        !sandwichedTextStart = max( decorationThicknessHorz , (decorationWidth-decorationTextLen)/2 )
        !sandwichedTextEnd   = min( decorationWidth , sandwichedTextStart + decorationTextLen - 1 )
        !decorationTextStart = 1 
        !decorationTextEnd   = sandwichedTextEnd - sandwichedTextStart + 1
        !sandwichedText(sandwichedTextStart:sandwichedTextEnd) = decorationText(decorationTextStart:decorationTextEnd)

    end function sandwich

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine write ( outputUnit    &
                            , marginTop     &
                            , marginBot     &
                            , count         &
                            , string        &
#if defined MEXPRINT_ENABLED
                            , advance       &
#endif
                            )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: write
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK, NLC
        implicit none
        integer(IK) , intent(in), optional  :: outputUnit
        integer(IK) , intent(in), optional  :: marginTop, marginBot, count
        character(*), intent(in), optional  :: string
        integer(IK)                         :: i, logFileUnit, thisManyTimes
#if defined MEXPRINT_ENABLED
        logical     , intent(in), optional  :: advance
        logical                             :: isStdout, advanceEnabled
        advanceEnabled = .true.; if (present(advance)) advanceEnabled = advance
#endif

        if (present(outputUnit)) then
#if defined MEXPRINT_ENABLED
            isStdout = output_unit == outputUnit
#endif
            logFileUnit = outputUnit
        else
#if defined MEXPRINT_ENABLED
            isStdout = .true.
#endif
            logFileUnit = output_unit
        end if

        if (present(marginTop)) then
            do i = 1, marginTop
#if defined MEXPRINT_ENABLED
                if (isStdout) then
                    call mexPrintf(NLC)
                else
                    write(logFileUnit,*)
                end if
#else
                write(logFileUnit,*)
#endif
            end do
        end if
    
        if (present(count)) then
            thisManyTimes = count
        else
            thisManyTimes = 1
        end if

        if (present(string)) then
            do i = 1, thisManyTimes
#if defined MEXPRINT_ENABLED
                if (isStdout) then
                    if (advanceEnabled) then
                        call mexPrintf(string//NLC)
                    else
                        call mexPrintf(string)
                    end if
                else
                    write(logFileUnit,"(g0)") string
                end if
#else
                write(logFileUnit,"(g0)") string
#endif
            end do
        elseif (.not. ( present(marginBot) .and. present(marginTop) ) ) then
            do i = 1, thisManyTimes
#if defined MEXPRINT_ENABLED
                if (isStdout) then
                    call mexPrintf(NLC)
                else
                    write(logFileUnit,*)
                end if
#else
                write(logFileUnit,*)
#endif
            end do
        end if
    
        if (present(marginBot)) then
            do i = 1, marginBot
#if defined MEXPRINT_ENABLED
                if (isStdout) then
                    call mexPrintf(NLC)
                else
                    write(logFileUnit,*)
                end if
#else
                write(logFileUnit,*)
#endif
            end do
        end if
  
    end subroutine write

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function wrapText(string,width,split,pad) result(ListOfLines)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: wrapText
#endif

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK

        implicit none

        character(*), parameter :: PROCEDURE_NAME = "@wrapText()"

        character(*), intent(in)            :: string
        integer(IK) , intent(in)            :: width
        character(*), intent(in), optional  :: split, pad
        type(CharVec_type), allocatable     :: ListOfLines(:)
        integer(IK)                         :: stringLen, splitLen, padLen, padLength, padLengthDynamic, newLineLen, oldLineLen
        integer(IK)                         :: istart, iend, numSplitEndLoc, counter, lineCount, indx, indxOld
        integer(IK), allocatable            :: IsEndOfSplitLoc(:), EndOfSplitLoc(:), EndOfLineLoc(:)
        logical                             :: isPadZone

        padLen = len(pad)
        splitLen = len(split)
        stringLen = len(string)
        if (stringLen==0) then
            allocate( ListOfLines(1) )
            ListOfLines(1)%record = ""
            return
        elseif (splitLen>=stringLen) then
            allocate( ListOfLines(1) )
            ListOfLines(1)%record = string
            return
        elseif (splitLen==0) then ! enforce wrapping at any character as necessary
            lineCount = stringLen / width + 1
            allocate( ListOfLines(lineCount) )
            do indx = 1, lineCount-1
                ListOfLines(indx)%record = string(width*(indx-1)+1:width*indx)
            end do
            ListOfLines(lineCount)%record = string(width*(lineCount-1)+1:stringLen)
            return
        end if

        ! get the initial pad size, and the locations of split ends.

        allocate(IsEndOfSplitLoc(stringLen))
        IsEndOfSplitLoc = 0_IK
        istart = 1_IK
        iend = istart + splitLen - 1_IK
        padLength = 0_IK
        isPadZone = .true.
        if (padLen==0) isPadZone = .false.
        blockFindSplit: do
            if (iend==stringLen) then
                IsEndOfSplitLoc(stringLen) = 1_IK
                exit blockFindSplit
            end if
            if ( isPadZone .and. mod(iend,padLen)==0 .and. string(istart:iend)==pad ) then
                padLength = iend
            else
                isPadZone = .false.
            end if
            if (string(istart:iend)==split) then
                IsEndOfSplitLoc(iend) = 1_IK
            else
                IsEndOfSplitLoc(iend) = 0_IK
            end if
            istart = istart + 1_IK
            iend = iend + 1_IK
        end do blockFindSplit

        ! create a vector of split-end indices

        numSplitEndLoc = sum(IsEndOfSplitLoc)
        if (numSplitEndLoc==0_IK) then
            allocate( ListOfLines(1) )
            ListOfLines(1)%record = string
            return
        else
            ! xxx: here goes another GFortran 7.3 bug: EndOfSplitLoc is assumed already allocated, despite the first appearance here.
            if (allocated(EndOfSplitLoc)) deallocate(EndOfSplitLoc)
            allocate(EndOfSplitLoc(numSplitEndLoc))
            counter = 0_IK
            do indx = 1,stringLen
                if (IsEndOfSplitLoc(indx)==1_IK) then
                    counter = counter + 1_IK
                    EndOfSplitLoc(counter) = indx
                end if
            end do
        end if
        EndOfSplitLoc = EndOfSplitLoc(1:counter)
        deallocate(IsEndOfSplitLoc)

        ! compute the number wrappings to be done

        ! xxx: here goes another GFortran 7.3 bug: EndOfLineLoc is assumed already allocated, despite the first appearance here.
        if (allocated(EndOfLineLoc)) deallocate(EndOfLineLoc)
        allocate( EndOfLineLoc(0:numSplitEndLoc+1) ) ! consider the maximum possible number of lines
        EndOfLineLoc = 0_IK
        lineCount = 0_IK
        padLengthDynamic = 0_IK ! first wrap does not need padding
        indxOld = 1_IK
        indx = 0_IK
        oldLineLen = -huge(oldLineLen)
        blockFindLine: do
            indx = indx + 1_IK
            if (indx>numSplitEndLoc) exit blockFindLine
            newLineLen = padLengthDynamic+EndOfSplitLoc(indx)-EndOfLineLoc(lineCount)
            if (newLineLen<=width) then
                oldLineLen = newLineLen
                cycle blockFindLine
            else
                ! swap the commented block with the uncommented to switch from better to good wrapping style.
                lineCount = lineCount + 1_IK
                if (indx-1_IK>indxOld) then ! ensure there is at least one split before the wrapping point
                    ! comment the following line to keep the max line length, strictly less than width (if possible).
                    if (width-oldLineLen>newLineLen-width) indx = indx + 1_IK ! removing the last token would make the line more beautiful
                    EndOfLineLoc(lineCount) = EndOfSplitLoc(indx-1)
                else
                    EndOfLineLoc(lineCount) = EndOfSplitLoc(indx)
                end if
                indxOld = indx
                padLengthDynamic = padLength
            end if
        end do blockFindLine

        ! add the remaining end of the string as a separate line

        if (EndOfLineLoc(lineCount)<stringLen .or. lineCount==0_IK) then
            lineCount = lineCount + 1_IK
            EndOfLineLoc(lineCount) = stringLen
        end if
        EndOfLineLoc = pack(EndOfLineLoc, mask=EndOfLineLoc/=0_IK)

        ! ensure the line count makes sense

        if ( lineCount /= size(EndOfLineLoc) ) then
            write(output_unit,"(*(g0,:,' '))")  MODULE_NAME // PROCEDURE_NAME // &
                                                ": Internal error occurred. lineCount /= size(EndOfLineLoc):", &
                                                lineCount, "/=", size(EndOfLineLoc), EndOfLineLoc
            write(output_unit,"(*(g0,:,' '))")  EndOfSplitLoc
            error stop
        end if

        ! construct the wrappings
        
        allocate( ListOfLines(lineCount) )
        indx = 1_IK
        ListOfLines(indx)%record = string(1:EndOfLineLoc(indx))
        do indx = 2, lineCount
            if ( padLength==0 .and. EndOfLineLoc(indx-1)+1>EndOfLineLoc(indx) ) then
                write(output_unit,"(*(g0,:,' '))")  MODULE_NAME // PROCEDURE_NAME // &
                                                    ": Fatal error occurred. " // &
                                                    "padLength==0 .and. EndOfLineLoc(indx-1)+1>EndOfLineLoc(indx) " // &
                                                    "for string: "
                write(output_unit,"(A)")            string
                error stop
            end if
            ListOfLines(indx)%record = string(1:padLength) // string(EndOfLineLoc(indx-1)+1:EndOfLineLoc(indx))
        end do

    end function wrapText

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getListOfLines(string,delimiter) result(ListOfLines)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getListOfLines
#endif
        use Constants_mod, only: IK
        implicit none
        character(len=*)  , intent(in)              :: string
        character(len=*)  , intent(in), optional    :: delimiter
        type(CharVec_type), allocatable             :: ListOfLines(:)
        character(len=:)  , allocatable             :: dumstr
        integer(IK)                                 :: stringLen, delimLen, delimLenMinusOne
        integer(IK)                                 :: maxNumSplit, counterString, counterLine, counterRecord
        logical                                     :: delimIsCStyle

        if (.not.present(delimiter)) then
            allocate(ListOfLines(1))
            ListOfLines(1)%record = string
            return
        end if

        stringLen = len(string)
        delimLen  = len(delimiter)
        delimLenMinusOne = delimLen - 1

        if (delimLen==0 .or. stringLen==0 .or. stringLen<delimLen) then
            allocate(ListOfLines(1))
            ListOfLines(1)%record = string
            return
        end if

        delimIsCStyle = delimLen==2 .and. delimiter=="\n"

        maxNumSplit = 1 + stringLen / delimLen
        allocate(ListOfLines(maxNumSplit))
        allocate( character(len=stringLen) :: dumstr )
        counterLine = 0
        counterRecord = 0
        counterString = 1
        loopParseString: do
            if (counterString+delimLenMinusOne>stringLen) then
                counterLine = counterLine + 1
                if (counterRecord==0) then
                    ListOfLines(counterLine)%record = string(counterString:stringLen)
                else
                    ListOfLines(counterLine)%record = dumstr(1:counterRecord) // string(counterString:stringLen)
                end if
                exit loopParseString
            end if
            if (string(counterString:counterString+delimLenMinusOne)==delimiter) then
                counterLine = counterLine + 1
                if (counterRecord==0) then
                    ListOfLines(counterLine)%record = ""
                else
                    ListOfLines(counterLine)%record = dumstr(1:counterRecord)
                    counterRecord = 0
                end if
                counterString = counterString + delimLen
                if (counterString>stringLen) then
                    counterLine = counterLine + 1
                    ListOfLines(counterLine)%record = ""
                    exit loopParseString
                end if
            elseif (delimIsCStyle .and. string(counterString:counterString)=="\") then
                counterString = counterString + 1
                counterRecord = counterRecord + 1
                dumstr(counterRecord:counterRecord) = "\"
                if (string(counterString:counterString)=="\") counterString = counterString + 1
            else
                counterRecord = counterRecord + 1
                dumstr(counterRecord:counterRecord) = string(counterString:counterString)
                counterString = counterString + 1
            end if
        end do loopParseString

        ListOfLines = ListOfLines(1:counterLine)

    end function getListOfLines

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule Routines_mod