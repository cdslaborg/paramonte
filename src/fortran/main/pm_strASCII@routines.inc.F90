!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_strASCII](@ref pm_strASCII).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(1,SKG), parameter :: SPACE_SKG = achar(32, SKG)
#if     getLocSpace_ENABLED
        do locSpace = 1_IK, len(str, kind = IK)
            if (str(locSpace:locSpace) /= SPACE_SKG) cycle
            return
        end do
        locSpace = 0_IK
#elif   getLocNonSpace_ENABLED
        do locNonSpace = 1_IK, len(str, kind = IK)
            if (str(locNonSpace:locNonSpace) == SPACE_SKG) cycle
            return
        end do
        locNonSpace = 0_IK
#elif   isCharDigit_ENABLED
        !integer(IK) :: j
        !charIsDigit = .false._LK
        !loopOverDigit: do j = 1_IK,10_IK
        !    if (chr == DIGIT_VEC_SK(j)) then
        !        charIsDigit = .true._LK
        !        exit loopOverDigit
        !    end if
        !end do loopOverDigit
        charIsDigit = SKG_"0" <= chr .and. chr <= SKG_"9"
#elif   isStrDigitAll_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharDigit(str(i:i))) cycle loopOverStr
                strIsDigitAll = .false._LK
                return
            end do loopOverStr
            strIsDigitAll = .true._LK
        else
            strIsDigitAll = .false._LK
        end if
#elif   isStrDigitAny_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharDigit(str(i:i))) then
                    strIsDigitAny = .true._LK
                    return
                end if
            end do loopOverStr
        end if
        strIsDigitAny = .false._LK
#elif   isStrDigit_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            StrIsNumeric(i) = isCharDigit(str(i:i))
        end do
#elif   isStrInteger_ENABLED
        integer(IK) :: i, lenStr
        strIsInteger = .false._LK
        lenStr = len(str, kind = IK)
        if (lenStr == 0_IK) return
        i = getLocNonSpace(str)
        if (i == 0_IK) return
        if (str(i:i) == SKG_"+" .or. str(i:i) == SKG_"-") i = i + 1_IK
        if (i > lenStr) return
        do
            if (i > lenStr) then
                strIsInteger = .true._LK
                return
            end if
            if (.not. isCharDigit(str(i:i))) exit
            i = i + 1_IK
        end do
        if (str(i:lenStr) == SPACE_SKG) strIsInteger = .true._LK
#elif   isStrComplex_ENABLED
        integer(IK)                     :: rebeg, refin
        integer(IK)                     :: imbeg, imfin
        integer(IK)                     :: i, lenStr
        strIsComplex = .false._LK
        lenStr = len(str, kind = IK)
        if (lenStr == 0_IK) return
        i = 0_IK
        do
            i = i + 1_IK
            if (i > lenStr) return
            if (str(i:i) == SPACE_SKG) cycle
            if (str(i:i) /= SKG_"(") return
            do
                i = i + 1_IK
                if (i > lenStr) return
                if (str(i:i) == SPACE_SKG) cycle
                rebeg = i
                do
                    i = i + 1_IK
                    if (i > lenStr) return
                    if (str(i:i) == SPACE_SKG .or. str(i:i) == SKG_",") exit
                    cycle
                end do
                refin = i - 1_IK
                if (str(i:i) == SPACE_SKG) then
                    do
                        i = i + 1_IK
                        if (i > lenStr) return
                        if (str(i:i) == SPACE_SKG) cycle
                        if (str(i:i) == SKG_",") exit
                    end do
                end if
                do
                    i = i + 1_IK
                    if (i > lenStr) return
                    if (str(i:i) == SPACE_SKG) cycle
                    exit
                end do
                imbeg = i
                do
                    i = i + 1_IK
                    if (i > lenStr) return
                    if (str(i:i) == SPACE_SKG .or. str(i:i) == SKG_")") exit
                end do
                imfin = i - 1_IK
                if (str(i:i) == SKG_")") then
                    do
                        i = i + 1_IK
                        if (i > lenStr) exit
                        if (str(i:i) == SPACE_SKG) cycle
                        return ! LCOV_EXCL_LINE
                    end do
                else ! str(i:i) == SPACE_SKG
                    do
                        i = i + 1_IK
                        if (i > lenStr) exit
                        if (str(i:i) == SPACE_SKG) cycle
                        if (str(i:i) /= SKG_")") return
                    end do
                end if
                strIsComplex = isStrReal(str(rebeg:refin)) .and. isStrReal(str(imbeg:imfin))
                return
            end do
        end do
#elif   isStrReal_ENABLED
        integer(IK) :: i, lenStr
        logical(LK) :: digitized
        strIsReal = .false._LK
        lenStr = len(str, kind = IK)
        if (lenStr == 0_IK) return
        i = getLocNonSpace(str)
        if (i == 0_IK) return ! str is all whitespace.
        !write(*,*) i, str(i:i)
        if (str(i:i) == SKG_"+" .or. str(i:i) == SKG_"-") then
            if (i == lenStr) return
        else
            i = i - 1_IK
        end if
        ! Skip any digits after sign
        digitized = .false._LK
        do
            i = i + 1_IK
            if (i > lenStr) then ! this never happens in the first round of loop.
                strIsReal = .true._LK
                return
            elseif (isCharDigit(str(i:i))) then
                digitized = .true._LK
                cycle
            end if
            exit
        end do
        !write(*,*) i, """"//str(i:i)//""""
        if (str(i:i) == SKG_".") then
            do
                i = i + 1_IK
                if (i > lenStr) then
                    strIsReal = digitized
                    !write(*,*) i, """"//str//""""
                    return
                elseif (isCharDigit(str(i:i))) then
                    digitized = .true._LK
                    cycle
                end if
                if (digitized) exit
                return
            end do
        end if
        !write(*,*) i, """"//str//""""
        if (str(i:i) == SKG_"e" .or. str(i:i) == SKG_"E" .or. str(i:i) == SKG_"d" .or. str(i:i) == SKG_"D") then
            i = i + 1_IK
            if (i > lenStr) return
            if (str(i:i) == SKG_"+" .or. str(i:i) == SKG_"-") then
                if (i == lenStr) return
                i = i + 1_IK
            end if
            if (.not. isCharDigit(str(i:i))) return
            do
                i = i + 1_IK
                if (i > lenStr) then ! This never happens on the first iteration.
                    strIsReal = .true._LK
                    return
                elseif (isCharDigit(str(i:i))) then
                    cycle
                end if
                exit ! LCOV_EXCL_LINE
            end do
        end if
        if (str(i:lenStr) == SPACE_SKG) strIsReal = .true._LK ! all the rest must be whitespace.
#elif   isStrNumber_ENABLED
        strIsNumber = isStrInteger(str) .or. isStrReal(str) .or. isStrComplex(str)
#elif   isCharUpper_ENABLED
        !charIsUpper = any(ALPHA_UPPER_VEC_SK == chr)
        charIsUpper = SKG_"A" <= chr .and. chr <= SKG_"Z"
#elif   isCharLower_ENABLED
        !charIsLower = any(ALPHA_LOWER_VEC_SK == chr)
        charIsLower = SKG_"a" <= chr .and. chr <= SKG_"z"
#elif   isStrUpperAll_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharUpper(str(i:i))) cycle
                strIsUpperAll = .false._LK
                return
            end do loopOverStr
            strIsUpperAll = .true._LK
        else
            strIsUpperAll = .false._LK
        end if
#elif   isStrLowerAll_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharLower(str(i:i))) cycle
                strIsLowerAll = .false._LK
                return
            end do loopOverStr
            strIsLowerAll = .true._LK
        else
            strIsLowerAll = .false._LK
        end if
#elif   isStrUpperAny_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharUpper(str(i:i))) then
                    strIsUpperAny = .true._LK
                    return
                end if
            end do loopOverStr
        end if
        strIsUpperAny = .false._LK
#elif   isStrLowerAny_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            loopOverStr: do i = 1_IK, len(str, kind = IK)
                if (isCharLower(str(i:i))) then
                    strIsLowerAny = .true._LK
                    return
                end if
            end do loopOverStr
        end if
        strIsLowerAny = .false._LK
#elif   isStrUpper_ENABLED
        integer(IK) :: i
        loopOverStr: do concurrent(i = 1_IK : len(str, kind = IK))
            StrIsUpper(i) = isCharUpper(str(i:i))
        end do loopOverStr
#elif   isStrLower_ENABLED
        integer(IK) :: i
        loopOverStr: do concurrent(i = 1_IK : len(str, kind = IK))
            StrIsLower(i) = isCharLower(str(i:i))
        end do loopOverStr
#elif   isCharAlphaNum_ENABLED
        charIsAlphaNum = (SKG_"0" <= chr .and. chr <= SKG_"9") .or. (SKG_"A" <= chr .and. chr <= SKG_"Z") .or. (SKG_"a" <= chr .and. chr <= SKG_"z")
#elif   isStrAlphaNumAll_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            do i = 1_IK, len(str, kind = IK)
                !if (any(ALPHANUM_VEC_SK == str(i:i))) cycle
                if ((SKG_"0" <= str(i:i) .and. str(i:i) <= SKG_"9") .or. (SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z")) cycle
                strIsAlphaNumAll = .false._LK
                return
            end do
            strIsAlphaNumAll = .true._LK
        else
            strIsAlphaNumAll = .false._LK
        end if
#elif   isStrAlphaNumAny_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            do i = 1_IK, len(str, kind = IK)
                if ((SKG_"0" <= str(i:i) .and. str(i:i) <= SKG_"9") .or. (SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z")) then
                    strIsAlphaNumAny = .true._LK
                    return
                end if
            end do
        end if
        strIsAlphaNumAny = .false._LK
#elif   isStrAlphaNum_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            StrIsAlphaNum(i) = logical((SKG_"0" <= str(i:i) .and. str(i:i) <= SKG_"9") .or. (SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z"), kind = LK)
        end do
#elif   isCharAlpha_ENABLED
        charIsAlpha = logical((SKG_"A" <= chr .and. chr <= SKG_"Z") .or. (SKG_"a" <= chr .and. chr <= SKG_"z"), LK)
#elif   isStrAlphaAll_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            do i = 1_IK, len(str, kind = IK)
                !if (any(ALPHANUM_VEC_SK == str(i:i))) cycle
                if ((SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z")) cycle
                strIsAlphaAll = .false._LK
                return
            end do
            strIsAlphaAll = .true._LK
        else
            strIsAlphaAll = .false._LK
        end if
#elif   isStrAlphaAny_ENABLED
        integer(IK) :: i
        if (len(str, kind = IK) > 0_IK) then
            do i = 1_IK, len(str, kind = IK)
                if ((SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z")) then
                    strIsAlphaAny = .true._LK
                    return
                end if
            end do
        end if
        strIsAlphaAny = .false._LK
#elif   isStrAlpha_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            StrIsAlpha(i) = logical((SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") .or. (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z"), kind = LK)
        end do
#elif   getStrUpper_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            if (SKG_"a" > str(i:i) .or. str(i:i) > SKG_"z") then
                strUpper(i:i) = str(i:i)
            else
                strUpper(i:i) = char(ichar(str(i:i), kind = IK) + UPPER_MINUS_LOWER_IK, kind = SKG)
            end if
        end do
#elif   getCharUpper_ENABLED
        if (SKG_"a" > chr .or. chr > SKG_"z") then
            chrUpper = chr
        else
            chrUpper = char(ichar(chr, kind = IK) + UPPER_MINUS_LOWER_IK, kind = SKG)
        end if
#elif   setCharUpper_ENABLED
        if (SKG_"a" <= chr .and. chr <= SKG_"z") chr = char(ichar(chr, kind = IK) + UPPER_MINUS_LOWER_IK, kind = SKG)
#elif   getCharLower_ENABLED
        if (SKG_"A" > chr .or. chr > SKG_"Z") then
            chrLower = chr
        else
            chrLower = char(ichar(chr, kind = IK) - UPPER_MINUS_LOWER_IK, kind = SKG)
        end if
#elif   setCharLower_ENABLED
        if (SKG_"A" <= chr .and. chr <= SKG_"Z") chr = char(ichar(chr, kind = IK) - UPPER_MINUS_LOWER_IK, kind = SKG)
#elif   setStrUpper_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            if (SKG_"a" <= str(i:i) .and. str(i:i) <= SKG_"z") str(i:i) = char(ichar(str(i:i), kind = IK) + UPPER_MINUS_LOWER_IK, kind = SKG)
        end do
#elif   getStrLower_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            if (SKG_"A" > str(i:i) .or. str(i:i) > SKG_"Z") then
                strLower(i:i) = str(i:i)
            else
                strLower(i:i) = char(ichar(str(i:i), kind = IK) - UPPER_MINUS_LOWER_IK, kind = SKG)
            end if
        end do
#elif   setStrLower_ENABLED
        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            if (SKG_"A" <= str(i:i) .and. str(i:i) <= SKG_"Z") str(i:i) = char(ichar(str(i:i), kind = IK) - UPPER_MINUS_LOWER_IK, kind = SKG)
        end do
#elif   getStrQuoted_ENABLED || setStrQuoted_ENABLED
        integer(IK) :: i, counter, pos, lenSeg, lenStr, lenStrQuoted, Loc(0:len(str))
        lenStr = len(str, IK)
        counter = 0_IK
        do i = 1_IK, lenStr
            if (str(i:i) == SKG_"""") then
                counter = counter + 1_IK
                Loc(counter) = i
            end if
        end do
        lenStrQuoted = lenStr + counter + 2_IK
        allocate(character(lenStrQuoted, SKG) :: strQuoted)
        strQuoted(1:1) = SKG_""""
        Loc(0) = 0_IK
        pos = 1_IK
        do i = 1_IK, counter
            lenSeg = Loc(i) - Loc(i-1_IK)
            strQuoted(pos + 1_IK : pos + lenSeg) = str(Loc(i-1_IK) + 1_IK : Loc(i))
            pos = pos + lenSeg + 1_IK
            strQuoted(pos : pos) = SKG_""""
        end do
        strQuoted(pos + 1_IK : lenStrQuoted - 1_IK) = str(Loc(counter) + 1 : lenStr)
        strQuoted(lenStrQuoted : lenStrQuoted) = SKG_""""
#elif   getAsciiFromEscaped_ENABLED || setAsciiFromEscaped_ENABLED
        integer(IK) :: i, j, lenStr, code
#if     getAsciiFromEscaped_ENABLED
        integer(IK) :: endloc
        lenStr = len(str, IK)
        allocate(character(lenStr,SKG) :: ascii)
#elif   setAsciiFromEscaped_ENABLED && Rep_ENABLED
        lenStr = len(str, IK)
#define ASCII str
#elif   setAsciiFromEscaped_ENABLED && New_ENABLED
        lenStr = len(str, IK)
        CHECK_ASSERTION(__LINE__, lenStr <= len(ascii,IK), SK_"@setAsciiFromEscaped(): The condition `len(str) <= len(ascii)` must hold. len(str), len(ascii) "//getStr([len(str,IK), len(ascii,IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        endloc = 0_IK
        i = 1_IK
        do
            if (i < lenStr) then
                code = -1_IK
                if (str(i:i) == SKG_"\") then
                    j = i + 1_IK
                    if (str(j:j) == SKG_"n") then
                        code = 10_IK
                    elseif (str(j:j) == SKG_"r") then
                        code = 13_IK
                    elseif (str(j:j) == SKG_"t") then
                        code = 9_IK
                    elseif (str(j:j) == SKG_"v") then
                        code = 11_IK
                    elseif (str(j:j) == SKG_"a") then
                        code = 7_IK
                    elseif (str(j:j) == SKG_"b") then
                        code = 8_IK
                    elseif (str(j:j) == SKG_"f") then
                        code = 12_IK
                    elseif (str(j:j) == SKG_"\") then
                        code = 92_IK
                    elseif (str(j:j) == SKG_"'") then
                        code = 39_IK
                    elseif (str(j:j) == SKG_'"') then
                        code = 34_IK
                    elseif (str(j:j) == SKG_"?") then
                        code = 63_IK
                    elseif (SKG_"0" <= str(j:j) .and. str(j:j) < SKG_"8") then ! is octal
                        do
                            if (j == lenStr) exit
                            j = j + 1_IK
                            if (SKG_"0" <= str(j:j) .and. str(j:j) < SKG_"8") cycle
                            j = j - 1_IK
                            exit
                        end do
                        code = getDecimal(str(i + 1_IK : j), 8_IK)
                    elseif (str(j:j) == SKG_"x") then ! is hex
                        do
                            if (j == lenStr) exit
                            j = j + 1_IK
                            if (isCharDigit(str(j:j)) .or. (SKG_"A" <= str(j:j) .and. str(j:j) < SKG_"G")) cycle
                            j = j - 1_IK
                            exit
                        end do
                        if (j > i + 1_IK) code = getDecimal(str(i + 2_IK : j), 16_IK)
                    elseif (str(j:j) == SKG_"u") then ! is UTF-8 four digit hex
#define SET_ASCII_CODE(STR_OFFSET)  \
do; \
if (j == lenStr) exit; \
j = j + 1_IK; \
if (isCharDigit(str(j:j)) .or. (SKG_"A" <= str(j:j) .and. str(j:j) < SKG_"G")) then; \
if (j < i + STR_OFFSET) cycle; \
exit; \
end if; \
j = j - 1_IK; \
exit; \
end do; \
if (j == i + STR_OFFSET) code = getDecimal(str(i + 2_IK : j), 16_IK);
                        SET_ASCII_CODE(5_IK) ! fpp
                    elseif (str(j:j) == SKG_"U") then ! is UTF-8 four digit hex
                        SET_ASCII_CODE(9_IK) ! fpp
#undef                  SET_ASCII_CODE
                    end if
                else
                    j = i
                end if
                endloc = endloc + 1_IK
                if (code < 0_IK .or. code > 127_IK) then
                    ASCII(endloc : endloc + j - i) = str(i:j) ! fpp
                    endloc = endloc + j - i
                else
                    ASCII(endloc : endloc) = achar(code, kind = SKG) ! fpp
                end if
                i = j + 1_IK
            else
                if (i == lenStr) then
                    endloc = endloc + 1_IK
                    ASCII(endloc : endloc) = str(i:i)
                end if
#if             getAsciiFromEscaped_ENABLED
                ASCII = ASCII(1:endloc) ! fpp
#endif
                CHECK_ASSERTION(__LINE__, endloc <= lenStr, SK_"The condition `endloc <= lenStr` must hold. endloc, lenStr = "//getStr([endloc, lenStr]))
                return
            end if
        end do
        error stop "Internal library error occurred. The procedure should not reach this line." ! LCOV_EXCL_LINE
#else
#error  "Unrecognized interface."
#endif
#undef  ASCII