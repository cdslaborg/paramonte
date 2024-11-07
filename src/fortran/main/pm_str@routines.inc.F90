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
!>  This file contains procedure implementations of [pm_str](@ref pm_str).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%
#if     alleq_ENABLED
        !%%%%%%%%%%%%

        integer(IK) :: i, lenStr1, lenStr2
        lenStr1 = len(str1, IK)
        lenStr2 = len(str2, IK)
        allEqual = logical(lenStr1 == lenStr2, LK)
        if (allEqual) then
            allEqual = logical(str1 == str2, LK)
        elseif (lenStr1 == 1_IK) then
            do i = 1, lenStr2
                allEqual = logical(str1(1:1) == str2(i:i), LK)
                if (.not. allEqual) return
            end do
        elseif (lenStr2 == 1_IK) then
            do i = 1, lenStr1
                allEqual = logical(str1(i:i) == str2(1:1), LK)
                if (.not. allEqual) return
            end do
        end if

        !%%%%%%%%%%%%%%%%%
#elif   getCharSeq_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        do concurrent(i = 1_IK : size(charVec, kind = IK))
            str(i:i) = charVec(i)
        end do

        !%%%%%%%%%%%%%%%%%
#elif   getCharVec_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        do concurrent(i = 1_IK : len(str, kind = IK))
            charVec(i) = str(i:i)
        end do

        !%%%%%%%%%%%%%%%%%%
#elif   isEndedWith_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     BSSK_ENABLED || PSSK_ENABLED
#define GET_VAL(X)X%val
#elif   SK_ENABLED
#define GET_VAL(X)X
#else
#error  "Unrecognized interface."
#endif
        endedWith = logical(GET_VAL(str)(max(1_IK, len(GET_VAL(str), IK) - len(GET_VAL(suffix), IK) + 1_IK) : len(GET_VAL(str), IK)) == GET_VAL(suffix), LK)
#undef  GET_VAL

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMinLoc_ENABLED || getMinVal_ENABLED || getMaxLoc_ENABLED || getMaxVal_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMinLoc_ENABLED || getMinVal_ENABLED
#define EXTREMUM_LOC minLoc
#define EXTREMUM_VAL minVal
#define COMPARES_WITH <
#elif   getMaxLoc_ENABLED || getMaxVal_ENABLED
#define EXTREMUM_LOC maxLoc
#define EXTREMUM_VAL maxVal
#define COMPARES_WITH >
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i
#if     getMinLoc_ENABLED || getMaxLoc_ENABLED
#define getLoc_ENABLED 1
        character(1,SKG) :: EXTREMUM_VAL
        EXTREMUM_LOC = 1_IK
#elif   !(getMinVal_ENABLED || getMaxVal_ENABLED)
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0_IK < len(str, IK), SK_": The condition `0 < len(str)` must hold. len(str) = "//getStr(len(str, IK)))
#if     Masked_ENABLED
        CHECK_ASSERTION(__LINE__, size(mask, 1, IK) == len(str, IK), SK_": The condition `size(mask) == len(str)` must hold. size(mask), len(str) = "//getStr([size(mask, 1, IK), len(str, IK)]))
#endif
#if     getLoc_ENABLED
        if (present(back)) then
            if (back) then
                EXTREMUM_VAL = str(len(str, kind = IK) : len(str, kind = IK))
                do i = len(str, kind = IK) - 1_IK, 1_IK, -1_IK
                    if  ( & ! LCOV_EXCL_LINE
#if                     Masked_ENABLED
                        mask(i) .and. & ! LCOV_EXCL_LINE
#endif
                        str(i:i) COMPARES_WITH EXTREMUM_VAL) then
                        EXTREMUM_VAL = str(i:i)
#if                     getLoc_ENABLED
                        EXTREMUM_LOC = i
#endif
                    end if
                end do
                return
            end if
        end if
#endif
        EXTREMUM_VAL = str(1:1)
        do i = 2_IK, len(str, kind = IK)
            if (& ! LCOV_EXCL_LINE
#if             Masked_ENABLED
                mask(i) .and. & ! LCOV_EXCL_LINE
#endif
                str(i:i) COMPARES_WITH EXTREMUM_VAL) then
                EXTREMUM_VAL = str(i:i)
#if             getLoc_ENABLED
                EXTREMUM_LOC = i
#endif
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%
#elif   getTrimmedTZ_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, lenstr, dotpos, tzstart
        logical(LK) :: isDigitSeq, leadingDigitIsMissing
        lenstr = len(str, IK)
        allocate(character(lenstr,SKG) :: strt)
        if (lenstr > 0_IK) then
            i = 1_IK
            j = 1_IK
            strt(j:j) = str(i:i)
            loopOverString: do
                if (str(i:i) == SKG_".") then
                    dotpos = i
                    loopOverDigitSeq: do
                        i = i + 1_IK
                        if (str(i:i) == SKG_"0") then
                            tzstart = i
                            loopOverTZ: do
                                if (i == lenstr) exit loopOverTZ
                                i = i + 1_IK
                                if (str(i:i) == SKG_"0") cycle loopOverTZ
                                isDigitSeq = logical(index(SKG_"123456789", str(i:i)) > 0, LK) ! isCharDigit(str(i:i))
                                if (isDigitSeq) then ! the zeros are significant, keep them.
                                    j = j + 1_IK
                                    strt(j:j+i-tzstart) = str(tzstart:i)
                                    j = j + i - tzstart
                                    if (i == lenstr) exit loopOverString
                                    cycle loopOverDigitSeq
                                end if
                                exit loopOverTZ
                            end do loopOverTZ
                            if (tzstart == dotpos + 1_IK) then
                                !   Ignore the all trailing zeros only if there is digit before the decimal point.
                                !   Otherwise, keep a zero as the leading digit before the decimal point.
                                leadingDigitIsMissing = logical(dotpos == 1_IK, LK)
                                if (.not. leadingDigitIsMissing) leadingDigitIsMissing = logical(index(SKG_"0123456789", str(dotpos-1:dotpos-1)) == 0, LK)
                                if (leadingDigitIsMissing) then
                                    strt(j : j + 1) = SKG_"0."
                                    j = j + 1_IK
                                    !j = j + 1_IK
                                    !strt(j:j+i-tzstart) = str(tzstart:i)
                                    !j = j + i - tzstart
                                end if
                            end if
                            if (i == lenstr) exit loopOverString
                        else
                            isDigitSeq = logical(index(SKG_"123456789", str(i:i)) > 0, LK) ! isCharDigit(str(i:i))
                        end if
                        j = j + 1_IK
                        strt(j:j) = str(i:i)
                        if (i == lenstr) exit loopOverString
                        if (isDigitSeq) cycle loopOverDigitSeq
                        cycle loopOverString
                    end do loopOverDigitSeq
                end if
                if (i == lenstr) exit loopOverString
                i = i + 1_IK
                j = j + 1_IK
                strt(j:j) = str(i:i)
            end do loopOverString
            strt = strt(1:j)
        else
            strt = SKG_""
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   getLenIndent_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenStr, lenPattern, ibeg, iend
        lenPattern = len(pattern, IK)
        lenStr = len(str, IK)
        if (lenStr == 0_IK .or. lenPattern == 0_IK .or. lenStr < lenPattern) then
            lenIndent = 0_IK
            return
        else
            ! Determine how many instances of pattern exist at the beginning of the string, and compute the indentation length.
            lenIndent = 0_IK
            iend = 0_IK
            do
                ibeg = iend
                iend = ibeg + lenPattern
                if (iend > lenStr) return
                if (str(ibeg + 1_IK : iend) /= pattern) return
                lenIndent = lenIndent + 1_IK
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%
#elif   getStrWrapped_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

       !logical(LK)                     :: isNewLine ! new line indicator
        logical(LK)                     :: isNewPar ! new paragraph indicator
        integer(IK)                     :: i, j, iwidth
        integer(IK)                     :: def_width, def_maxwidth
        integer(IK)                     :: widthCurrent, maxWidthCurrent
        integer(IK)                     :: lenStr, lenStrWrapped, lenStrWrappedTemp, countIndentPatternOld, countIndentPatternNew
        integer(IK)                     :: lenPrefix, lenIndent, lenIndentPattern, lenBreak, lenNewLine, lenLineFeed, lenPrefixIndent
        character(:,SKG), allocatable   :: strWrappedTemp, prefix_def, prefixIndent_def, indentPattern_def, break_def, newline_def, linefeed_def

        lenStr = len(str, kind = IK)
        if (lenStr == 0_IK) then
            strWrapped = SKG_""
            return
        end if

        ! Determine the prefix to each line.

        if (present(prefix)) then
            lenPrefix = len(prefix, IK)
            prefix_def = prefix
        else
            lenPrefix = 0_IK
            prefix_def = SKG_""
        end if

        ! Determine the indentation pattern.

        if (present(indent)) then
            lenIndentPattern = len(indent, IK)
            indentPattern_def = indent
        else
            lenIndentPattern = 1_IK
            indentPattern_def = SKG_" "
        end if

        ! Determine the allowed break characters.

        if (present(break)) then
            lenBreak = len(break, kind = IK)
            break_def = break
        else
            lenBreak = 1_IK
            break_def = SKG_" "
        end if

        ! Determine the newline.

        if (present(newline)) then
            lenNewLine = len(newline, IK)
            newline_def = newline
        else
            lenNewLine = 1_IK
            newline_def = new_line(SKG_"a")
        end if

        ! Determine the linefeed.

        if (present(linefeed)) then
            lenLineFeed = len(linefeed, IK)
            linefeed_def = linefeed
        else
            lenLineFeed = lenNewLine
            linefeed_def = newline_def
        end if

        ! Set the width of each line.

        if (present(width)) then
            def_width = max(width, 1_IK)
        else
            def_width = 132_IK
        end if

        ! Set the max-width of each line.

        if (present(maxwidth)) then
            def_maxwidth = maxwidth
        else
            def_maxwidth = (huge(0_IK) - 1_IK) / 2_IK
        end if
        def_maxwidth = max(def_maxwidth, def_width) ! `def_maxwidth > def_width` must always hold.

        ! Construct the wrapped array.

        lenStrWrapped = lenStr * 3_IK / 2_IK ! best guess
        allocate(character(lenStrWrapped,SKG) :: strWrappedTemp)

        if (lenIndentPattern == 0_IK) then
            lenPrefixIndent = lenPrefix
            prefixIndent_def = prefix_def
            maxWidthCurrent = def_maxwidth
            widthCurrent = def_width
            lenIndent = 0_IK
        end if

        i = 0_IK ! counter for `str`.
        j = 0_IK ! counter for `strWrapped`.
        isNewPar = .true._LK
       !isNewLine = .false._LK
        countIndentPatternOld = -1_IK
        loopOverLines: do

            ! Set the prefix + indentation for the new paragraph.

            if (isNewPar .and. lenIndentPattern > 0_IK) then
                countIndentPatternNew = getLenIndent(str(i + 1 : lenStr), indentPattern_def)
                if (countIndentPatternNew /= countIndentPatternOld) then ! This must be true at least the first time reaching here, which requires `countIndentPatternOld < 0_IK` for the first time.
                    lenIndent = lenIndentPattern * countIndentPatternNew
                    prefixIndent_def = prefix_def//repeat(indentPattern_def, countIndentPatternNew)
                    lenPrefixIndent = len(prefixIndent_def, IK)
                    countIndentPatternOld = countIndentPatternNew
                    if (lenIndentPattern * countIndentPatternNew > def_width) then
                        widthCurrent = lenIndentPattern * countIndentPatternNew + 1_IK
                        maxWidthCurrent = widthCurrent
                    else
                        maxWidthCurrent = def_maxwidth
                        widthCurrent = def_width
                    end if
                end if
            end if
            !write(*,*) "widthCurrent, def_width, maxWidthCurrent, def_maxwidth", widthCurrent, def_width, maxWidthCurrent, def_maxwidth

            ! Indent the line/paragraph.

            iwidth = lenIndent ! width of the current line.

            if (lenPrefixIndent > 0_IK) then
                if (j + lenPrefixIndent > lenStrWrapped) call resizeStrWrapped()
                strWrappedTemp(j + 1 : j + lenPrefixIndent) = prefixIndent_def
                j = j + lenPrefixIndent
                if (isNewPar) then
                    i = i + lenIndent ! The first indent is already taken into account in the above, so skip it.
                    if (i > lenStr) exit
                end if
            end if

            ! Wrap the string until reaching a newline.

            loopOverChars: do

                ! Check for the presence of a newline.

                if (i + lenNewLine <= lenStr) then
                    if (str(i + 1 : i + lenNewLine) == newline_def) then
                        if (j + lenLineFeed > lenStrWrapped) call resizeStrWrapped()
                        strWrappedTemp(j + 1 : j + lenLineFeed) = linefeed_def
                        j = j + lenLineFeed
                        i = i + lenNewLine
                        !isNewLine = .true._LK
                        if (i + lenNewLine <= lenStr) then ! Check for new paragraph
                            if (str(i + 1 : i + lenNewLine) == newline_def) then
                                if (j + lenPrefix + lenLineFeed > lenStrWrapped) call resizeStrWrapped()
                                strWrappedTemp(j + 1 : j + lenPrefix) = prefix_def
                                j = j + lenPrefix
                                strWrappedTemp(j + 1 : j + lenLineFeed) = linefeed_def
                                j = j + lenLineFeed
                                i = i + lenNewLine
                                isNewPar = .true._LK
                                cycle loopOverLines
                            end if
                            isNewPar = .false._LK
                        end if
                        cycle loopOverLines
                    end if
                end if

                ! Copy the next character to `strWrapped`.

                i = i + 1_IK; if (i > lenStr) exit loopOverLines
                j = j + 1_IK; if (j > lenStrWrapped) call resizeStrWrapped()
                strWrappedTemp(j : j) = str(i : i)
                iwidth = iwidth + 1_IK
                if (iwidth >= widthCurrent) then
                    if (lenBreak == 0_IK .or. index(break_def, str(i:i), kind = IK) > 0_IK .or. iwidth == maxWidthCurrent) then
                        if (j + lenLineFeed > lenStrWrapped) call resizeStrWrapped()
                        strWrappedTemp(j + 1 : j + lenLineFeed) = linefeed_def
                        j = j + lenLineFeed
                        !isNewLine = .false._LK
                        isNewPar = .false._LK
                        cycle loopOverLines
                    end if
                end if

            end do loopOverChars

        end do loopOverLines

        strWrapped = strWrappedTemp(1 : j)

    contains

        subroutine resizeStrWrapped()
            lenStrWrappedTemp = max(lenStrWrapped * 3_IK / 2_IK, lenStrWrapped + lenPrefixIndent)
            allocate(character(lenStrWrappedTemp,SKG) :: strWrapped)
            strWrapped(1_IK:lenStrWrapped) = strWrappedTemp(1_IK:lenStrWrapped)
            call move_alloc(strWrapped, strWrappedTemp)
#if         __GFORTRAN__
            if (allocated(strWrapped)) deallocate(strWrapped) ! Bypass gfortran <12 bug with Heap allocation.
#endif
            lenStrWrapped = lenStrWrappedTemp
        end subroutine
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  getLoc_ENABLED
#undef  COMPARES_WITH
#undef  EXTREMUM_LOC
#undef  EXTREMUM_VAL