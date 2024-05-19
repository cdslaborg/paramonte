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
!>  This file contains the implementation details of the routines for converting a logical or number of different types and kinds to char.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 9:44 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !> The maximum possible length of the string output from the functions under the generic interface `getStr`.
        !> Note that `STRLENMAX = 127` is very generous given that,
        !>  +   The maximum length of a complex of kind `real128` including signs, exponentiation, comma, and parentheses is `93` characters.
        !>      `(-1.189731495357231765085759326628007016E+4932,-1.189731495357231765085759326628007016E+4932)`
        !>  +   The maximum length of a real of kind `real128` including signs and exponentiation is `44` characters.
        !>      `-1.18973149535723176508575932662800702E+4932`
        !>  +   The maximum length of an integer of kind `int64` including signs `20` characters.
        !>      `-9223372036854775807`
        character(*), parameter :: SEP = ", "
        integer(IK) , parameter :: SEPLEN = len(SEP, kind = IK)
        integer(IK) , parameter :: STRLENMAX = 127_IK
#if     CK_ENABLED
        character(*, SK), parameter :: FORMAT_SIGNED = SK_"(*('(',sp,g0,'"//SEP//SK_"',g0,')',:,'"//SEP//SK_"'))"
        character(*, SK), parameter :: FORMAT_UNSIGNED = SK_"(*('(',g0,'"//SEP//SK_"',g0,')',:,'"//SEP//SK_"'))"
#elif   IK_ENABLED || RK_ENABLED
        character(*, SK), parameter :: FORMAT_SIGNED = SK_"(*(sp,g0,:,'"//SEP//SK_"'))"
        character(*, SK), parameter :: FORMAT_UNSIGNED = SK_"(*(g0,:,'"//SEP//SK_"'))"
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, allocated(val%val), SK_"@getStr(): The condition `allocated(val%val)` must hold.")
        if (present(length)) then
            call setResized(str, length)
            if (present(format)) then
                write(str, format) val%val
            else
                CHECK_ASSERTION(__LINE__, len(val%val, IK) <= length, SK_"@getStr(): The condition `len(val%val) <= length` must hold. len(val%val), length = "//getStr([len(val%val, IK), length]))
                str(1: len(val%val, IK)) = val%val
            end if
        elseif (present(format)) then
            block
                integer(IK) :: iostat
                character(127, SK) :: iomsg
                do
                    write(str, format, iostat = iostat, iomsg = iomsg) val%val
                    if (iostat == 0_IK) then
                        exit
                    elseif (is_iostat_eor(iostat)) then
                        call setResized(str)
                        cycle
                    else
                        error stop MODULE_NAME//SK_"@getStr(): "//trim(iomsg) ! LCOV_EXCL_LINE
                    end if
                end do
            end block
            str = trim(str)
        else
            str = val%val
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, allocated(val%val), SK_"@setStr(): The condition `allocated(val%val)` must hold.")
        if (present(format)) then
            write(str, format) val%val
        else
            CHECK_ASSERTION(__LINE__, len(val%val, IK) <= len(str, IK), SK_"@setStr(): The condition `len(val%val) <= len(str)` must hold. len(val%val), len(str) = "//getStr([len(val%val, IK), len(str, IK)]))
            str = val%val
        end if
        length = len_trim(str, IK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenstr
        CHECK_ASSERTION(__LINE__, all([(allocated(val(i)%val), i = 1, size(val, 1, IK))]), SK_"@getStr(): The condition `all([(allocated(val(i)%val), i = 1, size(val, 1, IK))])` must hold.")
        if (present(length)) then
            call setResized(str, length)
            if (present(format)) then
                write(str, format) (val(i)%val, i = 1, size(val, 1, IK))
            else
                if (0_IK < size(val, kind = IK)) then
                    str = val(1)%val
                    do i = 2, size(val, 1, IK)
                        str = str//SEP//val(i)%val
                    end do
                else
                    str = SKG_""
                end if
                CHECK_ASSERTION(__LINE__, len(str, IK) <= length, SK_"@getStr(): The condition `len(str) <= length` must hold. len(str), length = "//getStr([len(str, IK), length]))
            end if
        elseif (present(format)) then
            block
                integer :: iostat
                character(127, SK) :: iomsg
                lenstr = -SEPLEN
                do i = 1, size(val, 1, IK)
                    lenstr = lenstr + SEPLEN + len(val(i)%val, IK)
                end do
                call setResized(str, lenstr)
                do
                    write(str, format, iostat = iostat, iomsg = iomsg) (val(i)%val, i = 1, size(val, 1, IK))
                    if (iostat == 0_IK) then
                        exit
                    elseif (is_iostat_eor(iostat)) then
                        call setResized(str)
                        cycle
                    else
                        error stop MODULE_NAME//SK_"@getStr(): "//trim(iomsg) ! LCOV_EXCL_LINE
                    end if
                end do
            end block
            str = trim(str)
        else
            if (0_IK < size(val, kind = IK)) then
                str = val(1)%val
                do i = 2, size(val, 1, IK)
                    str = str//SEP//val(i)%val
                end do
            else
                str = SKG_""
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        CHECK_ASSERTION(__LINE__, all([(allocated(val(i)%val), i = 1, size(val, 1, IK))]), SK_"@getStr(): The condition `all([(allocated(val(i)%val), i = 1, size(val, 1, IK))])` must hold.")
        if (present(format)) then
            write(str, format) (val(i)%val, i = 1, size(val, 1, IK))
        else
            write(str, "(*(a,:,'"//SEP//"'))") (val(i)%val, i = 1, size(val, 1, IK))
        end if
        length = len_trim(str, IK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, lenstr
        CHECK_ASSERTION(__LINE__, all([((allocated(val(i,j)%val), i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))]), SK_"@getStr(): The condition `all([((allocated(val(i,j)%val), i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))])` must hold.")
        if (present(length)) then
            call setResized(str, length)
            if (present(format)) then
                write(str, format) ((val(i,j)%val, i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))
            else
                if (0_IK < size(val, kind = IK)) then
                    str = val(1,1)%val
                    do i = 2, size(val, 1, IK)
                        str = str//SEP//val(i,1)%val
                    end do
                    do j = 2, size(val, 2, IK)
                        do i = 2, size(val, 1, IK)
                            str = str//SEP//val(i,j)%val
                        end do
                    end do
                else
                    str = SKG_""
                end if
                CHECK_ASSERTION(__LINE__, len(str, IK) <= length, SK_"@getStr(): The condition `len(str) <= length` must hold. len(str), length = "//getStr([len(str, IK), length]))
            end if
        elseif (present(format)) then
            block
                integer :: iostat
                character(127, SK) :: iomsg
                lenstr = -SEPLEN
                do j = 1, size(val, 2, IK)
                    do i = 1, size(val, 1, IK)
                        lenstr = lenstr + SEPLEN + len(val(i,j)%val, IK)
                    end do
                end do
                call setResized(str, lenstr)
                do
                    write(str, format, iostat = iostat, iomsg = iomsg) ((val(i,j)%val, i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))
                    if (iostat == 0_IK) then
                        exit
                    elseif (is_iostat_eor(iostat)) then
                        call setResized(str)
                        cycle
                    else
                        error stop MODULE_NAME//SK_"@getStr(): "//trim(iomsg) ! LCOV_EXCL_LINE
                    end if
                end do
            end block
            str = trim(str)
        else
            if (0_IK < size(val, kind = IK)) then
                str = val(1,1)%val
                do i = 2, size(val, 1, IK)
                    str = str//SEP//val(i,1)%val
                end do
                do j = 2, size(val, 2, IK)
                    do i = 2, size(val, 1, IK)
                        str = str//SEP//val(i,j)%val
                    end do
                end do
            else
                str = SKG_""
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStr_ENABLED && (BSSK_ENABLED || PSSK_ENABLED) && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j
        CHECK_ASSERTION(__LINE__, all([((allocated(val(i,j)%val), i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))]), SK_"@getStr(): The condition `all([((allocated(val(i,j)%val), i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))])` must hold.")
        if (present(format)) then
            if (0_IK < size(val, kind = IK)) write(str, format) ((val(i, j)%val, i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))
        else
            if (0_IK < size(val, kind = IK)) write(str, "(*(a,:,'"//SEP//"'))") ((val(i, j)%val, i = 1, size(val, 1, IK)), j = 1, size(val, 2, IK))
        end if
        length = len_trim(str, IK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getStr_ENABLED
#define SET_LENGTH(i)
#define CHECK_STR_LEN(LINE)
#elif   setStr_ENABLED
#define SET_LENGTH(i) length = i
#define CHECK_STR_LEN(LINE) \
CHECK_ASSERTION(LINE, len(str, IK) >= length, SK_"@setStr(): The condition `len(str, IK) >= length` must hold. len(str, IK), length = "//getStr([len(str, IK), length]))
#else
#error  "Unrecognized interface."
#endif
        if (present(format)) then

#if         getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str)
                write(str, format) val
            else
#if             SK_ENABLED && D0_ENABLED
                ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                allocate(character(len(val, kind = IK),SKO) :: str)
#elif           SK_ENABLED && (D1_ENABLED || D2_ENABLED)
                ! extra 2 allows for possible separator.
                ! Fortran standard: Upon running the write statement,
                ! the untouched section of the record is padded with blanks.
                allocate(character(size(val, kind = IK) * (len(val, kind = IK) + 2_IK),SKO) :: str)
#elif           (IK_ENABLED || LK_ENABLED || RK_ENABLED || CK_ENABLED) && D0_ENABLED
                ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                allocate(character(STRLENMAX,SKO) :: str)
#elif           (IK_ENABLED || LK_ENABLED || RK_ENABLED || CK_ENABLED) && (D1_ENABLED || D2_ENABLED)
                ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                allocate(character(size(val, kind = IK) * STRLENMAX,SKO) :: str)
#else
#error          "Unrecognized interface."
#endif
                if (len(str, IK) > 0_IK) then
                    block
                        character(127, SK) :: iomsg
                        integer(IK) :: iostat
                        do
                            write(str, format, iostat = iostat, iomsg = iomsg) val
                            if (iostat == 0_IK) then
                                exit
                            elseif (is_iostat_eor(iostat)) then
                                call setResized(str)
                                cycle
                            else
                                error stop MODULE_NAME//SK_"@getStr(): "//trim(iomsg) ! LCOV_EXCL_LINE
                            end if
                        end do
                    end block
                    str = trim(str)
                end if
            end if
#elif       setStr_ENABLED
            length = len(str, IK)
            if (length > 0_IK) then
                write(str, format) val
                do
                    if (length == 0_IK) exit
                    if (str(length:length) /= SKO_" ") exit
                    length = length - 1_IK
                end do
            end if
#else
#error      "Unrecognized interface."
#endif
            return

        else

            !%%%%%%%%%
#if         SK_ENABLED
            !%%%%%%%%%

#if         D0_ENABLED && getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str)
                CHECK_ASSERTION(__LINE__, length >= len(val, IK), SK_"@getStr(): The condition `length >= len(val)` must hold. length, len(val) = "//getStr([length, len(val, IK)]))
                str(1:len(val, IK)) = val
                return
            end if
            str = trim(val)
#elif       D0_ENABLED && setStr_ENABLED
            SET_LENGTH(len_trim(val, IK)) ! fpp
            str(1:length) = val(1:length)
            CHECK_STR_LEN(__LINE__) ! fpp
#elif       D1_ENABLED || D2_ENABLED
            if (size(val, kind = IK) > 0_IK) then
                call setStrFromStr(val, str, length)
                CHECK_STR_LEN(__LINE__) ! fpp
            else
#if             getStr_ENABLED
                if (present(length)) then
                    str = repeat(SKO_" ", length)
                    return
                end if
#elif           !setStr_ENABLED
#error          "Unrecognized interface."
#endif
                SET_LENGTH(0_IK) ! fpp
                str = SKO_""
            end if
            CHECK_STR_LEN(__LINE__) ! fpp
#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%
#elif       LK_ENABLED
            !%%%%%%%%%

#if         D0_ENABLED
#if         getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str)
                if (val) then
                    CHECK_ASSERTION(__LINE__, length >= 4_IK, SK_"@getStr(): The condition `length >= 4_IK` must hold. length = "//getStr(length))
                    str = SKO_"TRUE"
                else
                    CHECK_ASSERTION(__LINE__, length >= 5_IK, SK_"@getStr(): The condition `length >= 5_IK` must hold. length = "//getStr(length))
                    str = SKO_"FALSE"
                end if
                return
            end if
#elif       !setStr_ENABLED
#error      "Unrecognized interface."
#endif
            if (val) then
                SET_LENGTH(4_IK) ! fpp
                str = SKO_"TRUE"
            else
                SET_LENGTH(5_IK) ! fpp
                str = SKO_"FALSE"
            end if
            CHECK_STR_LEN(__LINE__) ! fpp
#elif       D1_ENABLED || D2_ENABLED
            if (size(val, kind = IK) > 1_IK) then
                call setStrFromLogical(val, str, length)
            elseif (size(val, kind = IK) == 1_IK) then
                if (any(val)) then
                    SET_LENGTH(4_IK) ! fpp
                    str = SKO_"TRUE"
                else
                    SET_LENGTH(5_IK) ! fpp
                    str = SKO_"FALSE"
                end if
            else
                SET_LENGTH(0_IK) ! fpp
                str = SKO_""
            end if
            CHECK_STR_LEN(__LINE__) ! fpp
#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif       IK_ENABLED || CK_ENABLED || RK_ENABLED
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str)
            else
#if             D0_ENABLED
                ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                allocate(character(STRLENMAX,SKO) :: str)
#elif           D1_ENABLED || D2_ENABLED
                ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                allocate(character(STRLENMAX * size(val, kind = IK),SKO) :: str)
#else
#error          "Unrecognized interface."
#endif
            end if
#elif       !setStr_ENABLED
#error      "Unrecognized interface."
#endif
            if (present(signed)) then
                if (signed) then
                    write(str, FORMAT_SIGNED) val
#if                 getStr_ENABLED
                    if (.not. present(length)) str = trim(str) ! Fortran standard: char is by default left-adjusted when dealing with internal files.
#elif               setStr_ENABLED
                    length = len(str, IK)
                    if (length == 0_IK) return
                    do
                        if (str(length : length) /= SKO_" ") exit
                        length = length - 1_IK
                    end do
                    CHECK_STR_LEN(__LINE__) ! fpp
#endif
                    return
                end if
            end if
            if (len(str, IK) > 0_IK) then
                write(str, FORMAT_UNSIGNED) val
#if             getStr_ENABLED
                ! Fortran standard: char is by default left-adjusted when dealing with internal files.
                if (.not. present(length)) str = trim(str)
#elif           setStr_ENABLED
                length = len(str, IK)
                if (length == 0_IK) return
                do
                    if (str(length : length) /= SKO_" ") exit
                    length = length - 1_IK
                end do
                CHECK_STR_LEN(__LINE__) ! fpp
#endif
            else
                SET_LENGTH(0_IK) ! fpp
            end if
#else
#error      "Unrecognized interface."
#endif

        end if

    contains

#if     SK_ENABLED && (D1_ENABLED || D2_ENABLED)
        PURE subroutine setStrFromStr(ValVec, str, length)
#if         getStr_ENABLED
            integer(IK)                                     :: endpos
            integer(IK)     , intent(in)    , optional      :: length
            character(:,SKO), intent(out)   , allocatable   :: str
#elif       setStr_ENABLED
#define     endpos length
            character(*,SKO), intent(out)                   :: str
            integer(IK)     , intent(out)                   :: length
#elif       !setStr_ENABLED
#error      "Unrecognized interface."
#endif
            character(*,SKG), intent(in)                    :: ValVec(*)
            integer(IK)                                     :: i, iend, sizeVal, lenVal, startpos
            sizeVal = size(val, kind = IK)
            lenVal = len(val, kind = IK)
#if         getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str) ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
            else
                allocate(character(sizeVal * (lenVal + SEPLEN) - SEPLEN,SKO) :: str) ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
            end if
#endif
            do endpos = lenVal, 1_IK, -1_IK
                if (ValVec(1)(endpos:endpos) /= SKG_" ") exit
            end do
            str(1:endpos) = ValVec(1)(1:endpos)
            do i = 2_IK, sizeVal
                startpos = endpos + 1_IK
                endpos = endpos + SEPLEN
                str(startpos : endpos) = SEP
                do iend = lenVal, 1_IK, -1_IK
                    if (ValVec(i)(iend:iend) /= SKG_" ") exit
                end do
                startpos = endpos + 1_IK
                endpos = endpos + iend
                str(startpos:endpos) = ValVec(i)(1:iend)
            end do
!#if         getStr_ENABLED
!            ! This condition cannot be readily verified because the right blanks are trimmed, leading `endpos` values are smaller than the expected value.
!            ! Nevertheless, any length error is typically well captured by Intel and gfortran compilers.
!            check_assertion(__LINE__, endpos == sizeVal * (lenVal + SEPLEN) - SEPLEN .or. (present(length) .and. endpos >= sizeVal * (lenVal + SEPLEN) - SEPLEN), \
!            SK_"@getStr(): The condition `length >= sizeVal * (lenVal + SEPLEN) - SEPLEN` must hold. length, ... = "//getStr([endpos, sizeVal * (lenVal + SEPLEN) - SEPLEN]))
!            str = str(1:endpos)
!#endif
        end subroutine
#elif   LK_ENABLED && (D1_ENABLED || D2_ENABLED)
        PURE subroutine setStrFromLogical(ValVec, str, length)
#if         getStr_ENABLED
            integer(IK)                                     :: endpos
            integer(IK)     , intent(in)    , optional      :: length
            character(:,SKO), intent(out)   , allocatable   :: str
#elif       setStr_ENABLED
#define     endpos length
            character(*,SKO), intent(out)                   :: str
            integer(IK)     , intent(out)                   :: length
#elif       !setStr_ENABLED
#error      "Unrecognized interface."
#endif
            logical(LKG)    , intent(in)                    :: ValVec(*)
            integer(IK)                                     :: lenStr, sizeVal, i, startpos
            sizeVal = size(val, kind = IK)
            lenStr = count(val, kind = IK)
            lenStr = lenStr * (4_IK + SEPLEN) + (sizeVal - lenStr) * (5_IK + SEPLEN) - SEPLEN
#if         getStr_ENABLED
            if (present(length)) then
                allocate(character(length,SKO) :: str) ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
                CHECK_ASSERTION(__LINE__, length >= lenStr, SK_"@getStr(): The input `length` argument must be sufficiently large such that the input `val` fits within the string buffer. length, lenStr."//getStr([length, lenStr]))
            else
                allocate(character(lenStr,SKO) :: str) ! Fortran standard: Upon running the write statement, the untouched section of the record is padded with blanks.
            end if
#endif
            if (ValVec(1)) then
                endpos = 4_IK
                str(1:endpos) = SKO_"TRUE"
            else
                endpos = 5_IK
                str(1:endpos) = SKO_"FALSE"
            end if
            do i = 2_IK, sizeVal
                startpos = endpos + 1_IK
                if (ValVec(i)) then
                    endpos = endpos + SEPLEN + 4_IK
                    str(startpos:endpos) = SEP//"TRUE"
                else
                    endpos = endpos + SEPLEN + 5_IK
                    str(startpos:endpos) = SEP//"FALSE"
                end if
            end do
        end subroutine
#endif
#undef  CHECK_STR_LEN
#undef  SET_LENGTH
#undef  endpos
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
