#define SET_COPY \
call disp%show("From"); \
call disp%show( From ); \
call disp%show("To"); \
call disp%show( To ); \
call disp%show("call setCopyStrided(From, To, incf, inct)"); \
                call setCopyStrided(From, To, incf, inct); \
call disp%show("From"); \
call disp%show( From ); \
call disp%show("To"); \
call disp%show( To );

program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_arrayCopy, only: setCopyStrided
    use pm_distUnif, only: getUnifRand

    implicit none

    integer(IK) :: incf, inct

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy a string of characters of arbitrary kind in arbitrary stride.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_strASCII, only: ALPHA_STR_SK
        use pm_strASCII, only: ALPHA_LOWER_STR_SK
        use pm_strASCII, only: ALPHA_UPPER_STR_SK

        character(:, SK), allocatable :: From, To

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_STR_SK")
                        From = ALPHA_UPPER_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("incf = 0; inct = 2")
                        incf = 0; inct = 2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_STR_SK")
                        From = ALPHA_UPPER_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("incf = 2; inct = -2")
                        incf = 2; inct = -2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_STR_SK")
                        From = ALPHA_UPPER_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("incf = 2; inct = -2")
                        incf = 2; inct = -2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_STR_SK")
                        From = ALPHA_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("incf = 4; inct = 2")
                        incf = 4; inct = 2
        SET_COPY ! fpp
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy array of strings of the same length of arbitrary kind in arbitrary order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_strASCII, only: ALPHA_VEC_SK
        use pm_strASCII, only: ALPHA_LOWER_VEC_SK
        use pm_strASCII, only: ALPHA_UPPER_VEC_SK

        character(1, SK), allocatable :: From(:), To(:)

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_VEC_SK")
                        From = ALPHA_UPPER_VEC_SK
        call disp%show("To = ALPHA_LOWER_VEC_SK")
                        To = ALPHA_LOWER_VEC_SK
        call disp%show("incf = 2; inct = -2")
                        incf = 2; inct = -2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_LOWER_VEC_SK")
                        From = ALPHA_LOWER_VEC_SK
        call disp%show("To = ALPHA_LOWER_VEC_SK")
                        To = ALPHA_LOWER_VEC_SK
        call disp%show("incf = 2; inct = -2")
                        incf = 2; inct = -2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_VEC_SK")
                        From = ALPHA_VEC_SK
        call disp%show("To = ALPHA_LOWER_VEC_SK")
                        To = ALPHA_LOWER_VEC_SK
        call disp%show("incf = 4; inct = 2")
                        incf = 4; inct = 2
        SET_COPY ! fpp
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy array of integer values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer, allocatable :: From(:), To(:)

        call disp%skip()
        call disp%show("From = [13]")
                        From = [13]
        call disp%show("To = [(huge(0), incf = 1, 5)]")
                        To = [(huge(0), incf = 1, 5)]
        call disp%show("incf = 0; inct = 1")
                        incf = 0; inct = 1
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), incf = 1, 9)]")
                        To = [(huge(0), incf = 1, 9)]
        call disp%show("incf = 1; inct = 2")
                        incf = 1; inct = 2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), incf = 1, 9)]")
                        To = [(huge(0), incf = 1, 9)]
        call disp%show("incf = 1; inct = -2")
                        incf = 1; inct = -2
        SET_COPY ! fpp
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy array of logical values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical, allocatable :: From(:), To(:)

        call disp%skip()
        call disp%show("From = [.true.]")
                        From = [.true.]
        call disp%show("To = [(.false., incf = 1, 5)]")
                        To = [(.false., incf = 1, 5)]
        call disp%show("incf = 0; inct = 1")
                        incf = 0; inct = 1
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [.true., .true., .true., .true., .true.]")
                        From = [.true., .true., .true., .true., .true.]
        call disp%show("To = [(.false., incf = 1, 9)]")
                        To = [(.false., incf = 1, 9)]
        call disp%show("incf = 1; inct = 2")
                        incf = 1; inct = 2
        SET_COPY ! fpp
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy array of complex values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real, allocatable :: From(:), To(:)

        call disp%skip()
        call disp%show("From = cmplx([13], -[13])")
                        From = cmplx([13], -[13])
        call disp%show("To = cmplx([(huge(0), incf = 1, 5)], -[(huge(0), incf = 1, 5)])")
                        To = cmplx([(huge(0), incf = 1, 5)], -[(huge(0), incf = 1, 5)])
        call disp%show("incf = 0; inct = 1")
                        incf = 0; inct = 1
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])")
                        From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])
        call disp%show("To = cmplx([(huge(0), incf = 1, 9)], -[(huge(0), incf = 1, 9)])")
                        To = cmplx([(huge(0), incf = 1, 9)], -[(huge(0), incf = 1, 9)])
        call disp%show("incf = 1; inct = 2")
                        incf = 1; inct = 2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])")
                        From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])
        call disp%show("To = cmplx([(huge(0), incf = 1, 9)], -[(huge(0), incf = 1, 9)])")
                        To = cmplx([(huge(0), incf = 1, 9)], -[(huge(0), incf = 1, 9)])
        call disp%show("incf = 1; inct = -2")
                        incf = 1; inct = -2
        SET_COPY ! fpp
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy array of real values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real, allocatable :: From(:), To(:)

        call disp%skip()
        call disp%show("From = [13]")
                        From = [13]
        call disp%show("To = [(huge(0), incf = 1, 5)]")
                        To = [(huge(0), incf = 1, 5)]
        call disp%show("incf = 0; inct = 1")
                        incf = 0; inct = 1
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), incf = 1, 9)]")
                        To = [(huge(0), incf = 1, 9)]
        call disp%show("incf = 1; inct = 2")
                        incf = 1; inct = 2
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), incf = 1, 9)]")
                        To = [(huge(0), incf = 1, 9)]
        call disp%show("incf = 1; inct = -2")
                        incf = 1; inct = -2
        SET_COPY ! fpp
        call disp%skip()

    end block

end program example