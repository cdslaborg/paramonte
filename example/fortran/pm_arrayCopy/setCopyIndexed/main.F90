#define SET_COPY \
call disp%show("indexF"); \
call disp%show( indexF ); \
call disp%show("indexT"); \
call disp%show( indexT ); \
call disp%show("From"); \
call disp%show( From ); \
call disp%show("To"); \
call disp%show( To ); \
call disp%show("call setCopyIndexed(From, To, indexF, indexT)"); \
                call setCopyIndexed(From, To, indexF, indexT); \
call disp%show("From"); \
call disp%show( From ); \
call disp%show("To"); \
call disp%show( To );

program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_arrayCopy, only: setCopyIndexed
    use pm_distUnif, only: getUnifRand
    use pm_arrayReverse, only: getReversed
    use pm_arrayRange, only: getRange

    implicit none

    integer(IK), allocatable :: i, indexF(:), indexT(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Copy a string of characters of arbitrary kind in arbitrary orders.")
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
        call disp%show("indexF = getRange(1, len(From)); indexT = getRange(len(To), 1, -1)")
                        indexF = getRange(1, len(From)); indexT = getRange(len(To), 1, -1)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_STR_SK")
                        From = ALPHA_UPPER_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("indexF = getRange(1, len(From), 2); indexT = getRange(len(To), 1, -2)")
                        indexF = getRange(1, len(From), 2); indexT = getRange(len(To), 1, -2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_UPPER_STR_SK")
                        From = ALPHA_UPPER_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("indexF = getRange(1, len(From), 2); indexT = getRange(len(To), 1, -2)")
                        indexF = getRange(1, len(From), 2); indexT = getRange(len(To), 1, -2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_STR_SK")
                        From = ALPHA_STR_SK
        call disp%show("To = ALPHA_LOWER_STR_SK")
                        To = ALPHA_LOWER_STR_SK
        call disp%show("indexF = getRange(1, len(From), 4); indexT = getRange(1, len(To), 2)")
                        indexF = getRange(1, len(From), 4); indexT = getRange(1, len(To), 2)
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
        call disp%show("indexF = getRange(1, size(From), 2); indexT = getRange(size(To), 1, -2)")
                        indexF = getRange(1, size(From), 2); indexT = getRange(size(To), 1, -2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_LOWER_VEC_SK")
                        From = ALPHA_LOWER_VEC_SK
        call disp%show("To = ALPHA_LOWER_VEC_SK")
                        To = ALPHA_LOWER_VEC_SK
        call disp%show("indexF = getRange(1, size(From), 2); indexT = getRange(size(To), 1, -2)")
                        indexF = getRange(1, size(From), 2); indexT = getRange(size(To), 1, -2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = ALPHA_VEC_SK")
                        From = ALPHA_VEC_SK
        call disp%show("To = ALPHA_LOWER_VEC_SK")
                        To = ALPHA_LOWER_VEC_SK
        call disp%show("indexF = getRange(1, size(From), 4); indexT = getRange(1, size(To), 2)")
                        indexF = getRange(1, size(From), 4); indexT = getRange(1, size(To), 2)
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
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), i = 1, 9)]")
                        To = [(huge(0), i = 1, 9)]
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)")
                        indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), i = 1, 9)]")
                        To = [(huge(0), i = 1, 9)]
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)")
                        indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)
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
        call disp%show("From = [.true., .true., .true., .true., .true.]")
                        From = [.true., .true., .true., .true., .true.]
        call disp%show("To = [(.false., i = 1, 9)]")
                        To = [(.false., i = 1, 9)]
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)")
                        indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)
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
        call disp%show("From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])")
                        From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])
        call disp%show("To = cmplx([(huge(0), i = 1, 9)], -[(huge(0), i = 1, 9)])")
                        To = cmplx([(huge(0), i = 1, 9)], -[(huge(0), i = 1, 9)])
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)")
                        indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])")
                        From = cmplx([1, 2, 3, 4, 5], -[1, 2, 3, 4, 5])
        call disp%show("To = cmplx([(huge(0), i = 1, 9)], -[(huge(0), i = 1, 9)])")
                        To = cmplx([(huge(0), i = 1, 9)], -[(huge(0), i = 1, 9)])
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)")
                        indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)
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
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), i = 1, 9)]")
                        To = [(huge(0), i = 1, 9)]
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)")
                        indexF = getRange(1, size(From)); indexT = getRange(1, size(To), 2)
        SET_COPY ! fpp
        call disp%skip()

        call disp%skip()
        call disp%show("From = [1, 2, 3, 4, 5]")
                        From = [1, 2, 3, 4, 5]
        call disp%show("To = [(huge(0), i = 1, 9)]")
                        To = [(huge(0), i = 1, 9)]
        call disp%show("indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)")
                        indexF = getRange(1, size(From)); indexT = getRange(size(To), 1, -2)
        SET_COPY ! fpp
        call disp%skip()

    end block

end program example