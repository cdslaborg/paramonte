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

program getCompilerVersion

    use iso_fortran_env, only: compiler_version
    implicit none

    !integer    , allocatable   :: MinVersion(:)
    !integer    , parameter     :: GNU_VERSION(3) = [10,3,0]
    !integer    , parameter     :: INTEL_VERSION(3) = [21,0,0]
   !logical                     :: isIntel = .false., isGNU = .false.
   !integer                     :: fileunit, startindex, endindex, status, ipart, isep
   !character(:), allocatable   :: isParaMonteCompatibleCompiler, isGfortran10
    integer                     :: fileunit, ipart, isep, iostat
    character(:), allocatable   :: comver, outfile
    logical                     :: tripletFound

    type :: css_type
        character(:), allocatable :: val
    end type

    type(css_type), allocatable :: verlist(:)
    type(css_type), allocatable :: seplist(:)

    ! Set the output file name.

    if (command_argument_count() > 0) then
        allocate(character(1000) :: outfile)
        call get_command_argument(1, value = outfile, status = iostat)
        if (iostat == -1) error stop "Failed to fetch the output file name from the command line argument."
        outfile = trim(adjustl(outfile))
    else
        outfile = "getCompilerVersion.F90.tmp"
    end if
    comver = trim(adjustl(compiler_version()))

    ! Split the compiler version information to find the compiler version triplet.

    verlist = [css_type(comver)]
    seplist =   [ css_type(" ") &
                , css_type(".") &
                , css_type(":") &
                , css_type("(") &
                , css_type(")") &
                , css_type("v") &
                , css_type("V") &
                ]

    ipart = 0
    loopOverParts: do
        ipart = ipart + 1
        if (size(verlist) < ipart) exit loopOverParts
        loopOverSeps: do isep = 1, size(seplist)
            if (index(verlist(ipart)%val, seplist(isep)%val) == 0) cycle loopOverSeps
            verlist =   [ verlist(1 : ipart - 1) &
                        , getSplit(verlist(ipart)%val, seplist(isep)%val) &
                        , verlist(ipart + 1 :) &
                        ]
            ipart = ipart - 1
            exit loopOverSeps
        end do loopOverSeps
    end do loopOverParts

    ! The three consecutive all-digit parts are the compiler version.

    open(newunit = fileunit, file = outfile)
    loopTripletSearch: do ipart = 2, size(verlist) - 1
        tripletFound = isStrDigitAll(verlist(ipart)%val)
        tripletFound = tripletFound .and. isStrDigitAll(verlist(ipart - 1)%val)
        tripletFound = tripletFound .and. isStrDigitAll(verlist(ipart + 1)%val)
        if (.not. tripletFound) cycle loopTripletSearch
        write(fileunit, "(A)") verlist(ipart - 1)%val//"."//verlist(ipart)%val//"."//verlist(ipart + 1)%val
        exit loopTripletSearch
    end do loopTripletSearch
    close(fileunit)

    ! check for ParaMonteCompatibleCompiler

    !isParaMonteCompatibleCompiler = "false"

    !!if (lge(comver(1:3),"7.3") ) isParaMonteCompatibleCompiler = "true"
    !!if (lge(comver(1:6),"18.0.0") ) isParaMonteCompatibleCompiler = "true"

    !verlist = getSplit(comver, ".")
    !isGfortran10 = "false"
    !if (isIntel) then
    !    MinVersion = INTEL_VERSION
    !elseif (isGNU) then
    !    MinVersion = GNU_VERSION
    !    if (str2int(verlist(1)%val) > 9 ) isGfortran10 = "true"
    !end if
    !if  ( str2int(verlist(1)%val) > MinVersion(1) .or. &
    !    ( str2int(verlist(1)%val)== MinVersion(1) .and. str2int(verlist(2)%val) >= MinVersion(2) ) &
    !    ) isParaMonteCompatibleCompiler = "true"
    !
    !open(newunit=fileunit,file=outdir//"isParaMonteCompatibleCompiler.tmp")
    !write(fileunit,"(A)") isParaMonteCompatibleCompiler
    !close(fileunit)
    !
    !open(newunit=fileunit,file=outdir//"isGfortran10.tmp")
    !write(fileunit,"(A)") isGfortran10
    !close(fileunit)

contains

    pure elemental function isCharDigit(chr) result(charIsDigit)
        character(1), intent(in) :: chr
        logical :: charIsDigit
        charIsDigit = "0" <= chr .and. chr <= "9"
    end function

    pure elemental function isStrDigitAll(str) result(strIsDigitAll)
        character(1), intent(in) :: str
        logical :: strIsDigitAll
        integer :: i
        if (len(str) > 0) then
            loopOverStr: do i = 1, len(str)
                if (isCharDigit(str(i : i))) cycle loopOverStr
                strIsDigitAll = .false.
                return
            end do loopOverStr
            strIsDigitAll = .true.
        else
            strIsDigitAll = .false.
        end if
    end function

    pure function getSplit(comver, sep) result(parts)
        character(*)    , intent(in)    :: comver, sep
        character(:)    , allocatable   :: dummyStr
        type(css_type)  , allocatable   :: parts(:)
        integer                         :: maxNumSplit, stringLen, sepLen, splitCounter, currentPos
        dummyStr = comver
        sepLen = len(sep)
        stringLen = len(dummyStr)
        if (sepLen == 0) then
            allocate(parts(1))
            parts(1)%val = comver
            return
        end if
        maxNumSplit = 1 + stringLen / sepLen
        allocate(parts(maxNumSplit))
        splitCounter = 1
        loopParseString: do
            if (stringLen < sepLen) then
                parts(splitCounter)%val = dummyStr
                exit loopParseString
            elseif (stringLen == sepLen) then
                if (dummyStr == sep) then
                    parts(splitCounter)%val = ""
                else
                    parts(splitCounter)%val = dummyStr
                end if
                exit loopParseString
            elseif (dummyStr(1 : sepLen) == sep) then
                dummyStr = dummyStr(sepLen + 1 : stringLen)
                stringLen = len(dummyStr)
                cycle loopParseString
            else
                currentPos = 2
                loopSearchString: do
                    if (dummyStr(currentPos:currentPos + sepLen - 1) == sep) then
                        parts(splitCounter)%val = dummyStr(1 : currentPos - 1)
                        if (currentPos + sepLen > stringLen) then
                            exit loopParseString
                        else
                            splitCounter = splitCounter + 1
                            dummyStr = dummyStr(currentPos+sepLen:stringLen)
                            stringLen = len(dummyStr)
                            cycle loopParseString
                        end if
                    else
                        currentPos = currentPos + 1
                        if (stringLen < currentPos + sepLen - 1) then
                            parts(splitCounter)%val = dummyStr
                            exit loopParseString
                        end if
                        cycle loopSearchString
                    end if
                end do loopSearchString
            end if
        end do loopParseString
        parts = parts(1 : splitCounter)
    end function

    function str2int(str, iostat)
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        integer                         :: str2int
        if (present(iostat)) then
            iostat = 0
            read(str, *, iostat = iostat) str2int
        else
            read(str, *) str2int
        endif
    end function str2int

end