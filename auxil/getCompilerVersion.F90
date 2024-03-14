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
    !integer(IK) , parameter    :: INTEL_VERSION(3) = [21_IK,0_IK,0_IK]
    !integer(IK) , parameter    :: GNU_VERSION(3) = [10_IK,3_IK,0_IK]
    !integer(IK) , allocatable  :: MinVersion(:)
    integer, parameter          :: IK = kind(0)
    integer(IK)                 :: fileunit, startindex, endindex, status
    logical                     :: isIntel = .false.
    logical                     :: isGNU = .false.
    character(:), allocatable   :: string, outfile
    character(:), allocatable   :: isParaMonteCompatibleCompiler, isGfortran10

    type :: css_type
        character(:), allocatable :: record
    end type css_type
    type(css_type), allocatable :: VersionParts(:)

    if (command_argument_count() > 0) then
        allocate(character(1000) :: outfile)
        call get_command_argument(1, value = outfile, status = status)
        if (status == -1) error stop "Failed to fetch the output file name from the command line argument."
        outfile = trim(adjustl(outfile))
    else
        outfile = "getCompilerVersion.F90.tmp"
    end if
    string = trim(adjustl(compiler_version()))

    if (index(string,"GCC") > 0) then ! it is Gfortran
        startindex = index(string,"version") + 8
        isGNU = .true.
    elseif (index(string,"Intel") > 0) then ! it is Intel
        startindex = index(string,"Version") + 8
        endindex = index(string,"Build") - 2
        isIntel = .true.
    else ! skip
        startindex = len(string, kind = IK) + 1_IK
    end if

    ! find the end of the version number
    endindex = startindex
    do
        if (endindex >= len(string, kind = IK)) exit
        if (string(endindex + 1_IK : endindex + 1_IK) == " ") exit
        endindex = endindex + 1_IK
    end do

    ! this is hopefully the pure compiler version.

    string = trim(adjustl(string(startindex : endindex)))

    open(newunit = fileunit, file = outfile)
    write(fileunit, "(A)") string
    close(fileunit)

    ! check for ParaMonteCompatibleCompiler

    !isParaMonteCompatibleCompiler = "false"

    !!if (lge(string(1:3),"7.3") ) isParaMonteCompatibleCompiler = "true"
    !!if (lge(string(1:6),"18.0.0") ) isParaMonteCompatibleCompiler = "true"

    !VersionParts = splitStr(string, ".")
    !isGfortran10 = "false"
    !if (isIntel) then
    !    MinVersion = INTEL_VERSION
    !elseif (isGNU) then
    !    MinVersion = GNU_VERSION
    !    if (str2int(VersionParts(1)%record) > 9 ) isGfortran10 = "true"
    !end if
    !if  ( str2int(VersionParts(1)%record) > MinVersion(1) .or. &
    !    ( str2int(VersionParts(1)%record)== MinVersion(1) .and. str2int(VersionParts(2)%record) >= MinVersion(2) ) &
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

    function splitStr(string,sep,nPart) result(field)

        implicit none
        character(len=*)    , intent(in)            :: string, sep
        integer(IK)         , intent(out), optional :: nPart
        character(len=:)    , allocatable           :: dummyStr
        type(css_type)  , allocatable           :: field(:)
        integer(IK)                                 :: maxNumSplit
        integer(IK)                                 :: stringLen, sepLen, splitCounter, currentPos

        dummyStr  = string
        sepLen  = len(sep)
        stringLen = len(dummyStr)

        if (sepLen==0) then
            allocate(field(1))
            field(1)%record = string
            return
        end if

        maxNumSplit = 1 + stringLen / sepLen
        allocate(field(maxNumSplit))
        splitCounter = 1
        loopParseString: do
            if (stringLen<sepLen) then
                field(splitCounter)%record = dummyStr
                exit loopParseString
            elseif (stringLen==sepLen) then
                if (dummyStr==sep) then
                    field(splitCounter)%record = ""
                else
                    field(splitCounter)%record = dummyStr
                end if
                exit loopParseString
            elseif (dummyStr(1:sepLen)==sep) then
                dummyStr = dummyStr(sepLen+1:stringLen)
                stringLen = len(dummyStr)
                cycle loopParseString
            else
                currentPos = 2
                loopSearchString: do
                    if (dummyStr(currentPos:currentPos+sepLen-1)==sep) then
                        field(splitCounter)%record = dummyStr(1:currentPos-1)
                        if (currentPos+sepLen>stringLen) then
                            exit loopParseString
                        else
                            splitCounter = splitCounter + 1
                            dummyStr = dummyStr(currentPos+sepLen:stringLen)
                            stringLen = len(dummyStr)
                            cycle loopParseString
                        end if
                    else
                        currentPos = currentPos + 1
                        if (stringLen<currentPos+sepLen-1) then
                            field(splitCounter)%record = dummyStr
                            exit loopParseString
                        end if
                        cycle loopSearchString
                    end if
                end do loopSearchString
            end if
        end do loopParseString
        field = field(1:splitCounter)
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
