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

!>  \brief This module contains classes and procedures for decorating and outputting text.
!>  @author Amir Shahmoradi

module Decoration_mod

    use Constants_mod, only: IK
    use JaggedArray_mod, only: CharVec_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@Decoration_mod"

    integer(IK) , parameter :: DECORATION_WIDTH = 132
    integer(IK) , parameter :: DECORATION_THICKNESS_HORZ = 4
    integer(IK) , parameter :: DECORATION_THICKNESS_VERT = 1
    character(*), parameter :: STAR = "*"
    character(*), parameter :: TAB = "    "
    character(*), parameter :: INDENT = TAB // TAB

    character(*), parameter :: GENERIC_OUTPUT_FORMAT = "(*(g0,:,' '))"
    character(*), parameter :: GENERIC_TABBED_FORMAT = "('" // TAB // TAB // "',*(g0,:,' '))"

    ! ANSI string style

    integer(IK) , parameter             :: ANSI_ATTRIBUTE_LIST_LEN = 8
    integer(IK) , parameter             :: ANSI_COLOR_LIST_LEN = 16

    type, extends(CharVec_type)         :: Ansi_type
        character(:), allocatable       :: code
    end type

    type(Ansi_type), allocatable        :: mc_AnsiAttributeList(:)
    type(Ansi_type), allocatable        :: mc_AnsiForegroundColorList(:)
    type(Ansi_type), allocatable        :: mc_AnsiBackgroundColorList(:)

    !> The decoration class
    type :: decoration_type
        character(:), allocatable       :: tab
        character(:), allocatable       :: text
        character(:), allocatable       :: symbol
        type(CharVec_type), allocatable :: List(:)
    contains
        procedure, nopass :: write, writeDecoratedText, writeDecoratedList, wrapText, style
        !generic           :: write => writeDecoratedText, writeDecoratedList
    end type decoration_type

    interface decoration_type
        module procedure :: constructDecoration
    end interface decoration_type

    type :: wrapper_type
        type(CharVec_type), allocatable :: Line(:)
    contains
        procedure, nopass :: wrap => wrapText
    end type wrapper_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function constructDecoration(tabStr,symbol,text,List) result(Decoration)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDecoration
#endif
        use JaggedArray_mod, only: CharVec_type
        implicit none
        character(*), intent(in), optional          :: tabStr, symbol, text
        type(CharVec_type), intent(in), optional    :: List
        type(Decoration_type)                       :: Decoration
    end function constructDecoration
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
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
    end subroutine writeDecoratedText
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
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
    end subroutine writeDecoratedList
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    pure module function drawLine(symbol,width) result(line)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: drawLine
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in), optional  :: symbol
        integer(IK), intent(in) , optional  :: width
        character(:), allocatable           :: line
    end function drawLine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    pure module function sandwich(text,symbol,width,thicknessHorz) result(sandwichedText)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sandwich
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in), optional  :: text, symbol
        integer(IK), intent(in) , optional  :: width,thicknessHorz
        character(:), allocatable           :: sandwichedText
    end function sandwich
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
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
        use Constants_mod, only: IK
        implicit none
        integer(IK) , intent(in), optional  :: outputUnit
        integer(IK) , intent(in), optional  :: marginTop, marginBot, count
        character(*), intent(in), optional  :: string
#if defined MEXPRINT_ENABLED
        logical     , intent(in), optional  :: advance
#endif
    end subroutine write
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function getListOfLines(string,delimiter) result(ListOfLines)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getListOfLines
#endif
        implicit none
        character(len=*)  , intent(in)              :: string
        character(len=*)  , intent(in), optional    :: delimiter
        type(CharVec_type), allocatable             :: ListOfLines(:)
    end function getListOfLines
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function wrapText(string,width,split, pad) result(ListOfLines)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: wrapText
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)            :: string
        integer(IK) , intent(in)            :: width
        character(*), intent(in), optional  :: split, pad
        type(CharVec_type), allocatable     :: ListOfLines(:)
    end function wrapText
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    interface
!    module function style(string, attr, clbg, clfg, isUnixShell) result(modifiedString)
!#if defined DLL_ENABLED && !defined CFI_ENABLED
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGenericFormat
!#endif
!        implicit none
!        character(*), intent(in)            :: string
!        character(*), intent(in), optional  :: attr, clbg, clfg
!        logical, intent(in), optional       :: isUnixShell
!        character(:), allocatable           :: modifiedString
!    end function style
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Wrap the input text with ANSI/VT100 Control sequences to stylize and color the string.
    !>
    !> @param[in]   string      : The input string to be stylized.
    !> @param[in]   attr        : The requested attribute (optional). It can be:
    !>                            + `"off"`         : All attributes off (0).
    !>                            + `"bold"`        : Boldface text (1)
    !>                            + `"bright"`      : Bright text (1)
    !>                            + `"dim"`         : Dimmed text (2)
    !>                            + `"underlined"`  : Underlined text (4)
    !>                            + `"blinking"`    : Blinking text (5)
    !>                            + `"reverse"`     : Reversed attributes text (7)
    !>                            + `"hidden"`      : Hidden text (8)
    !>                            The **default** value is "off".
    !> @param[in]   clfg        : The Foreground color of the text (optional, see below for possible colors).
    !> @param[in]   clbg        : The Background color of the text (optional, see below for possible colors).
    !>
    !> \return
    !> `modifiedString` : The output string wrapped with the requested style and coloring.
    !>
    !> \remark
    !> Possible Foreground / Background colors are the following:
    !> + `"black"`          : (30) / (40)
    !> + `"red"`            : (31) / (41)
    !> + `"green"`          : (32) / (42)
    !> + `"yellow"`         : (33) / (43)
    !> + `"blue"`           : (34) / (44)
    !> + `"magenta"`        : (35) / (45)
    !> + `"cyan"`           : (36) / (46)
    !> + `"light gray"`     : (37) / (47)
    !> + `"dark gray"`      : (90) / (100)
    !> + `"light red"`      : (91) / (101)
    !> + `"light green"`    : (92) / (102)
    !> + `"light yellow"`   : (93) / (103)
    !> + `"light blue"`     : (94) / (104)
    !> + `"light magenta"`  : (95) / (105)
    !> + `"light cyan"`     : (96) / (106)
    !> + `"white"`          : (97) / (107)
    !
    !  For more information, see: https://misc.flogisoft.com/bash/tip_colors_and_formatting
    function style(string, attr, clfg, clbg) result(modifiedString)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenericFormat
#endif
        use Constants_mod, only: ESC
        implicit none
        character(*), intent(in)            :: string
        character(*), intent(in), optional  :: attr, clfg, clbg
        character(:), allocatable           :: modifiedString
        integer                             :: i

        ! cache the attribute / color list codes

        if (.not.allocated(mc_AnsiAttributeList)) then
            allocate(mc_AnsiAttributeList(ANSI_ATTRIBUTE_LIST_LEN))
            mc_AnsiAttributeList(1)%record = "off"          ; mc_AnsiAttributeList(1)%code = "0"
            mc_AnsiAttributeList(2)%record = "bold"         ; mc_AnsiAttributeList(2)%code = "1"
            mc_AnsiAttributeList(3)%record = "bright"       ; mc_AnsiAttributeList(3)%code = "1"
            mc_AnsiAttributeList(4)%record = "dim"          ; mc_AnsiAttributeList(4)%code = "2"
            mc_AnsiAttributeList(5)%record = "underlined"   ; mc_AnsiAttributeList(5)%code = "4"
            mc_AnsiAttributeList(6)%record = "blinking"     ; mc_AnsiAttributeList(6)%code = "5"
            mc_AnsiAttributeList(7)%record = "reverse"      ; mc_AnsiAttributeList(7)%code = "7"
            mc_AnsiAttributeList(8)%record = "hidden"       ; mc_AnsiAttributeList(8)%code = "8"
        end if

        if (.not.allocated(mc_AnsiForegroundColorList)) then
            allocate(mc_AnsiForegroundColorList(ANSI_COLOR_LIST_LEN))
            mc_AnsiForegroundColorList(1) %record = "black"         ; mc_AnsiForegroundColorList(1) %code = "30"
            mc_AnsiForegroundColorList(2) %record = "red"           ; mc_AnsiForegroundColorList(2) %code = "31"
            mc_AnsiForegroundColorList(3) %record = "green"         ; mc_AnsiForegroundColorList(3) %code = "32"
            mc_AnsiForegroundColorList(4) %record = "yellow"        ; mc_AnsiForegroundColorList(4) %code = "33"
            mc_AnsiForegroundColorList(5) %record = "blue"          ; mc_AnsiForegroundColorList(5) %code = "34"
            mc_AnsiForegroundColorList(6) %record = "magenta"       ; mc_AnsiForegroundColorList(6) %code = "35"
            mc_AnsiForegroundColorList(7) %record = "cyan"          ; mc_AnsiForegroundColorList(7) %code = "36"
            mc_AnsiForegroundColorList(8) %record = "light gray"    ; mc_AnsiForegroundColorList(8) %code = "37"
            mc_AnsiForegroundColorList(9) %record = "dark gray"     ; mc_AnsiForegroundColorList(9) %code = "90"
            mc_AnsiForegroundColorList(10)%record = "light red"     ; mc_AnsiForegroundColorList(10)%code = "91"
            mc_AnsiForegroundColorList(11)%record = "light green"   ; mc_AnsiForegroundColorList(11)%code = "92"
            mc_AnsiForegroundColorList(12)%record = "light yellow"  ; mc_AnsiForegroundColorList(12)%code = "93"
            mc_AnsiForegroundColorList(13)%record = "light blue"    ; mc_AnsiForegroundColorList(13)%code = "94"
            mc_AnsiForegroundColorList(14)%record = "light magenta" ; mc_AnsiForegroundColorList(14)%code = "95"
            mc_AnsiForegroundColorList(15)%record = "light cyan"    ; mc_AnsiForegroundColorList(15)%code = "96"
            mc_AnsiForegroundColorList(16)%record = "white"         ; mc_AnsiForegroundColorList(16)%code = "97"
        end if

        if (.not.allocated(mc_AnsiBackgroundColorList)) then
            allocate(mc_AnsiBackgroundColorList(ANSI_COLOR_LIST_LEN))
            mc_AnsiBackgroundColorList(1) %record = "black"         ; mc_AnsiBackgroundColorList(1) %code = "40"
            mc_AnsiBackgroundColorList(2) %record = "red"           ; mc_AnsiBackgroundColorList(2) %code = "41"
            mc_AnsiBackgroundColorList(3) %record = "green"         ; mc_AnsiBackgroundColorList(3) %code = "42"
            mc_AnsiBackgroundColorList(4) %record = "yellow"        ; mc_AnsiBackgroundColorList(4) %code = "43"
            mc_AnsiBackgroundColorList(5) %record = "blue"          ; mc_AnsiBackgroundColorList(5) %code = "44"
            mc_AnsiBackgroundColorList(6) %record = "magenta"       ; mc_AnsiBackgroundColorList(6) %code = "45"
            mc_AnsiBackgroundColorList(7) %record = "cyan"          ; mc_AnsiBackgroundColorList(7) %code = "46"
            mc_AnsiBackgroundColorList(8) %record = "light gray"    ; mc_AnsiBackgroundColorList(8) %code = "47"
            mc_AnsiBackgroundColorList(9) %record = "dark gray"     ; mc_AnsiBackgroundColorList(9) %code = "100"
            mc_AnsiBackgroundColorList(10)%record = "light red"     ; mc_AnsiBackgroundColorList(10)%code = "101"
            mc_AnsiBackgroundColorList(11)%record = "light green"   ; mc_AnsiBackgroundColorList(11)%code = "102"
            mc_AnsiBackgroundColorList(12)%record = "light yellow"  ; mc_AnsiBackgroundColorList(12)%code = "103"
            mc_AnsiBackgroundColorList(13)%record = "light blue"    ; mc_AnsiBackgroundColorList(13)%code = "104"
            mc_AnsiBackgroundColorList(14)%record = "light magenta" ; mc_AnsiBackgroundColorList(14)%code = "105"
            mc_AnsiBackgroundColorList(15)%record = "light cyan"    ; mc_AnsiBackgroundColorList(15)%code = "106"
            mc_AnsiBackgroundColorList(16)%record = "white"         ; mc_AnsiBackgroundColorList(16)%code = "107"
        end if

        ! construct the escape sequence

        modifiedString = ESC//"[0"

        if (present(attr)) then
            do i = 1, ANSI_ATTRIBUTE_LIST_LEN
                if (attr==mc_AnsiAttributeList(i)%record) then
                    modifiedString = modifiedString // ";" // mc_AnsiAttributeList(i)%code
                    exit
                end if
            end do
        end if

        if (present(clfg)) then
            do i = 1, ANSI_COLOR_LIST_LEN
                if (clfg==mc_AnsiForegroundColorList(i)%record) then
                    modifiedString = modifiedString // ";" // mc_AnsiForegroundColorList(i)%code
                    exit
                end if
            end do
        end if

        if (present(clbg)) then
            do i = 1, ANSI_COLOR_LIST_LEN
                if (clbg==mc_AnsiBackgroundColorList(i)%record) then
                    modifiedString = modifiedString // ";" // mc_AnsiBackgroundColorList(i)%code
                    exit
                end if
            end do
        end if

        modifiedString = modifiedString // "m" // string // ESC // "[0m"

    end function style

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an IO Fortran format given the input characteristics.
    !> @param[in]   width       :   The width of the target IO record (optional, default = dynamically set).
    !> @param[in]   precision   :   The precision of the target IO record if it happens to be a real number (optional, default = dynamically set).
    !> @param[in]   delim       :   The delimiter of the target IO record if it happens multiple entries (optional, default = "").
    !> @param[in]   prefix      :   The prefix of the target IO record (optional, default = "").
    !>
    !> \return
    !> `formatStr` : The output format string to be used in IO.
    pure function getGenericFormat(width,precision,delim,prefix) result(formatStr)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGenericFormat
#endif
        ! generates IO format strings, primarily for use in the output report files of ParaMonte
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        integer(IK) , intent(in), optional  :: width
        integer(IK) , intent(in), optional  :: precision
        character(*), intent(in), optional  :: delim
        character(*), intent(in), optional  :: prefix
        character(:), allocatable           :: widthStr
        character(:), allocatable           :: formatStr
        character(:), allocatable           :: precisionStr
        character(:), allocatable           :: delimDefault

        widthStr = "0"; if (present(width)) widthStr = num2str(width)
        precisionStr = ""; if (present(precision)) precisionStr = "."//num2str(precision)
        delimDefault = ""; if (present(delim)) delimDefault = ",:,'"//delim//"'"
        formatStr = "*(g"//widthStr//precisionStr//delimDefault//"))"
        if (present(prefix)) then
            formatStr = "('" // prefix // "'," // formatStr
        else
            formatStr = "(" // formatStr
        end if

    end function getGenericFormat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Decoration_mod ! LCOV_EXCL_LINE