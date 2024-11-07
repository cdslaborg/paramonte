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

!>  \bug
!>  ESC should be ideally and portably defined as `achar(27, kind)`, once the gfortran 11 bug is resolved.
#define ESC achar(27)
!#define ESC ""

!>  \brief
!>  This module contains procedures and generic interfaces for styling strings according for display on DEC VT100 or compatible terminals.
!>
!>  \details
!>  The current implementation includes only a subset of the ANSI standard text attributes and colors that
!>  are recognized by virtually all terminals and terminal emulators. These styles are particularly recognized
!>  by the Bash and Windows CMD (>2016) terminals.
!>
!>  **Notational conventions used in this module**<br>
!>  -#  Any object name that begins with an **f** is related to **foreground** effects.
!>  -#  Any object name that begins with an **b** is related to **background** effects.
!>  -#  Any object name that begins with an **end** contains the string that ends the corresponding ANSI style.
!>
!>  \test
!>  [test_pm_strANSI](@ref test_pm_strANSI)
!>
!>  \bug
!>  \status \unresolved
!>  \source \gfortran{10-12}
!>  \desc
!>  The character `ESC` should be ideally and portably defined as `achar(27, kind)`.<br>
!>  However, \gfortran as of version 11 cannot handle this in component initializations.<br>
!>  \remedy
!>  The character `ESC` is defined a direct copy of the actual character via a preprocessor macro.<br>
!>
!>  \todo
!>  \pvhigh The gfortran bug must be resolved as soon as possible.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
module pm_strANSI

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_strANSI"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !   \warning
        !   Any modifications to the following constants must be also reflected in the derived type styleSeq_type below.

        character(4 , SK)   :: bold             = ESC//"[1m"    !<  \public The ANSI escape code that makes subsequent texts appear **bold/bright**.
        character(4 , SK)   :: bright           = ESC//"[1m"    !<  \public The ANSI escape code that makes subsequent texts appear **bold/bright**.
        character(4 , SK)   :: dim              = ESC//"[2m"    !<  \public The ANSI escape code that makes subsequent texts appear **dim**.
        character(4 , SK)   :: italic           = ESC//"[3m"    !<  \public The ANSI escape code that makes subsequent texts appear **italic**.
        character(4 , SK)   :: underlined       = ESC//"[4m"    !<  \public The ANSI escape code that makes subsequent texts appear **underlined**.
        character(4 , SK)   :: blinking         = ESC//"[5m"    !<  \public The ANSI escape code that makes subsequent texts appear **blinking** (not universally supported).
        character(4 , SK)   :: reverse          = ESC//"[7m"    !<  \public The ANSI escape code that makes subsequent texts appear **in reverse**.
        character(4 , SK)   :: hidden           = ESC//"[8m"    !<  \public The ANSI escape code that makes subsequent texts **disappear** (hidden but not removed!).
        character(4 , SK)   :: strike           = ESC//"[9m"    !<  \public The ANSI escape code that makes subsequent texts to have **strike-through** (not universally supported).

        character(4 , SK)   :: endbold          = ESC//"[1m"    !<  \public The ANSI escape code that **ends** the ANSI style **bold/bright**.
        character(4 , SK)   :: endbright        = ESC//"[1m"    !<  \public The ANSI escape code that **ends** the ANSI style **bold/bright**.
        character(4 , SK)   :: enddim           = ESC//"[2m"    !<  \public The ANSI escape code that **ends** the ANSI style **dim**.
        character(4 , SK)   :: enditalic        = ESC//"[3m"    !<  \public The ANSI escape code that **ends** the ANSI style **italic**.
        character(4 , SK)   :: endunderlined    = ESC//"[4m"    !<  \public The ANSI escape code that **ends** the ANSI style **underlined**.
        character(4 , SK)   :: endblinking      = ESC//"[5m"    !<  \public The ANSI escape code that **ends** the ANSI style **blinking**.
        character(4 , SK)   :: endreverse       = ESC//"[7m"    !<  \public The ANSI escape code that **ends** the ANSI style **reverse**.
        character(4 , SK)   :: endhidden        = ESC//"[8m"    !<  \public The ANSI escape code that **ends** the ANSI style **hidden**.
        character(4 , SK)   :: endstrike        = ESC//"[9m"    !<  \public The ANSI escape code that **ends** the ANSI style **strike-through**.
        character(4 , SK)   :: endall           = ESC//"[0m"    !<  \public The ANSI escape code that **resets** all previously specified ANSI text styles (same as [](@ref pm_strANSI::reset)).
        character(4 , SK)   :: reset            = ESC//"[0m"    !<  \public The ANSI escape code that **resets** all previously specified ANSI text styles (same as [](@ref pm_strANSI::endall)).

        character(5 , SK)   :: fblack           = ESC//"[30m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       black               color.
        character(5 , SK)   :: fred             = ESC//"[31m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       red                 color.
        character(5 , SK)   :: fgreen           = ESC//"[32m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       green               color.
        character(5 , SK)   :: fyellow          = ESC//"[33m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       yellow              color.
        character(5 , SK)   :: fblue            = ESC//"[34m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       blue                color.
        character(5 , SK)   :: fmagenta         = ESC//"[35m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       magenta             color.
        character(5 , SK)   :: fcyan            = ESC//"[36m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       cyan                color.
        character(5 , SK)   :: fwhite           = ESC//"[37m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       white (light gray)  color.
        character(5 , SK)   :: fdefault         = ESC//"[39m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       default             color.
        character(6 , SK)   :: flblack          = ESC//"[90m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light black (dark gray)   color.
        character(6 , SK)   :: flred            = ESC//"[91m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light red                 color.
        character(6 , SK)   :: flgreen          = ESC//"[92m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light green               color.
        character(6 , SK)   :: flyellow         = ESC//"[93m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light yellow              color.
        character(6 , SK)   :: flblue           = ESC//"[94m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light blue                color.
        character(6 , SK)   :: flmagenta        = ESC//"[95m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light magenta             color.
        character(6 , SK)   :: flcyan           = ESC//"[96m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light cyan                color.
        character(6 , SK)   :: flwhite          = ESC//"[97m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light gray                color.

        character(5 , SK)   :: bblack           = ESC//"[40m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       black               color.
        character(5 , SK)   :: bred             = ESC//"[41m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       red                 color.
        character(5 , SK)   :: bgreen           = ESC//"[42m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       green               color.
        character(5 , SK)   :: byellow          = ESC//"[43m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       yellow              color.
        character(5 , SK)   :: bblue            = ESC//"[44m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       blue                color.
        character(5 , SK)   :: bmagenta         = ESC//"[45m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       magenta             color.
        character(5 , SK)   :: bcyan            = ESC//"[46m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       cyan                color.
        character(5 , SK)   :: bwhite           = ESC//"[47m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       white (light gray)  color.
        character(5 , SK)   :: bdefault         = ESC//"[49m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       default             color.
        character(6 , SK)   :: blblack          = ESC//"[100m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light black (dark gray)   color.
        character(6 , SK)   :: blred            = ESC//"[101m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light red                 color.
        character(6 , SK)   :: blgreen          = ESC//"[102m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light green               color.
        character(6 , SK)   :: blyellow         = ESC//"[103m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light yellow              color.
        character(6 , SK)   :: blblue           = ESC//"[104m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light blue                color.
        character(6 , SK)   :: blmagenta        = ESC//"[105m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light magenta             color.
        character(6 , SK)   :: blcyan           = ESC//"[106m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light cyan                color.
        character(6 , SK)   :: blwhite          = ESC//"[107m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light white               color.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [styleSeq_type](@ref styleSeq_type) Parameterized Derived Type (PDT) for styling Fortran scalar strings of arbitrary kinds
    !>  according to the conventions set by the ANSI. This derived type contains only a minimal subset of the standard recognized by
    !>  the terminals and terminal emulators, in particular, DEC VT100.
    !>
    !>  \details
    !>  This type is not meant to be used directly outside its host module.
    !>  Rather it is intended for the composition of the [styleSeq_type](@ref pm_strANSI::styleSeq_type) type.
    !>
    !>  \param[in]  kind    :   The kind type parameter of the PDT (**optionall**, default = \SK)
    !>
    !>  \warning
    !>  The ANSI **blink** escape code does not function as expected in most terminal emulators,
    !>  although it is known to work in the **tty** and **XTerm**.
    !>
    !>  \note
    !>  See [this page](https://misc.flogisoft.com/bash/tip_colors_and_formatting)
    !>  and [this page](https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797)
    !>  for more information on ANSI/VT100 standard subset.
    !>
    !>  \final{styleSeq_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                    :: styleSeq_type(kind)

        integer     , kind  :: kind             = SK

        !   \warning
        !   Any modifications to the following components must be also reflected in the list of constants above.

        character(4 , kind) :: bold             = ESC//"[1m"    !<  \public The ANSI escape code that makes subsequent texts appear **bold/bright**.
        character(4 , kind) :: bright           = ESC//"[1m"    !<  \public The ANSI escape code that makes subsequent texts appear **bold/bright**.
        character(4 , kind) :: dim              = ESC//"[2m"    !<  \public The ANSI escape code that makes subsequent texts appear **dim**.
        character(4 , kind) :: italic           = ESC//"[3m"    !<  \public The ANSI escape code that makes subsequent texts appear **italic**.
        character(4 , kind) :: underlined       = ESC//"[4m"    !<  \public The ANSI escape code that makes subsequent texts appear **underlined**.
        character(4 , kind) :: blinking         = ESC//"[5m"    !<  \public The ANSI escape code that makes subsequent texts appear **blinking** (not universally supported).
        character(4 , kind) :: reverse          = ESC//"[7m"    !<  \public The ANSI escape code that makes subsequent texts appear **in reverse**.
        character(4 , kind) :: hidden           = ESC//"[8m"    !<  \public The ANSI escape code that makes subsequent texts **disappear** (hidden but not removed!).
        character(4 , kind) :: strike           = ESC//"[9m"    !<  \public The ANSI escape code that makes subsequent texts to have **strike-through** (not universally supported).

        character(4 , kind) :: endbold          = ESC//"[1m"    !<  \public The ANSI escape code that **ends** the ANSI style **bold/bright**.
        character(4 , kind) :: endbright        = ESC//"[1m"    !<  \public The ANSI escape code that **ends** the ANSI style **bold/bright**.
        character(4 , kind) :: enddim           = ESC//"[2m"    !<  \public The ANSI escape code that **ends** the ANSI style **dim**.
        character(4 , kind) :: enditalic        = ESC//"[3m"    !<  \public The ANSI escape code that **ends** the ANSI style **italic**.
        character(4 , kind) :: endunderlined    = ESC//"[4m"    !<  \public The ANSI escape code that **ends** the ANSI style **underlined**.
        character(4 , kind) :: endblinking      = ESC//"[5m"    !<  \public The ANSI escape code that **ends** the ANSI style **blinking**.
        character(4 , kind) :: endreverse       = ESC//"[7m"    !<  \public The ANSI escape code that **ends** the ANSI style **reverse**.
        character(4 , kind) :: endhidden        = ESC//"[8m"    !<  \public The ANSI escape code that **ends** the ANSI style **hidden**.
        character(4 , kind) :: endstrike        = ESC//"[9m"    !<  \public The ANSI escape code that **ends** the ANSI style **strike-through**.
        character(4 , kind) :: endall           = ESC//"[0m"    !<  \public The ANSI escape code that **resets** all previously specified ANSI text styles.
        character(4 , SK)   :: reset            = ESC//"[0m"    !<  \public The ANSI escape code that **resets** all previously specified ANSI text styles (same as [](@ref pm_strANSI::endall)).

        character(5 , kind) :: fblack           = ESC//"[30m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       black               color.
        character(5 , kind) :: fred             = ESC//"[31m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       red                 color.
        character(5 , kind) :: fgreen           = ESC//"[32m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       green               color.
        character(5 , kind) :: fyellow          = ESC//"[33m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       yellow              color.
        character(5 , kind) :: fblue            = ESC//"[34m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       blue                color.
        character(5 , kind) :: fmagenta         = ESC//"[35m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       magenta             color.
        character(5 , kind) :: fcyan            = ESC//"[36m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       cyan                color.
        character(5 , kind) :: fwhite           = ESC//"[37m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       white (light gray)  color.
        character(5 , kind) :: fdefault         = ESC//"[39m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground       default             color.
        character(6 , kind) :: flblack          = ESC//"[90m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light black (dark gray)   color.
        character(6 , kind) :: flred            = ESC//"[91m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light red                 color.
        character(6 , kind) :: flgreen          = ESC//"[92m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light green               color.
        character(6 , kind) :: flyellow         = ESC//"[93m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light yellow              color.
        character(6 , kind) :: flblue           = ESC//"[94m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light blue                color.
        character(6 , kind) :: flmagenta        = ESC//"[95m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light magenta             color.
        character(6 , kind) :: flcyan           = ESC//"[96m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light cyan                color.
        character(6 , kind) :: flwhite          = ESC//"[97m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>f</b>oreground light gray                color.

        character(5 , kind) :: bblack           = ESC//"[40m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       black               color.
        character(5 , kind) :: bred             = ESC//"[41m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       red                 color.
        character(5 , kind) :: bgreen           = ESC//"[42m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       green               color.
        character(5 , kind) :: byellow          = ESC//"[43m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       yellow              color.
        character(5 , kind) :: bblue            = ESC//"[44m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       blue                color.
        character(5 , kind) :: bmagenta         = ESC//"[45m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       magenta             color.
        character(5 , kind) :: bcyan            = ESC//"[46m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       cyan                color.
        character(5 , kind) :: bwhite           = ESC//"[47m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       white (light gray)  color.
        character(5 , kind) :: bdefault         = ESC//"[49m"   !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground       default             color.
        character(6 , kind) :: blblack          = ESC//"[100m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light black (dark gray)   color.
        character(6 , kind) :: blred            = ESC//"[101m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light red                 color.
        character(6 , kind) :: blgreen          = ESC//"[102m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light green               color.
        character(6 , kind) :: blyellow         = ESC//"[103m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light yellow              color.
        character(6 , kind) :: blblue           = ESC//"[104m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light blue                color.
        character(6 , kind) :: blmagenta        = ESC//"[105m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light magenta             color.
        character(6 , kind) :: blcyan           = ESC//"[106m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light cyan                color.
        character(6 , kind) :: blwhite          = ESC//"[107m"  !<  \public The scalar `character` of kind \SKALL containing the ANSI code sequence for <b>b</b>ackground light white               color.

    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Generate and return a string containing the requested ANSI-standard text style and color or optionally,
!    !>  decorate the input string with the requested ANSI style.
!    !>
!    !>  \details
!    !>  The ANSI texts styles are useful for displaying colored/styled text in terminals that support the ANSI standard,
!    !>  most prominently VT100. All major terminals support ANSI escape sequences as of 2017.
!    !>
!    !>  \param[in]  attr    :   The input scalar, or array of the same rank and shape as other array-like arguments,of the same type
!    !>                          and kind as the output `strout`, containing the ANSI escape code for the desired text attribute.<br>
!    !>                          &nbsp; See the components of [styleSeq_type](@ref pm_strANSI::styleSeq_type) class for possible input values.<br>
!    !>                          &nbsp; (**optional**, default = [styleSeq_type%reset](@ref pm_strANSI::styleSeq_type))
!    !>  \param[in]  fg      :   The input scalar, or array of the same rank and shape as other array-like arguments, of the same type and
!    !>                          kind as the output `strout`, containing the ANSI escape code for the desired text **foreground** color.<br>
!    !>                          &nbsp; See the components of [styleSeq_type](@ref pm_strANSI::styleSeq_type) class for possible input values.<br>
!    !>                          &nbsp; (**optional**, default = [styleSeq_type%fg%default](@ref pm_strANSI::styleSeq_type))
!    !>  \param[in]  bg      :   The input scalar, or array of the same rank and shape as other array-like arguments, of the same type and
!    !>                          kind as the output `strout`, containing the ANSI escape code for the desired text **background** color.<br>
!    !>                          &nbsp; See the components of [styleSeq_type](@ref pm_strANSI::styleSeq_type) class for possible input values.<br>
!    !>                          &nbsp; (**optional**, default = [styleSeq_type%bg%default](@ref pm_strANSI::styleSeq_type))
!    !>  \param[in]  str     :   The input scalar, or array of the same rank and shape as other array-like arguments, of the same type and
!    !>                          kind as the output `strout`, containing the text to be styled.<br>
!    !>                          &nbsp; (**optional**, default = `""`)
!    !>
!    !>  \return
!    !>  `strout`            :   The output scalar, or array of the same rank and shape as the array-like input arguments
!    !>                          of type `character` of kind \SKALL containing the desired ANSI text style.<br>
!    !>                          If , the input argument `str` is present, then `strout` will be set as
!    !>                          `strout = style(attr, fg, bg) // str // endstyle` where the suffix
!    !>                          [endstyle](@ref pm_strANSI::endstyle) is a constant equivalent to
!    !>                          the `reset` component of [styleSeq_type](@ref pm_strANSI::styleSeq_type).
!    !>                          If all arguments are missing, the output is `
!    !>
!    !>  \interface
!    !>  \code{.F90}
!    !>
!    !>      use pm_strANSI, only: getStyle
!    !>
!    !>      style = getStyle(attr = attr, fg = fg, bg = bg, str = str) ! all arguments are optional.
!    !>      !
!    !>  \endcode
!    !>
!    !>  \warning
!    !>  The ANSI **blink** and **strike-through** escape codes do not function as expected in most terminal emulators,
!    !>  although it is known to work in the **tty** and **XTerm**.
!    !>
!    !>  \pure
!    !>
!    !>  \elemental
!    !>
!    !>  \note
!    !>  The sole purpose of this generic interface is to provide a **convenient** but **fast** method of generating a resized and centered
!    !>  copy of an array. See [pm_arrayResize](@ref pm_arrayResize) and [pm_strANSI](@ref pm_strANSI) for the relevant benchmarks.
!    !>
!    !>  \see
!    !>  [styleSeq_type](@ref pm_strANSI::styleSeq_type)<br>
!    !>
!    !>  \example
!    !>  \include{lineno} example/pm_strANSI/getStyle/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_strANSI/getStyle/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_strANSI](@ref test_pm_strANSI)
!    !>
!    !>  \todo
!    !>  \pmed Two new optional input scalar `lbcold` and `ubcold` arguments can be added to procedures to specify
!    !>  a subset of the contents of the original array that has to be kept in the newly allocated centered array.
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface getStyle
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    pure elemental module function style_D0_SK5(attr, fg, bg, str) result(strout)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: style_D0_SK5
!#endif
!        use pm_kind, only: IK, LK, SKG => SK5
!        character(1,SKG)            , intent(in)    , optional      :: attr
!        character(2,SKG)            , intent(in)    , optional      :: fg
!        character(*,SKG)            , intent(in)    , optional      :: bg
!        character(*,SKG)            , intent(in)    , optional      :: str
!    end function
!#endif
!
!#if SK4_ENABLED
!    pure elemental module function style_D0_SK4(attr, fg, bg, str) result(strout)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: style_D0_SK4
!#endif
!        use pm_kind, only: IK, LK, SKG => SK4
!        character(1,SKG)            , intent(in)    , optional      :: attr
!        character(2,SKG)            , intent(in)    , optional      :: fg
!        character(*,SKG)            , intent(in)    , optional      :: bg
!        character(*,SKG)            , intent(in)    , optional      :: str
!    end function
!#endif
!
!#if SK3_ENABLED
!    pure elemental module function style_D0_SK3(attr, fg, bg, str) result(strout)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: style_D0_SK3
!#endif
!        use pm_kind, only: IK, LK, SKG => SK3
!        character(1,SKG)            , intent(in)    , optional      :: attr
!        character(2,SKG)            , intent(in)    , optional      :: fg
!        character(*,SKG)            , intent(in)    , optional      :: bg
!        character(*,SKG)            , intent(in)    , optional      :: str
!    end function
!#endif
!
!#if SK2_ENABLED
!    pure elemental module function style_D0_SK2(attr, fg, bg, str) result(strout)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: style_D0_SK2
!#endif
!        use pm_kind, only: IK, LK, SKG => SK2
!        character(1,SKG)            , intent(in)    , optional      :: attr
!        character(2,SKG)            , intent(in)    , optional      :: fg
!        character(*,SKG)            , intent(in)    , optional      :: bg
!        character(*,SKG)            , intent(in)    , optional      :: str
!    end function
!#endif
!
!#if SK1_ENABLED
!    pure elemental module function style_D0_SK1(attr, fg, bg, str) result(strout)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: style_D0_SK1
!#endif
!        use pm_kind, only: IK, LK, SKG => SK1
!        character(1,SKG)            , intent(in)    , optional      :: attr
!        character(2,SKG)            , intent(in)    , optional      :: fg
!        character(*,SKG)            , intent(in)    , optional      :: bg
!        character(*,SKG)            , intent(in)    , optional      :: str
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_strANSI ! LCOV_EXCL_LINE