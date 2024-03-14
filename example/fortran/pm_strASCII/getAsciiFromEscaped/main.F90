program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_strASCII, only: getAsciiFromEscaped

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getAsciiFromEscaped('')")
    call disp%show( getAsciiFromEscaped('') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getAsciiFromEscaped('A')")
    call disp%show( getAsciiFromEscaped('A') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getAsciiFromEscaped('\paramonte')")
    call disp%show( getAsciiFromEscaped('\paramonte') , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("getAsciiFromEscaped('\\nn\004\t')")
    call disp%show( getAsciiFromEscaped('\\nn\004\t') , deliml = """" )
    call disp%show("'\nn'//achar(4)//achar(9) == getAsciiFromEscaped('\\nn\004\t')")
    call disp%show( '\nn'//achar(4)//achar(9) == getAsciiFromEscaped('\\nn\004\t') )
    call disp%skip()

    call disp%skip()
    call disp%show("getAsciiFromEscaped('\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007E\U+2661')")
    call disp%show( getAsciiFromEscaped('\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007E\U+2661') , deliml = """" )
    call disp%show("achar(10)//'n'//achar(12)//'n'//achar(4)//'? '//achar(13)//achar(8)//'\a'//achar(39)//achar(33)//achar(9)//'f\xZFa'//char(44)//char(126)//'A\UFBDA' == getAsciiFromEscaped('\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007EA\UFBDA')")
    call disp%show( achar(10)//'n'//achar(12)//'n'//achar(4)//'? '//achar(13)//achar(8)//'\a'//achar(39)//achar(33)//achar(9)//'f\xZFa'//char(44)//char(126)//'A\UFBDA' == getAsciiFromEscaped('\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007EA\UFBDA') )
    call disp%skip()

end program example