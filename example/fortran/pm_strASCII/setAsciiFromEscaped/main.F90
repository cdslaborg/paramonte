program example

    use pm_kind, only: IK
    use pm_kind, only: SK ! all other processor kinds are also supported.
    use pm_io, only: display_type
    use pm_strASCII, only: setAsciiFromEscaped

    implicit none

    character(:), allocatable   :: strascii
    character(255)              :: str, ascii
    integer(IK)                 :: endloc

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("str = ''")
                    str = ''
    call disp%show("call setAsciiFromEscaped(trim(str), ascii, endloc)")
                    call setAsciiFromEscaped(trim(str), ascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("ascii(1:endloc)")
    call disp%show( ascii(1:endloc) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("strascii = trim(str)")
                    strascii = trim(str)
    call disp%show("call setAsciiFromEscaped(strascii, endloc)")
                    call setAsciiFromEscaped(strascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("strascii(1:endloc)")
    call disp%show( strascii(1:endloc) , deliml = """" )
    call disp%show("strascii(1:endloc) == ascii(1:endloc)")
    call disp%show( strascii(1:endloc) == ascii(1:endloc) )
    call disp%skip()

    call disp%skip()
    call disp%show("str = 'A'")
                    str = 'A'
    call disp%show("call setAsciiFromEscaped(trim(str), ascii, endloc)")
                    call setAsciiFromEscaped(trim(str), ascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("ascii(1:endloc)")
    call disp%show( ascii(1:endloc) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("strascii = trim(str)")
                    strascii = trim(str)
    call disp%show("call setAsciiFromEscaped(strascii, endloc)")
                    call setAsciiFromEscaped(strascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("strascii(1:endloc)")
    call disp%show( strascii(1:endloc) , deliml = """" )
    call disp%show("strascii(1:endloc) == ascii(1:endloc)")
    call disp%show( strascii(1:endloc) == ascii(1:endloc) )
    call disp%skip()

    call disp%skip()
    call disp%show("str = '\paramonte'")
                    str = '\paramonte'
    call disp%show("call setAsciiFromEscaped(trim(str), ascii, endloc)")
                    call setAsciiFromEscaped(trim(str), ascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("ascii(1:endloc)")
    call disp%show( ascii(1:endloc) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("strascii = trim(str)")
                    strascii = trim(str)
    call disp%show("call setAsciiFromEscaped(strascii, endloc)")
                    call setAsciiFromEscaped(strascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("strascii(1:endloc)")
    call disp%show( strascii(1:endloc) , deliml = """" )
    call disp%show("strascii(1:endloc) == ascii(1:endloc)")
    call disp%show( strascii(1:endloc) == ascii(1:endloc) )
    call disp%skip()

    call disp%skip()
    call disp%show("str = '\\nn\004\t'")
                    str = '\\nn\004\t'
    call disp%show("call setAsciiFromEscaped(trim(str), ascii, endloc)")
                    call setAsciiFromEscaped(trim(str), ascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("ascii(1:endloc)")
    call disp%show( ascii(1:endloc) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("strascii = trim(str)")
                    strascii = trim(str)
    call disp%show("call setAsciiFromEscaped(strascii, endloc)")
                    call setAsciiFromEscaped(strascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("strascii(1:endloc)")
    call disp%show( strascii(1:endloc) , deliml = """" )
    call disp%show("strascii(1:endloc) == ascii(1:endloc)")
    call disp%show( strascii(1:endloc) == ascii(1:endloc) )
    call disp%skip()

    call disp%skip()
    call disp%show("str = '\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007E\U+2661'")
                    str = '\nn\fn\04\? \r\b\\a\47\41\x009f\xZFa\u002C\U0000007E\U+2661'
    call disp%show("call setAsciiFromEscaped(trim(str), ascii, endloc)")
                    call setAsciiFromEscaped(trim(str), ascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("ascii(1:endloc)")
    call disp%show( ascii(1:endloc) , deliml = """" )
    call disp%skip()

    call disp%skip()
    call disp%show("strascii = trim(str)")
                    strascii = trim(str)
    call disp%show("call setAsciiFromEscaped(strascii, endloc)")
                    call setAsciiFromEscaped(strascii, endloc)
    call disp%show("endloc")
    call disp%show( endloc )
    call disp%show("strascii(1:endloc)")
    call disp%show( strascii(1:endloc) , deliml = """" )
    call disp%show("strascii(1:endloc) == ascii(1:endloc)")
    call disp%show( strascii(1:endloc) == ascii(1:endloc) )
    call disp%skip()

end program example