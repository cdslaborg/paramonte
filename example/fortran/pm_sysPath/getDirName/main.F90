program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_sysPath, only: verbatim
    use pm_sysPath, only: getDirName

    implicit none

    character(:, SK), allocatable :: path

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the dirname of a POSIX-compliant path.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_''")
                    path = SK_''
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.'")
                    path = SK_'.'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'a'")
                    path = SK_'a'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'..'")
                    path = SK_'..'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/.'")
                    path = SK_'/.'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/..'")
                    path = SK_'/..'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./'")
                    path = SK_'./'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.//'")
                    path = SK_'.//'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'../'")
                    path = SK_'../'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'..///'")
                    path = SK_'..///'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/'")
                    path = SK_'/'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//'")
                    path = SK_'//'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//./'")
                    path = SK_'//./'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//.//'")
                    path = SK_'//.//'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//.//a'")
                    path = SK_'//.//a'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'paramonte'")
                    path = SK_'paramonte'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/paramonte'")
                    path = SK_'/paramonte'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'paramonte/'")
                    path = SK_'paramonte/'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte'")
                    path = SK_'./paramonte'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte/'")
                    path = SK_'./paramonte/'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte/parallel/library.txt'")
                    path = SK_'./paramonte/parallel/library.txt'
    call disp%show("getDirName(path, SK_'/')")
    call disp%show( getDirName(path, SK_'/') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/', verbatim)")
    call disp%show( getDirName(path, SK_'/', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the dirname of a Windows-compliant path.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_''")
                    path = SK_''
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.'")
                    path = SK_'.'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'..'")
                    path = SK_'..'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/.'")
                    path = SK_'/.'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'\.'")
                    path = SK_'\.'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/\..'")
                    path = SK_'/\..'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./'")
                    path = SK_'./'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'../'")
                    path = SK_'../'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/'")
                    path = SK_'/'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/\'")
                    path = SK_'/\'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//'")
                    path = SK_'//'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//.\'")
                    path = SK_'//.\'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'//.\/'")
                    path = SK_'//.\/'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'\\.//a'")
                    path = SK_'\\.//a'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'paramonte'")
                    path = SK_'paramonte'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'/paramonte\'")
                    path = SK_'/paramonte\'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'paramonte/\'")
                    path = SK_'paramonte/\'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte'")
                    path = SK_'./paramonte'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'./paramonte/'")
                    path = SK_'./paramonte/'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("path = SK_'.\paramonte/parallel\library.txt'")
                    path = SK_'.\paramonte/parallel\library.txt'
    call disp%show("getDirName(path, SK_'/\')")
    call disp%show( getDirName(path, SK_'/\') , deliml = SK_"""" )
    call disp%show("getDirName(path, SK_'/\', verbatim)")
    call disp%show( getDirName(path, SK_'/\', verbatim) , deliml = SK_"""" )
    call disp%skip()

end program example