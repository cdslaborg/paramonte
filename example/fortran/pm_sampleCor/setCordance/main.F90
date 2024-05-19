program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_sampleCor, only: setCordance
    use pm_arrayShuffle, only: getShuffled
    use pm_arrayRange, only: getRange
    use pm_sampleCor, only: css_type
    use pm_sampleCor, only: uppDia
    use pm_sampleCor, only: lowDia

    implicit none

    integer(IK) :: itry, ntry = 15
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cordance of two random string containers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) :: nsam, cordance, concordance, discordance, tiedx, tiedy
        type(css_type), allocatable :: x(:), y(:)
        do itry = 1, 2
            call disp%skip()
            call disp%show("nsam = getUnifRand(0, 20)")
                            nsam = getUnifRand(0, 20)
            call disp%show("nsam")
            call disp%show( nsam )
            call disp%show("x = css_type([character(10) :: 'ParaMonte', 'is', 'a', 'Monte', 'Carlo', 'and', 'Machine', 'Learning', 'library', '.'])")
                            x = css_type([character(10) :: 'ParaMonte', 'is', 'a', 'Monte', 'Carlo', 'and', 'Machine', 'Learning', 'library', '.'])
            call disp%show("x")
            call disp%show( x , deliml = SK_'''' )
            call disp%show("y = getShuffled(x)")
                            y = getShuffled(x)
            call disp%show("y")
            call disp%show( y , deliml = SK_'''' )
            call disp%show("call setCordance(cordance, tiedx, tiedy, x, y)")
                            call setCordance(cordance, tiedx, tiedy, x, y)
            call disp%show("[cordance, tiedx, tiedy]")
            call disp%show( [cordance, tiedx, tiedy] )
            call disp%show("call setCordance(concordance, discordance, tiedx, tiedy, x, y)")
                            call setCordance(concordance, discordance, tiedx, tiedy, x, y)
            call disp%show("[concordance, discordance, tiedx, tiedy]")
            call disp%show( [concordance, discordance, tiedx, tiedy] )
            call disp%show("if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'")
                            if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cordance of two random character strings.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) :: nsam, cordance, concordance, discordance, tiedx, tiedy
        character(:), allocatable :: x, y
        do itry = 1, ntry
            call disp%skip()
            call disp%show("nsam = getUnifRand(0, 20)")
                            nsam = getUnifRand(0, 20)
            call disp%show("nsam")
            call disp%show( nsam )
            call disp%show("x = getUnifRand(repeat('A', nsam), repeat('Z', nsam))")
                            x = getUnifRand(repeat('A', nsam), repeat('Z', nsam))
            call disp%show("x")
            call disp%show( x , deliml = SK_'''' )
            call disp%show("y = getUnifRand(repeat('A', nsam), repeat('Z', nsam))")
                            y = getUnifRand(repeat('A', nsam), repeat('Z', nsam))
            call disp%show("y")
            call disp%show( y , deliml = SK_'''' )
            call disp%show("call setCordance(cordance, tiedx, tiedy, x, y)")
                            call setCordance(cordance, tiedx, tiedy, x, y)
            call disp%show("[cordance, tiedx, tiedy]")
            call disp%show( [cordance, tiedx, tiedy] )
            call disp%show("call setCordance(concordance, discordance, tiedx, tiedy, x, y)")
                            call setCordance(concordance, discordance, tiedx, tiedy, x, y)
            call disp%show("[concordance, discordance, tiedx, tiedy]")
            call disp%show( [concordance, discordance, tiedx, tiedy] )
            call disp%show("if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'")
                            if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cordance of two random character vectors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) :: nsam, cordance, concordance, discordance, tiedx, tiedy
        character(2), allocatable :: x(:), y(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("nsam = getUnifRand(0, 20)")
                            nsam = getUnifRand(0, 20)
            call disp%show("nsam")
            call disp%show( nsam )
            call disp%show("x = getUnifRand('AA', 'ZZ', nsam)")
                            x = getUnifRand('AA', 'ZZ', nsam)
            call disp%show("x")
            call disp%show( x , deliml = SK_'''' )
            call disp%show("y = getUnifRand('AA', 'ZZ', nsam)")
                            y = getUnifRand('AA', 'ZZ', nsam)
            call disp%show("y")
            call disp%show( y , deliml = SK_'''' )
            call disp%show("call setCordance(cordance, tiedx, tiedy, x, y)")
                            call setCordance(cordance, tiedx, tiedy, x, y)
            call disp%show("[cordance, tiedx, tiedy]")
            call disp%show( [cordance, tiedx, tiedy] )
            call disp%show("call setCordance(concordance, discordance, tiedx, tiedy, x, y)")
                            call setCordance(concordance, discordance, tiedx, tiedy, x, y)
            call disp%show("[concordance, discordance, tiedx, tiedy]")
            call disp%show( [concordance, discordance, tiedx, tiedy] )
            call disp%show("if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'")
                            if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cordance of two random integer vectors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => IK ! All processor kinds are supported.
        integer(IK) :: nsam, cordance, concordance, discordance, tiedx, tiedy
        integer(TKG), allocatable :: x(:), y(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("nsam = getUnifRand(0, 20)")
                            nsam = getUnifRand(0, 20)
            call disp%show("nsam")
            call disp%show( nsam )
            call disp%show("x = getUnifRand(1, nsam, nsam)")
                            x = getUnifRand(1, nsam, nsam)
            call disp%show("x")
            call disp%show( x )
            call disp%show("y = getUnifRand(1, nsam, nsam)")
                            y = getUnifRand(1, nsam, nsam)
            call disp%show("y")
            call disp%show( y )
            call disp%show("call setCordance(cordance, tiedx, tiedy, x, y)")
                            call setCordance(cordance, tiedx, tiedy, x, y)
            call disp%show("[cordance, tiedx, tiedy]")
            call disp%show( [cordance, tiedx, tiedy] )
            call disp%show("call setCordance(concordance, discordance, tiedx, tiedy, x, y)")
                            call setCordance(concordance, discordance, tiedx, tiedy, x, y)
            call disp%show("[concordance, discordance, tiedx, tiedy]")
            call disp%show( [concordance, discordance, tiedx, tiedy] )
            call disp%show("if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'")
                            if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cordance of two random    real vectors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => IK ! All processor kinds are supported.
        integer(IK) :: nsam, cordance, concordance, discordance, tiedx, tiedy
        real(TKG), allocatable :: x(:), y(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("nsam = getUnifRand(0, 20)")
                            nsam = getUnifRand(0, 20)
            call disp%show("nsam")
            call disp%show( nsam )
            call disp%show("x = getUnifRand(1, nsam, nsam)")
                            x = getUnifRand(1, nsam, nsam)
            call disp%show("x")
            call disp%show( x )
            call disp%show("y = getUnifRand(1, nsam, nsam)")
                            y = getUnifRand(1, nsam, nsam)
            call disp%show("y")
            call disp%show( y )
            call disp%show("call setCordance(cordance, tiedx, tiedy, x, y)")
                            call setCordance(cordance, tiedx, tiedy, x, y)
            call disp%show("[cordance, tiedx, tiedy]")
            call disp%show( [cordance, tiedx, tiedy] )
            call disp%show("call setCordance(concordance, discordance, tiedx, tiedy, x, y)")
                            call setCordance(concordance, discordance, tiedx, tiedy, x, y)
            call disp%show("[concordance, discordance, tiedx, tiedy]")
            call disp%show( [concordance, discordance, tiedx, tiedy] )
            call disp%show("if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'")
                            if (concordance - discordance /= cordance) error stop 'Internal library error detected. This never happens.'
            call disp%skip()
        end do
    end block

end program example