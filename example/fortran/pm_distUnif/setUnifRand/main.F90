program example

    use pm_kind, only: SK, LK
    use pm_kind, only: IK => IKS ! All other kinds are also supported.
    use pm_kind, only: CK => CKS ! All other kinds are also supported.
    use pm_kind, only: RK => RKS ! All other kinds are also supported.
    use pm_io, only: display_type
    use pm_distUnif, only: setUnifRand
    use pm_io, only: getErrTableWrite

    implicit none

    integer(IK)     , parameter :: NP = 5_IK

    character(5,SK) :: rand_D0_SK
    character(2,SK) ::             rand_D1_SK(NP), rand_D2_SK(NP,NP)
    logical(LK)     :: rand_D0_LK, rand_D1_LK(NP), rand_D2_LK(NP,NP)
    integer(IK)     :: rand_D0_IK, rand_D1_IK(NP), rand_D2_IK(NP,NP)
    complex(CK)     :: rand_D0_CK, rand_D1_CK(NP), rand_D2_CK(NP,NP)
    real(RK)        :: rand_D0_RK, rand_D1_RK(NP), rand_D2_RK(NP,NP)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `logical`-valued random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar logical random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_LK) ! random binary logical in the default range [.false._LK, .true._LK].")
                    call setUnifRand(rand_D0_LK)
    call disp%show("rand_D0_LK")
    call disp%show( rand_D0_LK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of logical random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_LK) ! random binary logical in the default range [.false._LK, .true._LK].")
                    call setUnifRand(rand_D1_LK)
    call disp%show("rand_D1_LK")
    call disp%show( rand_D1_LK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of logical random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_LK) ! random binary integer in the default range [.false._LK, .true._LK].")
                    call setUnifRand(rand_D2_LK)
    call disp%show("rand_D2_LK")
    call disp%show( rand_D2_LK )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random `character` values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar character random value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_SK) ! random ASCII character in the default range [char(1), char(127)].")
                    call setUnifRand(rand_D0_SK)
    call disp%show("rand_D0_SK")
    call disp%show( rand_D0_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_SK, lb = 'AA0aa', ub = 'ZZ9zz') ! random ASCII character in the range [""A"", ""z""].")
                    call setUnifRand(rand_D0_SK, lb = 'AA0aa', ub = 'ZZ9zz')
    call disp%show("rand_D0_SK")
    call disp%show( rand_D0_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of character random values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_SK) ! random ASCII character in the default range [char(1), char(127)].")
                    call setUnifRand(rand_D1_SK)
    call disp%show("rand_D1_SK")
    call disp%show( rand_D1_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_SK, lb = 'Aa', ub = 'Zz') ! random ASCII character in the range [""Aa"", ""Zz""] for each character in string separately.")
                    call setUnifRand(rand_D1_SK, lb = 'Aa', ub = 'Zz')
    call disp%show("rand_D1_SK")
    call disp%show( rand_D1_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of character random values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_SK) ! random ASCII character in the default range [char(1), char(127)].")
                    call setUnifRand(rand_D2_SK)
    call disp%show("rand_D2_SK")
    call disp%show( rand_D2_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_SK, lb = 'aA', ub = 'zA') ! random ASCII character in the range [""aA"", ""zZ""] for each character in string separately.")
                    call setUnifRand(rand_D2_SK, lb = 'aA', ub = 'zA')
    call disp%show("rand_D2_SK")
    call disp%show( rand_D2_SK , deliml = SK_"""" )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `integer`-valued random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar integer random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_IK) ! random binary integer in the default range [0, 1].")
                    call setUnifRand(rand_D0_IK)
    call disp%show("rand_D0_IK")
    call disp%show( rand_D0_IK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_IK, -3_IK, 2_IK) ! random integer in range [-3 2].")
                    call setUnifRand(rand_D0_IK, -3_IK, 2_IK)
    call disp%show("rand_D0_IK")
    call disp%show( rand_D0_IK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of integer random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_IK) ! random binary integer in the default range [0, 1].")
                    call setUnifRand(rand_D1_IK)
    call disp%show("rand_D1_IK")
    call disp%show( rand_D1_IK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_IK, -3_IK, 2_IK) ! random integer in range [-3 2].")
                    call setUnifRand(rand_D1_IK, -3_IK, 2_IK)
    call disp%show("rand_D1_IK")
    call disp%show( rand_D1_IK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of integer random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_IK) ! random binary integer in the default range [0, 1].")
                    call setUnifRand(rand_D2_IK)
    call disp%show("rand_D2_IK")
    call disp%show( rand_D2_IK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_IK, -3_IK, 2_IK) ! random integer in range [-3 2].")
                    call setUnifRand(rand_D2_IK, -3_IK, 2_IK)
    call disp%show("rand_D2_IK")
    call disp%show( rand_D2_IK )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `complex`-valued random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar complex random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_CK) ! random binary complex in the default range [(0,0), (1,1)].")
                    call setUnifRand(rand_D0_CK)
    call disp%show("rand_D0_CK")
    call disp%show( rand_D0_CK )

    call disp%show("call setUnifRand(rand_D0_CK, (-3._CK, 3._CK), (2._CK, 5._CK))")
                    call setUnifRand(rand_D0_CK, (-3._CK, 3._CK), (2._CK, 5._CK))
    call disp%show("rand_D0_CK")
    call disp%show( rand_D0_CK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of complex random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_CK) ! random binary complex in the default range [(0,0), (1,1)).")
                    call setUnifRand(rand_D1_CK)
    call disp%show("rand_D1_CK")
    call disp%show( rand_D1_CK )

    call disp%show("call setUnifRand(rand_D1_CK, (-3._CK, 3._CK), (2._CK, 5._CK)).")
                    call setUnifRand(rand_D1_CK, (-3._CK, 3._CK), (2._CK, 5._CK))
    call disp%show("rand_D1_CK")
    call disp%show( rand_D1_CK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of complex random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_CK) ! random binary complex in the default range [(0,0), (1,1)).")
                    call setUnifRand(rand_D2_CK)
    call disp%show("rand_D2_CK")
    call disp%show( rand_D2_CK )

    call disp%show("call setUnifRand(rand_D2_CK, (-3._CK, 3._CK), (2._CK, 5._CK)) ! random complex in range [(-3._CK, 3._CK), (2._CK, -2._CK)).")
                    call setUnifRand(rand_D2_CK, (-3._CK, 3._CK), (2._CK, 5._CK))
    call disp%show("rand_D2_CK")
    call disp%show( rand_D2_CK )
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `real`-valued random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar real random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_RK) ! random binary real in the default range [0, 1).")
                    call setUnifRand(rand_D0_RK)
    call disp%show("rand_D0_RK")
    call disp%show( rand_D0_RK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D0_RK, -3._RK, 2._RK) ! random real in range [-3 2).")
                    call setUnifRand(rand_D0_RK, -3._RK, 2._RK)
    call disp%show("rand_D0_RK")
    call disp%show( rand_D0_RK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_RK) ! random binary real in the default range [0, 1).")
                    call setUnifRand(rand_D1_RK)
    call disp%show("rand_D1_RK")
    call disp%show( rand_D1_RK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D1_RK, -3._RK, 2._RK) ! random real in range [-3 2).")
                    call setUnifRand(rand_D1_RK, -3._RK, 2._RK)
    call disp%show("rand_D1_RK")
    call disp%show( rand_D1_RK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_RK) ! random binary real in the default range [0, 1).")
                    call setUnifRand(rand_D2_RK)
    call disp%show("rand_D2_RK")
    call disp%show( rand_D2_RK )
    call disp%skip()

    call disp%show("call setUnifRand(rand_D2_RK, -3._RK, 2._RK) ! random real in range [-3 2).")
                    call setUnifRand(rand_D2_RK, -3._RK, 2._RK)
    call disp%show("rand_D2_RK")
    call disp%show( rand_D2_RK )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: rand(5000)
        call setUnifRand(rand, -2, 3)
        if (0 /= getErrTableWrite(SK_"setUnifRand.IK.txt", rand)) error stop "Table writing failed."
    end block

    block
        complex :: rand(5000)
        call setUnifRand(rand, (-2., +2.), (3., 5.))
        if (0 /= getErrTableWrite(SK_"setUnifRand.CK.txt", rand)) error stop "Table writing failed."
    end block

    block
        real :: rand(5000)
        call setUnifRand(rand)
        if (0 /= getErrTableWrite(SK_"setUnifRand.RK.txt", rand)) error stop "Table writing failed."
    end block

end program example