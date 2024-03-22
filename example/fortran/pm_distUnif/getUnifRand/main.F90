program example

    use pm_kind, only: SK, LK
    use pm_kind, only: IK => IKS ! All other kinds are also supported.
    use pm_kind, only: CK => CKS ! All other kinds are also supported.
    use pm_kind, only: RK => RKS ! All other kinds are also supported.
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand

    implicit none

    integer(IK)     , parameter :: NP = 5_IK

    character(5,SK) :: rand_D0_SK
    character(2,SK) ::             rand_D1_SK(NP)
    integer(IK)     :: rand_D0_IK, rand_D1_IK(NP)
    complex(CK)     :: rand_D0_CK, rand_D1_CK(NP)
    real(RK)        :: rand_D0_RK, rand_D1_RK(NP)

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

    call disp%show("getUnifRand() ! single random coin flipping `logical` result in the default range [.false._LK, .true._LK].")
    call disp%show( getUnifRand() )
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

    call disp%show("rand_D0_SK = getUnifRand(lb = SK_'AA0aa', ub = SK_'ZZ9zz') ! random ASCII character in the specified range per single character.")
                    rand_D0_SK = getUnifRand(lb = SK_'AA0aa', ub = SK_'ZZ9zz')
    call disp%show("rand_D0_SK")
    call disp%show( rand_D0_SK , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of character random values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_D1_SK = getUnifRand(lb = SK_'Aa', ub = SK_'Zz', s1 = size(rand_D1_SK, kind = IK)) ! random ASCII character in the range [""Aa"", ""Zz""] per character in string separately.")
                    rand_D1_SK = getUnifRand(lb = SK_'Aa', ub = SK_'Zz', s1 = size(rand_D1_SK, kind = IK))
    call disp%show("rand_D1_SK")
    call disp%show( rand_D1_SK , deliml = SK_"""" )
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

    call disp%show("rand_D0_IK = getUnifRand(-3_IK, 2_IK) ! random integer in range [-3 2].")
                    rand_D0_IK = getUnifRand(-3_IK, 2_IK)
    call disp%show("rand_D0_IK")
    call disp%show( rand_D0_IK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of integer random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_D1_IK = getUnifRand(-3_IK, 2_IK, s1 = size(rand_D1_IK, kind = IK)) ! random integer in range [-3 2].")
                    rand_D1_IK = getUnifRand(-3_IK, 2_IK, s1 = size(rand_D1_IK, kind = IK))
    call disp%show("rand_D1_IK")
    call disp%show( rand_D1_IK )
    call disp%skip()

    call disp%show("getUnifRand([0_IK, 5_IK, 10_IK], 15_IK)")
    call disp%show( getUnifRand([0_IK, 5_IK, 10_IK], 15_IK) )
    call disp%skip()

    call disp%show("getUnifRand([0_IK, 6_IK, 11_IK], [5_IK, 10_IK, 15_IK])")
    call disp%show( getUnifRand([0_IK, 6_IK, 11_IK], [5_IK, 10_IK, 15_IK]) )
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

    call disp%show("rand_D0_CK = getUnifRand(lb = (-3._CK, -2._CK), ub = (2._CK, 3._CK)).")
                    rand_D0_CK = getUnifRand(lb = (-3._CK, -2._CK), ub = (2._CK, 3._CK))
    call disp%show("rand_D0_CK")
    call disp%show( rand_D0_CK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of complex random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()
    call disp%show("rand_D1_CK = getUnifRand((-3._CK, 2._CK), (2._CK, 4._CK), s1 = size(rand_D1_CK, kind = IK)).")
                    rand_D1_CK = getUnifRand((-3._CK, 2._CK), (2._CK, 4._CK), s1 = size(rand_D1_CK, kind = IK))
    call disp%show("rand_D1_CK")
    call disp%show( rand_D1_CK )
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

    call disp%show("rand_D0_RK = getUnifRand(-3._RK, 2._RK) ! random real in range [-3 2].")
                    rand_D0_RK = getUnifRand(-3._RK, 2._RK)
    call disp%show("rand_D0_RK")
    call disp%show( rand_D0_RK )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_D1_RK = getUnifRand(-3._RK, 2._RK, s1 = size(rand_D1_RK, kind = IK)) ! random real in range [-3 2].")
                    rand_D1_RK = getUnifRand(-3._RK, 2._RK, s1 = size(rand_D1_RK, kind = IK))
    call disp%show("rand_D1_RK")
    call disp%show( rand_D1_RK )
    call disp%skip()

    ! Write an example complex random array to the output file for plotting.

    block
        integer :: fileUnit,i
        open(newunit = fileUnit, file = "main.unif.rand.txt")
        write(fileUnit,"(2(g0,:,','))") getUnifRand((0._CK,100._CK), (100._CK,150._CK), s1 = 500_IK)
    end block

end program example