program example

    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: setUnifCDF
    use pm_arrayRange, only: getRange

    implicit none

    integer(IK) , parameter :: NP = 1000_IK

    ! 1-dimensional array of x values.

    integer(IK) , allocatable :: Point_IK(:)
    real(RK)    , allocatable :: Point_RK(:), CDF_RK(:), CDF_IK(:)
    complex(CK) , allocatable :: Point_CK(:), CDF_CK(:)

    integer(IK) :: range_IK, lower_IK = -3_IK           , upper_IK = +4_IK
    real(RK)    :: range_RK, lower_RK = -4._RK          , upper_RK = +4._RK
    complex(CK) :: range_CK, lower_CK = (-4._CK,-1._CK) , upper_CK = (+4._CK,+3._CK)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    range_IK = upper_IK - lower_IK
    range_RK = upper_RK - lower_RK
    range_CK = upper_CK - lower_CK

    Point_IK =    getRange(lower_IK - range_IK/2, upper_IK + range_IK/2)
    Point_RK = getLinSpace(lower_RK - range_RK/2, upper_RK + range_RK/2, count = NP)
    Point_CK = getLinSpace(lower_CK - range_CK/2, upper_CK + range_CK/2, count = NP)
    allocate(CDF_IK(size(Point_IK)))
    allocate(CDF_RK(size(Point_RK)))
    allocate(CDF_CK(size(Point_CK)))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the Uniform distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the discrete Uniform CDF at an input scalar integer value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setUnifCDF(CDF_IK(1), x = 0_IK)")
                    call setUnifCDF(CDF_IK(1), x = 0_IK)
    call disp%show("CDF_IK(1)")
    call disp%show( CDF_IK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("lower_IK")
    call disp%show( lower_IK )
    call disp%show("upper_IK")
    call disp%show( upper_IK )
    call disp%show("call setUnifCDF(CDF_IK(1), x = 0_IK, lower = lower_IK, upper = upper_IK)")
                    call setUnifCDF(CDF_IK(1), x = 0_IK, lower = lower_IK, upper = upper_IK)
    call disp%show("CDF_IK(1)")
    call disp%show( CDF_IK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the continuous Uniform CDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setUnifCDF(CDF_RK(1), x = 0._RK)")
                    call setUnifCDF(CDF_RK(1), x = 0._RK)
    call disp%show("CDF_RK(1)")
    call disp%show( CDF_RK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("lower_RK")
    call disp%show( lower_RK )
    call disp%show("upper_RK")
    call disp%show( upper_RK )
    call disp%show("call setUnifCDF(CDF_RK(1), x = 0._RK, lower = lower_RK, upper = upper_RK)")
                    call setUnifCDF(CDF_RK(1), x = 0._RK, lower = lower_RK, upper = upper_RK)
    call disp%show("CDF_RK(1)")
    call disp%show( CDF_RK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the continuous Uniform CDF at an input scalar complex value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setUnifCDF(CDF_CK(1), x = (0._CK,0._CK))")
                    call setUnifCDF(CDF_CK(1), x = (0._CK,0._CK))
    call disp%show("CDF_CK(1)")
    call disp%show( CDF_CK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("lower_CK")
    call disp%show( lower_CK )
    call disp%show("upper_CK")
    call disp%show( upper_CK )
    call disp%show("call setUnifCDF(CDF_CK(1), x = (0._CK,0._CK), lower = lower_CK, upper = upper_CK)")
                    call setUnifCDF(CDF_CK(1), x = (0._CK,0._CK), lower = lower_CK, upper = upper_CK)
    call disp%show("CDF_CK(1)")
    call disp%show( CDF_CK(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call setUnifCDF(CDF_IK, Point_IK, lower_IK, upper_IK)
    call setUnifCDF(CDF_RK, Point_RK, lower_RK, upper_RK)
    call setUnifCDF(CDF_CK, Point_CK, lower_CK, upper_CK)

    block

        integer :: fileUnit, i

        open(newunit = fileUnit, file = "main.unif.cdf.IK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Point_IK(i), CDF_IK(i), i = 1, size(Point_IK))
        close(fileUnit)

        open(newunit = fileUnit, file = "main.unif.cdf.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Point_RK(i), CDF_RK(i), i = 1, size(Point_RK))
        close(fileUnit)

        open(newunit = fileUnit, file = "main.unif.cdf.CK.txt")
        write(fileUnit,"(4(g0,:,' '))") (Point_CK(i), CDF_CK(i), i = 1, size(Point_CK))
        close(fileUnit)

    end block

end program example