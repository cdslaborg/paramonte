program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_matrixInit, only: getMatInit, dia, uppDia, lowDia, uppLow, uppLowDia
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange

    implicit none

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize diagonal matrices.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct       diagonal string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = 'OO')")
                        mat = getMatInit(shape(mat, IK), dia, vdia = 'OO')
        call disp%show("where(mat /= 'OO'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'OO'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = 'AA', ndia = 5_IK, roff = 3_IK)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = 'AA', ndia = 5_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = ['QQ', 'YY', 'WW', 'XX'], roff = 3_IK, coff = 2_IK)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = ['QQ', 'YY', 'WW', 'XX'], roff = 3_IK, coff = 2_IK)
        call disp%show("where(mat /= 'QQ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'QQ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct       diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = -2)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = -2)
        call disp%show("where(mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = getRange(1, 7), roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = getRange(1, 7), roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct       diagonal logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = .true.)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = .true.)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = [.true., .false., .true., .false.], roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = [.true., .false., .true., .false.], roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct       diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = -2., ndia = size(mat,1,IK))")
                        mat = getMatInit(shape(mat, IK), dia, vdia = -2., ndia = size(mat,1,IK))
        call disp%show("where(mat /= -2.); mat = 0.; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -2.); mat = 0.; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(shape(mat, IK), dia, vdia = -4., ndia = 7_IK, coff = 3_IK)")
                        mat = getMatInit(shape(mat, IK), dia, vdia = -4., ndia = 7_IK, coff = 3_IK)
        call disp%show("where(mat /= -4.); mat = 0.; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -4.); mat = 0.; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize upper-diagonal matrices.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'AA', vdia = 'ZZ', nrow = 4_IK, ncol = 7_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'AA', vdia = 'ZZ', nrow = 4_IK, ncol = 7_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'ZZ'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'ZZ'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'ZZ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'ZZ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -1_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -1_IK)
        call disp%show("where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -4_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -4_IK)
        call disp%show("where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -1, vdia = -2)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -1, vdia = -2)
        call disp%show("where(mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -3, vdia = -4, nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -3, vdia = -4, nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= -3 .and. mat /= -4); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -3 .and. mat /= -4); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        mat = .false.
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = .true., vdia = .true.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = .true., vdia = .true.)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

        mat = .false.
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = .true., vdia = .true., nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = .true., vdia = .true., nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -1., vdia = -2.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -1., vdia = -2.)
        call disp%show("where(mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -3., vdia = -4., nrow = 3_IK, ncol = 7_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppDia, vupp = -3., vdia = -4., nrow = 3_IK, ncol = 7_IK, coff = 3_IK)
        call disp%show("where(mat /= -3. .and. mat /= -4.); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -3. .and. mat /= -4.); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize lower-diagonal matrices.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'AA', vdia = 'ZZ', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'AA', vdia = 'ZZ', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'ZZ'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'ZZ'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'ZZ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'ZZ' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)
        call disp%show("where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
        call disp%show("where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()
        mat = ""

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -1, vdia = -2)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -1, vdia = -2)
        call disp%show("where(mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= -3 .and. mat /= -4); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -3 .and. mat /= -4); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        mat = .false.
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = .true., vdia = .true.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = .true., vdia = .true.)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

        mat = .false.
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = .true., vdia = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = .true., vdia = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -1., vdia = -2.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -1., vdia = -2.)
        call disp%show("where(mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = lowDia, vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= -3. .and. mat /= -4.); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= -3. .and. mat /= -4.); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize upper-lower matrices.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 4_IK, ncol = 7_IK, roff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'HH', nrow = 4_IK, ncol = 7_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'BB', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = 'AA', vlow = 'BB', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit(vupp = 'AA', vlow = 'BB', Shape = [10_IK, 10_IK], subset = uppLow, nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        mat = getMatInit(vupp = 'AA', vlow = 'BB', Shape = [10_IK, 10_IK], subset = uppLow, nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1, vlow = -1)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1, vlow = -1)
        call disp%show("where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1, vlow = -3, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1, vlow = -3, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= +1 .and. mat /= -3 .and. mat /= -4); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -3 .and. mat /= -4); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = .true., vlow = .true.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = .true., vlow = .true.)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = .true., vlow = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK, doff = +1_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = .true., vlow = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK, doff = +1_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1., vlow = -1.)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +1., vlow = -1.)
        call disp%show("where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +3., vlow = -3., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)")
                        mat = getMatInit(Shape = [10_IK, 10_IK], subset = uppLow, vupp = +3., vlow = -3., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= +3. .and. mat /= -3. .and. mat /= -4.); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +3. .and. mat /= -3. .and. mat /= -4.); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize upper-lower-diagonal matrices.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = 'OO', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = 'OO', nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'HH' .and. mat /= 'OO' .and. mat /= 'YY' .and. mat /= 'WW' .and. mat /= 'XX'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 6_IK, roff = 2_IK, coff = 3_IK, doff = +2_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
        call disp%show("where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= 'AA' .and. mat /= 'BB' .and. mat /= 'YY'); mat = '  '; endwhere
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1, vlow = -1, vdia = -2)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1, vlow = -1, vdia = -2)
        call disp%show("where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= +1 .and. mat /= -3 .and. mat /= -4); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -3 .and. mat /= -4); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = .true., vlow = .true., vdia = .false.)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = .true., vlow = .true., vdia = .false.)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = .true., vlow = .true., vdia = .false., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK, doff = +1_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = .true., vlow = .true., vdia = .false., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK, doff = +1_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = .false.

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1., vlow = -1., vdia = -2.)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +1., vlow = -1., vdia = -2.)
        call disp%show("where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +1 .and. mat /= -1 .and. mat /= -2); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

        call disp%skip()
        call disp%show("mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +3., vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)")
                        mat = getMatInit([10_IK, 10_IK], uppLowDia, vupp = +3., vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, coff = 3_IK)
        call disp%show("where(mat /= +3. .and. mat /= -3. .and. mat /= -4.); mat = 0; endwhere ! get rid of uninitialized values for better illustration.")
                        where(mat /= +3. .and. mat /= -3. .and. mat /= -4.); mat = 0; endwhere
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()
        mat = 0

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end program example