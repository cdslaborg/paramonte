program example

    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_matrixSubset, only: dia, uppDia, lowDia, uppLow, uppLowDia
    use pm_matrixInit, only: setMatInit
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
        call disp%show("mat = '' ! preset the matrix for better illustration.")
                        mat = ''
        call disp%show("call setMatInit(mat, dia, vdia = 'OO', ndia = minval(shape(mat, IK), 1))")
                        call setMatInit(mat, dia, vdia = 'OO', ndia = minval(shape(mat, IK), 1))
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = '' ! preset the matrix for better illustration.")
                        mat = ''
        call disp%show("call setMatInit(mat, dia, vdia = 'OO', ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)")
                        call setMatInit(mat, dia, vdia = 'OO', ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = '' ! preset the matrix for better illustration.")
                        mat = ''
        call disp%show("call setMatInit(mat, dia, vdia = 'OO', ndia = 5_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, dia, vdia = 'OO', ndia = 5_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = '' ! preset the matrix for better illustration.")
                        mat = ''
        call disp%show("call setMatInit(mat, dia, vdia = ['OO', 'YY', 'WW', 'XX'], roff = 3_IK, coff = 2_IK)")
                        call setMatInit(mat, dia, vdia = ['OO', 'YY', 'WW', 'XX'], roff = 3_IK, coff = 2_IK)
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
        call disp%show("mat = 0 ! preset the matrix for better illustration.")
                        mat = 0
        call disp%show("call setMatInit(mat, dia, vdia = -2, ndia = minval(shape(mat, IK), 1))")
                        call setMatInit(mat, dia, vdia = -2, ndia = minval(shape(mat, IK), 1))
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = 0 ! preset the matrix for better illustration.")
                        mat = 0
        call disp%show("call setMatInit(mat, dia, vdia = -2, ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)")
                        call setMatInit(mat, dia, vdia = -2, ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = 0 ! preset the matrix for better illustration.")
                        mat = 0
        call disp%show("call setMatInit(mat, dia, vdia = getRange(1, 7), roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, dia, vdia = getRange(1, 7), roff = 3_IK, coff = 3_IK)
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
        call disp%show("mat = .false. ! preset the matrix for better illustration.")
                        mat = .false.
        call disp%show("call setMatInit(mat, dia, vdia = .true., ndia = minval(shape(mat, IK), 1))")
                        call setMatInit(mat, dia, vdia = .true., ndia = minval(shape(mat, IK), 1))
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = .false. ! preset the matrix for better illustration.")
                        mat = .false.
        call disp%show("call setMatInit(mat, dia, vdia = .true., ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)")
                        call setMatInit(mat, dia, vdia = .true., ndia = minval(shape(mat, IK), 1), roff = 0_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = .false. ! preset the matrix for better illustration.")
                        mat = .false.
        call disp%show("call setMatInit(mat, dia, vdia = [.true., .false., .true., .false.], roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, dia, vdia = [.true., .false., .true., .false.], roff = 3_IK, coff = 3_IK)
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
        call disp%show("mat = 0. ! preset the matrix for better illustration.")
                        mat = 0.
        call disp%show("call setMatInit(mat, dia, vdia = -2., ndia = size(mat,1,IK))")
                        call setMatInit(mat, dia, vdia = -2., ndia = size(mat,1,IK))
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = 0. ! preset the matrix for better illustration.")
                        mat = 0.
        call disp%show("call setMatInit(mat, dia, vdia = -2., ndia = size(mat,1,IK), roff = 0_IK, coff = 0_IK)")
                        call setMatInit(mat, dia, vdia = -2., ndia = size(mat,1,IK), roff = 0_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

        call disp%skip()
        call disp%show("mat = 0. ! preset the matrix for better illustration.")
                        mat = 0.
        call disp%show("call setMatInit(mat, dia, vdia = -4., ndia = 7_IK, roff = 0_IK, coff = 3_IK)")
                        call setMatInit(mat, dia, vdia = -4., ndia = 7_IK, roff = 0_IK, coff = 3_IK)
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

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = 'AA', vdia = 'ZZ', nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppDia, vupp = 'AA', vdia = 'ZZ', nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppDia, vupp = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK)")
                        call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 1_IK, coff = 3_IK, doff = -1_IK)")
                        call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 1_IK, coff = 3_IK, doff = -1_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -4_IK)")
                        call setMatInit(mat, uppDia, vupp = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = -4_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = -3, vdia = -4, nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppDia, vupp = -3, vdia = -4, nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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
        call disp%show("call setMatInit(mat, uppDia, vupp = .true., vdia = .true., nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppDia, vupp = .true., vdia = .true., nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppDia, vupp = -3., vdia = -4., nrow = 3_IK, ncol = 7_IK, roff = 0_IK, coff = 3_IK)")
                        call setMatInit(mat, uppDia, vupp = -3., vdia = -4., nrow = 3_IK, ncol = 7_IK, roff = 0_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = 'AA', vdia = 'ZZ', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, lowDia, vlow = 'AA', vdia = 'ZZ', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, lowDia, vlow = 'AA', vdia = ['ZZ', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)")
                        call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)")
                        call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        call setMatInit(mat, lowDia, vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, lowDia, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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
        call disp%show("call setMatInit(mat, lowDia, vlow = .true., vdia = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, lowDia, vlow = .true., vdia = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, lowDia, vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)")
                        call setMatInit(mat, lowDia, vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal string matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        character(2) :: mat(10,10)

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK, doff = -3_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'HH', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK, doff = -3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        call setMatInit(mat, uppLow, vupp = 'AA', vlow = 'BB', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal integer matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = +1, vlow = -3, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLow, vupp = +1, vlow = -3, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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
        call disp%show("call setMatInit(mat, uppLow, vupp = .true., vlow = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLow, vupp = .true., vlow = .true., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct lower-diagonal real    matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLow, vupp = +1., vlow = -3., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLow, vupp = +1., vlow = -3., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

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
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = 'OO', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = 'OO', nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 4_IK, ncol = 7_IK, roff = 3_IK, coff = 0_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK, doff = -3_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'HH', vdia = ['OO', 'YY', 'WW', 'XX'], nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 0_IK, doff = -3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 7_IK, ncol = 4_IK, roff = 2_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 5_IK, ncol = 6_IK, roff = 1_IK, coff = 3_IK, doff = +1_IK)
        call disp%show("mat")
        call disp%show( mat , deliml = SK_"""" )
        call disp%skip()

        mat = ""
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)")
                        call setMatInit(mat, uppLowDia, vupp = 'AA', vlow = 'BB', vdia = 'YY', nrow = 4_IK, ncol = 7_IK, roff = 2_IK, coff = 3_IK, doff = +4_IK)
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

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = +1, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLowDia, vupp = +1, vlow = -3, vdia = -4, nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal logical matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        logical :: mat(10,10)

        mat = .false.
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = .true., vlow = .true., vdia = .false., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLowDia, vupp = .true., vlow = .true., vdia = .false., nrow = 7_IK, ncol = 4_IK, roff = 3_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct upper-lower-diagonal    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real :: mat(10,10)

        mat = 0
        call disp%skip()
        call disp%show("call setMatInit(mat, uppLowDia, vupp = +1., vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)")
                        call setMatInit(mat, uppLowDia, vupp = +1., vlow = -3., vdia = -4., nrow = 7_IK, ncol = 3_IK, roff = 0_IK, coff = 3_IK)
        call disp%show("mat")
        call disp%show( mat )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end program example