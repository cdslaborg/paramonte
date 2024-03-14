! Define an expansion procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REBIND_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('lb'); \
    call disp%show( LB ); \
    call disp%show('ub'); \
    call disp%show( UB ); \
    call disp%show('call setRebound(array, lb, ub)'); \
                    call setRebound(array, LB, UB)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

! Define an expansion with contents shifting procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REBIND_SHIFT_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('lb'); \
    call disp%show( LB ); \
    call disp%show('ub'); \
    call disp%show( UB ); \
    call disp%show('lbc'); \
    call disp%show( LBC ); \
    call disp%show('call setRebound(array, lb, ub, lbc)'); \
                    call setRebound(array, LB, UB, LBC)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

! Define an expansion/shrinkage with contents shifting and subsetting procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REBIND_SHIFT_SUBSET_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('lb'); \
    call disp%show( LB ); \
    call disp%show('ub'); \
    call disp%show( UB ); \
    call disp%show('lbc'); \
    call disp%show( LBC ); \
    call disp%show('lbcold'); \
    call disp%show( LBCOLD ); \
    call disp%show('ubcold'); \
    call disp%show( UBCOLD ); \
    call disp%show('call setRebound(array, lb, ub, lbc, lbcold, ubcold)'); \
                    call setRebound(array, LB, UB, LBC, LBCOLD, UBCOLD)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

program example

    use pm_kind, only: SK, IK
    use pm_kind, only: SKC => SK ! All kinds are supported.
    use pm_kind, only: LKC => LK ! All kinds are supported.
    use pm_kind, only: IKC => IK ! All kinds are supported.
    use pm_kind, only: CKC => CK ! All kinds are supported.
    use pm_kind, only: RKC => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRebind, only: setRebound

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand an array with specific lower and upper bounds.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB 2_IK
#define UB 10_IK
#define DECLARE character(2,SKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `integer` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB 2_IK
#define UB 10_IK
#define DECLARE integer(IKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [1, 2, 3, 4, 5, 6]
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `logical` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB 2_IK
#define UB 10_IK
#define DECLARE logical(LKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [.true., .true., .true., .true., .true., .true.]
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `complex` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB 2_IK
#define UB 10_IK
#define DECLARE complex(CKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [(1., -1.), (2., -2.), (3., -3.), (4., -4.), (5., -5.), (6., -6.)]
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `real` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB 2_IK
#define UB 10_IK
#define DECLARE real(RKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [1., 2., 3., 4., 5., 6.]
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [-1_IK, -1_IK]
#define UB [+7_IK, +7_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [-1_IK, -1_IK, 1_IK]
#define UB [+7_IK, +7_IK, 3_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(2:3,3:5,2:2)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [3, 2, 1])
REBIND_ARRAY
#undef  LB
#undef  UB
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand an array and shift its contents.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB -2_IK
#define UB 10_IK
#define LBC -2_IK
#define DECLARE character(2,SKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REBIND_SHIFT_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [+5_IK, +5_IK]
#define UB [10_IK, 10_IK]
#define LBC [7_IK, 7_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REBIND_SHIFT_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [+5_IK, +5_IK, 1_IK]
#define UB [10_IK, 10_IK, 2_IK]
#define LBC [7_IK, 7_IK, 2_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(1:2,1:3,1:1)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [3, 2, 1])
REBIND_SHIFT_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand an array and shift a subset of its contents.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB -2_IK
#define UB 10_IK
#define LBC -2_IK
#define LBCOLD 5_IK
#define UBCOLD 7_IK
#define DECLARE character(2,SKC), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REBIND_SHIFT_SUBSET_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  LBCOLD
#undef  UBCOLD
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [1_IK, 3_IK]
#define UB [2_IK, 4_IK]
#define LBC [1_IK, 3_IK]
#define LBCOLD [2_IK, 4_IK]
#define UBCOLD [3_IK, 5_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REBIND_SHIFT_SUBSET_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  LBCOLD
#undef  UBCOLD
#undef  DECLARE
#undef  CONSTRUCT

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define LB [1_IK, 2_IK, 3_IK]
#define UB [2_IK, 3_IK, 5_IK]
#define LBC [1_IK, 2_IK, 4_IK]
#define LBCOLD [1_IK, 2_IK, 2_IK]
#define UBCOLD [2_IK, 3_IK, 2_IK]
#define DECLARE character(2,SKC), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(1:2,1:3,1:2)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL"], shape = shape(array), order = [3, 2, 1])
REBIND_SHIFT_SUBSET_ARRAY
#undef  LB
#undef  UB
#undef  LBC
#undef  LBCOLD
#undef  UBCOLD
#undef  DECLARE
#undef  CONSTRUCT

end program example