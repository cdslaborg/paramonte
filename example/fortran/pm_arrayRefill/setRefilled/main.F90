! Define an expansion procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REFILL_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('fill'); \
    call disp%show( FILL , deliml = SK_"""" ); \
    call disp%show('size'); \
    call disp%show( SIZE ); \
    call disp%show('call setRefilled(array, fill, size)'); \
                    call setRefilled(array, FILL, SIZE)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

! Define an expansion with contents shifting procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REFILL_SHIFT_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('fill'); \
    call disp%show( FILL , deliml = SK_"""" ); \
    call disp%show('size'); \
    call disp%show( SIZE ); \
    call disp%show('lbc'); \
    call disp%show( LBC ); \
    call disp%show('call setRefilled(array, fill, size, lbc)'); \
                    call setRefilled(array, FILL, SIZE, LBC)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

! Define an expansion/shrinkage with contents shifting and subsetting procedure call macro to avoid duplications in the example source file. Check the output file for usage.
#define REFILL_SHIFT_SUBSET_ARRAY \
block; \
    DECLARE; \
    CONSTRUCT; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
    call disp%show('fill'); \
    call disp%show( FILL , deliml = SK_"""" ); \
    call disp%show('size'); \
    call disp%show( SIZE ); \
    call disp%show('lbc'); \
    call disp%show( LBC ); \
    call disp%show('lbcold'); \
    call disp%show( LBCOLD ); \
    call disp%show('ubcold'); \
    call disp%show( UBCOLD ); \
    call disp%show('call setRefilled(array, fill, size, lbc, lbcold, ubcold)'); \
                    call setRefilled(array, FILL, SIZE, LBC, LBCOLD, UBCOLD)  ; \
    call disp%show('array'); \
    call disp%show( array , deliml = SK_"""" ); \
    call disp%show('lbound(array)'); \
    call disp%show( lbound(array) ); \
    call disp%show('ubound(array)'); \
    call disp%show( ubound(array) ); \
end block;

program example

    use pm_kind, only: SK, IK
    use pm_kind, only: SKG => SK ! All kinds are supported.
    use pm_kind, only: LKG => LK ! All kinds are supported.
    use pm_kind, only: IKG => IK ! All kinds are supported.
    use pm_kind, only: CKG => CK ! All kinds are supported.
    use pm_kind, only: RKG => RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRefill, only: setRefilled

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand an array with specific size.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE 9_IK
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `integer` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE 9_IK
#define FILL -3_IKG
#define DECLARE integer(IKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [1, 2, 3, 4, 5, 6]
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `logical` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE 9_IK
#define FILL .false._LKG
#define DECLARE logical(LKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [.true., .true., .true., .true., .true., .true.]
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `complex` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE 9_IK
#define FILL cmplx(-9.,-9.,CKG)
#define DECLARE complex(CKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [(1., -1.), (2., -2.), (3., -3.), (4., -4.), (5., -5.), (6., -6.)]
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `real` vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE 9_IK
#define FILL -1._RKG
#define DECLARE real(RKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = [1., 2., 3., 4., 5., 6.]
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [9_IK, 9_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REFILL_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [9_IK, 9_IK, 3_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(2:3,3:5,2:2)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [3, 2, 1])
REFILL_ARRAY

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

#define SIZE 13_IK
#define LBC 5_IK
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REFILL_SHIFT_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [6_IK, 6_IK]
#define LBC [4_IK, 4_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REFILL_SHIFT_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [6_IK, 6_IK, 3_IK]
#define LBC [3_IK, 4_IK, 2_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(1:2,1:3,1:1)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [3, 2, 1])
REFILL_SHIFT_ARRAY

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

#define SIZE 13_IK
#define LBC 5_IK
#define LBCOLD 5_IK
#define UBCOLD 7_IK
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:)
#define CONSTRUCT allocate(array(3:8)); array(:) = ["AA", "BB", "CC", "DD", "EE", "FF"]
REFILL_SHIFT_SUBSET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [2_IK, 3_IK]
#define LBC [2_IK, 4_IK]
#define LBCOLD [2_IK, 4_IK]
#define UBCOLD [3_IK, 5_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:)
#define CONSTRUCT allocate(array(2:3,3:5)); array(:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF"], shape = shape(array), order = [2, 1])
REFILL_SHIFT_SUBSET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expand and shift `character` cube.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define SIZE [2_IK, 2_IK, 4_IK]
#define LBC [1_IK, 1_IK, 3_IK]
#define LBCOLD [1_IK, 2_IK, 2_IK]
#define UBCOLD [2_IK, 3_IK, 2_IK]
#define FILL SKG_"--"
#define DECLARE character(2,SKG), allocatable :: array(:,:,:)
#define CONSTRUCT allocate(array(1:2,1:3,1:2)); array(:,:,:) = reshape(["AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL"], shape = shape(array), order = [3, 2, 1])
REFILL_SHIFT_SUBSET_ARRAY

end program example