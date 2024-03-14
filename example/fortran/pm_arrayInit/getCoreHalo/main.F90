program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayInit, only: getCoreHalo

    implicit none

    integer(IK) :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#define GET_ARRAY \
    block; \
        DECLARATION; \
        call disp%show("size"); \
        call disp%show( size, deliml = SK_"""" ); \
        call disp%show("core"); \
        call disp%show( core, deliml = SK_"""" ); \
        call disp%show("halo"); \
        call disp%show( halo, deliml = SK_"""" ); \
        call disp%show("coffset"); \
        call disp%show( coffset, deliml = SK_"""" ); \
        call disp%show("array = getCoreHalo(size, core, halo, coffset)"); \
                        array = getCoreHalo(size, core, halo, coffset); \
        call disp%show("array"); \
        call disp%show( array, deliml = SK_"""" ); \
        call disp%skip(); \
    end block;

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize scalar with core and halo values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define DECLARATION \
    integer(IK) :: coffset = 3, size = 21; \
    character(:), allocatable   :: array, core, halo; \
    core = "In God We Trust"; halo = "*";
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset = 2, size = 41; \
    character(:), allocatable   :: array, core, halo; \
    core = "In Science We Invest"; halo = "~";
    GET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize vector with core and halo values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define DECLARATION \
    integer(IK) :: coffset = 3, size = 10; \
    character(5), allocatable :: array(:), core(:), halo; \
    core = [character(5) :: "In", "God", "We", "Trust"]; \
    halo = "*****";
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset = 3, size = 10; \
    character(5), allocatable :: array(:), core(:), halo; \
    character(7) :: array(10) = "", core(4) = [character(7) :: "In", "Science", "We", "Invest"]; \
    halo = "~~~~~~~"; \
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset = 2, size = 10; \
    integer, allocatable :: array(:), core(:), halo; \
    core = [1, 2, 3, 4]; \
    halo = -1;
    GET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize matrix with core and halo values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define DECLARATION \
    integer(IK) :: coffset(2) = [5, 3], size(2) = [10, 10]; \
    integer, allocatable :: array(:,:), core(:,:), halo; \
    core = reshape([1, 2, 3, 4, 5, 6], [2, 3]); \
    halo = -1;
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset(2) = [3, 5], size(2) = [10, 10]; \
    logical, allocatable :: array(:,:), core(:,:), halo; \
    core = reshape([(.false., i = 1, 9)], [3, 3]); \
    halo = .true.;
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset(2) = [3, 5], size(2) = [10, 10]; \
    complex, allocatable :: array(:,:), core(:,:), halo; \
    core = reshape([(i, i = 1, 9)], [3, 3]); \
    halo = -1;
    GET_ARRAY

#define DECLARATION \
    integer(IK) :: coffset(2) = [3, 5], size(2) = [10, 10]; \
    real, allocatable :: array(:,:), Core(:,:), halo; \
    core = reshape([(i, i = 1, 9)], [3, 3]); \
    halo = -1;
    GET_ARRAY

end program example