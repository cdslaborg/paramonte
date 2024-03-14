program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayInit, only: setCoreHalo

    implicit none

    integer(IK) :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#define SET_ARRAY \
    block; \
        DECLARATION; \
        call disp%show("array"); \
        call disp%show( array, deliml = SK_"""" ); \
        call disp%show("core"); \
        call disp%show( core, deliml = SK_"""" ); \
        call disp%show("halo"); \
        call disp%show( halo, deliml = SK_"""" ); \
        call disp%show("coffset"); \
        call disp%show( coffset, deliml = SK_"""" ); \
        call disp%show("call setCoreHalo(array, core, halo, coffset)"); \
                        call setCoreHalo(array, core, halo, coffset); \
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
    character(21)   :: array = ""; \
    character(1 )   :: halo = "*"; \
    character(15)   :: core = "In God We Trust"; \
    integer(IK)     :: coffset = 3;
    SET_ARRAY

#define DECLARATION \
    character(41)   :: array = ""; \
    character(1 )   :: halo = "~"; \
    character(20)   :: core = "In Science We Invest"; \
    integer(IK)     :: coffset = 2;
    SET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize vector with core and halo values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define DECLARATION \
    character(5)    :: array(10) = "", core(4) = [character(5) :: "In", "God", "We", "Trust"], halo = "*****"; \
    integer(IK)     :: coffset = 3;
    SET_ARRAY

#define DECLARATION \
    character(7)    :: array(10) = "", core(4) = [character(7) :: "In", "Science", "We", "Invest"], halo = "~~~~~~~"; \
    integer(IK)     :: coffset = 3;
    SET_ARRAY

#define DECLARATION \
    integer         :: array(10) = 0, core(4) = [1, 2, 3, 4], halo = -1; \
    integer(IK)     :: coffset = 2;
    SET_ARRAY

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Initialize matrix with core and halo values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define DECLARATION \
    integer         :: array(10,10) = 0, core(2, 3) = reshape([1, 2, 3, 4, 5, 6], [2, 3]), halo = -1; \
    integer(IK)     :: coffset(rank(array)) = [5, 3];
    SET_ARRAY

#define DECLARATION \
    logical         :: array(10,10) = .false., core(3, 3) = reshape([(.false., i = 1, 9)], [3, 3]), halo = .true.; \
    integer(IK)     :: coffset(rank(array)) = [3, 5];
    SET_ARRAY

#define DECLARATION \
    complex         :: array(10,10) = 0, core(3, 3) = reshape([(i, i = 1, 9)], [3, 3]), halo = -1; \
    integer(IK)     :: coffset(rank(array)) = [3, 5];
    SET_ARRAY

#define DECLARATION \
    real            :: array(10,10) = 0, core(3, 3) = reshape([(i, i = 1, 9)], [3, 3]), halo = -1; \
    integer(IK)     :: coffset(rank(array)) = [3, 5];
    SET_ARRAY

end program example