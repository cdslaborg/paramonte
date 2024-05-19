! Define example template macro to avoid duplications. See the example output for actual usage.
#define SET_POLY_ADD \
block; \
    TYPE(RKG), allocatable :: Lhs(:), Rhs(:), Add(:); \
    Lhs = LHS; \
    Rhs = RHS; \
    call disp%skip(); \
    call disp%show("getPolyStr(Lhs)"); \
    call disp%show( getPolyStr(Lhs) ); \
    call disp%show("getPolyStr(Rhs)"); \
    call disp%show( getPolyStr(Rhs) ); \
    call disp%show("allocate(Add(1 : max(size(Lhs), size(Rhs))))"); \
                    allocate(Add(1 : max(size(Lhs), size(Rhs)))); \
    call disp%show("call setPolyAdd(Add, Lhs, Rhs)"); \
                    call setPolyAdd(Add, Lhs, Rhs); \
    call disp%show("getPolyStr(Add)"); \
    call disp%show( getPolyStr(Add) ); \
    call disp%skip(); \
end block;

program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all processor real and complex kinds are supported.
    use pm_io, only: display_type
    use pm_polynomial, only: setPolyAdd
    use pm_polynomial, only: getPolyStr

    implicit none

    integer(IK) :: lenQuo
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#define LHS [real(RKG) :: +1., +1.]
#define RHS [real(RKG) :: -1., +1.]
#define TYPE real
SET_POLY_ADD

#define LHS [real(RKG) :: 2., 3., 1.]
#define RHS [real(RKG) :: 1., 1.]
#define TYPE real
SET_POLY_ADD

#define LHS [real(RKG) :: -42., 0., -12., 1.]
#define RHS [real(RKG) :: 1., -2., 1.]
#define TYPE real
SET_POLY_ADD

#define LHS [real(RKG) :: -42., 0., -12., 1.]
#define RHS [real(RKG) :: -2., 1.]
#define TYPE real
SET_POLY_ADD


#define LHS cmplx([real(RKG) :: -4., 0., -2., 1.], -[real(RKG) :: -4., 0., -2., 1.], RKG)
#define RHS cmplx([real(RKG) :: -3., 1.], -[real(RKG) :: -3., 1.], RKG)
#define TYPE complex
SET_POLY_ADD

#define LHS cmplx([real(RKG) :: 2., 3., 1.], -[real(RKG) :: 2., 3., 1.], RKG)
#define RHS cmplx([real(RKG) :: 1., 1.], -[real(RKG) :: 1., 1.], RKG)
#define TYPE complex
SET_POLY_ADD

#define LHS cmplx([real(RKG) :: -42., 0., -12., 1.], -[real(RKG) :: -42., 0., -12., 1.], RKG)
#define RHS cmplx([real(RKG) :: 1., -2., 1.], -[real(RKG) :: 1., -2., 1.], RKG)
#define TYPE complex
SET_POLY_ADD

#define LHS cmplx([real(RKG) :: -42., 0., -12., 1.], -[real(RKG) :: -42., 0., -12., 1.], RKG)
#define RHS cmplx([real(RKG) :: -2., 1.], -[real(RKG) :: -2., 1.])
#define TYPE complex
SET_POLY_ADD

end program example