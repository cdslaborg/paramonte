! Define example template macro to avoid duplications. See the example output for actual usage.
#define SET_POLY_SUB \
block; \
    TYPE(RKG), allocatable :: Lhs(:), Rhs(:), Sub(:); \
    Lhs = LHS; \
    Rhs = RHS; \
    call disp%skip(); \
    call disp%show("getPolyStr(Lhs)"); \
    call disp%show( getPolyStr(Lhs) ); \
    call disp%show("getPolyStr(Rhs)"); \
    call disp%show( getPolyStr(Rhs) ); \
    call disp%show("Sub = getPolySub(Lhs, Rhs)"); \
                    Sub = getPolySub(Lhs, Rhs); \
    call disp%show("getPolyStr(Sub)"); \
    call disp%show( getPolyStr(Sub) ); \
    call disp%skip(); \
end block;

program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all processor real and complex kinds are supported.
    use pm_io, only: display_type
    use pm_polynomial, only: getPolySub
    use pm_polynomial, only: getPolyStr

    implicit none

    integer(IK) :: lenQuo
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#define LHS [real(RKG) :: ]
#define RHS [real(RKG) :: ]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: ]
#define RHS [real(RKG) :: -1., +1.]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: -1., +1.]
#define RHS [real(RKG) :: ]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: +1., +1.]
#define RHS [real(RKG) :: -1., +1.]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: 2., 3., 1.]
#define RHS [real(RKG) :: 1., 1.]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: -42., 0., -12., 1.]
#define RHS [real(RKG) :: 1., -2., 1.]
#define TYPE real
SET_POLY_SUB

#define LHS [real(RKG) :: -42., 0., -12., 1.]
#define RHS [real(RKG) :: -2., 1.]
#define TYPE real
SET_POLY_SUB


#define LHS [complex(RKG) :: ]
#define RHS [complex(RKG) :: ]
#define TYPE complex
SET_POLY_SUB

#define LHS [complex(RKG) :: ]
#define RHS cmplx([real(RKG) :: -1., +1.], -[real(RKG) :: -1., +1.], RKG)
#define TYPE complex
SET_POLY_SUB

#define LHS cmplx([real(RKG) :: -1., +1.], -[real(RKG) :: -1., +1.], RKG)
#define RHS [complex(RKG) :: ]
#define TYPE complex
SET_POLY_SUB

#define LHS cmplx([real(RKG) :: -4., 0., -2., 1.], -[real(RKG) :: -4., 0., -2., 1.], RKG)
#define RHS cmplx([real(RKG) :: -3., 1.], -[real(RKG) :: -3., 1.], RKG)
#define TYPE complex
SET_POLY_SUB

#define LHS cmplx([real(RKG) :: 2., 3., 1.], -[real(RKG) :: 2., 3., 1.], RKG)
#define RHS cmplx([real(RKG) :: 1., 1.], -[real(RKG) :: 1., 1.], RKG)
#define TYPE complex
SET_POLY_SUB

#define LHS cmplx([real(RKG) :: -42., 0., -12., 1.], -[real(RKG) :: -42., 0., -12., 1.], RKG)
#define RHS cmplx([real(RKG) :: 1., -2., 1.], -[real(RKG) :: 1., -2., 1.], RKG)
#define TYPE complex
SET_POLY_SUB

#define LHS cmplx([real(RKG) :: -42., 0., -12., 1.], -[real(RKG) :: -42., 0., -12., 1.], RKG)
#define RHS cmplx([real(RKG) :: -2., 1.], -[real(RKG) :: -2., 1.])
#define TYPE complex
SET_POLY_SUB

end program example