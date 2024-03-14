! Define example template macro to avoid duplications. See the example output for actual usage.
#define SET_POLY_DIV \
block; \
    TYPE(RKC), allocatable :: dividend(:), divisor(:), quorem(:); \
    dividend = DIVIDEND; \
    divisor = DIVISOR; \
    call disp%skip(); \
    call disp%show("getPolyStr(dividend)"); \
    call disp%show( getPolyStr(dividend) ); \
    call disp%show("getPolyStr(divisor)"); \
    call disp%show( getPolyStr(divisor) ); \
    call disp%show("allocate(quorem, mold = dividend)"); \
                    allocate(quorem, mold = dividend); \
    call disp%show("call setPolyDiv(dividend, divisor, quorem, lenQuo)"); \
                    call setPolyDiv(dividend, divisor, quorem, lenQuo); \
    call disp%show("lenQuo - 1 ! Degree of quotient."); \
    call disp%show( lenQuo - 1 ); \
    call disp%show("getPolyStr(quorem(1:lenQuo)) ! Quotient."); \
    call disp%show( getPolyStr(quorem(1:lenQuo)) ); \
    call disp%show("getPolyStr(quorem(lenQuo + 1 :)) ! Remainder."); \
    call disp%show( getPolyStr(quorem(lenQuo + 1 :)) ); \
    call disp%show("getPolyStr(getPolyMul(divisor, quorem(1:lenQuo))) ! Reconstruct the dividend from the divisor, Quotient, and the Remainder."); \
    call disp%show( getPolyStr(getPolyMul(divisor, quorem(1:lenQuo))) ); \
    call disp%show("getPolyStr(getPolyAdd(quorem(lenQuo + 1 :), getPolyMul(divisor, quorem(1:lenQuo)))) ! Reconstruct the dividend from the divisor, Quotient, and the Remainder."); \
    call disp%show( getPolyStr(getPolyAdd(quorem(lenQuo + 1 :), getPolyMul(divisor, quorem(1:lenQuo)))) ); \
    call disp%skip(); \
end block;


program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK32 ! all processor real and complex kinds are supported.
    use pm_io, only: display_type
    use pm_polynomial, only: getPolyMul
    use pm_polynomial, only: getPolyAdd
    use pm_polynomial, only: setPolyDiv
    use pm_polynomial, only: getPolyStr

    implicit none

    integer(IK) :: lenQuo
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#define DIVIDEND [real(RKC) :: -4., 0., -2., 1.]
#define DIVISOR [real(RKC) :: -3., 1.]
#define TYPE real
SET_POLY_DIV ! [3., 1., 1.], [5.]

#define DIVIDEND [real(RKC) :: 2., 3., 1.]
#define DIVISOR [real(RKC) :: 1., 1.]
#define TYPE real
SET_POLY_DIV

#define DIVIDEND [real(RKC) :: -42., 0., -12., 1.]
#define DIVISOR [real(RKC) :: 1., -2., 1.]
#define TYPE real
SET_POLY_DIV ! [-32, -21]

#define DIVIDEND [real(RKC) :: -42., 0., -12., 1.]
#define DIVISOR [real(RKC) :: -2., 1.]
#define TYPE real
SET_POLY_DIV

#define DIVIDEND cmplx([real(RKC) :: -4., 0., -2., 1.], -[real(RKC) :: -4., 0., -2., 1.], RKC)
#define DIVISOR cmplx([real(RKC) :: -3., 1.], -[real(RKC) :: -3., 1.], RKC)
#define TYPE complex
SET_POLY_DIV

#define DIVIDEND cmplx([real(RKC) :: 2., 3., 1.], -[real(RKC) :: 2., 3., 1.], RKC)
#define DIVISOR cmplx([real(RKC) :: 1., 1.], -[real(RKC) :: 1., 1.], RKC)
#define TYPE complex
SET_POLY_DIV

#define DIVIDEND cmplx([real(RKC) :: -42., 0., -12., 1.], -[real(RKC) :: -42., 0., -12., 1.], RKC)
#define DIVISOR cmplx([real(RKC) :: 1., -2., 1.], -[real(RKC) :: 1., -2., 1.], RKC)
#define TYPE complex
SET_POLY_DIV

#define DIVIDEND cmplx([real(RKC) :: -42., 0., -12., 1.], -[real(RKC) :: -42., 0., -12., 1.], RKC)
#define DIVISOR cmplx([real(RKC) :: -2., 1.], -[real(RKC) :: -2., 1.])
#define TYPE complex
SET_POLY_DIV

end program example