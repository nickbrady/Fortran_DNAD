module number_of_variables
  implicit none
  integer, parameter :: N = 3
end module number_of_variables
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module dnadmod

    ! use user_input, only: N
  use number_of_variables
    implicit none

    integer, PARAMETER :: ndv = N*2   ! cprev_vars and dcdx_vars

    private

    real :: negative_one = -1.0
    type,public :: dual  ! make this private will create difficulty to use the
                        ! original write/read commands, hence x and dx are
                        ! variables which can be accessed using D%x and D%dx in
                        ! other units using this module in which D is defined
                        ! as type(dual).
        sequence
        real :: x  ! functional value
        real :: dx(ndv)  ! derivative
    end type dual


!******** Interfaces for operator overloading
    public assignment (=)
    interface assignment (=)
        module procedure assign_di  ! dual=integer, elemental
        module procedure assign_dr  ! dual=real, elemental
        module procedure assign_id  ! integer=dual, elemental
    end interface


    public operator (+)
    interface operator (+)
        module procedure add_d   ! +dual number, elemental
        module procedure add_dd  ! dual + dual, elemental
        module procedure add_di  ! dual + integer, elemental
        module procedure add_dr  ! dual + real, elemental
        module procedure add_id  ! integer + dual, elemental
        module procedure add_rd  ! real + dual, elemental
    end interface

    public operator (-)
    interface operator (-)
        module procedure minus_d   ! negate a dual number,elemental
        module procedure minus_dd  ! dual -dual,elemental
        module procedure minus_di  ! dual-integer,elemental
        module procedure minus_dr  ! dual-real,elemental
        module procedure minus_id  ! integer-dual,elemental
        module procedure minus_rd  ! real-dual,elemental
    end interface

    public operator (*)
    interface operator (*)
        module procedure mult_dd    ! dual*dual, elemental
        module procedure mult_di    ! dual*integer,elemental
        module procedure mult_dr    ! dual*real,elemental
        module procedure mult_id    ! integer*dual,elemental
        module procedure mult_rd    ! real*dual,elemental
    end interface

    public operator (/)
    interface operator (/)
        module procedure div_dd ! dual/dual,elemental
        module procedure div_di ! dual/integer, elemental
        module procedure div_dr ! dual/real,emental
        module procedure div_id ! integer/dual, elemental
        module procedure div_rd ! real/dual, elemental
    end interface

    public operator (**)
    interface operator (**)
        module procedure pow_i ! dual number to an integer power,elemental
        module procedure pow_r ! dual number to a real power, elemental
        module procedure pow_d ! dual number to a dual power, elemental
    end interface

    public operator (==)
    interface operator (==)
        module procedure eq_dd ! compare two dual numbers, elemental
        module procedure eq_di ! compare a dual and an integer, elemental
        module procedure eq_dr ! compare a dual and a real, elemental
        module procedure eq_id ! compare integer with a dual number, elemental
        module procedure eq_rd ! compare a real with a dual number, elemental
    end interface

    public operator (<=)
    interface operator (<=)
        module procedure le_dd  ! compare two dual numbers, elemental
        module procedure le_di  ! compare a dual and an integer, elemental
        module procedure le_dr  ! compare a dual and a real,elemental
        module procedure le_id ! compare integer with a dual number, elemental
        module procedure le_rd ! compare a real with a dual number, elemental
    end interface

    public operator (<)
    interface operator (<)
        module procedure lt_dd  !compare two dual numbers, elemental
        module procedure lt_di  !compare a dual and an integer, elemental
        module procedure lt_dr  !compare dual with a real, elemental
        module procedure lt_id ! compare integer with a dual number, elemental
        module procedure lt_rd ! compare a real with a dual number, elemental
    end interface

    public operator (>=)
    interface operator (>=)
        module procedure ge_dd ! compare two dual numbers, elemental
        module procedure ge_di ! compare dual with integer, elemental
        module procedure ge_dr ! compare dual with a real number, elemental
        module procedure ge_id ! compare integer with a dual number, elemental
        module procedure ge_rd ! compare a real with a dual number, elemental
    end interface

    public operator (>)
    interface operator (>)
        module procedure gt_dd  !compare two dual numbers, elemental
        module procedure gt_di  !compare a dual and an integer, elemental
        module procedure gt_dr  !compare dual with a real, elemental
        module procedure gt_id ! compare integer with a dual number, elemental
        module procedure gt_rd ! compare a real with a dual number, elemental
    end interface

    public operator (/=)
    interface operator (/=)
        module procedure ne_dd  !compare two dual numbers, elemental
        module procedure ne_di  !compare a dual and an integer, elemental
        module procedure ne_dr  !compare dual with a real, elemental
        module procedure ne_id ! compare integer with a dual number, elemental
        module procedure ne_rd ! compare a real with a dual number, elemental
    end interface


!------------------------------------------------
! Interfaces for intrinsic functions overloading
!------------------------------------------------
    public abs
    interface abs
        module procedure abs_d  ! absolute value of a dual number, elemental
    end interface

    public dabs
    interface dabs
        module procedure abs_d ! same as abs, used for some old fortran commands
    end interface

    public acos
    interface acos
        module procedure acos_d ! arccosine of a dual number, elemental
    end interface

    public asin
    interface asin
        module procedure asin_d ! arcsine of a dual number, elemental
    end interface

    public atan
    interface atan
        module procedure atan_d ! arctan of a dual number, elemental
    end interface

    public atan2
    interface atan2
        module procedure atan2_d ! arctan of a dual number, elemental
    end interface

    public cos
    interface cos
        module procedure cos_d ! cosine of a dual number, elemental
    end interface

    public dcos
    interface dcos
        module procedure cos_d ! cosine of a dual number, elemental
    end interface

    public dot_product
    interface dot_product
        module procedure dot_product_dd ! dot product two dual number vectors
    end interface

    public exp
    interface exp
        module procedure exp_d ! exponential of a dual number, elemental
    end interface

    public int
    interface int
        module procedure int_d ! integer part of a dual number, elemental
    end interface

    public log
    interface log
        module procedure log_d ! log of a dual number, elemental
    end interface

    public log10
    interface log10
        module procedure log10_d ! log of a dual number, elemental
    end interface

    public matmul
    interface matmul
        module procedure matmul_dd ! multiply two dual matrices
        module procedure matmul_dv ! multiply a dual matrix with a dual vector
        module procedure matmul_vd ! multiply a dual vector with a dual matrix
    end interface


    public max
    interface max
        module procedure max_dd ! max of from two to four dual numbers, elemental
        module procedure max_di ! max of a dual number and an integer, elemental
        module procedure max_dr ! max of a dual number and a real, elemental
        module procedure max_rd ! max of a real,and a dual number,  elemental
    end interface

    public dmax1
    interface dmax1
        module procedure max_dd ! max of from two to four dual numbers, elemental
    end interface

    public maxval
    interface maxval
        module procedure maxval_d ! maxval of a dual number vector
    end interface

    public min
    interface min
        module procedure min_dd ! min of from two to four dual numbers, elemental
        module procedure min_dr ! min of a dual and a real, elemental
    end interface

    public dmin1
    interface dmin1
        module procedure min_dd ! min of from two to four dual numbers, elemental
    end interface

    public minval
    interface minval
        module procedure minval_d ! obtain the maxval  of a dual number vectgor
    end interface

    public nint
    interface nint
        module procedure nint_d ! nearest integer to the argument, elemental
    end interface

    public sign
    interface  sign
      module procedure  sign_dd ! sign(a,b) with two dual numbers, elemental
      module procedure  sign_rd ! sign(a,b) with a real and a dual, elemental
    end interface

    public sin
    interface sin
        module procedure sin_d ! obtain sine of a dual number, elemental
    end interface

    public dsin
    interface dsin
        module procedure sin_d ! obtain sine of a dual number, elemental
    end interface

    public tan
    interface tan
        module procedure tan_d ! obtain sine of a dual number, elemental
    end interface

    public dtan
    interface dtan
        module procedure tan_d ! obtain sine of a dual number, elemental
    end interface

    public sqrt
    interface sqrt
        module procedure sqrt_d ! obtain the sqrt of a dual number, elemental
    end interface

    public sum
    interface sum
        module procedure sum_d ! sum a dual array
    end interface

    public maxloc
    interface maxloc
        module procedure maxloc_d ! location of max in a dual array
    end interface

    public sinh
    interface sinh
        module procedure sinh_d ! obtain sinh of a dual number, elemental
    end interface

    public cosh
    interface cosh
        module procedure cosh_d ! obtain cosh of a dual number, elemental
    end interface

    public tanh
    interface tanh
        module procedure tanh_d ! obtain tanh of a dual number, elemental
    end interface

    public asinh
    interface asinh
        module procedure asinh_d ! obtain asinh of a dual number, elemental
    end interface

    public acosh
    interface acosh
        module procedure acosh_d ! obtain acosh of a dual number, elemental
    end interface

    public atanh
    interface atanh
        module procedure atanh_d ! obtain atanh of a dual number, elemental
    end interface

contains

!*********Begin: functions/subroutines for overloading operators

!******* Begin: (=)
!---------------------

    !-----------------------------------------
    ! dual = integer
    ! <u, du> = <i, 0>
    !-----------------------------------------
    elemental subroutine assign_di(u, i)
         type(dual), intent(out) :: u
         integer, intent(in) :: i

         u%x = real(i)  ! This is faster than direct assignment
         u%dx = 0.0

    end subroutine assign_di


    !-----------------------------------------
    ! dual = real(double)
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dr(u, r)
        type(dual), intent(out) :: u
        real, intent(in) :: r

        u%x = r
        u%dx = 0.0

    end subroutine assign_dr


    !-----------------------------------------
    ! integer = dual
    ! i = <u, du>
    !-----------------------------------------
    elemental subroutine assign_id(i, v)
         type(dual), intent(in) :: v
         integer, intent(out) :: i

         i = int(v%x)

    end subroutine assign_id

!******* end: (=)
!---------------------


!******* Begin: (+)
!---------------------

    !-----------------------------------------
    ! Unary positive
    ! <res, dres> = +<u, du>
    !-----------------------------------------
    elemental function add_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res

         res = u  ! Faster than assigning component wise

    end function add_d


    !-----------------------------------------
    ! dual + dual
    ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
    !-----------------------------------------
    elemental function add_dd(u, v) result(res)
         type(dual), intent(in) :: u, v
         type(dual) :: res

         res%x = u%x + v%x
         res%dx = u%dx + v%dx

    end function add_dd


    !-----------------------------------------
    ! dual + integer
    ! <res, dres> = <u, du> + i = <u + i, du>
    !-----------------------------------------
    elemental function add_di(u, i) result(res)
         type(dual), intent(in) :: u
         integer, intent(in) :: i
         type(dual) :: res

         res%x = real(i) + u%x
         res%dx = u%dx

    end function add_di


    !-----------------------------------------
    ! dual + double
    ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
    !-----------------------------------------
    elemental function add_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        res%x = r + u%x
        res%dx = u%dx

    end function add_dr


    !-----------------------------------------
    ! integer + dual
    ! <res, dres> = <i, 0> + <v, dv> = <i + v, dv>
    !-----------------------------------------
    elemental function add_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(i) + v%x
        res%dx = v%dx

    end function add_id


    !-----------------------------------------
    ! double + dual
    ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
    !-----------------------------------------
    elemental function add_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r + v%x
        res%dx = v%dx

    end function add_rd

!******* end: (+)
!---------------------


!******* Begin: (-)
!---------------------

    !-------------------------------------------------
    ! negate a dual
    ! <res, dres> = -<u, du>
    !-------------------------------------------------
    elemental function minus_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = -u%x
        res%dx = -u%dx

    end function minus_d


    !-------------------------------------------------
    ! dual - dual
    ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
    !-------------------------------------------------
    elemental function minus_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x - v%x
        res%dx = u%dx - v%dx

    end function minus_dd

    !-------------------------------------------------
    ! dual - integer
    ! <res, dres> = <u, du> - i = <u - i, du>
    !-------------------------------------------------
    elemental function minus_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        res%x = u%x - real(i)
        res%dx = u%dx

    end function minus_di


    !-------------------------------------------------
    ! dual - double
    ! <res, dres> = <u, du> - r = <u - r, du>
    !-------------------------------------------------
    elemental function minus_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real,intent(in) :: r
        type(dual) :: res

        res%x = u%x - r
        res%dx = u%dx

    end function minus_dr


    !-------------------------------------------------
    ! integer - dual
    ! <res, dres> = i - <v, dv> = <i - v, -dv>
    !-------------------------------------------------
    elemental function minus_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(i) - v%x
        res%dx = -v%dx

    end function minus_id


    !-------------------------------------------------
    ! double - dual
    ! <res, dres> = r - <v, dv> = <r - v, -dv>
    !-------------------------------------------------
    elemental function minus_rd(r, v) result(res)
         real, intent(in) :: r
         type(dual), intent(in) :: v
         type(dual) :: res

        res%x = r - v%x
        res%dx = -v%dx

    end function minus_rd

!******* end: (-)
!---------------------


!******* BEGIN: (*)
!---------------------

    !----------------------------------------
    ! dual * dual
    ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
    !----------------------------------------
    elemental function mult_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x * v%x
        res%dx = u%x * v%dx + v%x * u%dx

    end function mult_dd


    !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * i = <u * i, du * i>
    !-----------------------------------------
    elemental function mult_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: r

        r = real(i)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di

    !-----------------------------------------
    ! dual * double
    ! <res, dres> = <u, du> * r = <u * r, du * r>
    !----------------------------------------
    elemental function mult_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        res%x = u%x * r
        res%dx = u%dx * r

    end function mult_dr


    !-----------------------------------------
    ! integer * dual
    ! <res, dres> = i * <v, dv> = <i * v, i * dv>
    !-----------------------------------------
    elemental function mult_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: r

        r = real(i)
        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_id


    !-----------------------------------------
    ! double * dual
    ! <res, dres> = r * <v, dv> = <r * v, r * dv>
    !-----------------------------------------
    elemental function mult_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_rd

!******* end: (*)
!---------------------


!******* BEGIN: (/)
!---------------------

    !-----------------------------------------
    ! dual / dual
    ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
    !-----------------------------------------
    elemental function div_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = u%x * inv
        res%dx = (u%dx - res%x * v%dx) * inv

    end function div_dd


    !-----------------------------------------
    ! dual / integer
    ! <res, dres> = <u, du> / i = <u / i, du / i>
    !-----------------------------------------
    elemental function div_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: inv

        inv = 1.0 / real(i)
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_di


    !-----------------------------------------
    ! dual / double
    ! <res, dres> = <u, du> / r = <u / r, du / r>
    !----------------------------------------
    elemental function div_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual):: res

        real :: inv

        inv = 1.0 / r
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_dr


    !-----------------------------------------
    ! integer / dual
    ! <res, dres> = i / <v, dv> = <i / v, -i / v^2 * du>
    !-----------------------------------------
    elemental function div_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = real(i) * inv
        res%dx = -res%x * inv * v%dx

    end function div_id


    !-----------------------------------------
    ! double / dual
    ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
    !-----------------------------------------
    elemental function div_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = r * inv
        res%dx = -res%x * inv * v%dx

    end function div_rd

!******* end: (/)
!---------------------

!******* BEGIN: (**)
!---------------------

    !-----------------------------------------
    ! power(dual, integer)
    ! <res, dres> = <u, du> ^ i = <u ^ i, i * u ^ (i - 1) * du>
    !-----------------------------------------
    elemental function pow_i(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: pow_x

        pow_x = u%x ** (i - 1)
        res%x = u%x * pow_x
        res%dx = real(i) * pow_x * u%dx

    end function pow_i

    !-----------------------------------------
    ! power(dual, double)
    ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
    !-----------------------------------------
    elemental function pow_r(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        real :: pow_x

        pow_x = u%x ** (r - 1.0)
        res%x = u%x * pow_x
        res%dx = r * pow_x * u%dx

    end function pow_r

    !-----------------------------------------
    ! POWER dual numbers to a dual power
    ! <res, dres> = <u, du> ^ <v, dv>
    !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
    !-----------------------------------------
    elemental function pow_d(u, v) result(res)
        type(dual), intent(in)::u, v
        type(dual) :: res

        res%x = u%x ** v%x
        res%dx = res%x * (v%x / u%x * u%dx + log(u%x) * v%dx)

    end function pow_d

!******* end: (**)
!---------------------


!******* BEGIN: (==)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x == rhs%x)

    end function eq_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x == real(rhs))

    end function eq_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical::res

        res = (lhs%x == rhs)

    end function eq_dr


    !-----------------------------------------
    ! compare an integer with a dual,
    ! simply compare the functional value.
    !----------------------------------------
    elemental function eq_id(lhs, rhs) result(res)
         integer, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_rd

!******* end: (==)
!---------------------


!******* BEGIN: (<=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function le_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x <= rhs%x)

    end function le_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_dr(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         real, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_dr


    !-----------------------------------------
    ! compare a dual number with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_id(i, rhs) result(res)
         integer, intent(in) :: i
         type(dual), intent(in) :: rhs
         logical :: res

         res = (i <= rhs%x)

    end function le_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs <= rhs%x)

    end function le_rd

!******* end: (<=)
!---------------------

!******* BEGIN: (<)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function lt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x < rhs%x)

    end function lt_dd

    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function lt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function lt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function lt_id(i, rhs) result(res)
         integer, intent(in) :: i
         type(dual), intent(in) :: rhs
         logical :: res

         res = (i < rhs%x)

    end function lt_id


    !-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    elemental function lt_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs < rhs%x)

    end function lt_rd

!******* end: (<)
!---------------------

!******* BEGIN: (>=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function ge_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x >= rhs%x)

    end function ge_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ge_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ge_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ge_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i >= rhs%x)

    end function ge_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ge_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs >= rhs%x)

    end function ge_rd

!******* end: (>=)
!---------------------

!******* BEGIN: (>)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x > rhs%x)

    end function gt_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function gt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function gt_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i > rhs%x)

    end function gt_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function gt_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs > rhs%x)

    end function gt_rd

!******* end: (>)
!---------------------

!******* BEGIN: (/=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x /= rhs%x)

    end function ne_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ne_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ne_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i /= rhs%x)

    end function ne_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ne_rd(lhs, rhs) result(res)
        real, intent(in) :: lhs
        type(dual), intent(in) :: rhs
        logical :: res

        res = (lhs /= rhs%x)

    end function ne_rd

!******* end: (/=)
!---------------------

    !---------------------------------------------------
    ! Absolute value of dual numbers
    ! <res, dres> = abs(<u, du>) = <abs(u), du * sign(u)>
    !---------------------------------------------------
    elemental function abs_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res
         integer :: i

         if(u%x > 0) then
            res%x = u%x
            res%dx = u%dx
         else if (u%x < 0) then
            res%x = -u%x
            res%dx = -u%dx
         else
            res%x = 0.0
            do i = 1, ndv
                if (u%dx(i) .eq. 0.0) then
                    res%dx(i) = 0.0
                else
                    res%dx(i) = set_NaN()
                end if
            end do
         endif

    end function abs_d


    !-----------------------------------------
    ! ACOS of dual numbers
    ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function acos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acos(u%x)
        if (u%x == 1.0 .or. u%x == -1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = -u%dx / sqrt(1.0 - u%x**2)
        end if

    end function acos_d


    !-----------------------------------------
    ! ASIN of dual numbers
    ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function asin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asin(u%x)
        if (u%x == 1.0 .or. u%x == -1.0) then
            res%dx = set_NaN()  ! Undefined derivative
        else
            res%dx = u%dx / sqrt(1.0 - u%x**2)
        end if

    end function asin_d


    !-----------------------------------------
    ! ATAN of dual numbers
    ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
    !----------------------------------------
    elemental function atan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atan(u%x)
        res%dx = u%dx / (1.0 + u%x**2)

    end function atan_d


    !-----------------------------------------
    ! ATAN2 of dual numbers
    ! <res, dres> = atan2(<u, du>, <v, dv>)
    !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
    !----------------------------------------
    elemental function atan2_d(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real :: usq_plus_vsq

        res%x = atan2(u%x, v%x)

        usq_plus_vsq = u%x**2 + v%x**2
        res%dx = v%x / usq_plus_vsq * u%dx - u%x / usq_plus_vsq * v%dx

    end function atan2_d


    !-----------------------------------------
    ! COS of dual numbers
    ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
    !----------------------------------------
    elemental function cos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cos(u%x)
        res%dx = -sin(u%x) * u%dx

    end function cos_d


    !-----------------------------------------
    ! DOT PRODUCT two dual number vectors
    ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
    !-----------------------------------------
    function dot_product_dd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:)
        type(dual) :: res

        integer :: i

        res%x = dot_product(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = dot_product(u%x, v%dx(i)) + dot_product(v%x, u%dx(i))
        end do

    end function dot_product_dd


    !-----------------------------------------
    ! EXPONENTIAL OF dual numbers
    ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
    !-----------------------------------------
    elemental function exp_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: exp_x

        exp_x = exp(u%x)
        res%x = exp_x
        res%dx = u%dx * exp_x

    end function exp_d


    !-----------------------------------------
    ! Convert dual to integer
    ! i = int(<u, du>) = int(u)
    !----------------------------------------
    elemental function int_d(u) result(res)
         type(dual), intent(in) :: u
         integer :: res

         res = int(u%x)

    end function int_d


    !-----------------------------------------
    ! LOG OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log(<u, du>) = <log(u), du / u>
    !----------------------------------------
    elemental function log_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: inv

        inv = 1.0 / u%x
        res%x = log(u%x)
        res%dx = u%dx * inv

    end function log_d


    !-----------------------------------------
    ! LOG10 OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    elemental function log10_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: inv

        inv = 1.0 / (u%x * log(10.0))
        res%x = log10(u%x)
        res%dx = u%dx * inv

    end function log10_d


    !-----------------------------------------
    ! MULTIPLY two dual number matrices
    ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
    !----------------------------------------
    function matmul_dd(u,v) result(res)
        type(dual), intent(in) :: u(:,:), v(:,:)
        type(dual) :: res(size(u,1), size(v,2))

        integer :: i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_dd


    !-----------------------------------------
    ! MULTIPLY a dual number matrix with a dual number
    ! vector
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_dv(u, v) result(res)
        type(dual), intent(in) :: u(:,:), v(:)
        type(dual) :: res(size(u,1))
        integer :: i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_dv


    !-----------------------------------------
    ! MULTIPLY a dual vector with a  dual matrix
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_vd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:,:)
        type(dual) :: res(size(v, 2))
        integer::i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_vd

    !-----------------------------------------
    ! Obtain the max of 2 to 5 dual numbers
    !----------------------------------------
    elemental function max_dd(val1, val2, val3, val4,val5) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4,val5
        type(dual) :: res

        if (val1%x > val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x < val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x < val4%x) res = val4
        endif
        if(present(val5))then
           if(res%x < val5%x) res = val5
        endif

    end function max_dd


    !-----------------------------------------
    ! Obtain the max of a dual number and an integer
    !----------------------------------------
    elemental function max_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        if (u%x > i) then
            res = u
        else
            res = i
        endif

    end function max_di

    !-----------------------------------------
    ! Obtain the max of a dual number and a real number
    !----------------------------------------
    elemental function max_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        if (u%x > r) then
            res = u
        else
            res = r
        endif

    end function max_dr


    !---------------------------------------------------
    ! Obtain the max of a real and a dual
    !---------------------------------------------------
     elemental function max_rd(n, u) result(res)
        real, intent(in) :: n
        type(dual), intent(in) :: u
        type(dual) :: res

        if (u%x > n) then
            res = u
        else
            res = n
        endif

    end function max_rd


    !-----------------------------------------
    ! Obtain the max value of vector u
    !----------------------------------------
    function maxval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=maxloc(u%x)
        res=u(iloc(1))

    end function maxval_d


    !-----------------------------------------
    ! Obtain the min of 2 to 4 dual numbers
    !----------------------------------------
    elemental function min_dd(val1, val2, val3, val4) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4
        type(dual) :: res

        if (val1%x < val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x > val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x > val4%x) res = val4
        endif

    end function min_dd


    !-----------------------------------------
    ! Obtain the min of a dual and a double
    !----------------------------------------
    elemental function min_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        if (u%x < r) then
            res = u
        else
            res = r
        endif

    end function min_dr


  !-----------------------------------------
    ! Obtain the min value of vector u
    !----------------------------------------
    function minval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=minloc(u%x)
        res=u(iloc(1))

    end function minval_d


    !------------------------------------------------------
    !Returns the nearest integer to u%x, ELEMENTAL
    !------------------------------------------------------
    elemental function nint_d(u) result(res)
        type(dual), intent(in) :: u
        integer :: res

        res=nint(u%x)

    end function nint_d


    !----------------------------------------------------------------
    ! SIGN(a,b) with two dual numbers as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_dd(val1, val2) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual) :: res

        if (val2%x < 0.0) then
            res = -abs(val1)
        else
            res =  abs(val1)
        endif

     end function sign_dd


    !----------------------------------------------------------------
    ! SIGN(a,b) with one real and one dual number as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_rd(val1, val2) result(res)
        real, intent(in) :: val1
        type(dual), intent(in) :: val2
        type(dual) :: res

        if (val2%x < 0.0) then
            res = -abs(val1)
        else
            res = abs(val1)
        endif

     end function sign_rd


    !-----------------------------------------
    ! SIN of dual numbers
    ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
    !----------------------------------------
    elemental function sin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sin(u%x)
        res%dx = cos(u%x) * u%dx

    end function sin_d


    !-----------------------------------------
    ! TAN of dual numbers
    ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
    !----------------------------------------
    elemental function tan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tan(u%x)
        res%dx = u%dx / cos(u%x)**2

    end function tan_d


    !-----------------------------------------
    ! SQRT of dual numbers
    ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
    !----------------------------------------
    elemental function sqrt_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res
        integer :: i

        res%x = sqrt(u%x)

        if (res%x .ne. 0.0) then
            res%dx = 0.5 * u%dx / res%x
        else
            do i = 1, ndv
                if (u%dx(i) .eq. 0.0) then
                    res%dx(i) = 0.0
                else
                    res%dx(i) = set_NaN()
                end if
            end do
        end if

    end function sqrt_d


    !-----------------------------------------
    ! Sum of a dual array
    !-----------------------------------------
    function sum_d(u) result(res)
        type(dual), intent(in) :: u(:)
        type(dual) :: res
        integer :: i

        res%x = sum(u%x)
        do i = 1, ndv
            res%dx(i) = sum(u%dx(i))
        end do

    end function sum_d


    !-----------------------------------------
    ! Find the location of the max value in an
    ! array of dual numbers
    !-----------------------------------------
    function maxloc_d(array) result(ind)
        type(dual), intent(in) :: array(:)
        integer :: ind(1)

        ind = maxloc(array%x)

    end function maxloc_d


    elemental function set_NaN() result(res)
        real :: res

        res = sqrt(negative_one)

    end function set_NaN


    !-----------------------------------------
    ! Hyperbolic functions: sinh, cosh, tanh
    ! and their inverses: asinh, acosh, atanh
    !-----------------------------------------
    !-----------------------------------------
    ! SINH OF dual numbers
    ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
    !-----------------------------------------
    elemental function sinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sinh(u%x)
        res%dx = u%dx * cosh(u%x)

    end function sinh_d

    !-----------------------------------------
    ! COSH OF dual numbers
    ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
    !-----------------------------------------
    elemental function cosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cosh(u%x)
        res%dx = u%dx * sinh(u%x)

    end function cosh_d

    !-----------------------------------------
    ! TANH OF dual numbers
    ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0/cosh(u)**2 * du>
    !-----------------------------------------
    elemental function tanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tanh(u%x)
        res%dx = u%dx * 1.0/cosh(u%x)**2

    end function tanh_d

    !-----------------------------------------
    ! ASINH OF dual numbers
    ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
    !-----------------------------------------
    elemental function asinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asinh(u%x)
        res%dx = u%dx * 1.0/sqrt(u%x**2 + 1.0)

    end function asinh_d

    !-----------------------------------------
    ! ACOSH OF dual numbers
    ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
    !-----------------------------------------
    elemental function acosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acosh(u%x)
        if (u%x <= 1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0/sqrt(u%x**2 - 1.0)
        end if

    end function acosh_d

    !-----------------------------------------
    ! ATAHN OF dual numbers
    ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
    !-----------------------------------------
    elemental function atanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atanh(u%x)
        if (abs(u%x) >= 1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0/(1.0 - u%x**2)
        end if

    end function atanh_d


end module dnadmod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module user_input
  use number_of_variables
  use dnadmod
  implicit none

  integer, parameter :: NJ = 102                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 30 * 3.6e3     ! Number of time steps
  integer, parameter :: maxIterations = 1e6
  real               :: delT = 3.6e9                       ! size of timestep [s]
  real               :: time                             ! [s]
  real, parameter    :: xmax = 1e-0                      ! 1 cm
  real, parameter    :: R_1  = 0
  real, parameter    :: R_NJ = R_1 + xmax

  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real :: cbulk = 0.001                       ! concentration [mol/cm3]
  real :: Delta_V
  real :: kappa = 1e-1

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module variables
  use user_input
  implicit none

  real, dimension(N, NJ)  :: cprev, delC
  real, dimension(NJ)     :: xx, delX

  ! ABDGXY_VARS
  real, dimension(N, N)     :: A, B, X, Y
  real, dimension(N, 2*N+1) :: D
  real, dimension(N)        :: G

  ! BAND and MATINV variables
  real, dimension(N, N+1, NJ) :: E
  real, dimension(N)          :: ID
  integer :: NP1

end module variables
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module echem_transport_mod
  use dnadmod
  use user_input

  implicit none

  real :: alpha_a = 0.5
  real :: alpha_c = 0.5
  real :: k_rxn   = 1.0e-7

  interface Soln_Pot
    module procedure Soln_Pot_dual
    module procedure Soln_Pot_real
  end interface Soln_Pot

  interface exchange_current_density
    module procedure exchange_current_density_dual
    module procedure exchange_current_density_real
  end interface exchange_current_density

  interface OpenCircuitPotential
    module procedure OpenCircuitPotential_dual
    module procedure OpenCircuitPotential_real
  end interface OpenCircuitPotential

  interface Butler_Volmer_Rxn
    ! 2 variables: c0, Phi_1
    module procedure Butler_Volmer_Rxn_dual
    module procedure Butler_Volmer_Rxn_real

    ! 3 variables: c0, Phi_1, Phi_2
    module procedure Butler_Volmer_Rxn3_dual
    module procedure Butler_Volmer_Rxn3_real
  end interface Butler_Volmer_Rxn

  contains

    function Diffusion_Coef(c_vars_dual)    Result(diff_coef_)
      type(dual)                            :: diff_coef_
      type(dual), dimension(N), intent(in)  :: c_vars_dual

      type(dual) :: c0

      real :: D_Li_0 = 1e-5
      real :: D_Li_min = 1e-7
      real :: c0_0 = 9.0e-4

      if ('Constant' == 'Constant') then
        diff_coef_%dx = 0.0
        diff_coef_%x  = D_Li_0

      else if ('Constant' == 'Variable') then
        c0 = c_vars_dual(1)
        diff_coef_ = 1.5*D_Li_0 * exp(-c0 / c0_0) + D_Li_min
      end if

    end function Diffusion_Coef


!*******************************************************************************
    function Soln_Pot_dual(c0)                   result(SolutionPotential)
      type(dual)                :: SolutionPotential
      type(dual), intent(in)    :: c0
      real                      :: U_0

      U_0 =  0.0

      SolutionPotential = U_0 + Rigc*Temp/Fconst * LOG( c0/cbulk )

    end function Soln_Pot_dual


    function Soln_Pot_real(c0)                   result(SolutionPotential)
      real                :: SolutionPotential
      real, intent(in)    :: c0
      real                :: U_0

      U_0 =  0.0

      SolutionPotential = U_0 + Rigc*Temp/Fconst * LOG( c0/cbulk )

    end function Soln_Pot_real
!*******************************************************************************

!*******************************************************************************
    function exchange_current_density_dual(c0)  result(ex_curr)
      type(dual)              :: ex_curr
      type(dual), intent(in)  :: c0

      ex_curr = Fconst * k_rxn * c0**alpha_a

    end function exchange_current_density_dual

    function exchange_current_density_real(c0)  result(ex_curr)
      real              :: ex_curr
      real, intent(in)  :: c0

      ex_curr = Fconst * k_rxn * c0**alpha_a

    end function exchange_current_density_real
!*******************************************************************************

!*******************************************************************************
    function OpenCircuitPotential_dual(c0)     result(OCP)
      type(dual)              :: OCP
      type(dual), intent(in)  :: c0

      real, parameter :: U_ref = 3.0, c0_ref = 1e-3

      OCP = U_ref + Rigc*Temp/Fconst * log(c0 / c0_ref)

    end function OpenCircuitPotential_dual


    function OpenCircuitPotential_real(c0)     result(OCP)
      real              :: OCP
      real, intent(in)  :: c0

      real, parameter :: U_ref = 3.0, c0_ref = 1e-3

      OCP = U_ref + Rigc*Temp/Fconst * log(c0 / c0_ref)

    end function OpenCircuitPotential_real
!*******************************************************************************

!*******************************************************************************
! Current is positive when eta > 0; negative when eta < 0
! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
    function Butler_Volmer_Rxn_dual(c0, Phi_1)    result(i_BV)

      type(dual)              :: i_BV
      type(dual), intent(in)  :: c0, Phi_1
      type(dual)              :: eta
      type(dual)              :: i_ex_curr
      type(dual)              :: U_ocp

      U_ocp = OpenCircuitPotential(c0)
      i_ex_curr = exchange_current_density(c0)
      eta = Phi_1 - U_ocp

      i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                          - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

    end function Butler_Volmer_Rxn_dual

    function Butler_Volmer_Rxn_real(c0, Phi_1)    result(i_BV)
      real              :: i_BV
      real, intent(in)  :: c0, Phi_1
      real              :: eta
      real              :: i_ex_curr
      real              :: U_ocp

      U_ocp = OpenCircuitPotential(c0)
      i_ex_curr = exchange_current_density(c0)
      eta = Phi_1 - U_ocp

      i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                          - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

    end function Butler_Volmer_Rxn_real

    function Butler_Volmer_Rxn3_dual(c0, Phi_1, Phi_2)    result(i_BV)
      type(dual)              :: i_BV
      type(dual), intent(in)  :: c0, Phi_1, Phi_2
      type(dual)              :: eta
      type(dual)              :: i_ex_curr
      type(dual)              :: U_ocp

      U_ocp = OpenCircuitPotential(c0)
      i_ex_curr = exchange_current_density(c0)
      eta = Phi_1 - Phi_2 - U_ocp

      i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                          - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

    end function Butler_Volmer_Rxn3_dual

    function Butler_Volmer_Rxn3_real(c0, Phi_1, Phi_2)    result(i_BV)
      real              :: i_BV
      real, intent(in)  :: c0, Phi_1, Phi_2
      real              :: eta
      real              :: i_ex_curr
      real              :: U_ocp

      U_ocp = OpenCircuitPotential(c0)
      i_ex_curr = exchange_current_density(c0)
      eta = Phi_1 - Phi_2 - U_ocp

      i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                          - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

    end function Butler_Volmer_Rxn3_real
!*******************************************************************************

end module echem_transport_mod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use echem_transport_mod

  implicit none

  real :: c0, Phi_1, Phi_2
  real :: Phi_1_1, Phi_1_NJ
  real :: i_rxn, i_BV, i_BV_kappa

  character(len=65) :: header_fmt = '(1(A7,1X), 20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(I7,1X), 20(ES15.5,1X))'

  real :: t_write

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/10.0           ! 3600 s = 1 hour

    call write_to_screen(it)
    call write_positional_information_to_file(it)

  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it == 0) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Iters', 'DeltaV', 'Conc',  'Phi_1', 'Phi_2', 'DeltaV', 'i_rxn', 'i_kappa'
      write(*, header_fmt) '###',   'Volts' , 'mol/L', 'Volts', 'Volts', 'Volts',  'A/cm2', 'A/cm2'

    else
      c0       = cprev(1,1)*1e3
      Phi_1_1  = cprev(2,1)
      Phi_1_NJ = cprev(2,NJ)
      Phi_2    = cprev(3,1)
      i_rxn    = Butler_Volmer_Rxn(cprev(1,1), Phi_1_1, Phi_2)
      i_BV_kappa = Butler_Volmer_Rxn(cbulk, Phi_1_1, Phi_2)

      write(*, data_fmt) it, Delta_V, c0,  Phi_1_1, Phi_2, Phi_1_1 - Phi_1_NJ, i_rxn, i_BV_kappa

    end if

  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it == 0) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Iters', 'DeltaV', 'Position', 'Conc', 'Phi_1', 'Phi_2', 'i_rxn', 'i_BV', 'i_kappa'
      write(56, header_fmt) '###',   'Volts' , 'cm'      , 'mol/L', 'Volts', 'Volts', 'A/cm2', 'A/cm2', 'A/cm2'
                                                      !
    else                                              !

      Phi_1 = cprev(2,1)
      Phi_2 = cprev(3,1)
      i_rxn = Butler_Volmer_Rxn(cprev(1,1), Phi_1, Phi_2)                                                !
      i_BV  = Butler_Volmer_Rxn(cbulk, Phi_1, 0.0)           ! kinetic limitation
      i_BV_kappa = Butler_Volmer_Rxn(cbulk, Phi_1, Phi_2)

      do j = 1, NJ                                    !
        c0    = cprev(1,j) * 1e3                      !
                                                      !
        write(56, data_fmt) it, Delta_V, xx(j), c0, Phi_1, Phi_2, i_rxn, i_BV, i_BV_kappa

      end do
    end if

  end subroutine write_positional_information_to_file

end module write_data_mod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module GOV_EQNS
  use user_input
  use variables
  use echem_transport_mod
  use dnadmod

  implicit none

  ! initialize the fill_mat variables - shared with auto_fill and ABDGXY
  real                  :: alphaW, alphaE, betaW, betaE
  real, dimension(N,N)  :: dW, dE, fW, fE, rj
  real, dimension(N)    :: smG
  real, dimension(NJ)   :: Cntrl_Vol
  real, dimension(2, NJ):: Crx_Area        ! dimesion = (2,NJ) because east-west
                                           ! cross-sectional area
  type(dual) :: c0, Phi_1, Phi_2

contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   N_0 = -D dc/dx
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: diff_

    type(dual) :: dc0dx

    c0       = c_vars_dual(1)
    dc0dx    = dcdx_vars_dual(1)

    diff_ = Diffusion_Coef(c_vars_dual)

    Flux_(1) = -diff_ * dc0dx
    Flux_(2) = 0.0
    Flux_(3) = 0.0

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  no reaction
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    c0    = c_vars_dual(1)
    Phi_1 = c_vars_dual(2)
    Phi_2 = c_vars_dual(3)

    Rxn_(1) = 0.0
    Rxn_(2) = Phi_1 - 0.0
    Rxn_(3) = Phi_2 - 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) dc0/dt
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c0    = c_vars_dual(1)

    Accum_(1) = c0/delT
    Accum_(2) = 0.0
    Accum_(3) = 0.0

    Accum_%x = 0.0

  end function ACCUM
! ******************************************************************************


! ******************************************************************************
! **************************** BOUNDARY CONDITIONS *****************************
! ******************************************************************************
! boundary conditions need to be written as homogeneous equations
! i.e.    c = 1.0   -->      c - 1.0 = 0.0    -->   BC_WEST_(N) = c - 1.0
! i.e. flux = 1.0   -->   flux - 1.0 = 0.0    -->   BC_WEST_(N) = flux - 1.0
! ********** special case **********
! i.e. dc/dt = rxn  -->   dc/dt - rxn = 0.0
! dc_n   = c_vars_dual(n)
! dc_n%x = 0.0      -->   dc        -->           BC_WEST_(2) = dc_n/delT - rxn
  function Boundary_WEST (c_vars_dual, dcdx_vars_dual) result(BC_WEST_)
    ! (1)   c0 = cbulk
    type(dual), dimension(N)               :: BC_WEST_, flux_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    real       :: U_ocp

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    c0        = c_vars_dual(1)
    Phi_1     = c_vars_dual(2)
    Phi_2     = c_vars_dual(3)
    ! Phi_1_NJ = Phi_1_1 - Delta_V

    ! BC_WEST_(1) = flux_temp(1) - i_applied_cm2/Fconst
    ! BC_WEST_(1) = flux_temp(1) - Butler_Volmer_Rxn(c0, Phi_1)/Fconst
    ! conservation of current at each electrode
    ! Delta_V = Phi_1_NJ - Phi_1_1 OCP(cbulk) - IR_drop
    ! Phi_1_NJ = U_ocp
    ! Phi_1_1  = - IR_drop
    ! Phi_1_1  = U_ocp - Delta_V - IR_drop
    ! Phi_1_1 - U_ocp + Delta_V + IR_drop = 0.0
    U_ocp   = OpenCircuitPotential(cbulk)
    BC_WEST_(2) = Phi_1 - U_ocp - Delta_V! + kappa * Butler_Volmer_Rxn(c0, Phi_1)

    ! eta = Phi_1 - Phi_2 - U_ocp
    ! Phi_2 = i_BV * delX/kappa
    ! Delta_V = Phi_1 - OCP(cbulk)
    ! 0.0 = Phi_1 - OCP(cbulk) - Delta_V
    ! BC_WEST_(2) = Phi_1 - U_ocp - Delta_V
    ! BC_WEST_(3) = c_vars_dual(3) - Butler_Volmer_Rxn(c0, Phi_1) * xmax/kappa

    BC_WEST_(1) = flux_temp(1)   - Butler_Volmer_Rxn(c0, Phi_1, Phi_2)/Fconst
    BC_WEST_(3) = c_vars_dual(3) - Butler_Volmer_Rxn(c0, Phi_1, Phi_2) * xmax/kappa



  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   c0 = 0.1 * cbulk
    type(dual), dimension(N)               :: BC_EAST_, flux_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    real       :: U_ocp

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    c0        = c_vars_dual(1)
    Phi_1     = c_vars_dual(2)
    Phi_2     = c_vars_dual(3)
    U_ocp     = OpenCircuitPotential(cbulk)

    ! BC_EAST_(1) = flux_temp(1) - i_applied_cm2/Fconst
    BC_EAST_(1) = c0 - cbulk
    BC_EAST_(2) = Phi_1 - U_ocp
    BC_EAST_(3) = Phi_2 - 0.0

  end function Boundary_EAST
! ******************************************************************************


  ! ****************************************************************************
  ! ******* functions to convert concentration and dcdx variables to dual ******
  ! ----------------------------------------------------------------------------
  function c_to_dual(c_vars) result(c_dual)
    type(dual), dimension(N)       :: c_dual
    real, dimension(N), intent(in) :: c_vars
    real, dimension(N*2)            :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(ic) = 1.0

      c_dual(ic) = dual(c_vars(ic), dx_array)
    end do
  end function c_to_dual


  function dcdx_to_dual(dcdx) result(dcdx_dual)
    type(dual), dimension(N)       :: dcdx_dual
    real, dimension(N), intent(in) :: dcdx
    real, dimension(N*2)            :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(N+ic) = 1.0

      dcdx_dual(ic) = dual(dcdx(ic), dx_array)
    end do

  end function dcdx_to_dual
! ******************************************************************************


  ! ****************************************************************************
  ! ** functions to define the node control volumes and cross-sectional area ***
  ! ----------------------------------------------------------------------------
  function Control_Volume(Geometry) result(ctrl_vol)
    use variables, only: xx, delX

    real, dimension(NJ) :: ctrl_vol
    character(len=:), allocatable, intent(in) :: Geometry
    integer :: j

    print*, Geometry
    ctrl_vol(1) = 0.0
    ctrl_vol(NJ) = 0.0

    if (Geometry == 'Rectangular') then
      do j = 2, NJ-1
        ctrl_vol(j) = delX(j)
      end do

    else if (Geometry == 'Cylindrical') then
      do j = 2, NJ-1
        ctrl_vol(j) = PI*( (xx(j) + delX(j)/2.0)**2 - (xx(j) - delX(j)/2.0)**2 )
      end do

    else if (Geometry == 'Spherical') then
      do j = 2, NJ-1
        ctrl_vol(j) = 4.0/3.0*PI*( (xx(j) + delX(j)/2.0)**3 - (xx(j) - delX(j)/2.0)**3 )
      end do

    end if

  end function Control_Volume


  function Cross_Sectional_Area(Geometry) result(Cross_Area)
    use variables, only: xx, delX
    real, dimension(2, NJ) :: Cross_Area      ! 1 - west side, 2 - east side
    character(len=:), allocatable, intent(in) :: Geometry
    integer :: j

    if (Geometry == 'Rectangular') then
      Cross_Area = 1.0

    else if (Geometry == 'Cylindrical') then
      do j = 1, NJ
          Cross_Area(1,j) = 2.0 * PI * (xx(j) - delX(j)/2.0)
          Cross_Area(2,j) = 2.0 * PI * (xx(j) + delX(j)/2.0)
      end do

    else if (Geometry == 'Spherical') then
      do j = 1, NJ
        Cross_Area(1,j) = 4.0 * PI * (xx(j) - delX(j)/2.0)**2.0
        Cross_Area(2,j) = 4.0 * PI * (xx(j) + delX(j)/2.0)**2.0
      end do

    end if

  end function Cross_Sectional_Area
! ******************************************************************************

end module GOV_EQNS
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py





! ******************************************************************************










! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use echem_transport_mod
  use variables, only : cprev, delC
  use write_data_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j
  real :: convergence_error
  real :: step_size = 1.0, tol = 1e-7, Delta_V_step = 1e-3, &
  & Delta_V_upper_limit = 400e-3, Delta_V_lower_limit = -400e-3

  it = 0

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()
  call write_condition(it) ! writes the headers

  Delta_V = 0e-3

  do while (Delta_V <= Delta_V_upper_limit)

      convergence_error = 1e3
      it = 0
      do while ((convergence_error > tol).AND.(it <= maxIterations))

        ! **********************************************************************
        ! ******************************** BOUND VAL ***************************
        ! **********************************************************************
          do j = 1, NJ
            call auto_fill(j)
            call ABDGXY(j)
            call BAND(j)
          end do

          cprev = cprev + step_size*delC      ! Update the dependent variables
          convergence_error = MAXVAL(ABS(delC))
          it = it + 1
          ! ********************************************************************

      end do

      call write_condition(it)
      Delta_V = Delta_V + Delta_V_step

    end do

    call initial_condition()
    Delta_V = -Delta_V_step
    do while (Delta_V >= Delta_V_lower_limit)

        convergence_error = 1e3
        it = 0
        do while ((convergence_error > tol).AND.(it <= maxIterations))

          ! ********************************************************************
          ! ******************************** BOUND VAL *************************
          ! ********************************************************************
            do j = 1, NJ
              call auto_fill(j)
              call ABDGXY(j)
              call BAND(j)
            end do

            cprev = cprev + step_size*delC      ! Update the dependent variables
            convergence_error = MAXVAL(ABS(delC))
            it = it + 1
            ! ******************************************************************

        end do

        call write_condition(it)
        Delta_V = Delta_V - Delta_V_step
      end do

  call system_clock(t2,clock_rate,clock_max)
  write ( *, * ) 'Elapsed real time =', real(t2-t1)/real(clock_rate )

end program unsteady

! ------------------------------------------------------------------------------
! ***************************** end MAIN PROGRAM *******************************
! ------------------------------------------------------------------------------


subroutine initial_condition()
  use user_input
  use variables
  use echem_transport_mod
  use GOV_EQNS, only: Control_Volume, Cross_Sectional_Area,  &    ! functions
                      Cntrl_Vol, Crx_Area                         ! variables

  implicit none
  real    :: h
  integer :: j

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-2)
  do j=1,NJ

    if (j.EQ.1) then
      xx(j) = R_1
    else if (j.EQ.NJ) then
      xx(NJ) = R_1 + xmax
    else
      xx(j) = xx(1) + h*float(j-1) - h/2.0
    end if

  end do

  do j=2, NJ-1
     delx(j) = h
  end do
     delx(1) = 0.0                      ! it is common in the finite volume
     delx(NJ) = 0.0                     ! algorithm to set the control
                                        ! volumes at the boundaries to zero
  do j = 1, NJ
    cprev(1,j) = cbulk
    cprev(2,j) = OpenCircuitPotential(cbulk)
    cprev(3,J) = 0.0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition



!------------------------------------------------------------------------------!
!********************************* AUTO_FILL **********************************!
!------------------------------------------------------------------------------!
subroutine auto_fill(j)
  use user_input, only : N, NJ
  use GOV_EQNS
  use dnadmod
  use variables

  implicit none

  integer :: j
  integer :: ic

  ! variables and their derivatives at the control volume interfaces
  real, dimension(N) :: cW, cE, dcdxW, dcdxE

  ! dual variables
  type(dual), dimension(N) :: cW_dual, dcdxW_dual, flux_dualW
  type(dual), dimension(N) :: cE_dual, dcdxE_dual, flux_dualE
  type(dual), dimension(N) :: cj_dual, reaction_dual, accumulation_dual
  type(dual), dimension(N) :: boundary_conditionW, boundary_conditionE

  ! set all matrix variables (dW, dE, fW, fE, rj, smG) to 0.0
  dW = 0.0
  dE = 0.0
  fW = 0.0
  fE = 0.0
  rj = 0.0
  smG = 0.0

  !-----------------------------------------------------------------------------
  ! Calculate cW, cE, dcdxW, dcdxE
  !-----------------------------------------------------------------------------
  if (j /= 1) then                          !-----------------------------------
    alphaW = delx(j-1)/(delx(j-1)+delx(j))  ! West Side Interface Variables:
    betaW = 2.0/(delx(j-1)+delx(j))         ! cW, dcdxW
    do ic=1,N                               !-----------------------------------
      cW(ic)    = alphaW*cprev(ic,j) + (1.0 - alphaW)*cprev(ic,j-1)
      dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
    end do
    cW_dual = c_to_dual(cW)
    dcdxW_dual = dcdx_to_dual(dcdxW)
  end if
                                            !-----------------------------------
  if (j /= NJ) then                         ! East Side Interface Variables:
    alphaE = delx(j)/(delx(j+1)+delx(j))    ! cE, dcdxE
    betaE = 2.0/(delx(j)+delx(j+1))         !-----------------------------------
    do ic=1,N
      cE(ic)    = alphaE*cprev(ic,j+1) + (1.0 - alphaE)*cprev(ic,j)
      dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
    end do
    cE_dual = c_to_dual(cE)
    dcdxE_dual = dcdx_to_dual(dcdxE)
  end if

  !-----------------------------------------------------------------------------
  ! Boundary Conditions
  !-----------------------------------------------------------------------------
  if (j == 1) then                              ! West Side Boundary Conditions
    boundary_conditionW = Boundary_WEST(cE_dual, dcdxE_dual) !------------------
    do ic = 1,N
        fE(ic, :) = boundary_conditionW(ic)%dx(1:N)
        dE(ic, :) = boundary_conditionW(ic)%dx(N+1:2*N)

        smG(ic)   = -(-boundary_conditionW(ic)%x)
    end do
                                                !-------------------------------
  else if (j == NJ) then                        ! East Side Boundary Conditions
      boundary_conditionE = Boundary_EAST(cW_dual, dcdxW_dual) !----------------
    do ic = 1,N
        fW(ic, :) = boundary_conditionE(ic)%dx(1:N)
        dW(ic, :) = boundary_conditionE(ic)%dx(N+1:2*N)

        smG(ic)   = -(boundary_conditionE(ic)%x)
    end do

  !-----------------------------------------------------------------------------
  ! Governing Equations
  !-----------------------------------------------------------------------------
  else                                        !------------ Fluxes -------------
    flux_dualW = FLUX(cW_dual, dcdxW_dual)    ! West Side Flux
    flux_dualE = FLUX(cE_dual, dcdxE_dual)    ! East Side Flux

    cj_dual = c_to_dual(cprev(:,j))           !---------------------------------
    reaction_dual = RXN(cj_dual)              ! Control Volume - Rxn, Accum
    accumulation_dual = ACCUM(cj_dual)        !---------------------------------

    do ic = 1,N                               ! ic - equation number
        fW(ic,:) = flux_dualW(ic)%dx(1:N)     * Crx_Area(1,j)
        fE(ic,:) = flux_dualE(ic)%dx(1:N)     * Crx_Area(2,j)
        dW(ic,:) = flux_dualW(ic)%dx(N+1:2*N) * Crx_Area(1,j)
        dE(ic,:) = flux_dualE(ic)%dx(N+1:2*N) * Crx_Area(2,j)

        rj(ic,:) = ( reaction_dual(ic)%dx(1:N)                &
                 &  - accumulation_dual(ic)%dx(1:N) )*Cntrl_Vol(j)

        smG(ic)  = -(   flux_dualW(ic)%x * Crx_Area(1,j)      &
                 &    - flux_dualE(ic)%x * Crx_Area(2,j)      &
                 &    + reaction_dual(ic)%x * Cntrl_Vol(j) )
    end do
  end if

end subroutine auto_fill

!------------------------------------------------------------------------------!
!*********************************** ABDGXY ***********************************!
!------------------------------------------------------------------------------!
subroutine ABDGXY(j)
  ! ABDGXY equates the large coefficents based on the small coefficents.
  ! The coefficents A, B, D, G, X, Y can be found in Newman appendix C.
      use user_input, only: N, NJ
      use GOV_EQNS
      use variables, only: A, B, D, G, X, Y
      implicit none
      integer :: j


      if (j.eq.1) then
        X = 0.0
        B = rj - (1.0 - alphaE)*fE + betaE*dE
        D(1:N,1:N) = -alphaE*fE - betaE*dE
        G = smG

        return
      end if

      if (j.eq.NJ) then
        Y = 0.0
        A = (1.d0 - alphaW)*fW - betaW*dW
        B = rj + betaW*dW + alphaW*fW
        G = smG

        return
      end if

      A = (1.d0 - alphaW)*fW - betaW*dW
      B = rj + betaW*dW + alphaW*fW - (1.0 - alphaE)*fE + betaE*dE
      D(1:N,1:N) = -alphaE*fE - betaE*dE
      G = smG

      return
end subroutine ABDGXY

!------------------------------------------------------------------------------!
!*********************************** MATINV ***********************************!
!------------------------------------------------------------------------------!
SUBROUTINE MATINV(N, M, DETERM)
  use variables, only: A, B, delC, D, ID ! A imported but not used
  implicit double precision (A-H,O-Z) ! implicits are not good coding practice
 ! use variables, only: delC ! A imported but not used
 ! use ABDGXY_VARS, only: A, B, D
 ! implicit double precision (A-H,O-Z)
 ! real, dimension(N) :: ID

      DETERM=1.0
      ! ID = 0.0
      DO 1 I=1,N
1       ID(I)=0
      DO 18 NN=1,N
        BMAX=1.1
        DO 6 I=1,N
          IF (ID(I).NE.0) GOTO 6
          BNEXT=0.0
          BTRY=0.0
          DO 5 J=1,N
            IF (ID(J).NE.0) GOTO 5
            IF (ABS(B(I,J)).LE.BNEXT) GOTO 5
            BNEXT=ABS(B(I,J))
            IF (BNEXT.LE.BTRY) GOTO 5
            BNEXT=BTRY
            BTRY=ABS(B(I,J))
            JC=J
5           CONTINUE
          IF (BNEXT.GE.BMAX*BTRY) GOTO 6
          BMAX=BNEXT/BTRY
          IROW=I
          JCOL=JC
6       CONTINUE
        IF (ID(JC).EQ.0) GOTO 8
        DETERM=0.0
        RETURN
8       ID(JCOL)=1
        IF (JCOL.EQ.IROW) GOTO 12
9       DO 10 J=1,N
          SAVE=B(IROW,J)
          B(IROW,J)=B(JCOL,J)
10        B(JCOL,J)=SAVE
        DO 11 K=1,M
          SAVE=D(IROW,K)
          D(IROW,K)=D(JCOL,K)
11        D(JCOL,K)=SAVE
12      F=1.0/B(JCOL,JCOL)
        DO 13 J=1,N
13        B(JCOL,J)=B(JCOL,J)*F
        DO 14 K=1,M
14        D(JCOL,K)=D(JCOL,K)*F
        DO 18 I=1,N
          IF (I.EQ.JCOL) GOTO 18
          F=B(I,JCOL)
          DO 16 J=1,N
16          B(I,J)=B(I,J)-F*B(JCOL,J)
          DO 17 K=1,M
17          D(I,K)=D(I,K)-F*D(JCOL,K)
18    CONTINUE
      RETURN
      end

!------------------------------------------------------------------------------!
!************************************ BAND ************************************!
!------------------------------------------------------------------------------!
SUBROUTINE BAND(J)
! BAND(J) computes delC and calls MATINV to solve the problem using gaussian elimination.
use variables, only: A, B, delC, D, G, X, Y, NP1, E
! use variables, only: delC
use user_input, only: N, NJ
! use ABDGXY_VARS
! use BAND_J_VARS
implicit double precision (A-H,O-Z)
! integer :: NP1

101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
        D(I,2*N+1)=G(I)
        DO 2 L=1,N
          LPN=L+N
2         D(I,LPN)=X(I,L)
      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J
4     DO 5 K=1,N
        E(K,NP1,1)=D(K,2*N+1)
        DO 5 L=1,N
          E(K,L,1)=-D(K,L)
          LPN=L+N
5         X(K,L)=-D(K,LPN)
      RETURN
6     DO 7 I=1,N
        DO 7 K=1,N
          DO 7 L=1,N
7           D(I,K)=D(I,K)+A(I,L)*X(L,K)
8     IF (J-NJ) 11,9,9
9     DO 10 I=1,N
        DO 10 L=1,N
          G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
          DO 10 M=1,N
10          A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)
11    DO 12 I=1,N
        D(I,NP1)=-G(I)
        DO 12 L=1,N
          D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
          DO 12 K=1,N
12          B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)
      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N
        DO 15 M=1,NP1
15        E(K,M,J)=-D(K,M)
      IF (J-NJ) 20,16,16
16    DO 17 K=1,N
17      delC(K,J)=E(K,NP1,J)
      DO 18 JJ=2,NJ
        M=NJ-JJ+1
        DO 18 K=1,N
          delC(K,M)=E(K,NP1,M)
          DO 18 L=1,N
18          delC(K,M)=delC(K,M) +E(K,L,M)*delC(L,M+1)
      DO 19 L=1,N
        DO 19 K=1,N
19        delC(K,1)=delC(K,1)+X(K,L)*delC(L,3)
20    RETURN
      end



!******************************************************************************
!* dual Number Automatic Differentiation (DNAD) of Fortran Codes
!*-----------------------------------------------------------------------------
!* COPYRIGHT (c) Joshua Hodson, All rights reserved, you are free to copy,
!* modify, or translate this code to other languages such as c/c++. This is a
!* fork of the original Fortran DNAD module developed by Dr. Wenbin Yu. See
!* original copyright information below. You can download the original version
!* at https://cdmhub.org/resources/374
!*
!* COPYRIGHT (c) Wenbin Yu, All rights reserved, you are free to copy,
!* modify or translate this code to other languages such as c/c++. If
!* you find a bug please let me know through wenbinyu.heaven@gmail.com. If
!* you added new functions and want to share with others, please let me know
!* too. You are welcome to share your successful stories with us through
!* http://groups.google.com/group/hifi-comp.
!******************************************************************************
!* Acknowledgements
!*-----------------------------------------------------------------------------
!* The development of DNAD is supported, in part, by the Chief Scientist
!* Innovative Research Fund at AFRL/RB WPAFB, and by Department of Army
!* SBIR (Topic A08-022) through Advanced Dynamics Inc. The views and
!* conclusions contained herein are those of the authors and should not be
!* interpreted as necessarily representing the official policies or
!* endorsement, either expressed or implied, of the funding agency.
!*
!* Additional development of DNAD has been supported under a Department of
!* Energy (DOE) Nuclear Energy University Program (NEUP) Graduate Fellowship.
!* Any opinions, findings, conclusions or recommendations expressed in this
!* publication are those of the authors and do not necessarily reflect the
!* views of the Department of Energy Office of Nuclear Energy.
!******************************************************************************
!* Citation
!*-----------------------------------------------------------------------------
!* Your citation of the following two papers is appreciated:
!* Yu, W. and Blair, M.: "DNAD, a Simple Tool for Automatic Differentiation of
!* Fortran Codes Using dual Numbers," Computer Physics Communications, vol.
!* 184, 2013, pp. 1446-1452.
!*
!* Spall, R. and Yu, W.: "Imbedded dual-Number Automatic Differentiation for
!* CFD Sensitivity Analysis," Journal of Fluids Engineering, vol. 135, 2013,
!* 014501.
!******************************************************************************
!* Quick Start Guide
!*-----------------------------------------------------------------------------
!* To integrate DNAD into an existing Fortran program, do the following:
!*
!*   1. Include the DNAD module in the source files by adding "use dnadmod" to
!*      the beginning of all modules, global functions, and global subroutines
!*      that include definitions of floating-point variables.
!*   2. Redefine all floating-point variables as type(dual). This can be done
!*      using precompiler directives so that the integration can be turned on
!*      or off at compile-time, eliminating the need for maintaining two
!*      separate code bases for the same project.
!*   3. All I/O involving floating-point variables will need to be examined.
!*      A method will need to be determined for inputting and outputting
!*      derivative values. This customization is typically unique for each
!*      piece of software and needs to be determined on a case-by-case basis.
!*   4. When compiling DNAD, use the compiler option "-Dndv=#", where # is the
!*      number of design variables desired. This sizes the derivative array
!*      that is stored with each floating point number.
!*   5. When compiling DNAD, use compiler options to specify precision. If no
!*      compiler options are specified, DNAD will default to single-precision
!*      floating-point arithmetic. Most popular Fortran compilers provide
!*      options for specifying precision at compile-time so that it does not
!*      have to be hard-coded into the source code. For example, use the
!*      "-fdefault-real-8" compiler in gfortran or the "-r8" compiler option
!*      with Intel Fortran to compile DNAD as double-precision.
!*   6. Modify the compilation process for the target software to include the
!*      DNAD module in the resulting executable or library.
!******************************************************************************
!* Change Log
!*-----------------------------------------------------------------------------
!*
!*  2016-04-29  Joshua Hodson
!*  - Updated copyright, acknowledgments, and quick start guide.
!*  - Removed overloads for single-precision reals.
!*  - Added tan, dtan, atan, and atan2 intrinsic function overloads.
!*  - Removed macro for precision and defined all floating-point variables as
!*    default real. Compiler options can now be used to set precision.
!*  - Added checks for undefined derivatives when only constants are used in
!*    the calculation (i.e. all partial derivatives are zero). This limits the
!*    perpetuation of NaN values in the code.
!*  - Combined the header and source files into a single file.
!*
!*  2015-07-29  Joshua Hodson
!*  - Added maxloc intrinsic function overload.
!*  - Converted UPPERCASE to lowercase for readability.
!*  - Added macros for defining precision and number of design variables.
!*  - Renamed module from dual_Num_Auto_Diff to dnadmod
!*  - Renamed dual number type from dual_NUM to dual
!*  - Renamed components of dual number type from (xp_ad_, xp_ad_) to (x, dx)
!*
!*  2014-06-05  Wenbin Yu
!*  - Forked from original DNAD repository, see https://cdmhub.org/resources/374
!*
!******************************************************************************

! Number of design variables (default = 1)
! #ifndef ndv
! #define ndv 1
! #endif

! ******************************************************************************
! ******************************************************************************
