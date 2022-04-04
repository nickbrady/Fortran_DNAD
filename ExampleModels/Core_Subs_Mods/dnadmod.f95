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
!*  2022-03-24  Nicholas Brady
!*  - added type dual_complex (these have various physical applications)
!*  - extended the overloading of many intrinsic functions (but not all) to
!*      include type dual_complex
!*  - defined erf_c to calculate the error function of a complex number
!*
!*  2021-05-29  Nicholas Brady
!*  - overloaded hyperbolic, inverse hyperbolic, and error functions (intrinsic)
!*      sinh, cosh, tanh, asinh, acosh, atanh, erf
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

module dnadmod
    use user_input, only: N

    implicit none
    integer, PARAMETER :: ndv = N*2   ! cprev_vars and dcdx_vars

    private

    real :: negative_one = -1.0
    real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

    type,public :: dual ! make this private will create difficulty to use the
                        ! original write/read commands, hence x and dx are
                        ! variables which can be accessed using D%x and D%dx in
                        ! other units using this module in which D is defined
                        ! as type(dual).
        sequence
        real :: x        ! functional value
        real :: dx(ndv)  ! derivative
    end type dual

    type,public :: dual_complex ! make this private will create difficulty to use the
                        ! original write/read commands, hence z and dz are
                        ! variables which can be accessed using D%z and D%dz in
                        ! other units using this module in which D is defined
                        ! as type(dual_complex).
        sequence
        complex :: z        ! functional value
        complex :: dz(ndv)  ! derivative
    end type dual_complex


!******** Interfaces for operator overloading
    public assignment (=)
    interface assignment (=)
        module procedure assign_di ! dual=integer, elemental
        module procedure assign_dr ! dual=real, elemental
        module procedure assign_id ! integer=dual, elemental

        module procedure assign_dc_i    ! dual_complex = integer
        module procedure assign_dc_r    ! dual_complex = real
        module procedure assign_dc_c    ! dual_complex = complex
        module procedure assign_c_dc    ! complex = dual_complex
    end interface


    public operator (+)
    interface operator (+)
        module procedure add_d  ! +dual number, elemental
        module procedure add_dd ! dual + dual, elemental
        module procedure add_di ! dual + integer, elemental
        module procedure add_dr ! dual + real, elemental
        module procedure add_id ! integer + dual, elemental
        module procedure add_rd ! real + dual, elemental

        module procedure add_dc     ! +dual_complex number, elemental
        module procedure add_dc_dc  ! dual_complex + dual_complex
        module procedure add_c_dc   ! complex + dual_complex
        module procedure add_dc_c   ! dual_complex + complex
        module procedure add_i_dc   ! integer + dual_complex
        module procedure add_dc_i   ! dual_complex + integer
        module procedure add_r_dc   ! real + dual_complex
        module procedure add_dc_r   ! dual_complex + real
    end interface

    public operator (-)
    interface operator (-)
        module procedure minus_d  ! negate a dual number,elemental
        module procedure minus_dd ! dual -dual,elemental
        module procedure minus_di ! dual-integer,elemental
        module procedure minus_dr ! dual-real,elemental
        module procedure minus_id ! integer-dual,elemental
        module procedure minus_rd ! real-dual,elemental

        module procedure minus_dc     ! -dual_complex number, elemental
        module procedure minus_dc_dc  ! dual_complex - dual_complex
        module procedure minus_c_dc   ! complex - dual_complex
        module procedure minus_dc_c   ! dual_complex - complex
        module procedure minus_i_dc   ! integer - dual_complex
        module procedure minus_dc_i   ! dual_complex - integer
        module procedure minus_r_dc   ! real - dual_complex
        module procedure minus_dc_r   ! dual_complex - real
    end interface

    public operator (*)
    interface operator (*)
        module procedure mult_dd ! dual*dual, elemental
        module procedure mult_di ! dual*integer,elemental
        module procedure mult_dr ! dual*real,elemental
        module procedure mult_id ! integer*dual,elemental
        module procedure mult_rd ! real*dual,elemental

        module procedure mult_dc_dc ! dual_complex*dual_complex, elemental
        module procedure mult_dc_i ! dual_complex*integer,elemental
        module procedure mult_dc_r ! dual_complex*real,elemental
        module procedure mult_i_dc ! integer*dual,elemental
        module procedure mult_r_dc ! real*dual_complex,elemental
        module procedure mult_dc_c ! dual_complex*complex,elemental
        module procedure mult_c_dc ! complex*dual_complex,elemental
    end interface

    public operator (/)
    interface operator (/)
        module procedure div_dd ! dual/dual,elemental
        module procedure div_di ! dual/integer, elemental
        module procedure div_dr ! dual/real,emental
        module procedure div_id ! integer/dual, elemental
        module procedure div_rd ! real/dual, elemental

        module procedure div_dc_dc ! dual_complex/dual_complex,elemental
        module procedure div_dc_i ! dual_complex/integer, elemental
        module procedure div_dc_r ! dual_complex/real,emental
        module procedure div_i_dc ! integer/dual, elemental
        module procedure div_r_dc ! real/dual_complex, elemental
        module procedure div_dc_c ! dual_complex/complex,emental
        module procedure div_c_dc ! complex/dual_complex, elemental
    end interface

    public operator (**)
    interface operator (**)
        module procedure pow_i ! dual number to an integer power,elemental
        module procedure pow_r ! dual number to a real power, elemental
        module procedure pow_d ! dual number to a dual power, elemental

        module procedure pow_dc_i ! dual_complex number to an integer power,elemental
        module procedure pow_dc_r ! dual_complex number to a real power, elemental
        module procedure pow_dc_c ! dual_complex number to a complex power, elemental
        module procedure pow_dc_dc ! dual_complex number to a dual_complex power, elemental
    end interface

    public operator (==)
    interface operator (==)
        module procedure eq_dd ! compare two dual numbers, elemental
        module procedure eq_di ! compare a dual and an integer, elemental
        module procedure eq_dr ! compare a dual and a real, elemental
        module procedure eq_id ! compare integer with a dual number, elemental
        module procedure eq_rd ! compare a real with a dual number, elemental

        module procedure eq_dc_dc ! compare two dual_complex numbers, elemental
        module procedure eq_dc_i ! compare a dual_complex and an integer, elemental
        module procedure eq_dc_r ! compare a dual_complex and a real, elemental
        module procedure eq_i_dc ! compare integer with a dual_complex number, elemental
        module procedure eq_r_dc ! compare a real with a dual_complex number, elemental
    end interface

    public operator (<=)
    interface operator (<=)
        module procedure le_dd ! compare two dual numbers, elemental
        module procedure le_di ! compare a dual and an integer, elemental
        module procedure le_dr ! compare a dual and a real,elemental
        module procedure le_id ! compare integer with a dual number, elemental
        module procedure le_rd ! compare a real with a dual number, elemental

        ! <= cannot be done for complex numbers
    end interface

    public operator (<)
    interface operator (<)
        module procedure lt_dd ! compare two dual numbers, elemental
        module procedure lt_di ! compare a dual and an integer, elemental
        module procedure lt_dr ! compare dual with a real, elemental
        module procedure lt_id ! compare integer with a dual number, elemental
        module procedure lt_rd ! compare a real with a dual number, elemental

        ! < cannot be done for complex numbers
    end interface

    public operator (>=)
    interface operator (>=)
        module procedure ge_dd ! compare two dual numbers, elemental
        module procedure ge_di ! compare dual with integer, elemental
        module procedure ge_dr ! compare dual with a real number, elemental
        module procedure ge_id ! compare integer with a dual number, elemental
        module procedure ge_rd ! compare a real with a dual number, elemental

        ! >= cannot be done for complex numbers
    end interface

    public operator (>)
    interface operator (>)
        module procedure gt_dd ! compare two dual numbers, elemental
        module procedure gt_di ! compare a dual and an integer, elemental
        module procedure gt_dr ! compare dual with a real, elemental
        module procedure gt_id ! compare integer with a dual number, elemental
        module procedure gt_rd ! compare a real with a dual number, elemental

        ! > cannot be done for complex numbers
    end interface

    public operator (/=)
    interface operator (/=)
        module procedure ne_dd ! compare two dual numbers, elemental
        module procedure ne_di ! compare a dual and an integer, elemental
        module procedure ne_dr ! compare dual with a real, elemental
        module procedure ne_id ! compare integer with a dual number, elemental
        module procedure ne_rd ! compare a real with a dual number, elemental

        module procedure ne_dc_dc ! compare two dual_complex numbers, elemental
        module procedure ne_dc_i ! compare a dual_complex and an integer, elemental
        module procedure ne_dc_r ! compare dual_complex with a real, elemental
        module procedure ne_i_dc ! compare integer with a dual_complex number, elemental
        module procedure ne_r_dc ! compare a real with a dual_complex number, elemental
    end interface


!------------------------------------------------
! Interfaces for intrinsic functions overloading
!------------------------------------------------
    public abs
    interface abs
        module procedure abs_d  ! absolute value of a dual number, elemental
        module procedure abs_dc
    end interface

    public dabs
    interface dabs
        module procedure abs_d  ! same as abs, used for some old fortran commands
        module procedure abs_dc
    end interface

    public acos
    interface acos
        module procedure acos_d ! arccosine of a dual number, elemental
        module procedure acos_dc ! arccosine of a dual number, elemental
    end interface

    public asin
    interface asin
        module procedure asin_d ! arcsine of a dual number, elemental
        module procedure asin_dc ! arcsine of a dual number, elemental
    end interface

    public atan
    interface atan
        module procedure atan_d ! arctan of a dual number, elemental
        module procedure atan_dc ! arctan of a dual number, elemental
    end interface

    public atan2
    interface atan2
        module procedure atan2_d ! arctan of a dual number, elemental
    end interface

    public cos
    interface cos
        module procedure cos_d ! cosine of a dual number, elemental
        module procedure cos_dc ! cosine of a compulex_dual number, elemental
    end interface

    public dcos
    interface dcos
        module procedure cos_d ! cosine of a dual number, elemental
        module procedure cos_dc ! cosine of a compulex_dual number, elemental
    end interface

    public dot_product
    interface dot_product
        module procedure dot_product_dd ! dot product two dual number vectors
    end interface

    public exp
    interface exp
        module procedure exp_d ! exponential of a dual number, elemental
        module procedure exp_dc ! exponential of a dual number, elemental
    end interface

    public int
    interface int
        module procedure int_d ! integer part of a dual number, elemental
    end interface

    public log
    interface log
        module procedure log_d ! log of a dual number, elemental
        module procedure log_dc ! log of a dual number, elemental
    end interface

    public log10
    interface log10
        module procedure log10_d ! log10 of a dual number, elemental
        module procedure log10_dc ! log10 of a dual_complex number, elemental
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
        module procedure sin_dc ! obtain sine of a dual_complex number, elemental
    end interface

    public dsin
    interface dsin
        module procedure sin_d ! obtain sine of a dual number, elemental
        module procedure sin_dc ! obtain sine of a dual_complex number, elemental
    end interface

    public tan
    interface tan
        module procedure tan_d ! obtain sine of a dual number, elemental
        module procedure tan_dc ! obtain sine of a dual_complex number, elemental
    end interface

    public dtan
    interface dtan
        module procedure tan_d ! obtain sine of a dual number, elemental
        module procedure tan_dc ! obtain sine of a dual_complex number, elemental
    end interface

    public sqrt
    interface sqrt
        module procedure sqrt_d ! obtain the sqrt of a dual number, elemental
        module procedure sqrt_dc ! obtain the sqrt of a dual_complex number, elemental
    end interface

    public sum
    interface sum
        module procedure sum_d ! sum a dual array
        module procedure sum_dc ! sum a dual_complex array
    end interface

    public maxloc
    interface maxloc
        module procedure maxloc_d ! location of max in a dual array
    end interface

    public sinh
    interface sinh
        module procedure sinh_d ! obtain sinh of a dual number, elemental
        module procedure sinh_dc ! obtain sinh of a dual number, elemental
    end interface

    public cosh
    interface cosh
        module procedure cosh_d ! obtain cosh of a dual number, elemental
        module procedure cosh_dc ! obtain cosh of a dual number, elemental
    end interface

    public tanh
    interface tanh
        module procedure tanh_d ! obtain tanh of a dual number, elemental
        module procedure tanh_dc ! obtain tanh of a dual number, elemental
    end interface

    public asinh
    interface asinh
        module procedure asinh_d ! obtain asinh of a dual number, elemental
        module procedure asinh_dc ! obtain asinh of a dual number, elemental
    end interface

    public acosh
    interface acosh
        module procedure acosh_d ! obtain acosh of a dual number, elemental
        module procedure acosh_dc ! obtain acosh of a dual number, elemental
    end interface

    public atanh
    interface atanh
        module procedure atanh_d ! obtain atanh of a dual number, elemental
        module procedure atanh_dc ! obtain atanh of a dual number, elemental
    end interface

    public erf
    interface erf
        module procedure erf_d
        module procedure erf_c
        module procedure erf_dc
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

    !-----------------------------------------
    ! ERF OF dual numbers
    ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(PI)*exp(-u**2) * du>
    !-----------------------------------------
    elemental function erf_d(u) result(res)
      type(dual), intent(in) :: u
      type(dual) :: res

      res%x  = erf(u%x)
      res%dx = u%dx * 2.0/sqrt(PI) * exp(-u%x**2)

    end function erf_d







! ******************************************************************************
! COMPLEX DUAL NUMBERS
! ******************************************************************************

!*********Begin: functions/subroutines for overloading operators

!******* Begin: (=)
!---------------------

    !-----------------------------------------
    ! dual_complex = integer
    ! <u, du> = <i, 0>
    !-----------------------------------------
    elemental subroutine assign_dc_i(u, i)
         type(dual_complex), intent(out) :: u
         integer, intent(in) :: i

         u%z = complex(real(i), 0.0)  ! This is faster than direct assignment
         u%dz = complex(0.0, 0.0)

    end subroutine assign_dc_i


    !-----------------------------------------
    ! dual_complex = real(double)
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dc_r(u, r)
        type(dual_complex), intent(out) :: u
        real, intent(in) :: r

        u%z = complex(r, 0.0)
        u%dz = complex(0.0, 0.0)

    end subroutine assign_dc_r

    !-----------------------------------------
    ! dual_complex = complex
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dc_c(u, r)
        type(dual_complex), intent(out) :: u
        complex, intent(in) :: r

        u%z = r
        u%dz = complex(0.0, 0.0)

    end subroutine assign_dc_c


    !-----------------------------------------
    ! integer = dual        Is there a situation where complex --> integer makes sense?
    ! i = <u, du>
    !-----------------------------------------
    ! elemental subroutine assign_id(i, v)
    !      type(dual_complex), intent(in) :: v
    !      integer, intent(out) :: i
    !
    !      i = int(v%z)
    !
    ! end subroutine assign_id

    !-----------------------------------------
    ! complex = dual_complex        Is there a situation where complex --> integer makes sense?
    ! i = <u, du>
    !-----------------------------------------
    elemental subroutine assign_c_dc(u, v)
         type(dual_complex), intent(in) :: v
         complex, intent(out) :: u

         u = v%z

    end subroutine assign_c_dc


!******* end: (=)
!---------------------

!******* Begin: (+)
!---------------------

    !-----------------------------------------
    ! Unary positive
    ! <res, dres> = +<u, du>
    !-----------------------------------------
    elemental function add_dc(u) result(res)
         type(dual_complex), intent(in) :: u
         type(dual_complex) :: res

         res = u  ! Faster than assigning component wise

    end function add_dc

    !-----------------------------------------
    ! dual_complex + dual_complex
    ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
    !-----------------------------------------
    elemental function add_dc_dc(u, v) result(res)
       type(dual_complex), intent(in) :: u, v
       type(dual_complex) :: res

       res%z = u%z + v%z
       res%dz = u%dz + v%dz

    end function add_dc_dc

    !-----------------------------------------
    ! complex + dual_complex
    ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
    !-----------------------------------------
    elemental function add_c_dc(r, v) result(res)
        complex, intent(in) :: r
        type(dual_complex), intent(in) :: v
        type(dual_complex) :: res

        res%z = r + v%z
        res%dz = v%dz

    end function add_c_dc

    !-----------------------------------------
    ! dual_complex + complex
    ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
    !-----------------------------------------
    elemental function add_dc_c(u, r) result(res)
      type(dual_complex), intent(in) :: u
      complex, intent(in) :: r
      type(dual_complex) :: res

      res%z = r + u%z
      res%dz = u%dz

    end function add_dc_c


  !-----------------------------------------
  ! dual_complex + integer
  ! <res, dres> = <u, du> + i = <u + i, du>
  !-----------------------------------------
  elemental function add_dc_i(u, i) result(res)
       type(dual_complex), intent(in) :: u
       integer, intent(in) :: i
       type(dual_complex) :: res

       res%z = real(i) + u%z
       res%dz = u%dz

  end function add_dc_i

  !-----------------------------------------
  ! integer + dual_complex
  ! <res, dres> = <u, du> + i = <u + i, du>
  !-----------------------------------------
  elemental function add_i_dc(i, u) result(res)
       type(dual_complex), intent(in) :: u
       integer, intent(in) :: i
       type(dual_complex) :: res

       res%z = real(i) + u%z
       res%dz = u%dz

  end function add_i_dc

  !-----------------------------------------
  ! dual_complex + integer
  ! <res, dres> = <u, du> + i = <u + i, du>
  !-----------------------------------------
  elemental function add_dc_r(u, r) result(res)
       type(dual_complex), intent(in) :: u
       real, intent(in) :: r
       type(dual_complex) :: res

       res%z = r + u%z
       res%dz = u%dz

  end function add_dc_r

  !-----------------------------------------
  ! integer + dual_complex
  ! <res, dres> = <u, du> + i = <u + i, du>
  !-----------------------------------------
  elemental function add_r_dc(r, u) result(res)
       type(dual_complex), intent(in) :: u
       real, intent(in) :: r
       type(dual_complex) :: res

       res%z = r + u%z
       res%dz = u%dz

  end function add_r_dc
!
! !******* end: (+)
! !---------------------
!
!
!******* Begin: (-)
!---------------------

  !-------------------------------------------------
  ! negate a dual_complex
  ! <res, dres> = -<u, du>
  !-------------------------------------------------
  elemental function minus_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = -u%z
      res%dz = -u%dz

  end function minus_dc


  !-------------------------------------------------
  ! dual_complex - dual_complex
  ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
  !-------------------------------------------------
  elemental function minus_dc_dc(u, v) result(res)
      type(dual_complex), intent(in) :: u, v
      type(dual_complex) :: res

      res%z = u%z - v%z
      res%dz = u%dz - v%dz

  end function minus_dc_dc

  !-------------------------------------------------
  ! dual_complex - integer
  ! <res, dres> = <u, du> - i = <u - i, du>
  !-------------------------------------------------
  elemental function minus_dc_i(u, i) result(res)
      type(dual_complex), intent(in) :: u
      integer, intent(in) :: i
      type(dual_complex) :: res

      res%z = u%z - real(i)
      res%dz = u%dz

  end function minus_dc_i

  !-------------------------------------------------
  ! dual_complex - real
  ! <res, dres> = <u, du> - r = <u - r, du>
  !-------------------------------------------------
  elemental function minus_dc_r(u, r) result(res)
      type(dual_complex), intent(in) :: u
      real, intent(in) :: r
      type(dual_complex) :: res

      res%z = u%z - r
      res%dz = u%dz

  end function minus_dc_r


  !-------------------------------------------------
  ! dual_complex - complex
  ! <res, dres> = <u, du> - r = <u - r, du>
  !-------------------------------------------------
  elemental function minus_dc_c(u, c) result(res)
      type(dual_complex), intent(in) :: u
      complex, intent(in) :: c
      type(dual_complex) :: res

      res%z = u%z - c
      res%dz = u%dz

  end function minus_dc_c


  !-------------------------------------------------
  ! integer - dual_complex
  ! <res, dres> = i - <v, dv> = <i - v, -dv>
  !-------------------------------------------------
  elemental function minus_i_dc(i, v) result(res)
      integer, intent(in) :: i
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      res%z = real(i) - v%z
      res%dz = -v%dz

  end function minus_i_dc

  !-------------------------------------------------
  ! real - dual_complex
  ! <res, dres> = r - <v, dv> = <r - v, -dv>
  !-------------------------------------------------
  elemental function minus_r_dc(r, v) result(res)
      real, intent(in) :: r
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      res%z = r - v%z
      res%dz = -v%dz

  end function minus_r_dc


  !-------------------------------------------------
  ! complex - dual_complex
  ! <res, dres> = c - <v, dv> = <c - v, -dv>
  !-------------------------------------------------
  elemental function minus_c_dc(c, v) result(res)
       complex, intent(in) :: c
       type(dual_complex), intent(in) :: v
       type(dual_complex) :: res

      res%z = c - v%z
      res%dz = -v%dz

  end function minus_c_dc
!
! !******* end: (-)
! !---------------------
!
!
!******* BEGIN: (*)
!---------------------

  !----------------------------------------
  ! dual_complex * dual_complex
  ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
  !----------------------------------------
  elemental function mult_dc_dc(u, v) result(res)
      type(dual_complex), intent(in) :: u, v
      type(dual_complex) :: res

      res%z = u%z * v%z
      res%dz = u%z * v%dz + v%z * u%dz

  end function mult_dc_dc


  !-----------------------------------------
  ! dual_complex * integer
  ! <res, dres> = <u, du> * i = <u * i, du * i>
  !-----------------------------------------
  elemental function mult_dc_i(u, i) result(res)
      type(dual_complex), intent(in) :: u
      integer, intent(in) :: i
      type(dual_complex) :: res

      real :: r

      r = real(i)
      res%z = r * u%z
      res%dz = r * u%dz

  end function mult_dc_i

  !-----------------------------------------
  ! dual_complex * real
  ! <res, dres> = <u, du> * r = <u * r, du * r>
  !----------------------------------------
  elemental function mult_dc_r(u, r) result(res)
      type(dual_complex), intent(in) :: u
      real, intent(in) :: r
      type(dual_complex) :: res

      res%z = u%z * r
      res%dz = u%dz * r

  end function mult_dc_r


  !-----------------------------------------
  ! integer * dual_complex
  ! <res, dres> = i * <v, dv> = <i * v, i * dv>
  !-----------------------------------------
  elemental function mult_i_dc(i, v) result(res)
      integer, intent(in) :: i
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      real :: r

      r = real(i)
      res%z = r * v%z
      res%dz = r * v%dz

  end function mult_i_dc


  !-----------------------------------------
  ! double * dual_complex
  ! <res, dres> = r * <v, dv> = <r * v, r * dv>
  !-----------------------------------------
  elemental function mult_r_dc(r, v) result(res)
      real, intent(in) :: r
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      res%z = r * v%z
      res%dz = r * v%dz

  end function mult_r_dc


  !-----------------------------------------
  ! complex * dual_complex
  ! <res, dres> = c * <v, dv> = <c * v, c * dv>
  !-----------------------------------------
  elemental function mult_c_dc(c, v) result(res)
      complex, intent(in) :: c
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      res%z = c * v%z
      res%dz = c * v%dz

  end function mult_c_dc


  !-----------------------------------------
  ! dual_complex * complex
  ! <res, dres> = c * <v, dv> = <c * v, c * dv>
  !-----------------------------------------
  elemental function mult_dc_c(v, c) result(res)
      complex, intent(in) :: c
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      res%z = c * v%z
      res%dz = c * v%dz

  end function mult_dc_c
!
! !******* end: (*)
! !---------------------
!
!
! !******* BEGIN: (/)
! !---------------------
!
  !-----------------------------------------
  ! dual_complex / dual_complex
  ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
  !-----------------------------------------
  elemental function div_dc_dc(u, v) result(res)
      type(dual_complex), intent(in) :: u, v
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / v%z
      res%z = u%z * inv
      res%dz = (u%dz - res%z * v%dz) * inv

  end function div_dc_dc
!
!
  !-----------------------------------------
  ! dual_complex / integer
  ! <res, dres> = <u, du> / i = <u / i, du / i>
  !-----------------------------------------
  elemental function div_dc_i(u, i) result(res)
      type(dual_complex), intent(in) :: u
      integer, intent(in) :: i
      type(dual_complex) :: res

      real :: inv

      inv = 1.0 / real(i)
      res%z = u%z * inv
      res%dz = u%dz * inv

  end function div_dc_i
!
!
  !-----------------------------------------
  ! dual_complex / double
  ! <res, dres> = <u, du> / r = <u / r, du / r>
  !----------------------------------------
  elemental function div_dc_r(u, r) result(res)
      type(dual_complex), intent(in) :: u
      real, intent(in) :: r
      type(dual_complex):: res

      real :: inv

      inv = 1.0 / r
      res%z = u%z * inv
      res%dz = u%dz * inv

  end function div_dc_r
!
!
  !-----------------------------------------
  ! integer / dual_complex
  ! <res, dres> = i / <v, dv> = <i / v, -i / v^2 * du>
  !-----------------------------------------
  elemental function div_i_dc(i, v) result(res)
      integer, intent(in) :: i
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / v%z
      res%z = real(i) * inv
      res%dz = -res%z * inv * v%dz

  end function div_i_dc


  !-----------------------------------------
  ! double / dual_complex
  ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
  !-----------------------------------------
  elemental function div_r_dc(r, v) result(res)
      real, intent(in) :: r
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / v%z
      res%z = r * inv
      res%dz = -res%z * inv * v%dz

  end function div_r_dc


  !-----------------------------------------
  ! complex / dual_complex
  ! <res, dres> = c / <u, du> = <c / u, -r / u^2 * du>
  !-----------------------------------------
  elemental function div_c_dc(c, v) result(res)
      complex, intent(in) :: c
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / v%z
      res%z = c * inv
      res%dz = -res%z * inv * v%dz

  end function div_c_dc


  !-----------------------------------------
  ! dual_complex / complex
  ! <res, dres> = c / <u, du> = <c / u, -c / u^2 * du>
  !-----------------------------------------
  elemental function div_dc_c(v, c) result(res)
      complex, intent(in) :: c
      type(dual_complex), intent(in) :: v
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / c
      res%z = v%z * inv
      res%dz = v%dz * inv

  end function div_dc_c

!******* end: (/)
!---------------------

!******* BEGIN: (**)
!---------------------

  !-----------------------------------------
  ! power(dual_complex, integer)
  ! <res, dres> = <u, du> ^ i = <u ^ i, i * u ^ (i - 1) * du>
  !-----------------------------------------
  elemental function pow_dc_i(u, i) result(res)
      type(dual_complex), intent(in) :: u
      integer, intent(in) :: i
      type(dual_complex) :: res

      real :: pow_x

      pow_x = u%z ** (i - 1)
      res%z = u%z * pow_x
      res%dz = real(i) * pow_x * u%dz

  end function pow_dc_i

  !-----------------------------------------
  ! power(dual_complex, double)
  ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
  !-----------------------------------------
  elemental function pow_dc_r(u, r) result(res)
      type(dual_complex), intent(in) :: u
      real, intent(in) :: r
      type(dual_complex) :: res

      real :: pow_x

      pow_x = u%z ** (r - 1.0)
      res%z = u%z * pow_x
      res%dz = r * pow_x * u%dz

  end function pow_dc_r

  !-----------------------------------------
  ! power(dual_complex, complex)
  ! <res, dres> = <u, du> ^ c = <u ^ c, c * u ^ (c - 1) * du>
  !-----------------------------------------
  elemental function pow_dc_c(u, c) result(res)
      type(dual_complex), intent(in) :: u
      complex, intent(in) :: c
      type(dual_complex) :: res

      complex :: pow_x

      pow_x = u%z ** (c - 1.0)
      res%z = u%z * pow_x
      res%dz = c * pow_x * u%dz

  end function pow_dc_c

  !-----------------------------------------
  ! POWER dual_complex numbers to a dual_complex power
  ! <res, dres> = <u, du> ^ <v, dv>
  !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
  !-----------------------------------------
  elemental function pow_dc_dc(u, v) result(res)
      type(dual_complex), intent(in)::u, v
      type(dual_complex) :: res

      res%z = u%z ** v%z
      res%dz = res%z * (v%z / u%z * u%dz + log(u%z) * v%dz)

  end function pow_dc_dc

!******* end: (**)
!---------------------


!******* BEGIN: (==)
!---------------------
  !-----------------------------------------
  ! compare two dual_complex numbers,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_dc_dc(lhs, rhs) result(res)
       type(dual_complex), intent(in) :: lhs, rhs
       logical :: res

       res = (lhs%z == rhs%z)

  end function eq_dc_dc


  !-----------------------------------------
  ! compare a dual_complex with an integer,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_dc_i(lhs, rhs) result(res)
       type(dual_complex), intent(in) :: lhs
       integer, intent(in) :: rhs
       logical :: res

       res = (lhs%z == complex(real(rhs), 0.0))

  end function eq_dc_i


  !-----------------------------------------
  ! compare a dual_complex number with a real number,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_dc_r(lhs, rhs) result(res)
      type(dual_complex), intent(in) :: lhs
      real, intent(in) :: rhs
      logical::res

      res = (lhs%z == complex(rhs, 0.0))

  end function eq_dc_r


  !-----------------------------------------
  ! compare an integer with a dual_complex,
  ! simply compare the functional value.
  !----------------------------------------
  elemental function eq_i_dc(lhs, rhs) result(res)
       integer, intent(in) :: lhs
       type(dual_complex), intent(in) :: rhs
       logical :: res

       res = (complex(real(lhs), 0.0) == rhs%z)

  end function eq_i_dc


  !-----------------------------------------
  ! compare a real with a dual_complex,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_r_dc(lhs, rhs) result(res)
       real, intent(in) :: lhs
       type(dual_complex), intent(in) :: rhs
       logical :: res

       res = (complex(lhs, 0.0) == rhs%z)

  end function eq_r_dc


  !-----------------------------------------
  ! compare a complex number with a dual_complex,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_c_dc(lhs, rhs) result(res)
       real, intent(in) :: lhs
       type(dual_complex), intent(in) :: rhs
       logical :: res

       res = (lhs == rhs%z)

  end function eq_c_dc


  !-----------------------------------------
  ! compare a dual_complex with a comlex,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function eq_dc_c(lhs, rhs) result(res)
       type(dual_complex), intent(in) :: lhs
       complex, intent(in) :: rhs
       logical :: res

       res = (lhs%z == rhs)

  end function eq_dc_c

!******* end: (==)
!---------------------


!******* can't compare complex numbers
! (<=)
! (<)
! (>=)
! (>)
!
!******* BEGIN: (/=)
!---------------------
  !-----------------------------------------
  ! compare two dual_complex numbers, simply compare
  ! the functional value.
  !-----------------------------------------
  elemental function ne_dc_dc(lhs, rhs) result(res)
      type(dual_complex), intent(in) :: lhs, rhs
      logical :: res

      res = (lhs%z /= rhs%z)

  end function ne_dc_dc

  !-----------------------------------------
  ! compare a dual_complex number with a complex number,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function ne_dc_c(lhs, rhs) result(res)
      type(dual_complex), intent(in) :: lhs
      complex, intent(in) :: rhs
      logical :: res

      res = (lhs%z /= rhs)

  end function ne_dc_c

  !-----------------------------------------
  ! compare a complex with a dual_complex
  !-----------------------------------------
  elemental function ne_c_dc(lhs, rhs) result(res)
      complex, intent(in) :: lhs
      type(dual_complex), intent(in) :: rhs
      logical :: res

      res = (lhs /= rhs%z)

  end function ne_c_dc


  !-----------------------------------------
  ! compare a dual_complex with an integer,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function ne_dc_i(lhs, rhs) result(res)
       type(dual_complex), intent(in) :: lhs
       integer, intent(in) :: rhs
       logical :: res

       res = (lhs%z /= complex(real(rhs), 0.0))

  end function ne_dc_i


  !-----------------------------------------
  ! compare a dual_complex number with a real number,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function ne_dc_r(lhs, rhs) result(res)
      type(dual_complex), intent(in) :: lhs
      real, intent(in) :: rhs
      logical::res

      res = (lhs%z /= complex(rhs, 0.0))

  end function ne_dc_r


  !-----------------------------------------
  ! compare an integer with a dual_complex,
  ! simply compare the functional value.
  !----------------------------------------
  elemental function ne_i_dc(lhs, rhs) result(res)
       integer, intent(in) :: lhs
       type(dual_complex), intent(in) :: rhs
       logical :: res

       res = (complex(real(lhs), 0.0) /= rhs%z)

  end function ne_i_dc


  !-----------------------------------------
  ! compare a real with a dual_complex,
  ! simply compare the functional value.
  !-----------------------------------------
  elemental function ne_r_dc(lhs, rhs) result(res)
       real, intent(in) :: lhs
       type(dual_complex), intent(in) :: rhs
       logical :: res

       res = (complex(lhs, 0.0) /= rhs%z)

  end function ne_r_dc

!******* end: (/=)
!---------------------
!
  !---------------------------------------------------
  ! Absolute value of dual_complex numbers
  ! <res, dres> = abs(<u, du>) = <abs(u), du * u / |u|>
  !---------------------------------------------------
  elemental function abs_dc(u) result(res)
       type(dual_complex), intent(in) :: u
       type(dual_complex) :: res

       res%z = complex(abs(u%z), 0.)
       res%dz = u%dz * u%z / abs(u%z)

  end function abs_dc
!
!
  !-----------------------------------------
  ! ACOS of dual_complex numbers
  ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
  !----------------------------------------
  elemental function acos_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = acos(u%z)
      if (u%z == complex(1.0, 0.0) .or. u%z == complex(-1.0, 0.0)) then
          res%dz = set_Nan()  ! Undefined derivative
      else
          res%dz = -u%dz / sqrt(1.0 - u%z**2)
      end if

  end function acos_dc


  !-----------------------------------------
  ! ASIN of dual_complex numbers
  ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
  !----------------------------------------
  elemental function asin_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = asin(u%z)
      if (u%z == complex(1.0, 0.0) .or. u%z == complex(-1.0, 0.0)) then
          res%dz = set_NaN()  ! Undefined derivative
      else
          res%dz = u%dz / sqrt(1.0 - u%z**2)
      end if

  end function asin_dc


  !-----------------------------------------
  ! ATAN of dual_complex numbers
  ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
  !----------------------------------------
  elemental function atan_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = atan(u%z)
      res%dz = u%dz / (1.0 + u%z**2)

  end function atan_dc


  !-----------------------------------------
  ! ATAN2 of dual_complex numbers
  ! <res, dres> = atan2(<u, du>, <v, dv>)
  !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
  !----------------------------------------
  ! elemental function atan2_dc(u, v) result(res)
  !     type(dual_complex), intent(in) :: u, v
  !     type(dual_complex) :: res
  !
  !     complex :: usq_plus_vsq
  !
  !     res%z = atan2(u%z, v%z)
  !
  !     usq_plus_vsq = u%z**2 + v%z**2
  !     res%dz = v%z / usq_plus_vsq * u%dz - u%z / usq_plus_vsq * v%dz
  !
  ! end function atan2_dc
!
!
  !-----------------------------------------
  ! COS of dual_complex numbers
  ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
  !----------------------------------------
  elemental function cos_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = cos(u%z)
      res%dz = -sin(u%z) * u%dz

  end function cos_dc
!
!
!   !-----------------------------------------
!   ! DOT PRODUCT two dual_complex number vectors
!   ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
!   !-----------------------------------------
!   function dot_product_dc_d(u, v) result(res)
!       type(dual_complex), intent(in) :: u(:), v(:)
!       type(dual_complex) :: res
!
!       integer :: i
!
!       res%z = dot_product(u%z, v%z)
!       do i = 1, ndv
!           res%dz(i) = dot_product(u%z, v%dz(i)) + dot_product(v%z, u%dz(i))
!       end do
!
!   end function dot_product_dc_d
!
!
  !-----------------------------------------
  ! EXPONENTIAL OF dual_complex numbers
  ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
  !-----------------------------------------
  elemental function exp_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      complex :: exp_x

      exp_x = exp(u%z)
      res%z = exp_x
      res%dz = u%dz * exp_x

  end function exp_dc
!
!
!   !-----------------------------------------
!   ! Convert dual_complex to integer
!   ! i = int(<u, du>) = int(u)
!   !----------------------------------------
!   elemental function int_dc(u) result(res)
!        type(dual_complex), intent(in) :: u
!        integer :: res
!
!        res = int(u%z)
!
!   end function int_dc
!
!
  !-----------------------------------------
  ! LOG OF dual_complex numbers,defined for u%z>0 only
  ! the error control should be done in the original code
  ! in other words, if u%z<=0, it is not possible to obtain LOG.
  ! <res, dres> = log(<u, du>) = <log(u), du / u>
  !----------------------------------------
  elemental function log_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / u%z
      res%z = log(u%z)
      res%dz = u%dz * inv

  end function log_dc


  !-----------------------------------------
  ! LOG10 OF dual_complex numbers,defined for u%z>0 only
  ! the error control should be done in the original code
  ! in other words, if u%z<=0, it is not possible to obtain LOG.
  ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
  ! LOG<u,up>=<LOG(u),up/u>
  !----------------------------------------
  elemental function log10_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      complex :: inv

      inv = 1.0 / (u%z * log(10.0))
      res%z = log(u%z) / log(10.0)
      res%dz = u%dz * inv

  end function log10_dc
!
!
!   !-----------------------------------------
!   ! MULTIPLY two dual_complex number matrices
!   ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
!   !----------------------------------------
!   function matmul_dc_d(u,v) result(res)
!       type(dual_complex), intent(in) :: u(:,:), v(:,:)
!       type(dual_complex) :: res(size(u,1), size(v,2))
!
!       integer :: i
!
!       res%z = matmul(u%z, v%z)
!       do i = 1, ndv
!           res%dz(i) = matmul(u%dz(i), v%z) + matmul(u%z, v%dz(i))
!       end do
!
!   end function matmul_dc_d
!
!
!   !-----------------------------------------
!   ! MULTIPLY a dual_complex number matrix with a dual_complex number
!   ! vector
!   !
!   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
!   !----------------------------------------
!   function matmul_dc_v(u, v) result(res)
!       type(dual_complex), intent(in) :: u(:,:), v(:)
!       type(dual_complex) :: res(size(u,1))
!       integer :: i
!
!       res%z = matmul(u%z, v%z)
!       do i = 1, ndv
!           res%dz(i) = matmul(u%dz(i), v%z) + matmul(u%z, v%dz(i))
!       end do
!
!   end function matmul_dc_v
!
!
!   !-----------------------------------------
!   ! MULTIPLY a dual_complex vector with a  dual_complex matrix
!   !
!   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
!   !----------------------------------------
!   function matmul_vd(u, v) result(res)
!       type(dual_complex), intent(in) :: u(:), v(:,:)
!       type(dual_complex) :: res(size(v, 2))
!       integer::i
!
!       res%z = matmul(u%z, v%z)
!       do i = 1, ndv
!           res%dz(i) = matmul(u%dz(i), v%z) + matmul(u%z, v%dz(i))
!       end do
!
!   end function matmul_vd
!
!   !-----------------------------------------
!   ! Obtain the max of 2 to 5 dual_complex numbers
!   !----------------------------------------
!   elemental function max_dc_d(val1, val2, val3, val4,val5) result(res)
!       type(dual_complex), intent(in) :: val1, val2
!       type(dual_complex), intent(in), optional :: val3, val4,val5
!       type(dual_complex) :: res
!
!       if (val1%z > val2%z) then
!           res = val1
!       else
!           res = val2
!       endif
!       if(present(val3))then
!          if(res%z < val3%z) res = val3
!       endif
!       if(present(val4))then
!          if(res%z < val4%z) res = val4
!       endif
!       if(present(val5))then
!          if(res%z < val5%z) res = val5
!       endif
!
!   end function max_dc_d
!
!
!   !-----------------------------------------
!   ! Obtain the max of a dual_complex number and an integer
!   !----------------------------------------
!   elemental function max_dc_i(u, i) result(res)
!       type(dual_complex), intent(in) :: u
!       integer, intent(in) :: i
!       type(dual_complex) :: res
!
!       if (u%z > i) then
!           res = u
!       else
!           res = i
!       endif
!
!   end function max_dc_i
!
!   !-----------------------------------------
!   ! Obtain the max of a dual_complex number and a real number
!   !----------------------------------------
!   elemental function max_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real, intent(in) :: r
!       type(dual_complex) :: res
!
!       if (u%z > r) then
!           res = u
!       else
!           res = r
!       endif
!
!   end function max_dc_r
!
!
!   !---------------------------------------------------
!   ! Obtain the max of a real and a dual_complex
!   !---------------------------------------------------
!    elemental function max_rd(n, u) result(res)
!       real, intent(in) :: n
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res
!
!       if (u%z > n) then
!           res = u
!       else
!           res = n
!       endif
!
!   end function max_rd
!
!
!   !-----------------------------------------
!   ! Obtain the max value of vector u
!   !----------------------------------------
!   function maxval_dc(u) result(res)
!       type(dual_complex), intent(in) :: u(:)
!       integer :: iloc(1)
!       type(dual_complex) :: res
!
!       iloc=maxloc(u%z)
!       res=u(iloc(1))
!
!   end function maxval_dc
!
!
!   !-----------------------------------------
!   ! Obtain the min of 2 to 4 dual_complex numbers
!   !----------------------------------------
!   elemental function min_dc_d(val1, val2, val3, val4) result(res)
!       type(dual_complex), intent(in) :: val1, val2
!       type(dual_complex), intent(in), optional :: val3, val4
!       type(dual_complex) :: res
!
!       if (val1%z < val2%z) then
!           res = val1
!       else
!           res = val2
!       endif
!       if(present(val3))then
!          if(res%z > val3%z) res = val3
!       endif
!       if(present(val4))then
!          if(res%z > val4%z) res = val4
!       endif
!
!   end function min_dc_d
!
!
!   !-----------------------------------------
!   ! Obtain the min of a dual_complex and a double
!   !----------------------------------------
!   elemental function min_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real, intent(in) :: r
!       type(dual_complex) :: res
!
!       if (u%z < r) then
!           res = u
!       else
!           res = r
!       endif
!
!   end function min_dc_r
!
!
! !-----------------------------------------
!   ! Obtain the min value of vector u
!   !----------------------------------------
!   function minval_dc(u) result(res)
!       type(dual_complex), intent(in) :: u(:)
!       integer :: iloc(1)
!       type(dual_complex) :: res
!
!       iloc=minloc(u%z)
!       res=u(iloc(1))
!
!   end function minval_dc
!
!
!   !------------------------------------------------------
!   !Returns the nearest integer to u%z, ELEMENTAL
!   !------------------------------------------------------
!   elemental function nint_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       integer :: res
!
!       res=nint(u%z)
!
!   end function nint_dc
!
!
!   !----------------------------------------------------------------
!   ! SIGN(a,b) with two dual_complex numbers as inputs,
!   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
!   !----------------------------------------------------------------
!   elemental function sign_dc_d(val1, val2) result(res)
!       type(dual_complex), intent(in) :: val1, val2
!       type(dual_complex) :: res
!
!       if (val2%z < 0.0) then
!           res = -abs(val1)
!       else
!           res =  abs(val1)
!       endif
!
!    end function sign_dc_d
!
!
!   !----------------------------------------------------------------
!   ! SIGN(a,b) with one real and one dual_complex number as inputs,
!   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
!   !----------------------------------------------------------------
!   elemental function sign_rd(val1, val2) result(res)
!       real, intent(in) :: val1
!       type(dual_complex), intent(in) :: val2
!       type(dual_complex) :: res
!
!       if (val2%z < 0.0) then
!           res = -abs(val1)
!       else
!           res = abs(val1)
!       endif
!
!    end function sign_rd
!
!
!   !-----------------------------------------
!   ! SIN of dual_complex numbers
!   ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
!   !----------------------------------------
  elemental function sin_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = sin(u%z)
      res%dz = cos(u%z) * u%dz

  end function sin_dc
!
!
  !-----------------------------------------
  ! TAN of dual_complex numbers
  ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
  !----------------------------------------
  elemental function tan_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = tan(u%z)
      res%dz = u%dz / cos(u%z)**2

  end function tan_dc
!
!
  !-----------------------------------------
  ! SQRT of dual_complex numbers
  ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
  !----------------------------------------
  elemental function sqrt_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res
      integer :: i

      res%z = sqrt(u%z)

      if (res%z /= complex(0.0, 0.0)) then
          res%dz = 0.5 * u%dz / res%z
      else
          do i = 1, ndv
              if (u%dz(i) == complex(0.0, 0.0)) then
                  res%dz(i) = complex(0.0, 0.0)
              else
                  res%dz(i) = set_NaN()
              end if
          end do
      end if

  end function sqrt_dc


  !-----------------------------------------
  ! Sum of a dual_complex array
  !-----------------------------------------
  function sum_dc(u) result(res)
      type(dual_complex), intent(in) :: u(:)
      type(dual_complex) :: res
      integer :: i

      res%z = sum(u%z)
      do i = 1, ndv
          res%dz(i) = sum(u%dz(i))
      end do

  end function sum_dc
!
!
!   !-----------------------------------------
!   ! Find the location of the max value in an
!   ! array of dual_complex numbers
!   !-----------------------------------------
!   function maxloc_dc(array) result(ind)
!       type(dual_complex), intent(in) :: array(:)
!       integer :: ind(1)
!
!       ind = maxloc(array%z)
!
!   end function maxloc_dc
!
!
!   elemental function set_NaN() result(res)
!       real :: res
!
!       res = sqrt(negative_one)
!
!   end function set_NaN
!
!
  !-----------------------------------------
  ! Hyperbolic functions: sinh, cosh, tanh
  ! and their inverses: asinh, acosh, atanh
  !-----------------------------------------
  !-----------------------------------------
  ! SINH OF dual_complex numbers
  ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
  !-----------------------------------------
  elemental function sinh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = sinh(u%z)
      res%dz = u%dz * cosh(u%z)

  end function sinh_dc

  !-----------------------------------------
  ! COSH OF dual_complex numbers
  ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
  !-----------------------------------------
  elemental function cosh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = cosh(u%z)
      res%dz = u%dz * sinh(u%z)

  end function cosh_dc

  !-----------------------------------------
  ! TANH OF dual_complex numbers
  ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0/cosh(u)**2 * du>
  !-----------------------------------------
  elemental function tanh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = tanh(u%z)
      res%dz = u%dz * 1.0/cosh(u%z)**2

  end function tanh_dc

  !-----------------------------------------
  ! ASINH OF dual_complex numbers
  ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
  !-----------------------------------------
  elemental function asinh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = asinh(u%z)
      res%dz = u%dz * 1.0/sqrt(u%z**2 + 1.0)

  end function asinh_dc

  !-----------------------------------------
  ! ACOSH OF dual_complex numbers
  ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
  !-----------------------------------------
  elemental function acosh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = acosh(u%z)
      if (u%z == complex(1.0, 0.0)) then
          res%dz = set_Nan()  ! Undefined derivative (∞)
      else
          res%dz = u%dz * 1.0/sqrt(u%z**2 - 1.0)
      end if


  end function acosh_dc

  !-----------------------------------------
  ! ATAHN OF dual_complex numbers
  ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
  !-----------------------------------------
  elemental function atanh_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z = atanh(u%z)
      if (u%z == complex(1.0, 0.0) .or. u%z == complex(-1.0, 0.0)) then
          res%dz = set_Nan()  ! Undefined derivative
      else
          res%dz = u%dz * 1.0/(1.0 - u%z**2)
      end if

  end function atanh_dc


  ! -----------------------------------------
  ! ERF OF complex number (see Abramowitz & Stegun: Handbook of Mathematical Functions)
  ! https://math.stackexchange.com/questions/712434/erfaib-error-function-separate-into-real-and-imaginary-part
  ! https://personal.math.ubc.ca/~cbm/aands/
  !
  ! erf(x + iy) = erf(x) + e⁻ˣ² / (2πx) [(1 - cos(2xy)) + i sin(2xy)] +
  !                 2/π e⁻ˣ² ∑ₐ₌₁^(∞) [e^(-a²/4)/(a² + 4x²) ( fₐ(x,y) + i gₐ(x,y) ) ] + ϵ(x,y)
  ! fₐ(x,y) = 2x(1 - cos(2xy) cosh(ay)) + a sin(2xy) sinh(ay)
  ! gₐ(x,y) = 2x sin(2xy) cosh(ay) + a cos(2xy) sinh(ay)

  ! for improved numerical stability and accuracy, push e⁻ˣ² * e^(-a²/4) inside of fₐ(x,y) and gₐ(x,y)
  ! i.e. cosh(ay) * e⁻ˣ² * e^(-a²/4) = 0.5*(exp(ay - x² - a²/4) + exp(-ay - x² - a²/4))
  ! i.e. sinh(ay) * e⁻ˣ² * e^(-a²/4) = 0.5*(exp(ay - x² - a²/4) - exp(-ay - x² - a²/4))
  ! this is important because sinh(ay) --> ∞ as a --> ∞
  ! but e^(-a²/4) --> 0 as a --> ∞
  ! -----------------------------------------
    elemental function erf_c(c) result(res)
        complex, intent(in) :: c
        complex :: res

        real :: x, y
        complex, parameter :: i = complex(0, 1)
        integer :: k

        x = real(c)
        y = aimag(c)

        ! res = erf(x) + exp(-x**2) / (2.*PI*x)*((1. - cos(2*x*y)) + i*sin(2*x*y))
        res = erf(x) + 1. / (2.*PI*x)*((1. - cos(2*x*y)) + i*sin(2*x*y))

        do k = 1, 100   ! inefficiently defined
            ! res = res + (2./PI * exp(-x**2)) * (exp(-k**2/4.0)/(k**2 + 4*x**2)) * (f_k(k,x,y) + i*g_k(k,x,y))
            res = res + 2./PI * (1./(k**2 + 4*x**2)) * (f_k(k,x,y) + i*g_k(k,x,y))
        end do


        contains
            elemental function f_k(k, x, y) result(res)
                implicit none

                integer, intent(in) :: k
                real, intent(in) :: x, y
                real :: res
                ! 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) = cosh(k*y) * exp(-x**2 - k**2/4.0)
                ! 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) = sinh(k*y) * exp(-x**2 - k**2/4.0)

                ! res = 2.*x*(1. - cos(2*x*y)*cosh(k*y)) + k*sin(2*x*y)*sinh(k*y)
                ! res = ( 2.*x*(1. - cos(2*x*y)*cosh(k*y)) + k*sin(2*x*y)*sinh(k*y) ) * exp(-x**2 - k**2/4.0)
                ! res = ( 2.*x*(1. - cos(2*x*y)*( 0.5*(exp(k*y)+exp(-k*y)) )   ) + &
                !     & k*sin(2*x*y)*( 0.5*(exp(k*y)-exp(-k*y)) ) ) * exp(-x**2 - k**2/4.0)
                res = 2.*x*(exp(-x**2 - k**2/4.0) &
                    & - cos(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) ) ) + &
                    & k*sin(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) )


            end function f_k

            elemental function g_k(k, x, y) result(res)
                implicit none

                integer, intent(in) :: k
                real, intent(in) :: x, y
                real :: res
                ! 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) = cosh(k*y) * exp(-x**2 - k**2/4.0)
                ! 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) = sinh(k*y) * exp(-x**2 - k**2/4.0)

                ! res = 2*x*sin(2*x*y)*cosh(k*y) + k*cos(2*x*y)*sinh(k*y)
                ! res = ( 2*x*sin(2*x*y)*cosh(k*y) + k*cos(2*x*y)*sinh(k*y) ) * exp(-x**2 - k**2/4.0)
                ! res = ( 2*x*sin(2*x*y)*( 0.5*(exp(k*y)+exp(-k*y)) ) + &
                !     & k*cos(2*x*y)*( 0.5*(exp(k*y)-exp(-k*y)) ) ) * exp(-x**2 - k**2/4.0)
                res = 2*x*sin(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) ) + &
                    &   k*cos(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) )
            end function g_k

    end function erf_c

    ! -----------------------------------------
    ! ERF OF dual_complex numbers
    ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(PI)*exp(-u**2) * du>
    ! -----------------------------------------
    elemental function erf_dc(u) result(res)
      type(dual_complex), intent(in) :: u
      type(dual_complex) :: res

      res%z  = erf(u%z)
      res%dz = u%dz * 2.0/sqrt(PI) * exp(-u%z**2)

    end function erf_dc

end module dnadmod
