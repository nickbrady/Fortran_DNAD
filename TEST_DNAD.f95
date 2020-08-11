!

MODULE DNAD

  ! INTEGER :: NDV_AD
  ! NDV_AD = 1

  TYPE DUAL
    real :: x        ! functional value
    real :: dx(1)    ! derivative
  END TYPE DUAL

  interface operator (**)
      module procedure pow_i ! dual number to an integer power,elemental
      module procedure pow_r ! dual number to a real power, elemental
      module procedure pow_d ! dual number to a dual power, elemental
  end interface

  interface operator (*)
      module procedure mult_dd    ! dual*dual, elemental
      module procedure mult_di    ! dual*integer,elemental
      module procedure mult_dr    ! dual*real,elemental
      module procedure mult_id    ! integer*dual,elemental
      module procedure mult_rd    ! real*dual,elemental
  end interface

CONTAINS

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

  !******* END: (**)
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

  !******* END: (*)
  !---------------------

END MODULE DNAD


PROGRAM CircleArea
  use DNAD

  TYPE (DUAL) :: PI=DUAL(4.0D0*ATAN(1.0D0),(/0.D0/))
  TYPE (DUAL) :: radius, area

  radius = DUAL(3.0, (/1.D0/) )

	Area = PI * radius ** 2

  WRITE(*,*) "AREA=", Area
	WRITE(*,*) 'radius = ', radius
	write(*,*) 'PI =', PI

END PROGRAM CircleArea
