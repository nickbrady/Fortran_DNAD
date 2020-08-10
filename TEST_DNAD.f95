! TEST DNAD
MODULE DNAD

  ! INTEGER :: NDV_AD
  ! NDV_AD = 1

  TYPE DUAL_NUM
    REAL(8) :: x_ad_
    ! REAL(8) :: xp_ad_(NDV_AD)           ! need to figure out how to dynamically
    REAL(8) :: xp_ad_(1)                  ! change NDV_AD

  END TYPE DUAL_NUM

  INTERFACE OPERATOR(**)
  	MODULE PROCEDURE POW_D
  END INTERFACE

CONTAINS

  ELEMENTAL FUNCTION POW_D (u, v) RESULT (res)
    TYPE (DUAL_NUM)               :: res
    TYPE (DUAL_NUM), INTENT(IN)   :: u, v
    REAL(8)                       :: uf, vf

  	uf = u%x_ad_
  	vf = v%x_ad_
  	res%x_ad_ = uf**vf
  	res%xp_ad_ = res%x_ad_*(vf/uf*u%xp_ad_ + LOG(uf)*v%xp_ad_)

  END FUNCTION  POW_D

END MODULE DNAD


PROGRAM CircleArea
  use DNAD
  ! REAL(8) :: PI = 4.0D0*ATAN(1.0D0)
  ! REAL(8) :: radius, area
  ! USE dnadmod
	! use DUAL_NUM
  TYPE (DUAL_NUM) :: PI=DUAL_NUM(4.0D0*ATAN(1.0D0),(/0.D0/))
  ! TYPE (DUAL_NUM) :: radius, area

  ! READ(*,*) radius
	! Area = PI * radius ** 2
  ! Area=radius**2
	! Area = POW_D(radius,2)

  ! WRITE(*,*) "AREA=", Area
	! WRITE(*,*) 'radius = ', radius
	write(*,*) 'PI =', PI

END PROGRAM CircleArea
