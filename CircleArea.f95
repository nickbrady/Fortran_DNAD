! TYPE, PUBLIC :: DUAL_NUM
! 			REAL(DBL_AD) :: x_ad_
! 			REAL(DBL_AD) :: xp_ad_(NDV_AD)
!  END TYPE DUAL_NUM

INTEGER(2), PARAMETER :: NDV_AD=1


TYPE DUAL_NUM
			REAL(8) :: x_ad_
			REAL(8) :: xp_ad_(NDV_AD)
END TYPE DUAL_NUM

! INTERFACE OPERATOR(**)

! 	MODULE PROCEDURE POW_D
! END INTERFACE
!
! ELEMENTAL FUNCTION POW_D(u,v) RESULT(res)

!
! 	TYPE (DUAL_NUM), INTENT(IN) :: u
! 	TYPE (DUAL_NUM), INTENT(IN) :: v
! 	REAL(8) :: uf,vf
! 	TYPE (DUAL_NUM)::res
!
! 	uf = u%x_ad_
! 	vf = v%x_ad_
! 	res%x_ad_ = uf**vf
! 	res%xp_ad_ = res%x_ad_*(vf/uf*u%xp_ad_ + LOG(uf)*v%xp_ad_)
!
! END FUNCTION POW_D

! INTERFACE OPERATOR(**)
! 	FUNCTION POW_D(u,v) RESULT(res)
! 		import DUAL_NUM
!
! 		TYPE (DUAL_NUM), INTENT(IN) :: u
! 		TYPE (DUAL_NUM), INTENT(IN) :: v
! 		REAL(8) :: uf,vf
! 		TYPE (DUAL_NUM)::res
!
! 		uf = u%x_ad_
! 		vf = v%x_ad_
! 		res%x_ad_ = uf**vf
! 		res%xp_ad_ = res%x_ad_*(vf/uf*u%xp_ad_ + LOG(uf)*v%xp_ad_)
!
! 	END FUNCTION POW_D
! END INTERFACE

! FUNCTION POW_D(u,v) RESULT(res)
!
! 		uf = u%x_ad_
! 		vf = v%x_ad_
! 		res%x_ad_ = uf**vf
! 		res%xp_ad_ = res%x_ad_*(vf/uf*u%xp_ad_ + LOG(uf)*v%xp_ad_)
!
! END FUNCTION POW_D


INTERFACE OPERATOR(**)
	FUNCTION POW_D(u,v) RESULT(res)
		import DUAL_NUM
		TYPE (DUAL_NUM), INTENT(IN) :: u
		TYPE (DUAL_NUM), INTENT(IN) :: v
		REAL(8) :: uf,vf
		TYPE (DUAL_NUM)::res
	END FUNCTION POW_D
END INTERFACE




! INTERFACE SIN
! 	MODULE PROCEDURE SIN_D
! END INTERFACE
!
! FUNCTION SIN_D(u) RESULT(res)
! 	TYPE (DUAL_NUM), INTENT(IN) :: u
! 	TYPE (DUAL_NUM) :: res
! 	REAL(8) :: tmp
!
! 	res%x_ad_ = SIN(u%x_ad_)
! 	tmp = COS(u%x_ad_)
! 	res%xp_ad_ = u%xp_ad_*tmp
! END FUNCTION SIN_D
!
!
! END






! PROGRAM CircleArea
  ! REAL(8) :: PI = 4.0D0*ATAN(1.0D0)
  ! REAL(8) :: radius, area
  ! USE dnadmod
	! use DUAL_NUM
  TYPE (DUAL_NUM) :: PI=DUAL_NUM(4.0D0*ATAN(1.0D0),(/0.D0/))
  TYPE (DUAL_NUM) :: radius, area

  READ(*,*) radius
	! Area = PI * radius ** 2
  ! Area=radius**2
	! Area = POW_D(radius,2)

  ! WRITE(*,*) "AREA=", Area
	WRITE(*,*) 'radius = ', radius
	write(*,*) 'PI =', PI

END! PROGRAM CircleArea



! INTERFACE EXP
! 		MODULE PROCEDURE EXP_D
! END INTERFACE
!


! ELEMENTAL FUNCTION EXP_D(u) RESULT(res)
!       TYPE (DUAL_NUM), INTENT(IN) :: u
!       REAL(DBL_AD) :: tmp
!       TYPE (DUAL_NUM) :: res
!
!       tmp=EXP(u%x_ad_)
!       res%x_ad_ = tmp
!       res%xp_ad_ =u%xp_ad_* tmp
! END FUNCTION EXP_D



! END
