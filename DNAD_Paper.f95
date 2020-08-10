TYPE, PUBLIC :: DUAL_NUM
  REAL(DBL_AD) :: x_ad_
  REAL(DBL_AD) :: xp_ad_(NDV_AD)
END TYPE DUAL_NUM

INTERFACE SIN
  MODULE PROCEDURE SIN_D
END INTERFACE

ELEMENTAL FUNCTION SIN_D(u) RESULT(res)
  TYPE (DUAL_NUM), INTENT(IN) :: u
  TYPE (DUAL_NUM) :: res
  REAL (DBL_AD)   :: tmp

  res%x_ad_ = SIN(u%x_ad_)
  tmp = COS(u%x_ad_)
  res%xp_ad_= u%xp_ad_*tmp
END FUNCTION SIN_D

INTERFACE OPERATOR(**)
  MODULE PROCEDURE POW_D
END INTERFACE

ELEMENTAL FUNCTION POW_D(u,v) RESULT(res)
  TYPE (DUAL_NUM), INTENT(IN) :: u
  TYPE (DUAL_NUM), INTENT(IN) :: v
  REAL (DBL_AD) :: uf, vf
  TYPE (DUAL_NUM) :: res

  uf = u%x_ad_
  vf = v%x_ad_
  res%x_ad_ = uf**vf
  res%xp_ad_ = res%x_ad_*(vf/uf*u%xp_ad_ + LOG(uf)*v%xp_ad_)
END FUNCTION POW_D


! PROGRAM CALL_SIN
! SIN_D(0.0)
! END PROGRAM CALL_SIN
