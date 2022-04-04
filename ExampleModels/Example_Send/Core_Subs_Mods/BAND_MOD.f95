module BAND_mod
  ! BAND_mod contains the subroutines that make up the numberical engine which solves these coupled non-linear PDE's
  ! They are written so that they are general for all problems and therefore only rarely need to be modified

  ! use variables

  contains  ! these files rarely need to be modified
! needs to be after GOV_EQNS
! include statements need to remain on one line (cannot use line continuation)
include '../Core_Subs_Mods/sub_auto_fill.f95'
include '../Core_Subs_Mods/sub_ABDGXY.f95'
include '../Core_Subs_Mods/sub_MATINV.f95'
include '../Core_Subs_Mods/sub_BAND.f95'

end module BAND_mod
