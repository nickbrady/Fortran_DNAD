! BAND_mod contains the subroutines that make up the numerical engine which solves these coupled non-linear PDE's
! They are written so that they are general for all problems and therefore rarely need to be modified
!
! subroutine auto_fill needs to be "included" after module GOV_EQNS, because auto_fill has a "use GOV_EQNS" statement
! ABDGXY, BAND, MATINV only need to be included after module variables
!
! The philosophy that guides this partitioning of code is two-fold:
! 1. We want to minimize the amount of code accessible to the general user. This will hopefully improve readability,
!    and help with debugging. Code that rarely needs to be changed or should only be changed by a highly knowledgable
!    person (i.e. expert) should not generaly be available to common user.
! 2. Make it easier to maintain a master code. When changes or improvements are made, these can then be easily be
!    passed along to other pieces of code, makes it easier to implement tests for compatibility
! CORE_MODULES_PATH = __CORE__MODULES_PATH__
module BAND_mod

  contains  ! these files rarely need to be modified
    ! include statements need to remain on one line (cannot use line continuation)
    ! this limits the length of the path name that can be used
    include '__CORE__MODULES_PATH__/sub_auto_fill.f95'
    include '__CORE__MODULES_PATH__/sub_ABDGXY.f95'
    include '__CORE__MODULES_PATH__/sub_MATINV.f95'
    include '__CORE__MODULES_PATH__/sub_BAND.f95'

end module BAND_mod
