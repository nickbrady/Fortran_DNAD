
! Dependency Tree
! Example_DNAD.anc : Example_DNAD.f95 dnadmod.f95 variables.f95 BAND_MOD.f95
! dnadmod.f95 : Example_DNAD.f95
! variables.f95 : Example_DNAD.f95
! BAND_MOD.f95 : Example_DNAD.f95

! module dependencies
! number_of_variables :
! variables           : number_of_variables
! dnadmod             : number_of_variables
! user_input          : number_of_variables
! write_data_mod      : user_input variables
! GOV_EQNS            : user_input variables dnadmod
! BAND_MOD            : (number_of_variables GOV_EQNS variables dnadmod)
                      ! subroutine auto_fill  : number_of_variables GOV_EQNS variables dnadmod
                      ! subroutine ABDGXY     : number_of_variables variables
                      ! subroutine MATINV     : variables
                      ! subroutine BAND       : number_of_variables variables

! (prgoram) unsteady  : user_input, variables, write_data_mod, BAND_MOD


! Example_DNAD.f95 modules
        ! number_of_variables
        ! user_input
        ! write_data_mod
        ! GOV_EQNS
        ! (program) unsteady


! gfortran -fsyntax-only -c -J .mod Example_DNAD.f95                          --> fatal error
! gfortran -fsyntax-only -c -J .mod ../Core_Subs_Mods/variables.f95
! gfortran -fsyntax-only -c -J .mod ../Core_Subs_Mods/dnadmod.f95
! gfortran -fsyntax-only -c -J .mod Example_DNAD.f95                          --> fatal error
! gfortran -fsyntax-only -c -J .mod ../Core_Subs_Mods/BAND_MOD.f95
! make                                                                        --> success


module number_of_variables
  implicit none
  integer, parameter :: N = 1
  integer, parameter :: NJ = 22                          ! Number of mesh points
end module number_of_variables


! include '../../Core_Subs_Mods/dnadmod.f95'
! include '../../Core_Subs_Mods/variables.f95'


module user_input
  use number_of_variables
  implicit none

  integer, parameter :: Numbertimesteps = 1e4     ! Number of time steps

end module user_input


! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  real :: c0

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment

end module write_data_mod


module GOV_EQNS
  use user_input
  use variables, only: xx, delX
  use dnadmod

  implicit none

contains

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)

  end function FLUX

end module GOV_EQNS


! include '../../Core_Subs_Mods/BAND_MOD.f95'

! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use BAND_mod

  implicit none

    do j = 1, NJ
      call auto_fill(j)               ! These can be changed to functions
      call ABDGXY(j)                  ! because they each have a set of inputs
      call BAND(j)                    ! and a set of outputs
    end do

end program unsteady

! ------------------------------------------------------------------------------
! ***************************** end MAIN PROGRAM *******************************
! ------------------------------------------------------------------------------


subroutine initial_condition()
  use user_input
  use variables
  use GOV_EQNS, only: Control_Volume, Cross_Sectional_Area,  &    ! functions
                      Cntrl_Vol, Crx_Area                         ! variables

  implicit none
  !...
  return                                                  ! is this necessary?
end subroutine initial_condition
