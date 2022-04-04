! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module user_input
  implicit none

  integer, parameter :: N = 1
  integer, parameter :: NJ = 42                          ! Number of mesh points

  integer, parameter :: Numbertimesteps = 30 * 3.6e3     ! Number of time steps
  real               :: delT = 1.0                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = ''

  real               :: xmax = 1e-4                      ! 500 um is 500e-4 cm

  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real :: k_rxn = 1e-4
  real :: diff  = 1e-13                       ! diffusion coefficient [cm^2/s]
  real :: cbulk = 0.001                       ! concentration [mol/cm3]

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
! ******************************************************************************

! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  real :: c0

  character(len=65) :: header_fmt = '(2(A12,1X),   20(A15,1X)) '
  character(len=65) :: data_fmt   = '(2(F12.5,1X), 20(ES15.5,1X))'

  real :: t_write

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/10.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it == 1) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    else if ( (time - last_write_time).GE.write_every_x_sec ) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Time',  'Conc'
      write(*, header_fmt) 'hours', 'mol/L'
    end if

    c0 = cprev(1,NJ)*1e3

    write(*, data_fmt) t_write,    c0

  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time',  'Position', 'Conc'
      write(56, header_fmt) 'hours', 'um'      , 'mol/L'
                                                    !
    end if                                          !
                                                    !
    do j = 1, NJ                                    !
      c0 = cprev(1,j) * 1e3                         !
                                                    !
      write(56, data_fmt) t_write,    xx(j)*1e4,    c0

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  use dnadmod

  implicit none

contains
! General Equations
! Analytic:
!     ∂cᵢ/∂t = -∇⋅𝐍ᵢ + Rᵢ
! Finite Volume (Control Volume):
!     ΔV ∂c/∂t = (Aₓᵢ⋅𝐍ᵢ - Aₓₒ⋅𝐍ₒ) + ΔV ⋅ Rⱼ
!
! frequently, the control volume (ΔV) and cross-sectional area (Aₓ) are given
! by the system geometry
! 𝐍ᵢ is the flux of specie i
! Rᵢ is the rate of generation (reaction rate) of specie i

! (1) ∂c/∂t = D∇²c
! 𝐍 = -D ∇c
! R = -kᵣₓ c

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

    type(dual) :: c0
    type(dual) :: dc0dx

    c0       = c_vars_dual(1)
    dc0dx    = dcdx_vars_dual(1)

    Flux_(1) = -diff * dc0dx

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  -k * c0
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    type(dual) :: c0

    c0    = c_vars_dual(1)

    Rxn_(1) = -k_rxn * c0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) dc0/dt
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    type(dual) :: c0

    c0    = c_vars_dual(1)

    Accum_(1) = c0/delT

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
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    type(dual) :: c0

    c0      = c_vars_dual(1)

    BC_WEST_(1) = c0 - cbulk

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   N0 = 0.0
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp

    type(dual) :: c0

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    c0    = c_vars_dual(1)

    BC_EAST_(1) = flux_temp(1) - 0.0                            ! c0 = 0.1*cbulk
    ! print*, flux_temp(1), dcdx_vars_dual(1)

  end function Boundary_EAST

end module GOV_EQNS

! ******************************************************************************
! auto_fill, ABDGXY, BAND, MATINV
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/BAND_MOD.f95'
! ******************************************************************************

! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************
program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use BAND_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  do it = 1, Numbertimesteps

    call write_condition(it)
                                        ! function to dynamically and integlligently change delT
    ! **************************************************************************
    ! ******************************** BOUND VAL *******************************
    ! **************************************************************************
      do j = 1, NJ
        call auto_fill(j)               ! These can be changed to functions
        call ABDGXY(j)                  ! because they each have a set of inputs
        call BAND(j)                    ! and a set of outputs
      end do

      cprev = cprev + delC              ! Update the dependent variables
      time  = time  + delT              ! Update the time
  ! ****************************************************************************

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

  implicit none
  real    :: h
  integer :: j

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-2)
  do j=1,NJ

    if (j.EQ.1) then
      xx(j) = 0.0
    else if (j.EQ.NJ) then
      xx(NJ) = xmax
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
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition
