! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module user_input
  implicit none

  integer, parameter :: N = 2
  integer, parameter :: NJ = 102                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 30 * 3.6e3     ! Number of time steps
  integer, parameter :: maxIterations = 1e6
  real               :: delT = 1.0e-2                       ! size of timestep [s]
  real               :: time                             ! [s]
  real, parameter    :: xmax = 1e2                      ! 1 cm

  logical            :: UPWIND = .TRUE.
  character(len=65)  :: direction = 'WestToEast'

  real, parameter    :: R_1  = 0
  real, parameter    :: R_NJ = R_1 + xmax

  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant
  real, parameter :: h_heat_transfer = 1.0
  real, parameter :: c_heat_cap = 4.0
  real, parameter :: T_ambient = 25 + 273
  real, parameter :: T_in = 200 + 273
  real, parameter :: vel_in = 10

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
! ******************************************************************************


module fluid_properties_mod
  use dnadmod
  use user_input

  implicit none

  interface Density
    module procedure Density_dual
    module procedure Density_real
  end interface Density

  contains

    function Density_dual(Temperature)    Result(density_)
      type(dual)              :: density_
      type(dual), intent(in)  :: Temperature

      density_ = -0.5*(Temperature/300 - 1)**2 + 1

    end function Density_dual

    function Density_real(Temperature)    Result(density_)
      real              :: density_
      real, intent(in)  :: Temperature

      density_ = -0.5*(Temperature/300 - 1)**2 + 1

    end function Density_real

end module fluid_properties_mod



! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use fluid_properties_mod

  implicit none

  character(len=65) :: header_fmt = '(1(A12,1X),   20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(F12.5,1X), 20(ES15.5,1X))'

  real :: t_write

  real :: Temp_1, Temp_NJ, vel_1, vel_NJ
  real :: Temp, vel, density_

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 1.0e-1           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment

    if (it == 0) then
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

    real :: mass_in, mass_out

    if (it == 0) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Time', 'Temp_1', 'Temp_NJ', 'vel_1', 'vel_NJ', 'Density_in', 'Density_out', 'mass_in', 'mass_out'
      write(*, header_fmt) 'hours', 'Kelvin', 'Kelvin', 'cm/s' , 'cm/s'  , 'g/cm3',      'g/cm3',       'g/(cm2 * s)', 'g/(cm2 * s)'
    end if

      Temp_1   = cprev(1,1)
      vel_1    = cprev(2,1)
      Temp_NJ  = cprev(1,NJ)
      vel_NJ   = cprev(2,NJ)
      mass_in  = vel_1 * density(Temp_1)
      mass_out = vel_NJ * density(Temp_NJ)

      ! print*,
      write(*, data_fmt) t_write, Temp_1, Temp_NJ, vel_1, vel_NJ, density(Temp_1), density(Temp_NJ), mass_in, mass_out


  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it == 0) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time', 'Position', 'Temp',   'Velocity', 'Density'
      write(56, header_fmt) 'hours', 'cm',      'Kelvin', 'cm/s' ,    'g/cm3'
                                                      !
    end if                                              !


    do j = 1, NJ                                    !
      Temp    = cprev(1,j)
      vel     = cprev(2,j)
      density_ = Density(Temp)    !
                                                    !
      write(56, data_fmt) t_write, xx(j), Temp, vel, density_

    end do

  end subroutine write_positional_information_to_file

end module write_data_mod







module GOV_EQNS
  use user_input
  use fluid_properties_mod
  use dnadmod

  implicit none

  type(dual) :: Temp, dTdx, vel, dveldx
  type(dual) :: density_

contains

! (1) Energy Balance (Heat)
! ∂(cᵨρT)/∂t = -∇⋅(cᵨρT 𝐯) - h(T - Tₐ)
! (2) Continuity Equation
! ∂(ρ)/∂t = -∇⋅(ρ𝐯)

! ρ = f(T)

! BC - WEST
! T = Tᵢ
! 𝐯 = 𝐯ᵢ
! BC - EAST
! ∇T = 0
! ∇⋅𝐯 = 0

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual

    Temp     = c_vars_dual(1)
    dTdx     = dcdx_vars_dual(1)

    vel      = c_vars_dual(2)
    dveldx   = dcdx_vars_dual(2)

    density_ = Density(Temp)

    Flux_(1) = c_heat_cap * density_ * Temp * vel
    Flux_(2) = density_ * vel

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    Temp     = c_vars_dual(1)

    Rxn_(1) = -h_heat_transfer * (Temp - T_ambient)
    Rxn_(2) = 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    Temp     = c_vars_dual(1)
    vel      = c_vars_dual(2)

    density_ = Density(Temp)

    Accum_(1) = c_heat_cap * density_ * Temp/delT
    Accum_(2) = density_/delT
    ! Accum_(2) = density_ * vel / delT

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
    type(dual), dimension(N)               :: BC_WEST_, flux_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    Temp        = c_vars_dual(1)
    vel         = c_vars_dual(2)

    BC_WEST_(1) = Temp - T_in
    BC_WEST_(2) = vel - vel_in


  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    type(dual), dimension(N)               :: BC_EAST_, flux_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    Temp        = c_vars_dual(1)
    dTdx        = dcdx_vars_dual(1)

    vel         = c_vars_dual(2)
    dveldx      = dcdx_vars_dual(2)

    density_    = Density(Temp)

    BC_EAST_(1) = dTdx - 0.0    ! temperature does not change past the exit

    ! ∇⋅(ρ𝐯) = 0 = 𝐯⋅∇ρ + ρ∇𝐯
    ! ∇ρ = ∂ρ/∂T ⋅ ∇T
    BC_EAST_(2) = density_%dx(1)*dTdx * vel + density_*dveldx - 0.0


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
  use fluid_properties_mod
  use variables, only : cprev, delC
  use write_data_mod
  use BAND_MOD

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j

  it = 0

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()
  call write_condition(it) ! writes the headers
  it = 1
  do while (time <= xmax/vel_in*2.0) !it = 1, 360*2 !Numbertimesteps

      ! **********************************************************************
      ! ******************************** BOUND VAL ***************************
      ! **********************************************************************
        do j = 1, NJ
          call auto_fill(j)
          call ABDGXY(j)
          call BAND(j)
        end do

        cprev = cprev + delC      ! Update the dependent variables
        time  = time + delT
        ! ********************************************************************
        call write_condition(it)

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
  use fluid_properties_mod

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
    cprev(1,j) = T_in
    cprev(2,j) = vel_in
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  delT = h/vel_in/10.0
  print*, delT

  return                                                  ! is this necessary?
end subroutine initial_condition
