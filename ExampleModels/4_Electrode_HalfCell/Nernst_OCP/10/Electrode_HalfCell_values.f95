! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py
! Comments that describe equations frequently use UNICODE characters to do so

module user_input
  implicit none

  integer, parameter :: N = 4
  integer, parameter :: NJ = 102                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 3.6d3*1000     ! Number of time steps
  real               :: delT = 1e-0                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = ''

  ! **************************** Physical Constants ****************************
  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real, parameter :: density_active_material  = 5            ! [g/cm3]
  real, parameter :: max_mAhg       = 200.0
  real, parameter :: mAhg_to_conc   = 3.6 * density_active_material / Fconst  ! mAh/g (x C/mAh * g/cm3 * mol/C --> mol/cm3)
  real, parameter :: max_Li_conc = max_mAhg * mAhg_to_conc

  ! **************************** Physical Parameters ***************************
  real :: diff_Li  = 1e-7                       ! diffusion coefficient [cm^2/s]
  real :: diff_PF6 = 1e-7                       ! diff coeff PF6
  real :: z_1 = 1.0                             ! cation charge z_Li
  real :: z_2 = -1.0                            ! anion charge, z_PF6

  real :: cbulk = 0.001
  real :: porosity = 0.4
  real :: volumetric_surface_area

  real :: xmax          = 100e-4              ! 500 um is 500e-4 cm
  real :: crystal_xmax  = 500e-7               ! Crystal radius (cm) [10 nm]

  ! real :: c_0_init  = cbulk               ! initial c0 conc  [mol/cm3]
  real :: c_x_init   = 1.0e-7
  real :: Phi_1_init = 3.4
  real :: Phi_2_init = 0.0
  real :: sigma      = 3.0e-4

  real :: CCV_Voltage = 3.0
  real :: applied_specific_current             ! [mA/g]  926 mAh/g / xxx hours
  real :: applied_current_A
  real :: mAhg = 0.0                           ! Cummulative mAh/g

  integer :: cycle_number = 1                   ! cycle number counter


  character(len=3) :: state = 'D'               ! Discharge - D, Charge - C, Recovery - R

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
! ******************************************************************************

! ******************************************************************************
module echem_mod
  use user_input, only: N, Rigc, Fconst, Temp, max_Li_conc, cbulk
  use dnadmod

  implicit none

  real :: alpha_a = 0.5
  real :: alpha_c = 0.5
  real :: k_rxn   = 1.0e-7

contains


  function OCP(c0, c_x) result(OpenCircuitVoltage)
    type(dual)                :: OpenCircuitVoltage
    type(dual), intent(in)    :: c0, c_x
    type(dual)                :: theta
    real                      :: U_ref


    theta = c_x / max_Li_conc

    U_ref =  3.3

    OpenCircuitVoltage = U_ref + Rigc*Temp/Fconst * LOG( c0/cbulk * (1.0 - theta)/theta )

  end function OCP


  function exchange_current_density(c0, c_x)  result(ex_curr)
  type(dual)              :: ex_curr
  type(dual), intent(in)  :: c0, c_x

  ex_curr = Fconst * k_rxn * c0**alpha_a * c_x**alpha_c * (max_Li_conc - c_x)**alpha_a

  end function exchange_current_density


  function Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)    result(i_BV)
    ! positive when eta > 0; negative when eta < 0
    ! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
    ! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
    ! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
    type(dual)              :: i_BV
    type(dual), intent(in)  :: c0, c_x, Phi_1, Phi_2
    type(dual)              :: eta, U_ocp
    type(dual)              :: i_ex_curr

    i_ex_curr = exchange_current_density(c0, c_x)

    U_ocp = OCP(c0, c_x)
    eta = Phi_1 - Phi_2 - U_ocp

    i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                        - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

  end function Butler_Volmer_Rxn


end module echem_mod


! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use echem_mod

  implicit none

  real :: c_0_Li, c_x_Li, Phi_1, Phi_2

  character(len=65) :: header_fmt = '(1(A5,1X), 1(A5,1X), 3(A12,1X), 20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(I5,1X), 1(A5,1X), 3(F12.5,1X), 20(ES15.5,1X))'


  real :: t_write
  real :: Equiv
  real :: Li_balance


contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0
    Equiv       = mAhg / max_mAhg! * 8.

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/50.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it == 1) then
      ! call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    else if ( (time - last_write_time).GE.write_every_x_sec ) then
      ! call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Cycle', 'State', 'Time', 'Equiv',    'Voltage', 'Solution_Pot', 'c_Li_0', 'conc_x', 'Capacity'
      write(*, header_fmt) '#', 'CDR',  'hours', 'LixH', 'Volts'  ,      'Volts',        'mol/L',  'LixH', 'mAh/g'
    end if


    c_0_Li = cprev(1,NJ)*1e3
    c_x_Li = cprev(2,1) / max_Li_conc !* 8.0        ! Li_xFe_3O_4
    Phi_1  = cprev(3,NJ)
    Phi_2  = cprev(4,NJ)

    write(*, data_fmt) cycle_number, state, t_write, Equiv, Phi_1, Phi_2, c_0_Li, c_x_Li, mAhg

  end subroutine write_to_screen





  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Cycle', 'State', 'Time', 'Equivalence', &
      & 'Position', 'Voltage', 'Soln_Conc', 'Solid_Conc', 'Solution_Pot', &
      & 'Capacity'
      write(56, header_fmt) '#',     'CDR'  , 'hours', 'LixH' ,  &
      & 'um'      , 'Volts'  , 'mol/cm3'  , 'LixH', 'Volts', &
      & 'mAh/g'
                                                    !           !
    end if                                          !           !

    do j = 1, NJ                                    !           !
      c_0_Li = cprev(1,j) * 1e3
      c_x_Li = cprev(2,j) / max_Li_conc !* 8.0        ! Li_xFe_3O_4
      Phi_1  = cprev(3,j)
      Phi_2  = cprev(4,j)
                                                    !           !
      write(56, data_fmt) cycle_number, state,    t_write,    Equiv,       xx(j)*1e4,    Phi_1,      c_0_Li,  c_x_Li, Phi_2, mAhg

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  use dnadmod
  use echem_mod
  implicit none

  type(dual) :: c0, c_x, Phi_1, Phi_2
  type(dual) :: dc0dx, dc_xdx, dPhi_1dx, dPhi_2dx
  type(dual), dimension(N)               :: flux_temp, accum_temp, rxn_temp

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

  ! (1) ϵ ∂cₒ/∂t = D∇²cₒ + a iᵣₓ / F
  ! 𝐍ₒ = -ϵ * (D_Li * ∇cₒ + z_Li * u_Li * cₒ F ∇Φ₂)
  ! Rₒ =  a iᵣₓ / F
  ! BC-WEST : cₒ = cbulk
  ! BC-EAST : 𝐍ₒ = 0

  ! (2) (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F
  ! 𝐍ₓ = 0
  ! Rₓ = -a iᵣₓ / F
  ! BC-WEST : (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F
  ! BC-EAST : (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F

  ! (3) 0 = -∇⋅𝐢₁ - a iᵣₓ
  ! 𝐢₁ = -(1-ϵ) σ ∇Φ₁
  ! Rₓ =  a iᵣₓ
  ! BC-WEST : 𝐢₁ = 0                (Solid-State Current = 0)
  ! BC-EAST : Φ₂ = 0               (arbitrary ref Voltage)

  ! (4) 0 = -∇⋅𝐢₂ + a iᵣₓ
  ! 𝐢₂/F = ∑ᵢ (zᵢ 𝐍ᵢ)
  ! Rₓ = -a iᵣₓ
  ! BC-WEST : 𝐢₁ = i_applied     (All current carried in solid state)
  ! BC-EAST : 𝐢₂ = 0             (Solution Current = 0)

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   𝐍ₒ = -ϵ * (D_Li * ∇cₒ + z_Li * u_Li * cₒ F ∇Φ₂)
    ! (2)   𝐍ₓ = 0
    ! (3)   𝐢₁ = -(1-ϵ) σ ∇Φ₁
    ! (4)   𝐢₂/F = -ϵ * [ (z₁ D₁ + z₂ D₂) ∇cₒ
    !                       ((z₁)² * u₁ + (z₂)² * u₂) F cₒ ∇Φ₂]
    !       i.e. 𝐢₂/F = ∑ᵢ (zᵢ 𝐍ᵢ)
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: diff_
    real                                  :: u_1, u_2    ! mobility

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    dc0dx    = dcdx_vars_dual(1)
    dc_xdx   = dcdx_vars_dual(2)
    dPhi_1dx = dcdx_vars_dual(3)
    dPhi_2dx = dcdx_vars_dual(4)

    ! Nernst-Einstein Relationship
    u_1 = diff_Li / (Rigc * Temp)
    u_2 = diff_PF6 / (Rigc * Temp)

    Flux_(1) = -porosity**1.5 * (diff_Li * dc0dx + z_1 * u_1 * c0 * Fconst * dPhi_2dx)
    Flux_(2) = 0.0
    Flux_(3) = -(1 - porosity) * sigma * dPhi_1dx
    Flux_(4) = -porosity**1.5 * Fconst * ( (z_1 * diff_Li + z_2 * diff_PF6) * dc0dx &
        + (z_1**2 * u_1 + z_2**2 * u_2) * Fconst * c0 * dPhi_2dx )

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  a * iᵣₓ / F
    ! (2) -a * iᵣₓ / F
    ! (3) -a * iᵣₓ
    ! (4)  a * iᵣₓ
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual
    type(dual) :: i_rxn

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    i_rxn = Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)

    Rxn_(1) = +volumetric_surface_area * i_rxn / Fconst     ! c_0
    Rxn_(2) = -volumetric_surface_area * i_rxn / Fconst     ! c_x
    Rxn_(3) = -volumetric_surface_area * i_rxn              ! i_1
    Rxn_(4) = +volumetric_surface_area * i_rxn              ! i_2

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) ϵ * ∂cₒ/∂t
    ! (2) (1-ϵ) * ∂cₓ/∂t
    ! (3) 0                       (no accumulation of e-)
    ! (4) 0                       (electroneutrality - ions)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    Accum_(1) = (porosity) * c0/delT
    Accum_(2) = (1.0 - porosity) * c_x/delT
    Accum_(3) = 0.0                  ! no accumulation of e-
    Accum_(4) = 0.0                  ! electroneutrality

    Accum_%x  = 0.0

  end function ACCUM
! ******************************************************************************


! ******************************************************************************
! **************************** BOUNDARY CONDITIONS *****************************
! ******************************************************************************
! boundary conditions need to be written as homogeneous equations
! i.e.    c = 1.0   -->      c - 1.0 = 0.0    -->   BC_WEST_(N) = c - 1.0
! i.e. flux = 1.0   -->   flux - 1.0 = 0.0    -->   BC_WEST_(N) = flux - 1.0
! ********** special case **********
! no spatial gradients in governing equation - can repeat gov eqn at boundary
! BC_WEST_(N) = accum_temp(N) - rxn_temp(N)
  function Boundary_WEST (c_vars_dual, dcdx_vars_dual) result(BC_WEST_)
    ! (1)   cₒ = cbulk
    ! (2)   (1 - ϵ) ∂c/∂t = -a iᵣₓ/F
    ! (3)   𝐢₁ = 0                               (Solid-State Current = 0)
    ! (4)   Φ₂ = 0                             (arbitrary ref Voltage)
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    flux_temp   = FLUX(c_vars_dual, dcdx_vars_dual)
    accum_temp  = ACCUM(c_vars_dual)
    rxn_temp    = RXN(c_vars_dual)

    BC_WEST_(1) = c0 - cbulk                                      ! c0 = cbulk
    BC_WEST_(2) = accum_temp(2) - rxn_temp(2)        ! (1-eps)dc/dt = -ai_rxn/F
    BC_WEST_(3) = flux_temp(3) - 0.0                              !   i_1 = 0.0
    BC_WEST_(4) = Phi_2 - 0.0                                     ! Phi_2 = 0.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   𝐍ₒ = 0
    ! (2)   (1 - ϵ) ∂c/∂t = -a iᵣₓ/F
    ! (3)   𝐢₁ = i_applied
    ! (4)   𝐢₂ = 0                           (Solution Current = 0)
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    flux_temp   = FLUX(c_vars_dual, dcdx_vars_dual)
    accum_temp  = ACCUM(c_vars_dual)
    rxn_temp    = RXN(c_vars_dual)

    BC_EAST_(1) = flux_temp(1) - 0.0
    BC_EAST_(2) = accum_temp(2) - rxn_temp(2)        ! (1-eps)dc/dt = -ai_rxn/F
    BC_EAST_(3) = flux_temp(3) - applied_current_A              ! i_1 = i_app
    BC_EAST_(4) = flux_temp(4) - 0.0                            ! Phi_2 = 0.0

  end function Boundary_EAST

end module GOV_EQNS



module experiment_mod
  use user_input
  use variables, only: cprev, delX
  use echem_mod
  use write_data_mod
  implicit none

  real    :: c_rate
  real    :: discharge_time
  real    :: discharge_time_0
  real    :: charge_time, rest_time, CCV_time
  character(len=:), allocatable :: experimental_program

contains

  subroutine experimental_conditions(it)
  integer :: it

  experimental_program = 'C_200_Discharge_Rest_Charge'

! *****************************************************************************
  if (experimental_program == 'C_200_Discharge_Rest_Charge') then

    c_rate = 1.0/10
    state = 'D'
    applied_specific_current = +max_mAhg * c_rate

    ! discharge_time = 1.0/c_rate * ( 8.0 / 8.0) * 3600.0

    ! if (time <= discharge_time) then
    !   state = 'D'
    !   applied_specific_current = +max_mAhg * c_rate
    !
    ! else if (time <= (discharge_time + rest_time)) then
    !   state = 'R'
    !   applied_specific_current = 0.0
    !
    ! else if (time > (discharge_time + rest_time)) then
    !   state = 'C'
    !   applied_specific_current = -max_mAhg * c_rate
    !
    ! end if

! ******************************************************************************
!   else if (experimental_program == 'Other_Exp_Program') then
! ...
! ...
! ...
!
  end if

! ******************************************************************************
  mAhg = mAhg + applied_specific_current * delT/3600.0

  if (trim(geometry) == 'Rectangular') then
    applied_current_A = applied_specific_current/1000 * (1 - porosity) * density_active_material * xmax

  else if (trim(geometry) == 'Cylindrical') then
    applied_current_A = applied_specific_current/1000 * (1 - porosity) * density_active_material * xmax/2.0

  else if (trim(geometry) == 'Spherical') then
    applied_current_A = applied_specific_current/1000 * (1 - porosity) * density_active_material * xmax/3.0

  end if

  if (time > 100 * 3600) then
    print*, 'Time > 1000 hours'
    stop
  else if (cprev(1,NJ) < 0.0) then
    call t_write__Equiv__Li_Bal_Assignment
    ! call write_to_screen(0)
    print*, 'c0 < 0.0'
    stop
  else if (any(isnan(cprev))) then
    call t_write__Equiv__Li_Bal_Assignment
    ! call write_to_screen(0)
    print*, 'Dependent Variable is NaN'
    print*, 'iterations', it
    stop
  end if

  end subroutine experimental_conditions


end module experiment_mod

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
  use echem_mod
  use experiment_mod
  use BAND_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max, it
  integer :: i, k, j
  real :: crystal_volume, crystal_surface_area

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  crystal_surface_area    = 4.0 * PI * crystal_xmax**2
  crystal_volume          = 4/3 * PI * crystal_xmax**3
  volumetric_surface_area = (1 - porosity)/crystal_volume * crystal_surface_area

  do it = 1, Numbertimesteps

    call experimental_conditions(it)
    call write_condition(it)

    do j = 1, NJ
      call auto_fill(j)               ! These can be changed to functions
      call ABDGXY(j)                  ! because they each have a set of inputs
      call BAND(j)                    ! and a set of outputs
    end do

    cprev = cprev + delC              ! Update the dependent variables
    time  = time  + delT              ! Update the time

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
  ! use GOV_EQNS, only: Control_Volume, Cross_Sectional_Area,  &    ! functions
  !                     Cntrl_Vol, Crx_Area                         ! variables

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
      xx(j) = h*float(j-1) - h/2.0
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
    cprev(2,j) = c_x_init
    cprev(3,j) = Phi_1_init
    cprev(4,j) = Phi_2_init
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition
