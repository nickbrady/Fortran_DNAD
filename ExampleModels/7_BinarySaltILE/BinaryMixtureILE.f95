! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module user_input
  implicit none

  integer, parameter :: N = 2
  integer, parameter :: NJ = 202                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 1000*3600+1 !3.6d3*1000     ! Number of time steps
  real               :: delT = 1e-2                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .TRUE.
  character(len=65)  :: direction = 'EastToWest'

  real :: xmax          = 1.0              ! 500 um is 500e-4 cm

  ! **************************** Physical Constants ****************************
  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real, parameter :: z_1_Li    =  1.0           ! cation charge Li^+
  real, parameter :: z_2_EMI   =  1.0           ! cation charge, EMI^+
  real, parameter :: z_3_TFSI  = -1.0           ! anion charge of TFSI^-

  real, parameter :: nu_A_Li   = 1.0
  real, parameter :: nu_A_TFSI = 1.0
  real, parameter :: nu_A      = nu_A_Li + nu_A_TFSI
  real, parameter :: nu_B_EMI  = 1.0
  real, parameter :: nu_B_TFSI = 1.0
  real, parameter :: nu_B      = nu_B_EMI + nu_B_TFSI

  real, parameter :: MW_A_LiTFSI  = 287.09     ! molecular weight Li-TFSI  g/mol
  real, parameter :: MW_B_EMITFSI = 391.31     ! molecular weight EMI-TFSI g/mol
  real, parameter :: MW_Li = 6.941
  real, parameter :: MW_TFSI = MW_A_LiTFSI - MW_Li * nu_A_Li
  real, parameter :: MW_EMI = MW_B_EMITFSI - MW_TFSI * nu_B_TFSI
  ! **************************** Physical Parameters ***************************
  real :: diff_12 = 1e-7                       ! diffusion coefficient [cm^2/s]
  real :: diff_13 = 2e-7                       ! diff coeff PF6
  real :: diff_23 = 3e-7

  real :: cbulk = 1e-5

  real :: applied_current_A = -1.0e-5

  character(len=3) :: state = 'D'               ! Discharge - D, Charge - C, Recovery - R

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
! ******************************************************************************

module TRANSPORT_MODULE
  use user_input!, only: N, Rigc, Fconst, Temp
  use dnadmod

  implicit none

  ! type(dual)                              :: cA, cB


  interface density
    module procedure density_dual
    module procedure density_real
  end interface

contains

!*******************************************************************************
! Physical Properties that are functions of concentration
!   - mole fraction, mass fraction, chemical potential, density
!*******************************************************************************
  function c_A_to_x_A(c_A)                  result(x_A)
    type(dual)                :: x_A
    type(dual), intent(in)    :: c_A
    type(dual)                :: rho, c_B

    rho = density(c_A)
    c_B = (rho - MW_A_LiTFSI * c_A) / MW_B_EMITFSI
    x_A = c_A / (c_A + c_B)

  end function c_A_to_x_A


  function c_A_to_w_A(cA)                  result(wA)
    type(dual)                :: wA
    type(dual), intent(in)    :: cA
    type(dual)                :: rho

    rho = density(cA)
    wA = MW_A_LiTFSI * cA / rho

  end function c_A_to_w_A


  function density_dual(cA)       result(density)
    ! density data is empirically measured
    type(dual)                :: density
    type(dual), intent(in)    :: cA

    ! density%x  = 1.5
    ! density%dx = 0.0
    !
    density = 1.65 - (1.5 - 1.65)/2e-3 * (cA - 2e-3)

  end function density_dual

  function density_real(cA)       result(density)
    ! density data is empirically measured
    real                :: density
    real, intent(in)    :: cA

    ! density  = 1.5
    density = 1.65 - (1.5 - 1.65)/2e-3 * (cA - 2e-3)

  end function density_real


  function molar_volume(cA)       result(V_molar)
    type(dual)                :: V_molar
    type(dual), intent(in)    :: cA
    type(dual)                :: rho
    type(dual)                :: wA, wB

    rho = density(cA)
    wA  = c_A_to_w_A(cA)
    wB  = 1. - wA

    V_molar = 1./(rho * (1./(wA * MW_A_LiTFSI + wB * MW_B_EMITFSI) ) )

  end function molar_volume


  function partial_molar_volume_1(cA)     result(part_mol_vol)
    type(dual)                :: part_mol_vol
    type(dual), intent(in)    :: cA
    type(dual)                :: mol_vol
    type(dual)                :: xA

    mol_vol    = molar_volume(cA)
    xA         = c_A_to_x_A(cA)
    ! dV_mol_dxA =
    !
    ! part_mol_vol%x  = mol_vol *
    ! part_mol_vol%dx = 0.0

  end function partial_molar_volume_1


  function partial_molar_volume_2(cA)     result(part_mol_vol)
    type(dual)                :: part_mol_vol
    type(dual), intent(in)    :: cA

    part_mol_vol%x  = 1.0
    part_mol_vol%dx = 0.0

  end function partial_molar_volume_2


  function chemical_potential(x_A)        result(chem_pot_A)
    type(dual)                :: chem_pot_A
    type(dual), intent(in)    :: x_A          ! mol fraction of salt A

    chem_pot_A = Rigc * Temp * log(x_A ** nu_A_Li)

  end function chemical_potential
  ! du_A = u_A%dx * dxA_dx

  function activity_coefficient_gamma(x_A)   result(activity)
    type(dual)                :: activity
    type(dual), intent(in)    :: x_A          ! mol fraction of salt A

    activity%x  = 1.0
    activity%dx = 0.0

  end function activity_coefficient_gamma
!*******************************************************************************


!*******************************************************************************
! Transport Properties that are functions of concentration
!   - diffusion coefficient
!   - transference number
!*******************************************************************************
  function fundamental_diff(cA)               result(diff_fund)
    type(dual)                :: diff_fund
    type(dual), intent(in)    :: cA
    type(dual)                :: cB, c_Total
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB
    c_Total = c1 + c2 + c3

    diff_fund = z_3_TFSI**2 * c_Total / (nu_A_Li * nu_B_EMI) / &
             & (z_3_TFSI**2 * c3 / diff_12 &
             & + z_2_EMI**2 * c2 / diff_13 &
             & +  z_1_Li**2 * c1 / diff_23)

  end function fundamental_diff


  function practical_diff(cA)               result(diff_prac)
    type(dual)                :: diff_prac
    type(dual), intent(in)    :: cA
    type(dual)                :: cB, c3
    type(dual)                :: rho

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI
    c3  = nu_A_TFSI * cA + nu_B_TFSI * cB

    ! need to add activity relationship
    diff_prac = fundamental_diff(cA) * (c3 / (cA + cB)) * (nu_A_Li)

  end function practical_diff


  function transference_1_common_ion(cA)      result(t_1_c)
    type(dual)                :: t_1_c
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB

    ! cA  = c_vars_dual(1)
    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_1_c = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
    & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
    & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )

  end function transference_1_common_ion



  function transference_1_molar(cA)      result(t_1_molar)
    type(dual)                :: t_1_molar
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3, c_Total
    type(dual)                :: rho
    type(dual)                :: cB
    type(dual)                :: t_1_c

    ! cA  = c_vars_dual(1)
    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB
    c_Total = c1 + c2 + c3

    t_1_c = transference_1_common_ion(cA)

    t_1_molar = (nu_A_Li * nu_B * c3 * t_1_c - nu_B_EMI * nu_A_TFSI * c1) &
              & / (nu_A_Li * nu_B_TFSI * c_Total)

  end function transference_1_molar



  function Q_volume(cA)      result(t_1_volume)
    type(dual)                :: t_1_volume
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB

    ! cA  = c_vars_dual(1)
    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_1_volume = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
    & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
    & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )

  end function Q_volume


  function transference_1_mass(cA)      result(t_1_mass)
    type(dual)                :: t_1_mass
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB
    type(dual)                :: t_1_c

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_1_c    = transference_1_common_ion(cA)

    t_1_mass = (c3 * MW_B_EMITFSI * t_1_c - nu_A_TFSI * nu_B_EMI * cA * MW_EMI) / &
    & (nu_B_TFSI * rho)

  end function transference_1_mass


  function transference_2_mass(cA)      result(t_2_mass)
    type(dual)                :: t_2_mass
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB
    type(dual)                :: t_2_c

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_2_c    = 1.0 - transference_1_common_ion(cA)

    t_2_mass = (c3 * MW_A_LiTFSI * t_2_c - nu_B_TFSI * nu_A_Li * cB * MW_Li) / &
    & (nu_A_TFSI * rho)

  end function transference_2_mass


end module TRANSPORT_MODULE


! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx, delX
  use TRANSPORT_MODULE

  implicit none

  real :: cA, vel, cB, rho

  character(len=65) :: header_fmt = '( 1(A5,1X), 3(A12,1X), 20(A15,1X) ) '
  character(len=65) :: data_fmt   = '( 1(A5,1X), 3(F12.5,1X), 20(ES15.5,1X) )'


  real :: t_write
  real :: Li_balance


contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/100.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it <= 10*60) then
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
      write(*, header_fmt) 'State', 'Time', 'c_A_1', 'c_A_NJ', 'vel_1', 'vel_NJ'
      write(*, header_fmt) 'CDR',  'hours', 'mol/L',  'cm/s' , 'cm/s' , 'cm/s'
    end if


    cA  = cprev(1,1)*1e3
    vel = cprev(2,NJ)

    ! write(*, *) maxval(abs(cprev(2,:)))*density(cprev(1,NJ)) * delX(NJ/2),  diff_12
    write(*, data_fmt) state, t_write, cA, cprev(1,NJ)*1e3, cprev(2,1), cprev(2,NJ), cprev(2,1) * density(cprev(1,1)), &
    & cprev(2,NJ) * density(cprev(1,NJ))

  end subroutine write_to_screen





  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j
    real :: wA, xA


    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'State', 'Time', 'Position', 'Soln_Conc', 'Velocity', 'Density', 'c_B', 'Mass_Frac_A', 'Mol_Frac_A'
      write(56, header_fmt) 'CDR',  'hours',  'um'      , 'mol/L',    'nm/s',  'g/cm3', 'mol/L', '-', '-'
                                                    !           !
    end if                                          !           !

    do j = 1, NJ                                    !           !
      cA  = cprev(1,j)
      vel = cprev(2,j)  * 1e7
      rho = density(cA)
      cB  = (rho - MW_A_LiTFSI * cA)/MW_B_EMITFSI
      wA  = MW_A_LiTFSI * cA / rho
      xA  = cA / (cA + cB)
      cA  = cA * 1e3
      cB  = cB * 1e3
                                                    !           !
      write(56, data_fmt) state,    t_write,   xx(j)*1e4,   &
      & cA, vel, rho, cB, wA, xA

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod












module GOV_EQNS
  use user_input
  use dnadmod
  use TRANSPORT_MODULE
  implicit none

  type(dual) :: rho, Diff
  type(dual) :: cA, vel, cB
  type(dual) :: dcA_dx
contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   N_1 = \nu_1^A * \rho
    ! (2)   N_mass = rho * vel
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: diff, t_1_c, t_1_mass, t_2_mass
    type(dual)                            :: w_A, dwA_dx, w_B, dwB_dx
    type(dual)                            :: x_A, dxA_dx
    type(dual)                            :: u_A, duA_dx
    type(dual)                            :: vel
    type(dual)                            :: rho
    real                                  :: i_sol

    i_sol = applied_current_A

    cA     = c_vars_dual(1)
    vel    = c_vars_dual(2)

    dcA_dx = dcdx_vars_dual(1)

    ! x_A    = c_A_to_x_A(cA)              ! x_A, dxA/dcA
    ! dxA_dx = dcA_dx * x_A%dx(1)           ! dxA/dx = dcA/dx * dxA/dcA

    w_A    = c_A_to_w_A(cA)
    dwA_dx = dcA_dx * w_A%dx(1)
    w_B    = 1. - w_A
    dwB_dx = -dwA_dx

    ! u_A  = chemical_potential(x_A)        ! u_A, duA/dcA
    ! duA_dx = dcA_dx * u_A%dx(1)             ! duA/dx = dcA/dx * duA/dcA
                                          ! duA/dx = dcA/dx * dxA/dx * duA/dxA
    rho = density(cA)

    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    diff  = practical_diff(cA)
    ! t_1_c = transference_1_common_ion(cA)
    t_1_mass = transference_1_mass(cA)
    t_2_mass = transference_2_mass(cA)


    ! Flux_(1) = -nu_A_Li * nu_B_TFSI * diff /(Rigc * Temp) * cA * du_A  &
    !          & + t_1_c * i_sol / (z_1_Li * Fconst) + nu_A_Li * cA * vel
    ! Flux of Li^+
    ! 𝐍₁ = -νₐᴸ / MWₐ ρ D ∇ωₐ + t₁⋅𝐢₂ / (z₁ F) + νₐᴸ cₐ 𝐯
    Flux_(1) = -nu_A_Li / MW_A_LiTFSI * rho * diff * dwA_dx + &
              & (t_1_mass * i_sol) / (z_1_Li * Fconst) + nu_A_Li * cA * vel
    ! Flux_(2) = -nu_B_EMI / MW_B_EMITFSI * rho * diff * dwB_dx + &
    !           & t_2_mass / z_2_EMI * i_sol / Fconst + nu_B_EMI * cB * vel
    Flux_(2) = rho * vel

    ! 𝐍₁ = -νₐᴸ / MWₐ ρ D ∇ωₐ + t₁⋅𝐢₂ / (z₁ F) + νₐᴸ cₐ 𝐯
    ! 𝐍₂ = ρ𝐯
    ! print*, Flux_(1)%dx(1)*c_vars_dual(1)%x, Flux_(1)%dx(N+1), &
    ! & Flux_(1)%dx(1)*c_vars_dual(1)%x / abs(Flux_(1)%dx(N+1))
  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! no homogeneous reactions
    ! (1) rxn_1 = 0
    ! (2) rxn   = 0
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    cA     = c_vars_dual(1)
    vel    = c_vars_dual(2)

    Rxn_(1) = 0.0
    Rxn_(2) = 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) ∂(c₁)/∂t = ∂(νₐᴸ ρ ωₐ / MWₐ)/∂t = ∂(νₐᴸ * cA)/∂t
    ! (2) ∂ρ/∂t
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual
    ! type(dual)                            :: cA, vel
    ! type(dual)                            :: rho

    cA     = c_vars_dual(1)
    vel    = c_vars_dual(2)

    rho = density(cA)

    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    Accum_(1) = (nu_A_Li * cA)/delT
    ! Accum_(2) = (nu_B_EMI * cB)/delT

    Accum_(2) = rho/delT

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
    ! (1)   z_1 N_1 = i/F
    ! (2)   ρ𝐯 = MW_1 * i/(z1*F)
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp
    type(dual)                             :: cA, vel, rho

    cA  = c_vars_dual(1)
    vel = c_vars_dual(2)
    ! rho = density(cA)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    BC_WEST_(1) = z_1_Li * flux_temp(1) - applied_current_A/(Fconst)
    ! BC_WEST_(2) = flux_temp(2) - 0.0
    BC_WEST_(2) = flux_temp(2) - MW_Li * applied_current_A/(z_1_Li * Fconst)

    ! print*, BC_WEST_(2)
    ! print*, flux_temp(2)
    ! print*, MW_Li * applied_current_A/(z_1_Li * Fconst), cprev(2,1)*density(cprev(1,1)), cprev(2,2)*density(cprev(1,2))

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   z_1 N_1 = i/F
    ! (2)   ρ𝐯 = MW_1 * i/(z1*F)
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp
    type(dual)                             :: cA, vel, rho

    cA  = c_vars_dual(1)
    vel = c_vars_dual(2)
    ! rho = density(cA)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    BC_EAST_(1) = z_1_Li * flux_temp(1) - applied_current_A/Fconst
    ! BC_EAST_(2) = flux_temp(2) - 0.0
    BC_EAST_(2) = flux_temp(2) - MW_Li * applied_current_A/(z_1_Li * Fconst)


  end function Boundary_EAST

! ******************************************************************************
end module GOV_EQNS



module experiment_mod
  use user_input
  use variables, only: cprev, delX
  ! use echem_mod
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

  experimental_program = 'Restricted_Diffusion'
  ! experimental_program = 'Transference_Polarization'

! *****************************************************************************
  if (experimental_program == 'Restricted_Diffusion') then

    discharge_time = 0.15*3600!5 * 3600
    rest_time      = 5 * 3600
    delT = 1.0

    if (time <= discharge_time) then
      state = 'D'

    else if (time <= (discharge_time + rest_time)) then
      state = 'R'
      applied_current_A = 0.0

    end if

! ******************************************************************************
  else if (experimental_program == 'Transference_Polarization') then
    discharge_time = 1 * 60             ! 60 second polarization
    rest_time      = 1 * 3600           ! 1 hour recovery
    delT = 1e-1

    if (time <= discharge_time) then
      state = 'D'

    else if (time <= (discharge_time + rest_time)) then
      state = 'R'
      applied_current_A = 0.0

    end if
  end if

! ******************************************************************************

  ! if (trim(geometry) == 'Rectangular') then
  !
  !
  ! else if (trim(geometry) == 'Cylindrical') then
  !
  !
  ! else if (trim(geometry) == 'Spherical') then
  !
  !
  ! end if

  if (time > 100 * 3600) then
    print*, 'Time > 1000 hours'
    stop
  else if (time > discharge_time + rest_time) then
    stop
  else if (any(cprev(1,:) < 0.0)) then
    call t_write__Equiv__Li_Bal_Assignment
    call write_to_screen(0)
    print*, 'c0 < 0.0'
    stop
  else if (any(isnan(cprev))) then
    call t_write__Equiv__Li_Bal_Assignment
    call write_to_screen(0)
    print*, 'Dependent Variable is NaN'
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
  ! use echem_mod
  use dnadmod
  use TRANSPORT_MODULE
  use experiment_mod
  use BAND_mod


  implicit none
  integer :: t1, t2, clock_rate, clock_max      ! timing variables
  integer :: it, j                              ! loop indices
  type(dual), dimension(N) :: cdual_

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  ! cdual_ = c_to_dual(cprev(:,1))
  ! print*, cdual_(1)
  ! print*, cdual_(2)
  ! print*, density(cdual_(1))
  ! print*, (density(cdual_(1)) - MW_A_LiTFSI * cdual_(1)) / MW_B_EMITFSI


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

    ! print*, cprev(2,1)*density(cprev(1,1)), cprev(2,NJ)*density(cprev(1,NJ))

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
    cprev(2,j) = 0.0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition
