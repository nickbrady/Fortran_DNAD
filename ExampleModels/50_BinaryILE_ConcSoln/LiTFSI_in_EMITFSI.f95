module user_input
  implicit none

  integer, parameter :: N = 2   ! LiPF6 salt concentration, mass-average velocity
  integer, parameter :: NJ = 102                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 1e3*3600 !3.6d3*1000     ! Number of time steps
  real               :: delT = 1e-1/3                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .True.
  character(len=65)  :: direction = 'WestToEast'

  ! **************************** Physical Constants ****************************
  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant
  ! **************************** Physical Parameters ***************************

  real, parameter :: z_Li = +1.0, z_EMI = +1.0, z_TFSI = -1.0
  real, parameter :: nu_A_Li = 1.0, nu_A_TFSI = 1.0
  real, parameter :: nu_B_EMI = 1.0, nu_B_TFSI = 1.0
  real, parameter :: nu_A = nu_A_Li + nu_A_TFSI
  real, parameter :: nu_B = nu_B_EMI + nu_B_TFSI
  real, parameter :: MW_A_LiTFSI = 287.1                ! MW_LiTFSI
  real, parameter :: MW_B_EMITFSI = 391.31               
  real, parameter :: MW_Li    = 6.941                 ! MW_Li         6.941 g/mol
  real, parameter :: MW_TFSI  = (MW_A_LiTFSI - nu_A_Li*MW_Li)/nu_A_TFSI  
  real, parameter :: MW_EMI   = (MW_B_EMITFSI - nu_B_TFSI*MW_TFSI)/nu_B_EMI

  real, parameter :: t_Li_c  = 0.3
  real, parameter :: t_EMI_c = 1 - t_Li_c

  real :: Diff = 2e-6
  ! real :: kappa = 10e-3     ! 10 mS/cm

  real :: i_app_cm2 = 1e-4*3 !* Fconst / 10.
  ! electrochemical reactions
  real :: s_Li = 0.0, s_EMI = 0.0, s_TFSI = 0.0
  real :: n_rxn = 0.0

  real :: xmax          = 1.0              ! 500 um is 500e-4 cm

  ! Initial Conditions
  real :: c_LiTFSI_bulk = 1.0e-3
  real :: vel_0 = 0.0

  character(len=3) :: state = 'D'               ! Discharge - D, Charge - C, Recovery - R

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module TRANSPORT_MODULE
  use user_input!, only: N, Rigc, Fconst, Temp
  use dnadmod

  implicit none

  ! type(dual)                              :: cA, cB


  interface density
    module procedure density_dual
    module procedure density_real
  end interface

  interface c_A_to_x_A
    module procedure c_A_to_x_A_dual
    module procedure c_A_to_x_A_real
  end interface

  ! practical and fundament diffusion coefficient
  interface diffusion_coef_LiTFSI_in_EMI_TFSI
    module procedure diffusion_coef_LiTFSI_in_EMI_TFSI_dual
    module procedure diffusion_coef_LiTFSI_in_EMI_TFSI_real
  end interface

!   interface conductivity_LiTFSI_in_EMI_TFSI
!     module procedure conductivity_LiTFSI_in_EMI_TFSI_dual
!     module procedure conductivity_LiTFSI_in_EMI_TFSI_real
!   end interface

  ! tranference number with respect to different reference velocities
  ! mass-averaged
  ! common ion
  ! volume-averaged
  ! molar-averaged
  interface transference_Li_mass_avg
    module procedure transference_Li_mass_avg_dual
    module procedure transference_Li_mass_avg_real
  end interface

  interface transference_Li_common_ion
    module procedure transference_Li_common_ion_dual
    module procedure transference_Li_common_ion_real
  end interface

contains

!*******************************************************************************
! Physical Properties that are functions of concentration
!   - mole fraction, mass fraction, chemical potential, density
!*******************************************************************************
! cₐ 
  function c_A_to_x_A_dual(c_A)                  result(x_A)
    type(dual)                :: x_A
    type(dual), intent(in)    :: c_A
    type(dual)                :: rho, c_B

    rho = density(c_A)
    c_B = (rho - MW_A_LiTFSI * c_A) / MW_B_EMITFSI
    x_A = c_A / (c_A + c_B)

  end function c_A_to_x_A_dual

  function c_A_to_x_A_real(c_A)                  result(x_A)
    real                :: x_A
    real, intent(in)    :: c_A
    real                :: rho, c_B

    rho = density(c_A)
    c_B = (rho - MW_A_LiTFSI * c_A) / MW_B_EMITFSI
    x_A = c_A / (c_A + c_B)

  end function c_A_to_x_A_real


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
    type(dual)                :: c

    c = cA*1000
    ! density = 0.3883889166768137*w_A + 1.5139132335788112
    density = 0.0682173430940733*cA + 1.5149898756792253

  end function density_dual

  function density_real(cA)       result(density)
    ! density data is empirically measured
    real                :: density
    real, intent(in)    :: cA
    real                :: c

    c = cA*1000
    ! density = 0.3883889166768137*w_A + 1.5139132335788112
    density = 0.0682173430940733*cA + 1.5149898756792253

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
  ! function fundamental_diff(cA)               result(diff_fund)
  !   type(dual)                :: diff_fund
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: cB, c_Total
  !   type(dual)                :: c1, c2, c3
  !   type(dual)                :: rho

  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

  !   c1      = nu_A_Li * cA
  !   c2      = nu_B_EMI * cB
  !   c3      = nu_A_TFSI * cA + nu_B_TFSI * cB
  !   c_Total = c1 + c2 + c3

  !   diff_fund = z_3_TFSI**2 * c_Total / (nu_A_Li * nu_B_EMI) / &
  !            & (z_3_TFSI**2 * c3 / diff_12 &
  !            & + z_2_EMI**2 * c2 / diff_13 &
  !            & +  z_1_Li**2 * c1 / diff_23)

  ! end function fundamental_diff


  ! function practical_diff(cA)               result(diff_prac)
  !   type(dual)                :: diff_prac
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: cB, c3
  !   type(dual)                :: rho

  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI
  !   c3  = nu_A_TFSI * cA + nu_B_TFSI * cB

  !   ! need to add activity relationship
  !   diff_prac = fundamental_diff(cA) * (c3 / (cA + cB)) * (nu_A_Li)

  ! end function practical_diff


  ! function transference_1_common_ion(cA)      result(t_1_c)
  !   type(dual)                :: t_1_c
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: c1, c2, c3
  !   type(dual)                :: rho
  !   type(dual)                :: cB

  !   ! cA  = c_vars_dual(1)
  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

  !   c1      = nu_A_Li * cA
  !   c2      = nu_B_EMI * cB
  !   c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

  !   t_1_c = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
  !   & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
  !   & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )

  ! end function transference_1_common_ion



  ! function transference_1_molar(cA)      result(t_1_molar)
  !   type(dual)                :: t_1_molar
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: c1, c2, c3, c_Total
  !   type(dual)                :: rho
  !   type(dual)                :: cB
  !   type(dual)                :: t_1_c

  !   ! cA  = c_vars_dual(1)
  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

  !   c1      = nu_A_Li * cA
  !   c2      = nu_B_EMI * cB
  !   c3      = nu_A_TFSI * cA + nu_B_TFSI * cB
  !   c_Total = c1 + c2 + c3

  !   ! t_1_c = transference_1_common_ion(cA)

  !   t_1_molar = (nu_A_Li * nu_B * c3 * t_Li_c - nu_B_EMI * nu_A_TFSI * c1) &
  !             & / (nu_A_Li * nu_B_TFSI * c_Total)

  ! end function transference_1_molar



  ! function Q_volume(cA)      result(t_1_volume)
  !   type(dual)                :: t_1_volume
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: c1, c2, c3
  !   type(dual)                :: rho
  !   type(dual)                :: cB

  !   ! cA  = c_vars_dual(1)
  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

  !   c1      = nu_A_Li * cA
  !   c2      = nu_B_EMI * cB
  !   c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

  !   t_1_volume = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
  !   & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
  !   & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )

  ! end function Q_volume


  function transference_Li_common_ion_dual(cA)      result(t_Li_c)
    type(dual)                :: t_Li_c
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB

    ! cA  = c_vars_dual(1)
    ! rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    ! cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    ! c1      = nu_A_Li * cA
    ! c2      = nu_B_EMI * cB
    ! c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    ! t_Li_c = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
    ! & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
    ! & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )
    t_Li_c = -1.2

  end function transference_Li_common_ion_dual


  function transference_Li_common_ion_real(cA)      result(t_Li_c)
    real                :: t_Li_c
    real, intent(in)    :: cA
    real                :: c1, c2, c3
    real                :: rho
    real                :: cB

    ! cA  = c_vars_dual(1)
    ! rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    ! cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    ! c1      = nu_A_Li * cA
    ! c2      = nu_B_EMI * cB
    ! c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    ! t_Li_c = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) / &
    ! & ( (z_2_EMI / diff_13 - z_3_TFSI / diff_12) &
    ! & + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23 - z_3_TFSI / diff_12) )
    t_Li_c = -1.2

  end function transference_Li_common_ion_real


  function transference_Li_mass_avg_dual(cA)      result(t_1_mass)
    implicit none
    type(dual)                :: t_1_mass
    type(dual), intent(in)    :: cA
    type(dual)                :: c1, c2, c3
    type(dual)                :: rho
    type(dual)                :: cB
    type(dual)                :: t_Li_c

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_Li_c    = transference_Li_common_ion(cA)

    t_1_mass = (c3 * MW_B_EMITFSI * t_Li_c - nu_A_TFSI * nu_B_EMI * cA * MW_EMI) / &
    & (nu_B_TFSI * rho)

  end function transference_Li_mass_avg_dual

  function transference_Li_mass_avg_real(cA)      result(t_1_mass)
    implicit none
    real                :: t_1_mass
    real, intent(in)    :: cA
    real                :: c1, c2, c3
    real                :: rho
    real                :: cB
    real                :: t_Li_c

    rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_Li_c    = transference_Li_common_ion(cA)

    t_1_mass = (c3 * MW_B_EMITFSI * t_Li_c - nu_A_TFSI * nu_B_EMI * cA * MW_EMI) / &
    & (nu_B_TFSI * rho)

  end function transference_Li_mass_avg_real

  ! function transference_2_mass(cA)      result(t_2_mass)
  !   type(dual)                :: t_2_mass
  !   type(dual), intent(in)    :: cA
  !   type(dual)                :: c1, c2, c3
  !   type(dual)                :: rho
  !   type(dual)                :: cB
  !   type(dual)                :: t_2_c

  !   rho = density(cA)    ! rho = MW_A * cA + MW_B * cB
  !   cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

  !   c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

  !   t_2_c    = 1.0 - transference_1_common_ion(cA)

  !   t_2_mass = (c3 * MW_A_LiTFSI * t_2_c - nu_B_TFSI * nu_A_Li * cB * MW_Li) / &
  !   & (nu_A_TFSI * rho)

  ! end function transference_2_mass


  function diffusion_coef_LiTFSI_in_EMI_TFSI_dual(c_LiTFSI)    result(D_salt)
    type(dual), intent(in)    :: c_LiTFSI
    type(dual)                :: D_salt
    type(dual)                :: c
    real                      :: p1, p2, p3, p4
    real                      :: T

    p1 =  1.01e+03
    p2 =  1.01e+00
    p3 = -1.56e+03
    p4 = -4.87e+02

    c = c_LiTFSI * 1000
    T = Temp

    D_salt =  p1*exp(p2*c) * exp(p3/T) * exp(p4/T*c) * 1e-6
  end function diffusion_coef_LiTFSI_in_EMI_TFSI_dual
    
  function diffusion_coef_LiTFSI_in_EMI_TFSI_real(c_LiTFSI)    result(D_salt)
    real, intent(in)    :: c_LiTFSI
    real                :: D_salt
    real                :: c
    real                :: p1, p2, p3, p4
    real                :: T

    p1 =  1.01e+03
    p2 =  1.01e+00
    p3 = -1.56e+03
    p4 = -4.87e+02

    c = c_LiTFSI * 1000
    T = Temp

    D_salt =  p1*exp(p2*c) * exp(p3/T) * exp(p4/T*c) * 1e-6
  end function diffusion_coef_LiTFSI_in_EMI_TFSI_real


  function conductivity_LiTFSI_in_EMI_TFSI_dual(c_LiTFSI)    result(conductivity)
    type(dual), intent(in)    :: c_LiTFSI
    type(dual)                :: conductivity
    type(dual)                :: xA, p_val
    real                      :: p1, p2, p3
    real                      :: T

    p1 = -2.1628449260762674 
    p2 = -1.5938863797561     
    p3 = 0.9475473152610514

    xA = c_A_to_x_A(c_LiTFSI)
    T = Temp

    p_val = p1*xA**2 + p2*xA + p3

    conductivity = 10.**p_val
    conductivity = conductivity/1000

  end function conductivity_LiTFSI_in_EMI_TFSI_dual

  function conductivity_LiTFSI_in_EMI_TFSI_real(c_LiTFSI)    result(conductivity)
    implicit none

    real, intent(in)          :: c_LiTFSI
    real                      :: conductivity
    real                      :: xA, p_val
    real                      :: p1, p2, p3
    real                      :: T

    p1 = -2.1628449260762674 
    p2 = -1.5938863797561     
    p3 = 0.9475473152610514

    xA = c_A_to_x_A(c_LiTFSI)
    T = Temp

    p_val = p1*xA**2 + p2*xA + p3

    conductivity = 10**p_val
    conductivity = conductivity/1000

  end function conductivity_LiTFSI_in_EMI_TFSI_real

end module TRANSPORT_MODULE

!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx, delX
  use TRANSPORT_MODULE

  implicit none

  real :: c_LiTFSI, vel, Peclet_Number

  character(len=65) :: header_fmt_file    = '( 1(A5,1X), 2(A12,1X),   20(A15,1X) )'
  character(len=65) :: data_fmt_file      = '( 1(A5,1X), 2(F12.5,1X), 20(ES15.5E3,1X) )'

  character(len=65) :: header_fmt_screen  = '( 1(A5,1X), 3(A12,1X),   20(A15,1X) )'
  character(len=65) :: data_fmt_screen    = '( 1(A5,1X), 3(F12.5,1X), 20(ES15.5E3,1X) )'


  real :: t_write
  real :: Li_balance
  real :: write_every_x_sec = 3600.0


contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    ! real :: write_every_x_sec = 3600.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it <= 2) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      ! last_write_time = time
    else if ( (time - last_write_time).GE.write_every_x_sec ) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it, j
    real :: mass_balance_e, mass_balance_0, mass_balance_total
    real :: vel_1, vel_NJ, density_1, density_NJ, c_LiTFSI_1, c_LiTFSI_NJ, diff_1, diff_NJ, Phi_2_1, Phi_2_NJ

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt_screen) 'State', 'Time', 'c_A_1', 'c_A_NJ', 'ρv_1', 'ρv_NJ', 'Φ_1', 'Φ_NJ'
      write(*, header_fmt_screen) 'CDR',  'hours', 'mol/L',  'cm/s' , 'cm/s' , 'cm/s',   'V',    'V'
    end if

    c_LiTFSI_1  = cprev(1,1)
    c_LiTFSI_NJ = cprev(1,NJ-1)

    vel_1  = cprev(2,1)
    vel_NJ = cprev(2,NJ-1)

    ! Phi_2_1 = cprev(3,1)
    ! Phi_2_NJ = cprev(3,NJ)

    density_1   = density(c_LiTFSI_1)
    density_NJ  = density(c_LiTFSI_NJ)

    ! diff_1      = diffusion_coef_LiTFSI_in_EC_EMC(c_LiTFSI_1)
    ! diff_NJ     = diffusion_coef_LiTFSI_in_EC_EMC(c_LiTFSI_NJ)

    write(*, data_fmt_screen) state, t_write, c_LiTFSI_1*1e3, c_LiTFSI_NJ*1e3, vel_1*density_1, vel_NJ*density_NJ, &
    & Phi_2_1 - Phi_2_NJ, Phi_2_NJ - Phi_2_NJ

  end subroutine write_to_screen

  subroutine write_all_data_to_terminal
    integer :: j
    real    :: Peclet_Number

    do j = 1, NJ
      c_LiTFSI = cprev(1,j)
      vel     = cprev(2,j)

      Peclet_Number = vel*delX(j)/Diff

      write(*, data_fmt_screen) state, xx(j), c_LiTFSI, vel, Peclet_Number

    end do

  end subroutine write_all_data_to_terminal





  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j,jj, i, interval = 1
    real :: density_, w_a, Phi_2, Phi_2_NJ

    integer, parameter :: interval_max = NJ/100
    integer, dimension(100+2) :: write_node_list


    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
    write(56, header_fmt_file) 'State', 'Time', 'Position', 'c_LiPF6', 'velocity', 'Φ_2', 'Density', 'ω_LiPF6', 'Current'
    write(56, header_fmt_file) 'CDR',  'hours',  'um'     ,  'mol/mL',     'cm/s', 'V',    'g/mL',         '-', 'A/cm2'
                                                    !           !
    end if                                          !           !

    ! if (interval_max > 0) then
    !   interval = interval_max
    ! end if


    ! do i = 1, size(write_node_list)
    !   if (i == 1) then
    !     write_node_list(i) = 1
    !   else if (i < size(write_node_list)) then
    !     write_node_list(i) = (i-1)*interval
    !   else if (i == size(write_node_list)) then
    !     write_node_list(i) = NJ
    !   end if
    ! end do

    ! do jj = 1, size(write_node_list)!1, NJ, interval
    !   j = write_node_list(jj)

    do j = 1, NJ
      c_LiTFSI   = cprev(1,j)
      vel       = cprev(2,j)
      ! Phi_2     = cprev(3,j) - Phi_2_NJ
      density_  = density(c_LiTFSI)
      w_a       = MW_A_LiTFSI * c_LiTFSI / density_

      write(56, data_fmt_file) state,    t_write,   xx(j),  c_LiTFSI, vel, Phi_2, density_, w_a, i_app_cm2

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module experiment_mod
  use user_input
  use variables, only: cprev, delX
  ! use echem_mod
  use write_data_mod
  use TRANSPORT_MODULE
  implicit none

  real    :: c_rate
  real    :: discharge_time
  real    :: discharge_time_0
  real    :: charge_time, rest_time, CCV_time
  real    :: end_time_prev_cyc
  real    :: experiment_time = 100*3600.
  character(len=:), allocatable :: experimental_program

contains

  subroutine experimental_conditions(it)
    

    implicit none

    integer :: it
    integer :: cyc = 1, number_of_cycles

    ! experimental_program = 'Restricted_Diffusion'
    experimental_program = 'Transference_Polarization'

  ! *****************************************************************************
    if (experimental_program == 'Restricted_Diffusion') then

      discharge_time = 30 * 3600
      rest_time      = 20 * 3600

      if (time <= discharge_time) then
        state = 'D'

      else if (time <= (discharge_time + rest_time)) then
        state = 'R'
        i_app_cm2 = 0.0

      end if

  ! ******************************************************************************
    else if (experimental_program == 'Transference_Polarization') then
      discharge_time = 1 * 60             ! 60 second polarization
      rest_time      = 0.5 * 3600           ! 1 hour recovery
      number_of_cycles  = 7
      experiment_time   = number_of_cycles * (discharge_time + rest_time) + 2 * 3600.
      write_every_x_sec = 1.0 ! every second

      if ( ((time - end_time_prev_cyc) >= (rest_time + discharge_time)).AND.(cyc <= number_of_cycles) ) then
        end_time_prev_cyc = time
        cyc = cyc + 1
        call initial_condition()
        ! print*, cyc, cprev(1,NJ), cprev(1,1)
      end if
  
      if ((time - end_time_prev_cyc) <= 60) then  ! current ON
        state = 'D'
        i_app_cm2 = 1.0e-4 * cyc * 3./10
        
      else                                        ! turn current OFF
        state = 'R'
        i_app_cm2 = 0.0
  
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

    if (time > experiment_time) then
      print*, 'Experiment Finished'
      call write_all_data_to_terminal
      print*, diffusion_coef_LiTFSI_in_EMI_TFSI(cprev(1,NJ/2)), transference_Li_mass_avg(cprev(1,NJ/2))
      stop
    ! else if (time > discharge_time) then
    !   stop
    else if (any(cprev(1,:) < 0.0)) then
      call t_write__Equiv__Li_Bal_Assignment
      call write_all_data_to_terminal
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
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module GOV_EQNS
  use user_input
  use dnadmod
  use TRANSPORT_MODULE
  implicit none

  type(dual) :: w_a, dwdx_a
  type(dual) :: w_b, dwdx_b
  type(dual) :: vel, dveldx
  type(dual) :: vel_mass_avg
  type(dual) :: c_LiTFSI, dcdx_LiTFSI
  type(dual) :: density_
contains

! (1) ∂cₐ/∂t = -∇⋅(𝐍ₐ) + R
!     ∂(ρωₐ/MWₐ)/∂t = -∇⋅()
!     where cₐ = ρωₐ/MWₐ
!     𝐍ₐ = -D∇cₐ + tᴸⁱ⋅𝐢₂/(ν₁ᴬ z₊ F) + cₐ𝐯
! (2) ∂ρ/∂t = -∇⋅(ρ𝐯)      ! continuity equation
! (3) 0 = -∇⋅𝐢₂ + R        ! 𝐢₂ = -κ∇Φ₂ - κ/F(tᴸⁱ)/(ν₊z₊) ∇μₑ
!   
!   𝐢₂ = -κ∇Φ₂ - κ/F*[(s₁/(n ν₁ᴬ) + t₁ᶜ/z₁ ν₁ᴬ)(1 + ν₃ᴬ cₐ / (ν₃ᵇ c_B)) - s₃ cₐ / (n c_B ν₃ᵇ)] ∇μₐ
!   μ_A = RT ln (x\_A ^ (\nu_1^A) \gamma_A \lambda_A^0)
!   x_A = c_A / (c_A + c_B)


  ! function dPhi2_dx_function(c_vars_dual, dcdx_vars_dual)  Result(dPhi2_dx)
  !   type(dual)                            :: dPhi2_dx
  !   type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
  !   type(dual)                            :: dmudx_e
  !   type(dual)                            :: w_e, w_0
  !   type(dual)                            :: c_e, c_0
  !   𝐢₂ = -κ∇Φ₂ - κ/F*[(s₁/(n ν₁ᴬ) + t₁ᶜ/z₁ ν₁ᴬ)(1 + ν₃ᴬ cₐ / (ν₃ᵇ c_B)) - s₃ cₐ / (n c_B ν₃ᵇ)] ∇μₐ
  !   μ_A = RT ln (x\_A ^ (\nu_1^A) \gamma_A \lambda_A^0)
  !   x_A = c_A / (c_A + c_B)
  !
  !   dmudx_e = dmudx_e_function(c_vars_dual, dcdx_vars_dual)
  !
  !   w_e = c_vars_dual(1)
  !   w_0 = 1 - w_e
  !   c_e = w_e * density / MW_e
  !   c_0 = w_0 * density / MW_0
  !
  !   ! (s_cat/(n_rxn*nu_cat) - s_0*c_0/(n_rxn * c_e) + t_cat/(nu_cat*z_cat) )
  !   dPhi2_dx = -i_app_cm2/kappa &
  !         & - 1/Fconst*(t_cat/(nu_cat*z_cat)) * dmudx_e
  !
  ! end function dPhi2_dx_function

  ! function dmudx_e_function(c_vars_dual, dcdx_vars_dual)  Result(dmudx_e)
  !   type(dual)                            :: dmudx_e
  !   type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
  !   type(dual)                            :: w_e, w_0
  !   type(dual)                            :: density_0
  !   real :: activity = 1.0
  !   ! ∇μ_e = νRT * (1 - dln(c_0)/dln(c_e)) (1 + dln(γ)/dln(m)) ∇c/c
  !   ! ∇ω_e = ρ_0 MW_e/ρ * (1 - dln(c_0)/dln(c_e))∇c
  !   ! c = ρ ω_e / MW_e
  !   ! ∇μ_e = νRT * (1 + dln(γ)/dln(m)) /c * ρ*∇ω_e/(ρ_0 * MW_e)
  !   ! ∇μ_e = νRT * (1 + dln(γ)/dln(m)) 1/ρ_0 * ∇ω_e/ω_e
  !   ! ρ = ρ_e + ρ_0
  !   ! ρ_0 = c_0 * MW_0 = ρ(1 - ω_e)
  !   w_e = c_vars_dual(1)
  !   w_0 = 1.0 - w_e
  !
  !   dwdx_e = dcdx_vars_dual(1)
  !
  !   density_0 = density * (1 - w_e)
  !
  !   dmudx_e = nu_e * Rigc*Temp * (activity)/density_0 * dwdx_e/w_e
  !
  ! end function dmudx_e_function

  ! function solution_current_i2_function(c_vars_dual, dcdx_vars_dual)  Result(i2)
  !   type(dual)                            :: i2
  !   type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
  !   type(dual)                            :: c_LiTFSI, kappa, t_Li
  !   type(dual)                            :: dPhi2_dx
  !   type(dual)                            :: mu_e, dmudx_e
  !   real                                  :: dmu_dc, t_Li_1

  !   ! i = -κ∇Φ_2 - κ/F*(s_1/(n ν_1^A + t_1^c/(z_1 v_1^A)) - s_0 c_0 / (n c_e) + t^0_+/(ν_+ z_+)) ∇μ_e
  !   c_LiTFSI = c_vars_dual(1)
  !   ! kappa    = conductivity_LiPF6_in_EC_EMC(c_LiPF6)
  !   t_Li     = transference_1_mass_avg(c_LiTFSI)

  !   mu_a    = chemical_potential_LiPF6(c_LiTFSI)
  !   dmu_dc  = mu_e%dx(1)
  !   dmudx_a = dmu_dc * dcdx_vars_dual(1) ! dμ/dx = dμ/dc * dc/dx

  !   ! dPhi2_dx = dcdx_vars_dual(3)

  !   ! i2 = -kappa*dPhi2_dx - kappa/Fconst*(t_Li)/(nu_cat*z_cat) * dmudx_e

  !   ! i2 = -kappa*dPhi2_dx

  !   ! print*, kappa%x

  ! end function solution_current_i2_function

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual)  result(Flux_)
    ! (1)   N_e = N_+/ν_+ = -ρ/MW_e D dω_e/dx + ρ ω_e/MW_e * v
    ! (2)   massFlux = ρv
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual) :: Diff_, t_Li_mass, t_Li
    type(dual) :: i2_sol_curr  ! solution current

    c_LiTFSI     = c_vars_dual(1)
    dcdx_LiTFSI  = dcdx_vars_dual(1)
    density_    = density(c_LiTFSI)
    w_a         = MW_A_LiTFSI * c_LiTFSI / density_
    w_b         = 1.0 - w_a
    dwdx_a      = MW_A_LiTFSI * dcdx_LiTFSI / density_
    dwdx_b      = -dwdx_a

    vel_mass_avg = c_vars_dual(2)

    Diff_ = diffusion_coef_LiTFSI_in_EMI_TFSI(c_LiTFSI)
    t_Li_mass  = transference_Li_mass_avg(c_LiTFSI)

    ! i2_sol_curr = solution_current_i2_function(c_vars_dual, dcdx_vars_dual)

    ! FLUX_(1) = -Diff_pract * dcdx_LiTFSI + t_Li_mass * i2_sol_curr/(nu_A_Li*z_cat*Fconst) + c_LiTFSI * vel_mass_avg

    FLUX_(1) = -nu_A_Li*Diff_ * dcdx_LiTFSI + t_Li_mass * i_app_cm2/(nu_A_Li*z_Li*Fconst) + nu_A_Li*c_LiTFSI * vel_mass_avg
    FLUX_(2) =  density_*vel_mass_avg
    ! FLUX_(3) =  i2_sol_curr

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! no homogeneous reactions
    ! (1) rxn_1 = 0
    ! (2) rxn   = 0
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual
    type(dual) :: c_e, c_cat, c_an

    Rxn_(1) = 0.0
    Rxn_(2) = 0.0
    ! Rxn_(3) = 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) d(c_LiPF6)/dt = d(ρω_e6/MW_e)/dt
    ! (2) dρ/dt
    ! (3) 0   no accumulation of charge
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c_LiTFSI    = c_vars_dual(1)
    density_    = density(c_LiTFSI)
    w_a         = MW_A_LiTFSI * c_LiTFSI / density_
    w_b         = 1.0 - w_a

    ! Accum_(1) = density_/MW_A_LiTFSI * w_a/delT
    Accum_(1) = nu_A_Li * c_LiTFSI/delT
    Accum_(2) = density_/delT
    ! Accum_(3) = 0.0

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
    ! (2)   ρv = M_1 * i/(z1*F) = MW_e * N_e
    ! (3)   Φ_2 = 0   /   i_2 = i_app_cm2
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp, rxn_temp
    type(dual)                             :: w_e, dwdx_e
    type(dual)                             :: w_0, dwdx_0
    type(dual)                             :: Phi_2

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    ! Phi_2 = c_vars_dual(3)

    BC_WEST_(1) = flux_temp(1) - i_app_cm2/(z_Li*nu_A_Li*Fconst)
    BC_WEST_(2) = flux_temp(2) - flux_temp(1)*MW_Li
    ! BC_WEST_(3) = flux_temp(3) - i_app_cm2


  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   z_1 N_1 = i/F
    ! (2)   ρv = M_1 * i/(z1*F) = MW_e * N_e
    ! (3)   i_2 = i_app_cm2
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp, rxn_temp
    type(dual)                             :: w_e, dwdx_e
    type(dual)                             :: w_0, dwdx_0
    type(dual)                             :: i2, Phi_2
    type(dual)                             :: d_density_dx

    ! Phi_2 = c_vars_dual(3)

    vel          = c_vars_dual(2)
    dveldx       = dcdx_vars_dual(2)

    c_LiTFSI     = c_vars_dual(1)
    dcdx_LiTFSI  = dcdx_vars_dual(1)

    density_     = density(c_LiTFSI)
    d_density_dx = density_%dx(1)*dcdx_LiTFSI  ! ∇ρ = ∂ρ/∂c ⋅ ∇c

    flux_temp    = FLUX(c_vars_dual, dcdx_vars_dual)

    BC_EAST_(1)  = flux_temp(1) - i_app_cm2/(z_Li*nu_A_Li*Fconst)
    ! ∇⋅(ρ𝐯) = 0 = 𝐯⋅∇ρ + ρ∇𝐯
    BC_EAST_(2)  = d_density_dx*vel + density_*dveldx - 0.0
    ! BC_EAST_(3) = Phi_2 - 0.0

  end function Boundary_EAST
! ******************************************************************************
end module GOV_EQNS
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*


! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use experiment_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max      ! timing variables
  integer :: it, j                              ! loop indices
  type(dual), dimension(N) :: cdual_

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

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

  implicit none
  real    :: h
  integer :: j
  real    :: delX_init, delX_max

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-2)
  
  ! linear mesh
  delx(1) = 0.0               ! it is common in the finite volume
  delx(2:NJ-1) = h            ! algorithm to set the control
  delx(NJ) = 0.0              ! volumes at the boundaries to zero
  ! do j = 1, NJ
  !   print*, j, delX(j)
  ! end do
  
  ! calculate exponential / logarithmic mesh using secant-method
  delX_max = h*4
  ! use small values of x0 and x1 (1e-11, 1e-10)
  ! plotting the function reveals that convergence is more likely when
  ! the initial guess values are too small as opposed to too big
                            !  x0,    x1,             tol, max_it
  delX_init = secant_method(1e-11, 1e-10, delX_max, 1e-15, 100)
  delX(:NJ/2) = delX_mesh(delX_init, delX_max, NJ/2, xmax/2)
  delX(NJ/2) = xmax/2.0 - sum(delX(:NJ/2-1))
  delX(NJ/2+1:NJ-1) = delX(NJ/2:2:-1)

  ! xx positions are calculated from delX
  xx(1) = 0.0
  do j = 2, NJ
    xx(j) = xx(j-1) + (delX(j) + delX(j-1))/2
  end do

  do j = 1, NJ
    cprev(1,j) = c_LiTFSI_bulk
    cprev(2,j) = vel_0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?



contains

  function delX_mesh(initial_delX, delX_max, mesh_pts, xmax)    result(delX)
    ! defines an exponential grid mesh
    implicit none
    integer, intent(in) :: mesh_pts
    real, intent(in)    :: initial_delX, xmax, delX_max
    real, dimension(mesh_pts) :: delX
    integer             :: j
    real                :: h_log
    
    ! delX_max = 4 * xmax/(mesh_pts-1)

    h_log = (log10(delX_max) - log10(initial_delX)) / float(mesh_pts-2)

    do j = 1, mesh_pts
      delX(j+1) = log10(initial_delX) + h_log * (j-1)
    end do
    delX(1) = 0.0
    delX(2:) = 10**(delX(2:))
    
    return 
  end function delX_mesh

  function secant_method(x0_init, x1_init, delX_max, tolerance, max_iterations)    result(x2)
    ! """Return the root calculated using the secant method."""
    ! finds the exponential mesh size that fits the system boundaries
    use user_input, only: NJ, xmax
    implicit none
    integer, intent(in) :: max_iterations
    real, intent(in)    :: tolerance
    real, intent(in)    :: x0_init, x1_init, delX_max
    real                :: x0, x1, x2
    integer             :: i
    real                :: x1_x0_diff, f_x0, f_x1
    
    x0 = x0_init
    x1 = x1_init

    i = 0
    x1_x0_diff = 1e30

    print*, 'secant function'

    do while ((abs(x1_x0_diff) > tolerance) .AND. (i < max_iterations))
      f_x0 = sum(delX_mesh(x0, delX_max, NJ/2, xmax/2)) - xmax/2
      f_x1 = sum(delX_mesh(x1, delX_max, NJ/2, xmax/2)) - xmax/2
      x1_x0_diff = f_x1 - f_x0
      if (x1_x0_diff == 0.0) then
        x2 = x1
        return
      end if

      print*, i, x0, f_x0, x1, f_x1

      x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
      x0 = x1
      x1 = x2

      i = i + 1
    end do

    return

  end function secant_method


end subroutine initial_condition

include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_auto_fill.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_ABDGXY.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_MATINV.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_BAND.f95'
