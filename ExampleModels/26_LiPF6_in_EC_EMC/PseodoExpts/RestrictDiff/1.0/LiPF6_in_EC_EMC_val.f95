module user_input
  implicit none

  integer, parameter :: N = 3   ! LiPF6 salt concentration, mass-average velocity
  integer, parameter :: NJ = 42                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 1e3*3600 !3.6d3*1000     ! Number of time steps
  real               :: delT = 1e0                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .True.
  character(len=65)  :: direction = 'WestToEast'

  ! **************************** Physical Constants ****************************
  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant
  ! **************************** Physical Parameters ***************************

  real, parameter :: z_cat  = +1.0
  real, parameter :: z_an   = -1.0
  real, parameter :: nu_cat = 1.0
  real, parameter :: nu_an  = 1.0
  real, parameter :: nu_e   = nu_cat + nu_an
  real, parameter :: MW_cat = 6.941 ! Li^+          ! MW_Li         6.941 g/mol
  real, parameter :: MW_e   = 151.905               ! MW LiPF6    151.906 g/mol
  real, parameter :: MW_an  = MW_e - nu_cat*MW_cat  ! PF6^- !
  ! real, parameter :: MW_e   = nu_cat * MW_cat + nu_an * MW_an ! 151.905

  real, parameter :: t_cat = 0.3
  real, parameter :: t_an  = 1 - t_cat

  real, parameter :: MW_EC = 88.06    ! g/mol - ethylene carbonate
  real, parameter :: MW_EMC = 104.104 ! g/mol - ethyl methyl carbonate
  real, parameter :: mols_100g = 30 / MW_EC + 70 / MW_EMC ! EC:EMC (3:7)
  real, parameter :: MW_0 = 100 / mols_100g
  real :: Diff = 1e-6
  ! real :: kappa = 10e-3     ! 10 mS/cm

  real :: i_app_cm2 = 1e-4*3 !* Fconst / 10.
  ! electrochemical reactions
  real :: s_cat = 0.0, s_an = 0.0, s_0 = 0.0
  real :: n_rxn = 0.0

  real :: xmax          = 1.0              ! 500 um is 500e-4 cm

  ! Initial Conditions
  real :: c_LiPF6_bulk = 1.0e-3
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

  interface density
    module procedure density_dual
    module procedure density_real
  end interface

  interface diffusion_coef_LiPF6_in_EC_EMC
    module procedure diffusion_coef_LiPF6_in_EC_EMC_dual
    module procedure diffusion_coef_LiPF6_in_EC_EMC_real
  end interface

  interface conductivity_LiPF6_in_EC_EMC
    module procedure conductivity_LiPF6_in_EC_EMC_dual
    module procedure conductivity_LiPF6_in_EC_EMC_real
  end interface

  interface transference_number
    module procedure transference_number_dual
    module procedure transference_number_real
  end interface

  interface transference_number_mass_avg
    module procedure transference_number_mass_avg_dual
    module procedure transference_number_mass_avg_real
  end interface

contains

  ! functional forms and values of density, diffusion coefficient, transference number, and conductivity come from
  ! Temperature and Concentration Dependence of the Ionic TransportProperties of Lithium-Ion Battery Electrolytes
  ! Johannes Landesfeind and Hubert A. Gasteiger 2019 J. Electrochem. Soc. 166 A3079

  ! From supporting information
  ! density of LiPF6 in EC:EMC(3:7)
  ! equation
  ! ρ(c,T) = ρ1 + ρ2*c + ρ3*T       [g/mL]
  !
  ! conc [=] mol/cm3
  ! Temp [=] Kelvin
  function density_dual(c_LiPF6)     result(density_)
      type(dual), intent(in)    :: c_LiPF6
      type(dual)                :: density_
      real                      :: p1, p2, p3

      p1 =  1.41e+00
      p2 =  9.12e-02
      p3 = -1.09e-03

      density_ = p1 + p2*(c_LiPF6*1000) + p3*(Temp - 273)

  end function density_dual

  function density_real(c_LiPF6)     result(density_)
      real, intent(in)    :: c_LiPF6
      real                :: density_
      real                :: p1, p2, p3

      p1 =  1.41e+00
      p2 =  9.12e-02
      p3 = -1.09e-03

      density_ = p1 + p2*(c_LiPF6*1000) + p3*(Temp - 273)

  end function density_real

  ! Equation 15
  ! κ(c,T) = p1(1 + (T-p2)) * c * (1 + p3 √(c) + p4*(1 + p5*exp(1000/T))*c ) / (1 + c^4*(p6*exp(1000/T)))
  !
  ! κ    [=] mS/cm
  ! conc [=] mol/L (M)
  ! Temp [=] Kelvin
  function conductivity_LiPF6_in_EC_EMC_dual(c_LiPF6)    result(conductivity)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: conductivity
    type(dual)                :: c
    real                      :: p1, p2, p3, p4, p5, p6
    real                      :: T
    p1 =  5.21e-01
    p2 =  2.28e+02
    p3 = -1.06e+00
    p4 =  3.53e-01
    p5 = -3.59e-03
    p6 =  1.48e-03

    c = c_LiPF6 * 1000
    T = Temp

    conductivity = p1*(1 + (T-p2))*c*(1 + p3*sqrt(c) + p4*(1 + p5*exp(1000/T))*c) / (1 + c**4 * p6*exp(1000/T))
    conductivity = conductivity/1000

  end function conductivity_LiPF6_in_EC_EMC_dual

  function conductivity_LiPF6_in_EC_EMC_real(c_LiPF6)    result(conductivity)
    real, intent(in)          :: c_LiPF6
    real                      :: conductivity
    real                      :: c
    real                      :: p1, p2, p3, p4, p5, p6
    real                      :: T
    p1 =  5.21e-01
    p2 =  2.28e+02
    p3 = -1.06e+00
    p4 =  3.53e-01
    p5 = -3.59e-03
    p6 =  1.48e-03

    c = c_LiPF6 * 1000
    T = Temp

    conductivity = p1*(1 + (T-p2))*c*(1 + p3*sqrt(c) + p4*(1 + p5*exp(1000/T))*c) / (1 + c**4 * p6*exp(1000/T))
    conductivity = conductivity/1000

  end function conductivity_LiPF6_in_EC_EMC_real


  ! Equation 18
  ! D_±(c,T) = p1*exp(p2*c) * exp(p3/T) * exp(p4/T*c) * 1e-6
  !
  ! D_±  [=] cm2/s
  ! conc [=] mol/L (M)
  ! Temp [=] Kelvin
  function diffusion_coef_LiPF6_in_EC_EMC_dual(c_LiPF6)    result(D_salt)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: D_salt
    type(dual)                :: c
    real                      :: p1, p2, p3, p4
    real                      :: T

    p1 =  1.01e+03
    p2 =  1.01e+00
    p3 = -1.56e+03
    p4 = -4.87e+02

    c = c_LiPF6 * 1000
    T = Temp

    D_salt =  p1*exp(p2*c) * exp(p3/T) * exp(p4/T*c) * 1e-6
  end function diffusion_coef_LiPF6_in_EC_EMC_dual

  function diffusion_coef_LiPF6_in_EC_EMC_real(c_LiPF6)    result(D_salt)
    real, intent(in)    :: c_LiPF6
    real                :: D_salt
    real                :: c
    real                :: p1, p2, p3, p4
    real                :: T

    p1 =  1.01e+03
    p2 =  1.01e+00
    p3 = -1.56e+03
    p4 = -4.87e+02

    c = c_LiPF6 * 1000
    T = Temp

    D_salt =  p1*exp(p2*c) * exp(p3/T) * exp(p4/T*c) * 1e-6
  end function diffusion_coef_LiPF6_in_EC_EMC_real


  ! Equation 20
  ! t_+(c,T) = p1 + p2*c + p3*T + p4*c^2 + p5*c*T + p6*T^2 + p7*c^3 + p8*c^2*T + p9*c*T^2
  !
  ! conc [=] mol/L (M)
  ! Temp [=] Kelvin
  function transference_number_dual(c_LiPF6)           result(t_Li)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: t_Li
    type(dual)                :: c
    real                      :: p1, p2, p3, p4, p5, p6, p7, p8, p9
    real                      :: T

    p1 = -1.28e+01
    p2 = -6.12e+00
    p3 =  8.21e-02
    p4 =  9.04e-01
    p5 =  3.18e-02
    p6 = -1.27e-04
    p7 =  1.75e-02
    p8 = -3.12e-03
    p9 = -3.96e-05

    c = c_LiPF6 * 1000
    T = Temp

    t_Li = p1 + p2*c + p3*T + p4*c**2 + p5*c*T + p6*T**2 + p7*c**3 + p8*c**2*T + p9*c*T**2
  end function transference_number_dual

  function transference_number_real(c_LiPF6)           result(t_Li)
    real, intent(in)          :: c_LiPF6
    real                      :: t_Li
    real                      :: c
    real                      :: p1, p2, p3, p4, p5, p6, p7, p8, p9
    real                      :: T

    p1 = -1.28e+01
    p2 = -6.12e+00
    p3 =  8.21e-02
    p4 =  9.04e-01
    p5 =  3.18e-02
    p6 = -1.27e-04
    p7 =  1.75e-02
    p8 = -3.12e-03
    p9 = -3.96e-05

    c = c_LiPF6 * 1000
    T = Temp

    t_Li = p1 + p2*c + p3*T + p4*c**2 + p5*c*T + p6*T**2 + p7*c**3 + p8*c**2*T + p9*c*T**2
  end function transference_number_real


  ! t_Li_mass = (ρ_0 * t_Li - ρ_-)/ρ
  function transference_number_mass_avg_dual(c_LiPF6)           result(t_Li_mass)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: t_Li_0, t_Li_mass
    type(dual)                :: density_, density_0, density_an
    type(dual)                :: w_LiPF6, w_PF6, w_0

    density_ = density(c_LiPF6)

    w_LiPF6 = c_LiPF6*MW_e/density_
    w_0 = 1 - w_LiPF6
    w_PF6 = w_LiPF6 * MW_an/MW_e

    density_0 = density_ * w_0
    density_an = density_ * w_PF6

    ! print*, density_%x, density_0%x, density_an%x

    t_Li_0 = transference_number(c_LiPF6)
    t_Li_mass = (density_0 * t_Li_0 - density_an) / density_

  end function transference_number_mass_avg_dual

  function transference_number_mass_avg_real(c_LiPF6)           result(t_Li_mass)
    real, intent(in)    :: c_LiPF6
    real                :: t_Li_0, t_Li_mass
    real                :: density_, density_0, density_an
    real                :: w_LiPF6, w_PF6, w_0

    density_ = density(c_LiPF6)

    w_LiPF6 = c_LiPF6*MW_e/density_
    w_0 = 1 - w_LiPF6
    w_PF6 = w_LiPF6 * MW_an/MW_e

    density_0 = density_ * w_0
    density_an = density_ * w_PF6

    t_Li_0 = transference_number(c_LiPF6)
    t_Li_mass = (density_0 * t_Li_0 - density_an) / density_

  end function transference_number_mass_avg_real


  function ln_activity_LiPF6(c_LiPF6)      result(ln_a)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: ln_a
    type(dual)                :: c
    real                      :: p1, p2, p3, p4, p5, p6, p7!, p8, p9
    ! conc [=] mol/L (M)
    ! ln_a = p1*c**5 + p2*c**4 + p3*c**3 + p4*c**2 + p5*c + p6
    ! p1 =   0.4586279805977089
    ! p2 = - 4.040221139218544
    ! p3 =  14.030477222226859
    ! p4 = -19.687568269108052
    ! p5 =  17.974512970548624
    ! p6 = - 0.1228907319593211

    p1 =   0.2229543091079215
    p2 = - 2.1060293253384716
    p3 =   6.257954105092719
    p4 = - 4.586569226614202
    p5 = - 4.1364565292657645
    p6 =  12.533674173344446
    p7 =   0.047808099204433

    c = c_LiPF6 * 1000

    ! ln_a = p1*c**5 + p2*c**4 + p3*c**3 + p4*c**2 + p5*c**1 + p6
    ln_a = p1*c**6 + p2*c**5 + p3*c**4 + p4*c**3 + p5*c**2 + p6*c**1 + p7*c**0
    ! ln_a = ln_a / 1000.
    ! print*, ln_a%x, c_LiPF6%x*1e3
    ! ln_a = ln_a/1000 ! mV --> V

  end function ln_activity_LiPF6


  function chemical_potential_LiPF6(c_LiPF6)           result(mu_LiPF6)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: mu_LiPF6
    type(dual)                :: ln_a
    ! μ_e = ν_e RT ln(a)
    !
    ! conc [=] mol/L (M)
    ! Temp [=] Kelvin
    ln_a     = ln_activity_LiPF6(c_LiPF6)
    mu_LiPF6 = nu_e*Rigc*Temp * ln_a

    ! print*, mu_LiPF6%x/Fconst, c_LiPF6%x*1e3

  end function chemical_potential_LiPF6

  ! dμ/dx = dμ/dc * dc/dx


  function dchem_pot_dc(c_LiPF6)           result(dmu_dc)
    type(dual), intent(in)    :: c_LiPF6
    type(dual)                :: dmu_dc
    type(dual)                :: c
    type(dual)                :: ln_a, a_
    real                      :: p1, p2, p3, p4, p5, p6, p7!, p8, p9
    ! μ_e = ν_e RT ln(a)
    ! dμ_e/dc = ν_e RT / a * da/dc
    !
    ! conc [=] mol/L (M)
    ! Temp [=] Kelvin
    !
    ! ln_a = p1*c**5 + p2*c**4 + p3*c**3 + p4*c**2 + p5*c + p6
    ! p1 =    5.466845528724684
    ! p2 = - 48.15943597948498
    ! p3 =  167.24328848894388
    ! p4 = -234.67581376776752
    ! p5 =  214.25619460893918
    ! p6 = -  1.464857524954459

    p1 =   0.2229543091079215
    p2 = - 2.1060293253384716
    p3 =   6.257954105092719
    p4 = - 4.586569226614202
    p5 = - 4.1364565292657645
    p6 =  12.533674173344446
    p7 =   0.047808099204433

    c = c_LiPF6 * 1000

    ln_a = p1*c**6 + p2*c**5 + p3*c**4 + p4*c**3 + p5*c**2 + p6*c**1 + p7*c**0
    a_   = exp(ln_a)
    dmu_dc = nu_e*Rigc*Temp / a_

  end function dchem_pot_dc

end module TRANSPORT_MODULE
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx, delX
  use TRANSPORT_MODULE

  implicit none

  real :: c_LiPF6, vel, Peclet_Number

  character(len=65) :: header_fmt_file    = '( 1(A5,1X), 2(A12,1X),   20(A15,1X) )'
  character(len=65) :: data_fmt_file      = '( 1(A5,1X), 2(F12.5,1X), 20(ES15.5E3,1X) )'

  character(len=65) :: header_fmt_screen  = '( 1(A5,1X), 3(A12,1X),   20(A15,1X) )'
  character(len=65) :: data_fmt_screen    = '( 1(A5,1X), 3(F12.5,1X), 20(ES15.5E3,1X) )'


  real :: t_write
  real :: Li_balance


contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it <= 10*60) then
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
    real :: vel_1, vel_NJ, density_1, density_NJ, c_LiPF6_1, c_LiPF6_NJ, diff_1, diff_NJ, Phi_2_1, Phi_2_NJ

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt_screen) 'State', 'Time', 'c_A_1', 'c_A_NJ', 'ρv_1', 'ρv_NJ', 'Φ_1', 'Φ_NJ'
      write(*, header_fmt_screen) 'CDR',  'hours', 'mol/L',  'cm/s' , 'cm/s' , 'cm/s',   'V',    'V'
    end if

    c_LiPF6_1  = cprev(1,1)
    c_LiPF6_NJ = cprev(1,NJ-1)

    vel_1  = cprev(2,1)
    vel_NJ = cprev(2,NJ-1)

    Phi_2_1 = cprev(3,1)
    Phi_2_NJ = cprev(3,NJ)

    density_1   = density(c_LiPF6_1)
    density_NJ  = density(c_LiPF6_NJ)

    diff_1      = diffusion_coef_LiPF6_in_EC_EMC(c_LiPF6_1)
    diff_NJ     = diffusion_coef_LiPF6_in_EC_EMC(c_LiPF6_NJ)

    write(*, data_fmt_screen) state, t_write, c_LiPF6_1*1e3, c_LiPF6_NJ*1e3, vel_1*density_1, vel_NJ*density_NJ, &
    & Phi_2_1 - Phi_2_NJ, Phi_2_NJ - Phi_2_NJ

  end subroutine write_to_screen

  subroutine write_all_data_to_terminal
    integer :: j
    real    :: Peclet_Number

    do j = 1, NJ
      c_LiPF6 = cprev(1,j)
      vel     = cprev(2,j)

      Peclet_Number = vel*delX(j)/Diff

      write(*, data_fmt_screen) state, xx(j), c_LiPF6, vel, Peclet_Number

    end do

  end subroutine write_all_data_to_terminal





  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j
    real :: density_, w_e, Phi_2, Phi_2_NJ


    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
    write(56, header_fmt_file) 'State', 'Time', 'Position', 'c_LiPF6', 'velocity', 'Φ_2', 'Density', 'ω_LiPF6'
    write(56, header_fmt_file) 'CDR',  'hours',  'um'     ,  'mol/mL',     'cm/s', 'V',    'g/mL',       '-'
                                                    !           !
    end if                                          !           !

    Phi_2_NJ = cprev(3,NJ)
    do j = 1, NJ
      ! if (j == NJ) then
        ! dwdx_e = (cprev(1,NJ) - cprev(1,NJ-1))/(xx(NJ) - xx(NJ-1))
        ! vel = cprev(2,NJ-1)
      ! else
      !   dwdx_e = (cprev(1,j+1) - cprev(1,j))/(xx(j+1) - xx(j))

      ! end if

      c_LiPF6   = cprev(1,j)
      vel       = cprev(2,j)
      Phi_2     = cprev(3,j) - Phi_2_NJ
      density_  = density(c_LiPF6)
      w_e       = MW_e * c_LiPF6 / density_

      ! flux_e = -density/MW_e * Diff * dwdx_e + t_cat * i_app_cm2/(z_cat*Fconst) + density*w_e/MW_e * vel
      ! flux_0 = -density/MW_0 * Diff * dwdx_0 + density*w_0/MW_0 * vel + &
      !         & -1/MW_0 * (MW_cat/z_cat*t_cat + MW_an/z_an*t_an) * i_app_cm2/Fconst

      write(56, data_fmt_file) state,    t_write,   xx(j),  c_LiPF6, vel, Phi_2, density_, w_e!, flux_e, flux_0

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

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

    if (time <= discharge_time) then
      state = 'D'

    else if (time <= (discharge_time + rest_time)) then
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

  if (time > discharge_time + rest_time) then
    print*, 'Experiment Finished'
    call write_all_data_to_terminal
    stop
  ! else if (time > discharge_time) then
  !   stop
  else if (any(cprev(1,1:NJ-1) < 0.0)) then
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

  type(dual) :: w_e, dwdx_e
  type(dual) :: w_0, dwdx_0
  type(dual) :: vel, dveldx
  type(dual) :: vel_mass_avg
  type(dual) :: c_LiPF6, dcdx_LiPF6
  type(dual) :: density_
contains

! (1) ∂c/∂t = -∇⋅(𝐍) + R
!     ∂(ρωₑ/MWₑ)/∂t = -∇⋅()
!     where cₑ = ρωₑ/MWₑ
!     𝐍 = -D∇cₑ + tᴸⁱ⋅𝐢₂/(z₊F) + cₑ𝐯
! (2) ∂ρ/∂t = -∇⋅(ρ𝐯)      ! continuity equation
! (3) 0 = -∇⋅𝐢₂ + R        ! 𝐢₂ = -κ∇Φ₂ - κ/F(tᴸⁱ)/(ν₊z₊) ∇μₑ

  ! function dPhi2_dx_function(c_vars_dual, dcdx_vars_dual)  Result(dPhi2_dx)
  !   type(dual)                            :: dPhi2_dx
  !   type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
  !   type(dual)                            :: dmudx_e
  !   type(dual)                            :: w_e, w_0
  !   type(dual)                            :: c_e, c_0
  !   ! i = -κ∇Φ_2 - κ/F*(s_+/(n ν_+) - s_0 c_0 / (n c_e) + t_+/(ν_+ z_+)) ∇μ_e
  !   ! ∇Φ_2 = -i/κ - 1/F*(s_+/(n ν_+) - s_0 c_0 / (n c_e) + t_+/(ν_+ z_+)) ∇μ_e
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

  function solution_current_i2_function(c_vars_dual, dcdx_vars_dual)  Result(i2)
    type(dual)                            :: i2
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: c_LiPF6, kappa, t_Li
    type(dual)                            :: dPhi2_dx
    type(dual)                            :: mu_e, dmudx_e
    real                                  :: dmu_dc, t_Li_1

    ! i = -κ∇Φ_2 - κ/F*(s_+/(n ν_+) - s_0 c_0 / (n c_e) + t^0_+/(ν_+ z_+)) ∇μ_e
    c_LiPF6 = c_vars_dual(1)
    kappa   = conductivity_LiPF6_in_EC_EMC(c_LiPF6)
    t_Li    = transference_number(c_LiPF6)

    mu_e    = chemical_potential_LiPF6(c_LiPF6)
    dmu_dc  = mu_e%dx(1)
    dmudx_e = dmu_dc * dcdx_vars_dual(1) ! dμ/dx = dμ/dc * dc/dx

    dPhi2_dx = dcdx_vars_dual(3)

    i2 = -kappa*dPhi2_dx - kappa/Fconst*(t_Li)/(nu_cat*z_cat) * dmudx_e

    ! i2 = -kappa*dPhi2_dx

    ! print*, kappa%x

  end function solution_current_i2_function

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
    type(dual) :: Diff_LiPF6, t_Li_mass, t_Li
    type(dual) :: i2_sol_curr  ! solution current

    c_LiPF6     = c_vars_dual(1)
    dcdx_LiPF6  = dcdx_vars_dual(1)
    density_    = density(c_LiPF6)
    w_e         = MW_e * c_LiPF6 / density_
    w_0         = 1.0 - w_e
    dwdx_e      = MW_e * dcdx_LiPF6 / density_
    dwdx_0      = -dwdx_e

    vel_mass_avg = c_vars_dual(2)

    Diff_LiPF6 = diffusion_coef_LiPF6_in_EC_EMC(c_LiPF6)
    t_Li_mass  = transference_number_mass_avg(c_LiPF6)
    t_Li       = transference_number(c_LiPF6)
    ! print*, t_Li_mass%x, t_Li%x

    i2_sol_curr = solution_current_i2_function(c_vars_dual, dcdx_vars_dual)

    ! FLUX_(1) = -density_/MW_e * Diff_LiPF6 * dwdx_e + t_Li_mass*i2_sol_curr/(z_cat*Fconst) + w_e*density_/MW_e * vel_mass_avg
    ! FLUX_(1) = -density_/MW_e * Diff_LiPF6 * dwdx_e + t_Li_mass*i_app_cm2/(z_cat*Fconst) + w_e*density_/MW_e * vel_mass_avg
    ! FLUX_(1) = -Diff_LiPF6 * dcdx_LiPF6 + 0.2*i_app_cm2/(z_cat*Fconst) + c_LiPF6 * vel_mass_avg
    ! FLUX_(1) = -density_/MW_e * Diff_LiPF6 * dwdx_e + w_e*density_/MW_e * vel_mass_avg
    ! FLUX_(1) = -Diff_LiPF6 * dcdx_LiPF6 + c_LiPF6 * vel_mass_avg
    ! FLUX_(1) = -Diff_LiPF6 * dcdx_LiPF6 + i_app_cm2/(z_cat*Fconst) + c_LiPF6 * vel_mass_avg
    ! FLUX_(1) = -Diff_LiPF6 * dcdx_LiPF6 + 0.2 * i2_sol_curr/(z_cat*Fconst) + c_LiPF6 * vel_mass_avg
    FLUX_(1) = -Diff_LiPF6 * dcdx_LiPF6 + t_Li_mass * i2_sol_curr/(z_cat*Fconst) + c_LiPF6 * vel_mass_avg
    FLUX_(2) =  density_*vel_mass_avg
    FLUX_(3) =  i2_sol_curr

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
    Rxn_(3) = 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) d(c_LiPF6)/dt = d(ρω_e6/MW_e)/dt
    ! (2) dρ/dt
    ! (3) 0   no accumulation of charge
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c_LiPF6     = c_vars_dual(1)
    density_    = density(c_LiPF6)
    w_e         = MW_e * c_LiPF6 / density_
    w_0         = 1.0 - w_e

    Accum_(1) = density_/MW_e*w_e/delT
    Accum_(2) = density_/delT
    Accum_(3) = 0.0

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
    Phi_2 = c_vars_dual(3)

    BC_WEST_(1) = flux_temp(1) - i_app_cm2/(z_cat*nu_cat*Fconst)
    BC_WEST_(2) = flux_temp(2) - flux_temp(1)*MW_cat
    BC_WEST_(3) = flux_temp(3) - i_app_cm2


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

    Phi_2 = c_vars_dual(3)

    vel = c_vars_dual(2)
    dveldx = dcdx_vars_dual(2)

    c_LiPF6     = c_vars_dual(1)
    dcdx_LiPF6  = dcdx_vars_dual(1)

    density_    = density(c_LiPF6)
    d_density_dx = density_%dx(1)*dcdx_LiPF6  ! ∇ρ = ∂ρ/∂c ⋅ ∇c

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    BC_EAST_(1) = flux_temp(1) - i_app_cm2/(z_cat*nu_cat*Fconst)
    ! ∇⋅(ρ𝐯) = 0 = 𝐯⋅∇ρ + ρ∇𝐯
    BC_EAST_(2) = d_density_dx*vel + density_*dveldx - 0.0
    BC_EAST_(3) = Phi_2 - 0.0

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
  ! use echem_mod
  use dnadmod
  ! use TRANSPORT_MODULE
  use GOV_EQNS
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
    cprev(1,j) = c_LiPF6_bulk
    cprev(2,j) = vel_0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition

include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_auto_fill.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_ABDGXY.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_MATINV.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_BAND.f95'
