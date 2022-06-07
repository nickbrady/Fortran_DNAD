! auto_fill is a subroutine that uses dual numbers (automatic differentiation) to
! compute the matrices 𝐫, 𝐟ᴱ, 𝐟ᵂ, 𝐝ᴱ, 𝐝ᵂ, 𝐠 from the governing equations, which
! are given in the module GOV_EQNS.
! Interface variables are created from the dependent variables from the previous interation,
! cprev; these are the interpolated values of the dependent variables at the interfaces
! 𝐜ᵂ = αᵂ⋅𝐜ⱼ + (1 - αᵂ)⋅𝐜ⱼ₋₁
! 𝐜ᴱ = αᴱ⋅𝐜ⱼ₊₁ + (1 - αᴱ)⋅𝐜ⱼ
! ∇𝐜ᵂ = βᵂ⋅(𝐜ⱼ - 𝐜ⱼ₋₁)
! ∇𝐜ᴱ = βᴱ⋅(𝐜ⱼ₊₁ - 𝐜ⱼ)
! where:
! αᵂ = Δxⱼ / (Δxⱼ₋₁ + Δxⱼ)
! αᴱ = Δxⱼ / (Δxⱼ + Δxⱼ₊₁)
! βᵂ =   2 / (Δxⱼ₋₁ + Δxⱼ)
! βᴱ =   2 / (Δxⱼ + Δxⱼ₊₁)
! Dual variables are created from the dependent variables from the previous interation,
! cprev, (also called control volume variables) as well as the interface variables,
! 𝐜ᵂ, 𝐜ᴱ, ∇𝐜ᵂ, ∇𝐜ᴱ. The values of 𝐫, 𝐟ᴱ, 𝐟ᵂ, 𝐝ᴱ, 𝐝ᵂ, 𝐠 are then calculated as:
! 𝐫  = (∂R/∂𝐜ⱼ - ∂(Accum)/∂𝐜ⱼ)⋅ΔV
! 𝐟ᴱ = ∂𝐍/∂𝐜ᴱ  ⋅ Aₓᴱ
! 𝐟ᵂ = ∂𝐍/∂𝐜ᵂ  ⋅ Aₓᵂ
! 𝐝ᴱ = ∂𝐍/∂∇𝐜ᴱ ⋅ Aₓᴱ
! 𝐝ᵂ = ∂𝐍/∂∇𝐜ᵂ ⋅ Aₓᵂ
! 𝐠  = -( 𝐍(𝐜ᵂ,∇𝐜ᵂ)⋅Aₓᵂ - 𝐍(𝐜ᴱ,∇𝐜ᵂ)⋅Aₓᴱ + R(𝐜ⱼ))
!
! At the boundaries ...

! edits: June 3, 2022
! type(dual) --> type(dual_complex)
! real, dimension(N) :: cW, cE, dcdxW, dcdxE --> complex, dimension(N) :: cW, cE, dcdxW, dcdxE
!
! real, dimension(N), intent(in) :: c_vars   --> complex, dimension(N), intent(in) :: c_vars
! real, dimension(N*2)           :: dx_array --> complex, dimension(N*2)           :: dx_array
! real, dimension(N), intent(in) :: dcdx     --> complex, dimension(N), intent(in) :: dcdx
! real, dimension(N*2)           :: dx_array --> complex, dimension(N*2)           :: dx_array

subroutine auto_fill(j)
  use user_input, only: N, NJ, UPWIND, direction
  use GOV_EQNS
  use dnadmod
  use variables, only: alphaE, alphaW, betaE, betaW, &
                     & rj, fE, fW, dE, dW, smG, &
                     & cprev, xx, delX, Cntrl_Vol, Crx_Area

  implicit none

  ! integer, intent(in) :: j
  ! intent(in) :: cprev
  ! intent(out) :: dW, dE, fW, fE, rj, smG

  integer :: j
  integer :: ic

  ! variables and their derivatives at the control volume interfaces
  complex, dimension(N) :: cW, cE, dcdxW, dcdxE

  ! dual variables
  type(dual_complex), dimension(N) :: cW_dual, dcdxW_dual, flux_dualW
  type(dual_complex), dimension(N) :: cE_dual, dcdxE_dual, flux_dualE
  type(dual_complex), dimension(N) :: cj_dual, reaction_dual, accumulation_dual
  type(dual_complex), dimension(N) :: boundary_conditionW, boundary_conditionE

  ! set all matrix variables (dW, dE, fW, fE, rj, smG) to 0.0
  ! dW = 0.0
  ! dE = 0.0
  ! fW = 0.0
  ! fE = 0.0
  ! rj = 0.0
  ! smG = 0.0
  ! setting these matrices to 0 is redundant except for rj at the boundaries
  ! dW, dE, fW, fE, rj, smG are all set for the interior points
  ! at boundary_west ( j=1), only dE, fE are used in ABDGXY
  ! at boundary_east (j=NJ), only dW, fW are used in ABDGXY
  ! smG is set even at the boundaries
  ! special care only needs to be given to rj. But it seems this could be reworked
  ! because Cntrl_Vol = 0 at these interface regions, but for now just setting
  ! rj = 0 at the boundaries seems to work

  !-----------------------------------------------------------------------------
  ! Calculate cW, cE, dcdxW, dcdxE
  !-----------------------------------------------------------------------------
  if (UPWIND) then
    if (trim(direction) == 'EastToWest') then
      ! Upwind scheme - flow from East to West (negative current)
      alphaW = 1.0
      alphaE = 1.0
    else if (trim(direction) == 'WestToEast') then
      ! Upwind scheme - flow from West to East (positive current)
      alphaW = 0.0
      alphaE = 0.0
    end if

  else
    if (j /= 1) then                          !-----------------------------------
      alphaW = delx(j-1)/(delx(j-1)+delx(j))
    end if
    if (j /= NJ) then                         ! East Side Interface Variables:
      alphaE = delx(j)/(delx(j+1)+delx(j))
    end if
  end if

  if (j /= 1) then                          !-----------------------------------
    ! alphaW = delx(j-1)/(delx(j-1)+delx(j))  ! West Side Interface Variables:
    betaW = 2.0/(delx(j-1)+delx(j))         ! cW, dcdxW
    do ic=1,N                               !-----------------------------------
      cW(ic)    = alphaW*cprev(ic,j) + (1.0 - alphaW)*cprev(ic,j-1)
      dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
    end do
    cW_dual = c_to_dual(cW)
    dcdxW_dual = dcdx_to_dual(dcdxW)
  end if
                                            !-----------------------------------
  if (j /= NJ) then                         ! East Side Interface Variables:
    ! alphaE = delx(j)/(delx(j+1)+delx(j))    ! cE, dcdxE
    betaE = 2.0/(delx(j)+delx(j+1))         !-----------------------------------
    do ic=1,N
      cE(ic)    = alphaE*cprev(ic,j+1) + (1.0 - alphaE)*cprev(ic,j)
      dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
    end do
    cE_dual = c_to_dual(cE)
    dcdxE_dual = dcdx_to_dual(dcdxE)
  end if

  !-----------------------------------------------------------------------------
  ! Boundary Conditions
  !-----------------------------------------------------------------------------
  if (j == 1) then                              ! West Side Boundary Conditions
    boundary_conditionW = Boundary_WEST(cE_dual, dcdxE_dual) !------------------
    do ic = 1,N
        fE(ic, :) = boundary_conditionW(ic)%dz(1:N)
        dE(ic, :) = boundary_conditionW(ic)%dz(N+1:2*N)
        rj = 0.0

        smG(ic)   = -(-boundary_conditionW(ic)%z)
    end do
                                                !-------------------------------
  else if (j == NJ) then                        ! East Side Boundary Conditions
      boundary_conditionE = Boundary_EAST(cW_dual, dcdxW_dual) !----------------
    do ic = 1,N
        fW(ic, :) = boundary_conditionE(ic)%dz(1:N)
        dW(ic, :) = boundary_conditionE(ic)%dz(N+1:2*N)
        rj = 0.0

        smG(ic)   = -(boundary_conditionE(ic)%z)
    end do

  !-----------------------------------------------------------------------------
  ! Governing Equations
  !-----------------------------------------------------------------------------
  else                                        !------------ Fluxes -------------
    flux_dualW = FLUX(cW_dual, dcdxW_dual)    ! West Side Flux
    flux_dualE = FLUX(cE_dual, dcdxE_dual)    ! East Side Flux

    cj_dual = c_to_dual(cprev(:,j))           !---------------------------------
    reaction_dual = RXN(cj_dual)              ! Control Volume - Rxn, Accum
    accumulation_dual = ACCUM(cj_dual)        !---------------------------------

    do ic = 1,N                               ! ic - equation number
        fW(ic,:) = flux_dualW(ic)%dz(1:N)     * Crx_Area(1,j)
        fE(ic,:) = flux_dualE(ic)%dz(1:N)     * Crx_Area(2,j)
        dW(ic,:) = flux_dualW(ic)%dz(N+1:2*N) * Crx_Area(1,j)
        dE(ic,:) = flux_dualE(ic)%dz(N+1:2*N) * Crx_Area(2,j)

        rj(ic,:) = ( reaction_dual(ic)%dz(1:N)                &
                 &  - accumulation_dual(ic)%dz(1:N) )*Cntrl_Vol(j)

        smG(ic)  = -(   flux_dualW(ic)%z * Crx_Area(1,j)      &
                 &    - flux_dualE(ic)%z * Crx_Area(2,j)      &
                 &    + reaction_dual(ic)%z * Cntrl_Vol(j) )
    end do
  end if

contains

  ! ****************************************************************************
  ! ******* functions to convert concentration and dcdx variables to dual ******
  ! ----------------------------------------------------------------------------
  function c_to_dual(c_vars)          result(c_dual)
    type(dual_complex), dimension(N)  :: c_dual
    complex, dimension(N), intent(in) :: c_vars
    complex, dimension(N*2)           :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(ic) = 1.0

      c_dual(ic) = dual_complex(c_vars(ic), dx_array)
    end do
  end function c_to_dual


  function dcdx_to_dual(dcdx)         result(dcdx_dual)
    type(dual_complex), dimension(N)  :: dcdx_dual
    complex, dimension(N), intent(in) :: dcdx
    complex, dimension(N*2)           :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(N+ic) = 1.0

      dcdx_dual(ic) = dual_complex(dcdx(ic), dx_array)
    end do

  end function dcdx_to_dual
! ******************************************************************************

end subroutine auto_fill
