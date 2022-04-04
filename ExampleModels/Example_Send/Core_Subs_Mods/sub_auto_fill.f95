subroutine auto_fill(j)
  use number_of_variables! N, NJ
  use GOV_EQNS
  use dnadmod
  use variables, only: alphaE, alphaW, betaE, betaW, rj, fE, fW, dE, dW, smG, cprev

  implicit none

  integer :: j
  integer :: ic

  ! variables and their derivatives at the control volume interfaces
  real, dimension(N) :: cW, cE, dcdxW, dcdxE

  ! dual variables
  type(dual), dimension(N) :: cW_dual, dcdxW_dual, flux_dualW
  type(dual), dimension(N) :: cE_dual, dcdxE_dual, flux_dualE
  type(dual), dimension(N) :: cj_dual, reaction_dual, accumulation_dual
  type(dual), dimension(N) :: boundary_conditionW, boundary_conditionE

  ! set all matrix variables (dW, dE, fW, fE, rj, smG) to 0.0
  dW = 0.0
  dE = 0.0
  fW = 0.0
  fE = 0.0
  rj = 0.0
  smG = 0.0

  !-----------------------------------------------------------------------------
  ! Calculate cW, cE, dcdxW, dcdxE
  !-----------------------------------------------------------------------------
  if (j /= 1) then                          !-----------------------------------
    alphaW = delx(j-1)/(delx(j-1)+delx(j))  ! West Side Interface Variables:
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
    alphaE = delx(j)/(delx(j+1)+delx(j))    ! cE, dcdxE
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
        fE(ic, :) = boundary_conditionW(ic)%dx(1:N)
        dE(ic, :) = boundary_conditionW(ic)%dx(N+1:2*N)

        smG(ic)   = -(-boundary_conditionW(ic)%x)
    end do
                                                !-------------------------------
  else if (j == NJ) then                        ! East Side Boundary Conditions
      boundary_conditionE = Boundary_EAST(cW_dual, dcdxW_dual) !----------------
    do ic = 1,N
        fW(ic, :) = boundary_conditionE(ic)%dx(1:N)
        dW(ic, :) = boundary_conditionE(ic)%dx(N+1:2*N)

        smG(ic)   = -(boundary_conditionE(ic)%x)
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
        fW(ic,:) = flux_dualW(ic)%dx(1:N)     * Crx_Area(1,j)
        fE(ic,:) = flux_dualE(ic)%dx(1:N)     * Crx_Area(2,j)
        dW(ic,:) = flux_dualW(ic)%dx(N+1:2*N) * Crx_Area(1,j)
        dE(ic,:) = flux_dualE(ic)%dx(N+1:2*N) * Crx_Area(2,j)

        rj(ic,:) = ( reaction_dual(ic)%dx(1:N)                &
                 &  - accumulation_dual(ic)%dx(1:N) )*Cntrl_Vol(j)

        smG(ic)  = -(   flux_dualW(ic)%x * Crx_Area(1,j)      &
                 &    - flux_dualE(ic)%x * Crx_Area(2,j)      &
                 &    + reaction_dual(ic)%x * Cntrl_Vol(j) )
    end do
  end if

end subroutine auto_fill
