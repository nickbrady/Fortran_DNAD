! See Electrochemical Systems 3rd Edition by John Newman and Karen E. Thomas-Alyea, Appendix C
! Specifically C.2 Convergence Over Nonlinearities

! If 𝐠ⱼ is a matrix of the coupled-partial differential equations that describe a particular system
! 𝐠ⱼ = f(𝐜ⱼ₋₁, 𝐜ⱼ, 𝐜ⱼ₊₁)
! We can solve the equation
! 𝐠ⱼ(𝐜ⱼ₋₁, 𝐜ⱼ, 𝐜ⱼ₊₁) = 0                                                                    (C.8)
! using a Newton-Raphson method (at each node point)
! 0 = 𝐠ⱼᵒ + (∂𝐠ⱼ/∂𝐜ⱼ₋₁)ₒ Δ𝐜ⱼ₋₁ + (∂𝐠ⱼ/∂𝐜ⱼ)ₒ Δ𝐜ⱼ + (∂𝐠ⱼ/∂𝐜ⱼ₊₁)ₒ Δ𝐜ⱼ₊₁                         (C.12)
! where
! 𝐠ⱼᵒ = 𝐠ⱼ(𝐜ᵒⱼ₋₁, 𝐜ᵒⱼ, 𝐜ᵒⱼ₊₁)
! and (∂𝐠ⱼ/∂𝐜ⱼ)ₒ denotes the evaluation of ∂𝐠ⱼ/∂𝐜ⱼ at 𝐜ᵒⱼ
! where 𝐜ᵒⱼ are the values of the dependent variables from the previous time-step (or iteration)
! then the matrices 𝐀, 𝐁, 𝐃, 𝐗, 𝐘, 𝐆 can be defined as follows                                (C.14, C.15)
! 𝐀ⱼ = ∂𝐠ⱼ/∂𝐜ⱼ₋₁
! 𝐁ⱼ = ∂𝐠ⱼ/∂𝐜ⱼ
! 𝐃ⱼ = ∂𝐠ⱼ/∂𝐜ⱼ₊₁
! 𝐗ⱼ = ∂𝐠ⱼ/∂𝐜ⱼ₊₂
! 𝐘ⱼ = ∂𝐠ⱼ/∂𝐜ⱼ₋₂
! 𝐆ⱼ = 𝐠ⱼᵒ
! yielding
! 0 = 𝐠ⱼᵒ + 𝐀ⱼΔ𝐜ⱼ₋₁ + 𝐁ⱼΔ𝐜ⱼ + 𝐃ⱼΔ𝐜ⱼ₊₁                                                        (C.13)
! And at the boundaries:
! (j = 1)     0 = 𝐠ⱼᵒ                     + 𝐁ⱼΔ𝐜ⱼ + 𝐃ⱼΔ𝐜ⱼ₊₁ + 𝐗ⱼΔ𝐜ⱼ₊₂
! (j = NJ)    0 = 𝐠ⱼᵒ + 𝐘ⱼΔ𝐜ⱼ₋₂ + 𝐀ⱼΔ𝐜ⱼ₋₁ + 𝐁ⱼΔ𝐜ⱼ
!
! Equation C.13 and the boundary condtions are subsequently solved using subroutine BAND(J)

subroutine ABDGXY(j)
      use user_input, only: N, NJ
      use variables, only: A, B, D, G, X, Y, alphaE, alphaW, betaE, betaW, rj, fE, fW, dE, dW, smG
      implicit none
      integer :: j

      if (j.eq.1) then
        X = 0.0
        B = rj - (1.0 - alphaE)*fE + betaE*dE
        D(1:N,1:N) = -alphaE*fE - betaE*dE
        G = smG

        return
      end if

      if (j.eq.NJ) then
        Y = 0.0
        A = (1.d0 - alphaW)*fW - betaW*dW
        B = rj + betaW*dW + alphaW*fW
        G = smG

        return
      end if

      A = (1.d0 - alphaW)*fW - betaW*dW
      B = rj + betaW*dW + alphaW*fW - (1.0 - alphaE)*fE + betaE*dE
      D(1:N,1:N) = -alphaE*fE - betaE*dE
      G = smG

      return

end subroutine ABDGXY
