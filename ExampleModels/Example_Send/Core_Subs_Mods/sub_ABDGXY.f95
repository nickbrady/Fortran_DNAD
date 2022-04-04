subroutine ABDGXY(j)
  ! ABDGXY equates the large coefficents based on the small coefficents.
  ! The coefficents A, B, D, G, X, Y can be found in Newman appendix C.
      use number_of_variables, only: N, NJ
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
