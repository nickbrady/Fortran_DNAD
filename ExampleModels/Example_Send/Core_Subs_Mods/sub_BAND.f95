SUBROUTINE BAND(J)
! BAND(J) computes delC and calls MATINV to invert B matrix using gaussian elimination.
use user_input, only: N, NJ
use variables, only: A, B, delC, D, G, X, Y, E
implicit none

integer :: j, M
real :: determ

    ! '''
    !     Outputs the arrays 𝐄, 𝛏, and 𝐱
    !     From the block tridiagonal matrices 𝐀, 𝐁, 𝐃, 𝐗, 𝐘 and 𝐆
    !
    !     BAND at j = 1
    !     𝐁𝐜ⱼ + 𝐃𝐜ⱼ₊₁ + 𝐗𝐜ⱼ₊₂ = 𝐆
    !     𝐄ⱼ = 𝐁⁻¹ 𝐃
    !     𝛏ⱼ = 𝐁⁻¹ 𝐆
    !     𝐱  = 𝐁⁻¹ 𝐗
    !
    !     BAND at j = 2
    !     𝐀𝐜ⱼ₋₁ + 𝐁𝐜ⱼ + 𝐃𝐜ⱼ₊₁ = 𝐆
    !     𝐄ⱼ = -[𝐁 + 𝐀𝐄ⱼ₋₁]⁻¹ [𝐀𝐱 + 𝐃]
    !     𝛏ⱼ =  [𝐁 + 𝐀𝐄ⱼ₋₁]⁻¹ [𝐆 - 𝐀𝛏ⱼ₋₁]
    !
    !     BAND at 2 < j < NJ
    !     𝐀𝐜ⱼ₋₁ + 𝐁𝐜ⱼ + 𝐃𝐜ⱼ₊₁ = 𝐆
    !     𝐄ⱼ = -[𝐁 + 𝐀𝐄ⱼ₋₁]⁻¹ 𝐃
    !     𝛏ⱼ =  [𝐁 + 𝐀𝐄ⱼ₋₁]⁻¹ [𝐆 - 𝐀 𝛏ⱼ₋₁]
    !
    !     BAND at j = NJ
    !     𝐘𝐜ⱼ₋₂ + 𝐀𝐜ⱼ₋₁ + 𝐁𝐜ⱼ = 𝐆
    !     𝐜ⱼ = -[𝐀 + [𝐘𝐄ⱼ₋₂] 𝐄ⱼ₋₁ + 𝐁]⁻¹ [𝐆 - 𝐘𝛏ⱼ₋₂ - [𝐀 + 𝐘𝐄ⱼ₋₂] 𝛏ⱼ₋₁ ]
    !     𝐀 <- 𝐀 + 𝐘𝐄ⱼ₋₂
    !     𝐆 <- 𝐆 - 𝐘𝛏ⱼ₋₂
    !     𝐜ⱼ = -[𝐁 + 𝐀 𝐄ⱼ₋₁]⁻¹ [𝐆 - 𝐀 𝛏ⱼ₋₁ ]
    !
    !     Calculate 𝐜ⱼ
    !     𝐜ⱼ = 𝛏ⱼ
    !     𝐜ⱼ = 𝐜ⱼ + 𝐄ⱼ𝐜ⱼ₊₁  (j < NJ)
    !     𝐜₁ = 𝐜₁ + 𝐱𝐜₃
    ! '''

101   FORMAT(15H DETERM=0 AT J=,I4)

      if (j == 1) then
        D(:, 2*N+1)   = G
        D(:, N+1:2*N) = X


        ! CALL MATINV(N, 2*N+1, DETERM, B, D)
        CALL MATINV(N, 2*N+1, DETERM)
        if (DETERM == 0.0) PRINT 101, J


        E(:, N+1, j) =  D(:, 2*N+1)                                             ! C.18
        E(:, 1:N, j) = -D(:, 1:N)                                               ! C.19
        X = -D(:, N+1:2*N)                                                      ! C.20
        RETURN

      else if (j == 2) then
        D(:, 1:N) = D(:, 1:N) + matmul(A, X)

      else if (j == NJ) then
        G = G - matmul(Y, E(:, N+1, j-2))                                       ! C.27
        A = A + matmul(Y, E(:, 1:N, j-2))                                       ! C.28
      end if


      ! ############### if j /= 1
      D(:, N+1) = -G + matmul( A, E(:, N+1, j-1) )                              ! C.23
      B         =  B + matmul( A, E(:, 1:N, j-1) )                              ! C.25


      ! CALL MATINV(N, N+1, DETERM, B, D)     ! returns D = matmul(B^-1, D)
      CALL MATINV(N, N+1, DETERM)

      if (DETERM == 0.0) PRINT 101, J

      E(:,:,j) = -D(:,1:N+1)

      if (j == NJ) then
        delC(:,j) = E(:, N+1, j)      ! delC(:,:) = E(:, N+1, :)                !

        do M = NJ-1, 1, -1
          delC(:,M) = E(:, N+1, M) + matmul(E(:, 1:N, M), delC(:,M+1))          ! C.22
        end do

        delC(:,1) = delC(:,1) + matmul(X, delC(:,3))                            ! Appendix C.17
      end if

      RETURN
END SUBROUTINE BAND
