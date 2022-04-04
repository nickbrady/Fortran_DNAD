SUBROUTINE BAND(J)
! BAND(J) computes delC and calls MATINV to solve the problem using gaussian elimination.
use variables, only: A, B, delC, D, G, X, Y, NP1, E
use user_input, only: N, NJ
implicit real (A-H,O-Z)

101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
        D(I,2*N+1)=G(I)
        DO 2 L=1,N
          LPN=L+N
2         D(I,LPN)=X(I,L)

      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J

4     DO 5 K=1,N
        E(K,NP1,1)=D(K,2*N+1)
        DO 5 L=1,N
          E(K,L,1)=-D(K,L)
          LPN=L+N
5         X(K,L)=-D(K,LPN)

      RETURN

6     DO 7 I=1,N
        DO 7 K=1,N
          DO 7 L=1,N
7           D(I,K)=D(I,K)+A(I,L)*X(L,K)

8     IF (J-NJ) 11,9,9

9     DO 10 I=1,N
        DO 10 L=1,N
          G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
          DO 10 M=1,N
10          A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)

11    DO 12 I=1,N
        D(I,NP1)=-G(I)
        DO 12 L=1,N
          D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
          DO 12 K=1,N
12          B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)

      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J

14    DO 15 K=1,N
        DO 15 M=1,NP1
15        E(K,M,J)=-D(K,M)

      IF (J-NJ) 20,16,16

16    DO 17 K=1,N
17      delC(K,J)=E(K,NP1,J)

      DO 18 JJ=2,NJ
        M=NJ-JJ+1
        DO 18 K=1,N
          delC(K,M)=E(K,NP1,M)
          DO 18 L=1,N
18          delC(K,M)=delC(K,M) +E(K,L,M)*delC(L,M+1)

      DO 19 L=1,N
        DO 19 K=1,N
19        delC(K,1)=delC(K,1)+X(K,L)*delC(L,3)

20    RETURN

      end
