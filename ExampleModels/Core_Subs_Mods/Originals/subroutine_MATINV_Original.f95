SUBROUTINE MATINV(N, M, DETERM)
  use variables, only: A, B, delC, D, ID ! A imported but not used
  implicit real (A-H,O-Z) ! implicits are not good coding practice

 !     BTRY    - largest value in a row
 !     BNEXT   - second largest value in a row
 !     BMAX    - ratio of the two largest magnitude values in a row (BNEXT / BTRY)
 !     ID      - keeps track of available pivot rows
 !
 ! pivot row    - row with the smallest BMAX value (excluding already used pivot columns)
 ! pivot column - max magnitude value in pivot row

 ! indentations included for readability


      DETERM=1.0
      DO 1 I=1,N
1         ID(I)=0

      DO 18 NN=1,N
          BMAX=1.1

          DO 6 I=1,N
              IF (ID(I).NE.0) GOTO 6
                  BNEXT=0.0
                  BTRY=0.0

                  DO 5 J=1,N
                      IF (ID(J).NE.0) GOTO 5

                          IF (ABS(B(I,J)).LE.BNEXT) GOTO 5
                              BNEXT=ABS(B(I,J))

                              IF (BNEXT.LE.BTRY) GOTO 5
                                  BNEXT=BTRY
                                  BTRY=ABS(B(I,J))
                                  JC=J
5                 CONTINUE

                  IF (BNEXT.GE.BMAX*BTRY) GOTO 6
                    BMAX=BNEXT/BTRY
                    IROW=I
                    JCOL=JC
6         CONTINUE

          IF (ID(JC).EQ.0) GOTO 8
            DETERM=0.0
            RETURN

8         ID(JCOL)=1
          IF (JCOL.EQ.IROW) GOTO 12

9           DO 10 J=1,N
              SAVE=B(IROW,J)
              B(IROW,J)=B(JCOL,J)
10            B(JCOL,J)=SAVE

            DO 11 K=1,M
              SAVE=D(IROW,K)
              D(IROW,K)=D(JCOL,K)
11            D(JCOL,K)=SAVE

12        F=1.0/B(JCOL,JCOL)
          DO 13 J=1,N
13          B(JCOL,J)=B(JCOL,J)*F
          DO 14 K=1,M
14          D(JCOL,K)=D(JCOL,K)*F

          DO 18 I=1,N
            IF (I.EQ.JCOL) GOTO 18
              F=B(I,JCOL)
              DO 16 J=1,N
16              B(I,J)=B(I,J)-F*B(JCOL,J)
              DO 17 K=1,M
17              D(I,K)=D(I,K)-F*D(JCOL,K)

18    CONTINUE

      RETURN

      end
