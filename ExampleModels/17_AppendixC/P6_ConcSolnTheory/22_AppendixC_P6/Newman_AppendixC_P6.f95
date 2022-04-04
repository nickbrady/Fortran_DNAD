    ! Transient, multicomponent migration and diffusion.

program Main
    IMPLICIT REAL*8 (A-H,O-Z)
    COMMON A(11,11),B(11,11),C(11,402),D(11,23),G(11) ,X(11,11), Y(11,11), N ,NJ
    dimension cold(11, 402), dif(7,7), z(7), s(7), ref(7), u(7), cin(7)

    character(len=65) :: header_fmt
    character(len=65) :: data_fmt

101 format (i4,a6,2e15.5,2f6.1)
102 format (4f8.4,A6)

! 99  READ *, n,cr0
    N = 5
    cr0 = 0.0

    if(n.le.0) stop

    ! U(i) not used in program
    ! read 102, (U(i), dif(1,i), z(i), s(i), ref(i), i=2,n)
    ! ref(2) = ' OH'
    ! ref(3) = ' Cl'
    ! ref(4) = ' Na'
    ! ref(5) = ' O2'

    dif = 0.0
    dif(1,2) = 5.26000  ! OH
    dif(1,3) = 2.03200  ! Cl
    dif(1,4) = 1.33400  ! Na
    dif(1,5) = 2.00000  ! O2

    z = 0
    z(2) = -1
    z(3) = -1
    z(4) = +1
    z(5) = 0

    s = 0
    s(2) = -4.0
    s(5) = 1.0

    ! read *, mode,nj,Scm3
    mode = 3
    nj   = 101
    Scm3 = 0.032

! 96  read *, (cin(i),i=2,n)
    cin(2) = 0.0
    cin(3) = 1.0
    cin(4) = 1.0
    cin(5) = 2.e-4

    IF(CIN(2).LT.0.0) STOP

    cin(n+1) = 0.0
    nm = n                            ! number of species, including the solvent
    n = 2*n+1                         ! number of unknowns, including the potential
    ctot = 55.5d-3                    ! mol/cm3 c

    cin(1) = ctot*1000.
    do i = 2,nm
      Dif(1,i)  = dif(1,i)*ctot        ! mol/cm-s, Dli times total concentration
      Dif(i,1)  = Dif(1,i)
      cin(1)    = cin(1) - cin(i)      ! solvent concentration
    enddo

    do i = 1,nm
      do k = 1,nm
        if(dif(i,k).eq.0.0) dif(i,k) = 1.d6 ! large value for solutes
      enddo
    enddo

    h = 6.0/dble(nj-1)                  ! xmax = 6
    do j = 1,nj                         ! initial values
      do i = 1,nm
        c(i,j)    = 0.0                 ! flux densities
        c(nm+i,j) = cin(i)/ctot/1000    ! mole fractions
      enddo
      c(n,j)  = 0.0                     ! potential
    enddo

    ! ref(1) = ' H2O'
    z(1) = 0.0
    s(1) = 0.0
    do i = 1,nm
      print 101, i, ref(i), dif(1,i)/ctot, cin(i), z(i), s(i)
    enddo

    jcount = 0
5   jcount = jcount+1
    j = 0
    do i = 1,n
      do k = 1,n
        x(i,k) = 0.0
        y(i,k) = 0.0
      enddo
    enddo

7   j = j+1
    do i = 1,n
      g(i) = 0.0
      cold(i,j) = c(i,j)
      do k = 1,n
        a(i,k) = 0.0
        b(i,k) = 0.0
        d(i,k) = 0.0
      enddo
    enddo
! ----------------------------------------------------------------------------------------------------------------------------------
    g(n-1) = 1.0                                      ! sum of mole fractions = 1
    g(n)   = 0.0                                      ! electroneutrality
    do i = 1,nm
      g(n-1)      = g(n-1) - c(nm+i,j)
      b(n-1,nm+i) = 1.0
      g(n)        = g(n)-z(i)*c(nm+i,j)
      b(n,nm+i)   = z(i)
    enddo

    if(j.gt.1) then                                   ! material balances
      do i = 1,nm
        g(i)      = c(i,j) - c(i,j-1) - 2.0*h*dble(j-1+j-2)/2.0*(c(nm+i,j)-c(nm+i,j-1))
        b(i,i)    = -1.0
        a(i,i)    = 1.0
        b(i,nm+i) = 2.0*h*dble(j-1+j-2)/2.0
        a(i,nm+i) = -2.0*h*dble(j-1+j-2)/2.0
      enddo

    else                                              ! boundary conditions at electrode
      g(nm)         = 0.0 - c(nm+nm,j)                ! zero concentration of limiting reactant
      b(nm,nm+nm)   = 1.0
      do i = 1,nm-1
        g(i)        = s(i)/s(nm)*c(nm,j) - c(i,j)      ! relate flux density
        b(i,i)      = 1.0                             ! to that of limiting reactant
        b(i,nm)     = -s(i)/s(nm)
      enddo
    endif

    if(j.lt.nj) then                                  ! multicomponent-diffusion equations
      do i = 1,nm-1
        g(nm+i)       = -(c(nm+i,j+1)-c(nm+i,j))/h -z (i) * (c (nm+i, j + 1) +c (nm+i, j))/2.0*(c(n,j + 1)-c(n,j))/h
        d(nm+i,nm+i)  =  1.0/h+z(i)/2.0*(c(n,j+1) - c(n,j))/h
        b(nm+i,nm+i)  = -1.0/h+z(i)/2.0*(c(n,j+1) - c(n,j))/h
        b(nm+i,n)     = -z(i)*(c(nm+i,j +1) + c(nm+i,j))/2.0/h
        d(nm+i,n)     =  z(i)*(c(nm+i,j+1) + c(nm+i,j))/2.0/h

        do k = 1,nm
          if (k.ne.i) then
            g(nm+i) =g(nm+i) + (c(nm+i,j)*c(k,j) +c(nm+i,j+1)*c(k,j+1) -c(nm+k,j)*c(i,j) -c(nm+k,j+1)*c(i,j+1))*dif(1,3)/dif(i,k)/2
            b(nm+i,nm+i)  =  b(nm+i,nm+i) - c(k,j )*dif(1,3)/dif(i,k)/2.0
            d(nm+i,nm+i)  =  d(nm+i,nm+i) - c(k,j+1)*dif(1,3)/dif(i,k)/2.0
            b(nm+i,nm+k)  =  b(nm+i,nm+k) + c(i,j )*dif(1,3)/dif(i,k)/2.0
            d(nm+i,nm+k)  =  d(nm+i,nm+k) +c(i,j+1)*dif(1,3)/dif(i,k)/2.0

            b(nm+i,k)     = -c(nm+i, j )*dif(1,3)/dif(i,k)/2.0
            d(nm+i,k)     = -c(nm+i,j+1)*dif(1,3)/dif(i,k)/2.0
            b(nm+i,i)     =  b(nm+i,i) + c(nm+k,j )*dif(1,3)/dif(i,k)/2.0
            d(nm+i,i)     =  d(nm+i,i) + c(nm+k,j+1)*dif(1,3)/dif(i,k)/2.0
          endif
        enddo
      enddo

    else                                              ! boundary conditions at bulk
      do i = 3,nm+1
        g(nm+i-2)       = c(nm+i,j) - c(nm+i,j)       ! set concentration or potential
        b(nm+i-2,nm+i)  = 1.0
      enddo
    endif
!
    call band(j)

    if(j.lt.nj) go to 7

    nerr=0
    do j = 1,nj
      do i = 1,n
        if(dabs(c(i,j)).gt.1.d-10*dabs(cold (i,j)).and.dabs(c(i,j)).gt.1.d-18) nerr = nerr+1
        c(i,j) = cold(i,j)+c(i,j)
      enddo
    enddo
! ----------------------------------------------------------------------------------------------------------------------------------
    print *, nerr,'fluxes ',(c(i,1),i=1,nm)
    if ((nerr.gt.0).and.(jcount.lt.30)) go to 5
    ! print *, ' z    OH-    Cl-    Na+    O2'
    ! print 999, (h*dble(j-1), ctot*1000.*c(7,j), ctot*1000.*c(8,j), ctot*1000.*c(9,j), ctot*1000.*c(10,j), j=1,nj)
    print *, 'solution number ', jcount

    header_fmt = '(1(A12,2X),   20(A15,2X)) '
    data_fmt   = '(1(F12.5,2X), 20(ES15.8,2X))'


! ***** WRITE TO SCREEN **** !
    write(*, header_fmt) 'Position', 'C_O2', 'C_OH',  'C_Na',  'C_Cl',  'Phi',   'C_H2O'
    write(*, header_fmt) 'η',       'mol/L', 'mol/L', 'mol/L', 'mol/L', 'Volts', 'mol/L'
    do j = 1, NJ                                    !
      c_H2O   = ctot*1000.*c(6,j)
      c_OH    = ctot*1000.*c(7,j)
      c_Cl    = ctot*1000.*c(8,j)
      c_Na    = ctot*1000.*c(9,j)
      c_O2    = ctot*1000.*c(10,j)
      Phi     = c(11,j)
                                                    !
      write(*, data_fmt) h*dble(j-1), c_O2, c_OH, c_Na, c_Cl, Phi, c_H2O
    end do



! ***** WRITE TO FILE ****** !
    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    ! write the headers on the first entrance into write all voltage
    write(56, header_fmt) 'Position', 'C_O2', 'C_OH',  'C_Na',  'C_Cl',  'Phi',   'C_H2O'
    write(56, header_fmt) 'η',       'mol/L', 'mol/L', 'mol/L', 'mol/L', 'Volts', 'mol/L'

    do j = 1, NJ
      c_H2O   = ctot*1000.*c(6,j)
      c_OH    = ctot*1000.*c(7,j)
      c_Cl    = ctot*1000.*c(8,j)
      c_Na    = ctot*1000.*c(9,j)
      c_O2    = ctot*1000.*c(10,j)
      Phi     = c(11,j)
                                                    !
      write(56, data_fmt) h*dble(j-1), c_O2, c_OH, c_Na, c_Cl, Phi, c_H2O
    end do

    ! j = NJ
    ! print*, c(1,j), c(2,j), c(3,j), c(4,j), c(5,j), c(6,j), c(7,j), c(8,j), c(9,j), c(8,j), c(11,j)

999 format (f10.3, 4f15.6)
997 format (f10.3, 1p3e15.7)

end program Main


!***********************************MATINV*****************************************

SUBROUTINE MATINV(N,M,DETERM)
! use variables, only: A,B,delC,D,ID
implicit double precision (A-H,O-Z)
COMMON A(11,11), B(11,11), C(11,402), D(11,23)
DIMENSION ID(11)


      DETERM=1.0
      DO 1 I=1,N
1     ID(I)=0
      DO 18 NN=1,N
      BMAX=1.1
      DO 6 I=1,N
      IF (ID(I).NE.0) GOTO 6
      BNEXT=0.0
      BTRY=0.0
      DO 5 J=1,N
      IF (ID(J).NE.0) GOTO 5
      IF (DABS(B(I,J)).LE.BNEXT) GOTO 5
      BNEXT=DABS(B(I,J))
      IF (BNEXT.LE.BTRY) GOTO 5
      BNEXT=BTRY
      BTRY=DABS(B(I,J))
      JC=J
5     CONTINUE
      IF (BNEXT.GE.BMAX*BTRY) GOTO 6
      BMAX=BNEXT/BTRY
      IROW=I
      JCOL=JC
6     CONTINUE
      IF (ID(JC).EQ.0) GOTO 8
      DETERM=0.0
      RETURN
8     ID(JCOL)=1
      IF (JCOL.EQ.IROW) GOTO 12
9     DO 10 J=1,N
      SAVE=B(IROW,J)
      B(IROW,J)=B(JCOL,J)
10    B(JCOL,J)=SAVE
      DO 11 K=1,M
      SAVE=D(IROW,K)
      D(IROW,K)=D(JCOL,K)
11    D(JCOL,K)=SAVE
12    F=1.0/B(JCOL,JCOL)
      DO 13 J=1,N
13    B(JCOL,J)=B(JCOL,J)*F
      DO 14 K=1,M
14    D(JCOL,K)=D(JCOL,K)*F
      DO 18 I=1,N
      IF (I.EQ.JCOL) GOTO 18
      F=B(I,JCOL)
      DO 16 J=1,N
16    B(I,J)=B(I,J)-F*B(JCOL,J)
      DO 17 K=1,M
17    D(I,K)=D(I,K)-F*D(JCOL,K)
18    CONTINUE
      RETURN
END SUBROUTINE MATINV


!*************************************BAND******************************************

SUBROUTINE BAND(J)
! use variables, only: A,B,delC,D,G,X,Y,NP1,E
! use user_input, only: N,NJ
implicit double precision (A-H,O-Z)
COMMON A(11,11), B(11,11), C(11,402), D(11,23), G(11), X(11,11), Y(11,11), N ,NJ
DIMENSION E(11,12,402)
SAVE E, NP1


101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
      D(I,2*N+1)=G(I)
      DO 2 L=1,N
      LPN=L+N
2     D(I,LPN)=X(I,L)
      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J
4     DO 5 K=1,N
      E(K,NP1,1)=D(K,2*N+1)
      DO 5 L=1,N
      E(K,L,1)=-D(K,L)
      LPN=L+N
5     X(K,L)=-D(K,LPN)
      RETURN
6     DO 7 I=1,N
      DO 7 K=1,N
      DO 7 L=1,N
7     D(I,K)=D(I,K)+A(I,L)*X(L,K)
8     IF (J-NJ) 11,9,9
9     DO 10 I=1,N
      DO 10 L=1,N
      G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
      DO 10 M=1,N
10    A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)
11    DO 12 I=1,N
      D(I,NP1)=-G(I)
      DO 12 L=1,N
      D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
      DO 12 K=1,N
12    B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)
      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N
      DO 15 M=1,NP1
15    E(K,M,J)=-D(K,M)
      IF (J-NJ) 20,16,16
16    DO 17 K=1,N
17    C(K,J)=E(K,NP1,J)
      DO 18 JJ=2,NJ
      M=NJ-JJ+1
      DO 18 K=1,N
      C(K,M)=E(K,NP1,M)
      DO 18 L=1,N
18    C(K,M)=C(K,M) +E(K,L,M)*C(L,M+1)
      DO 19 L=1,N
      DO 19 K=1,N
19    C(K,1)=C(K,1)+X(K,L)*C(L,3)
20    RETURN
END SUBROUTINE BAND
