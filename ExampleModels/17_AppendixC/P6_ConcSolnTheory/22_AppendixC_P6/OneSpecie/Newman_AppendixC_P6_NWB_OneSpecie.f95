    ! Transient, multicomponent migration and diffusion.

program Main
    IMPLICIT REAL*8 (A-H,O-Z)
    COMMON A(11,11),B(11,11),C(11,402),D(11,23),G(11) ,X(11,11), Y(11,11), N ,NJ
    dimension cold(11, 402), dif(7,7), z(7), s(7), ref(7), u(7), cin(7)

    character(len=65) :: header_fmt
    character(len=65) :: data_fmt, matrix_fmt
    matrix_fmt = '(20(ES15.8,2X))'

    Fconst = 96485
    Rigc   = 8.314
    Temp   = 298

110 format (i4,2e15.5,2f6.1)
102 format (4f8.4,A6)

! 99  READ *, n,cr0
    N = 1
    cr0 = 0.0

    if(n.le.0) stop

    ! U(i) not used in program
    ! read 102, (U(i), dif(1,i), z(i), s(i), ref(i), i=2,n)
    ! ref(2) = ' OH'
    ! ref(3) = ' Cl'
    ! ref(4) = ' Na'
    ! ref(5) = ' O2'

    dif = 0.0
    ! dif(1,2) = 5.26000  ! OH
    ! dif(1,3) = 2.03200  ! Cl
    ! dif(1,4) = 1.33400  ! Na
    ! dif(1,5) = 2.00000  ! O2

    dif(1,2) = 2.00000

    ! z = 0
    ! z(2) = -1
    ! z(3) = -1
    ! z(4) = +1
    ! z(5) = 0
    z = 0.0
    z(1) = 0.0

    ! O2 + 4e- + 2 H2O --> 4 OH-
    ! s = 0
    ! s(1) = 2.0
    ! s(2) = -4.0
    ! s(5) = 1.0

    s(2) = 1.0

    ! read *, mode,nj,Scm3
    mode = 3
    nj   = 101
    Scm3 = 0.032

! 96  read *, (cin(i),i=2,n)
    ! cin(2) = 0.0
    ! cin(3) = 1.0
    ! cin(4) = 1.0
    ! cin(5) = 2.e-4

    cin(2) = 2.e-4
    cbulk_O2 = 2.e-4

    IF(CIN(2).LT.0.0) STOP

    cin(n+1) = 0.0
    nm = n                            ! number of species, including the solvent
    n = 2*n!+1                         ! number of unknowns, including the potential
    ctot = 55.5d-3                    ! mol/cm3 c

    cin(1) = ctot*1000.
    do i = 2,nm
      Dif(1,i)  = dif(1,i)*ctot        ! mol/cm-s, Dli times total concentration
      Dif(i,1)  = Dif(1,i)
      cin(1)    = cin(1) - cin(i)      ! solvent concentration
    enddo

    do i = 1,nm
      do k = 1,nm
        if(dif(i,k).eq.0.0) dif(i,k) = 1.d90 ! large value for solutes
      enddo
    enddo

    h = 6.0/dble(nj-1)                  ! xmax = 6
    do j = 1,nj                         ! initial values
      do i = 1,nm
        c(i,j)    = 0.0                 ! flux densities
        c(nm+i,j) = cin(i)/ctot/1000    ! mole fractions
      enddo
      ! c(n,j)  = 0.0                     ! potential
    enddo

    ! ref(1) = ' H2O'
    z(1) = 0.0
    s(1) = 0.0
    do i = 1,nm
      print 110, i, dif(1,i)/ctot, cin(i), z(i), s(i)
    enddo
    ! print*,

    iteration = 0
5   iteration = iteration+1
    j = 0

    X = 0.0
    Y = 0.0
    cold = c


    ! print*, n, nm, z
    ! print*, cold(1,1:5)
    ! print*, cold(2,1:5)
    ! print*, cold(3,1:5)*ctot*1000
    ! print*, cold(4,1:5)*ctot*1000

  do j = 1, NJ

    A = 0.0
    B = 0.0
    D(1:N, 1:N) = 0.0
    G = 0.0

    if (j > 1) then
    ! Material Balance:  dN/dη = 2η dc/dη
                      ! dN/dη = (c(1,j) - c(1,j-1))/h
                      ! dc/dη = (c(2,j) - c(2,j-1))/h
                      ! η_avg = [h*(j-1) + h*(j-2)]/2 = h*(j-1+j-2)/2
                      ! (c(1,j) - c(1,j-1))/h = 2h*(j-1+j-2)/2 * (c(2,j) - c(2,j-1))/h
                      ! (c(1,j) - c(1,j-1)) = 2h*(j-1+j-2)/2 * (c(2,j) - c(2,j-1))
                      ! (c(1,j) - c(1,j-1)) - 2h*(j-1+j-2)/2 * (c(2,j) - c(2,j-1)) = 0
      g(1) = -((c(1,j) - c(1,j-1)) - 2*h*dble(j-1+j-2)/2 * (c(2,j) - c(2,j-1)))
      a(1,1) = -1.0
      a(1,2) = + 2*h*dble(j-1+j-2)/2
      b(1,1) = 1.0
      b(1,2) = - 2*h*dble(j-1+j-2)/2
    else
      ! c(2,j) = 0.0
      g(1) = 0.0 - c(2,j)
      b(1,2) = 1.0
    end if

    if (j < NJ) then
      ! Diffusion equation:   N = -dc/dη
                          !   N_avg = [c(1,j+1) - c(1,j)]/2
                          !   dc/dη = (c(2,j) - c(2,j-1))/h
                          !   [c(1,j+1) - c(1,j)]/2 + (c(2,j) - c(2,j-1))/h = 0.0
      g(2) = -( (c(1,j+1) - c(1,j))/2 + (c(2,j) - c(2,j-1))/h )
      b(2,1) = -1.0/2
      b(2,2) = -1.0/h
      d(2,1) = +1.0/2
      d(2,2) = +1.0/h

    else
      g(2) = cbulk_O2/ctot - c(2,j)
      b(2,2) = 1.0

    end if
    ! do i = 1,n
    !   g(i) = 0.0
    !   cold(i,j) = c(i,j)
    !   do k = 1,n
    !     a(i,k) = 0.0
    !     b(i,k) = 0.0
    !     d(i,k) = 0.0
    !   enddo
    ! enddo
! ----------------------------------------------------------------------------------------------------------------------------------
    ! c(1:nm)       Flux_i  : H2O, OH, Cl, Na, O2
    ! c(nm+1:2*nm)  x_i's   : H2O, OH, Cl, Na, O2
    ! c(8)          Cl
    ! c(9)          Na
    ! c(10)         O2
    ! c(11)         Phi
    ! dc_i/dη + z_i*F*c_i/RT*dΦ/dη = Σ( D_OR/D_ik * (c_i*N_k - c_k*N_i) )
    ! dx_i/dη + z_i*F*x_i/RT*dΦ/dη = Σ( D_OR/D_ik * (x_i*N_k - x_k*N_i) )
    ! dN_i/dη = 2η/c_t * dc_i/dη      -->     ! dN_i/dη = 2η * dx_i/dη

    !       Governing Equations                                     BC - Electrode                        BC - Bulk
    ! (1)   dN_H2O/dη = 2η * dx_H2O/dη                              s_H2O / s_O2 * N_O2 - N_H2O = 0       dN_H2O/dη = 2η * dx_H2O/dη
    ! (2)   dN_OH/dη = 2η * dx_OH/dη                                s_OH / s_O2 * N_O2 - N_OH = 0         dN_OH/dη = 2η * dx_OH/dη
    ! (3)   dN_Cl/dη = 2η * dx_Cl/dη                                s_Cl / s_O2 * N_O2 - N_Cl = 0         dN_Cl/dη = 2η * dx_Cl/dη
    ! (4)   dN_Na/dη = 2η * dx_Na/dη                                s_Na / s_O2 * N_O2 - N_Na = 0         dN_Na/dη = 2η * dx_Na/dη
    ! (5)   dN_O2/dη = 2η * dx_O2/dη                                x_O2 = 0.0                            dN_O2/dη = 2η * dx_O2/dη
    ! (6)   dx_i/dη + z_i*F*x_i/RT*dΦ/dη = Σ(D_OR/D_ik * (x_i*N_k - x_k*N_i))                             x_Cl = x_bulk_Cl
    ! (7)   dx_i/dη + z_i*F*x_i/RT*dΦ/dη = Σ(D_OR/D_ik * (x_i*N_k - x_k*N_i))                             x_Na = x_bulk_Na
    ! (8)   dx_i/dη + z_i*F*x_i/RT*dΦ/dη = Σ(D_OR/D_ik * (x_i*N_k - x_k*N_i))                             x_O2 = x_bulk_O2
    ! (9)   dx_i/dη + z_i*F*x_i/RT*dΦ/dη = Σ(D_OR/D_ik * (x_i*N_k - x_k*N_i))                             Phi  = 0.0
    ! (10)  Σx_i = 1                                                Σx_i = 1                              Σx_i = 1
    ! (11)  Σ(z_i*x_i) = 0                                          Σ(z_i*x_i) = 0                        Σ(z_i*x_i) = 0


!     g(n) = 1.0                                      ! sum of mole fractions = 1
!     ! g(n)   = 0.0                                      ! electroneutrality
!     do i = 1,nm
!       g(n)      = g(n) - c(nm+i,j)
!       b(n,nm+i) = 1.0
!
!       ! g(n)        = g(n)-z(i)*c(nm+i,j)
!       ! b(n,nm+i)   = z(i)
!     enddo
!
!     if(j.gt.1) then                                   ! material balances - uses left dcdx --> ( c(j) - c(j-1) ) / h
!       do i = 1,nm                                     ! dN/dη = 2η/cT * dc/dη       -->   dN/dη = 2η * dx/dη
!                                                       ! h*dble(j-1+j-2)/2.0 <-- avg of two positions
!         g(i)      = c(i,j) - c(i,j-1) - 2.0*h*dble(j-1+j-2)/2.0*(c(nm+i,j)-c(nm+i,j-1))
!         ! g(i)      = c(i,j) - c(i,j-1) - 2.0*h*dble(j-1)*(c(nm+i,j)-c(nm+i,j-1))
!         b(i,i)    = -1.0
!         a(i,i)    = 1.0
!
!         b(i,nm+i) = 2.0*h*dble(j-1+j-2)/2.0
!         ! b(i,nm+i) = 2.0*h*dble(j-1)
!         a(i,nm+i) = -2.0*h*dble(j-1+j-2)/2.0
!         ! a(i,nm+i) = -2.0*h*dble(j-1)
!       enddo
!
!     else                                              ! boundary conditions at electrode
!       g(nm)         = 0.0 - c(nm+nm,j)                ! zero concentration of limiting reactant
!       b(nm,nm+nm)   = 1.0
!       do i = 1,nm-1       ! s_i / s_R * N_R - N_i = 0
!         g(i)        = s(i)/s(nm)*c(nm,j) - c(i,j)     ! relate flux density
!         b(i,i)      = 1.0                             ! to that of limiting reactant
!         b(i,nm)     = -s(i)/s(nm)
!       enddo
!     endif
!
!     if(j.lt.nj) then                                  ! multicomponent-diffusion equations
!                                                       ! uses right dcdx --> ( c(j+1) - c(j) ) / h
!       do i = 1,nm-1     ! -dx_i/dη - z_i*F*x_i/RT*dΦ/dη
!         g(nm+i)       = -(c(nm+i,j+1)-c(nm+i,j))/h !- z(i) * (c(nm+i,j+1) + c(nm+i,j))/2.0*(c(n,j+1)-c(n,j))/h
!         d(nm+i,nm+i)  =  1.0/h !+ z(i)/2.0*(c(n,j+1) - c(n,j))/h
!         b(nm+i,nm+i)  = -1.0/h !+ z(i)/2.0*(c(n,j+1) - c(n,j))/h
!         ! b(nm+i,n)     = -z(i)*(c(nm+i,j+1) + c(nm+i,j))/2.0/h
!         ! d(nm+i,n)     =  z(i)*(c(nm+i,j+1) + c(nm+i,j))/2.0/h
!
!         do k = 1,nm
!           if (k.ne.i) then    ! Σ(D_OR/D_ik * (x_i*N_k - x_k*N_i))
!             g(nm+i) =g(nm+i) + (c(nm+i,j)*c(k,j) +c(nm+i,j+1)*c(k,j+1) -c(nm+k,j)*c(i,j) -c(nm+k,j+1)*c(i,j+1))*dif(1,3)/dif(i,k)/2
!             b(nm+i,nm+i)  =  b(nm+i,nm+i) - c(k,j )*dif(1,3)/dif(i,k)/2.0
!             d(nm+i,nm+i)  =  d(nm+i,nm+i) - c(k,j+1)*dif(1,3)/dif(i,k)/2.0
!             b(nm+i,nm+k)  =  b(nm+i,nm+k) + c(i,j )*dif(1,3)/dif(i,k)/2.0
!             d(nm+i,nm+k)  =  d(nm+i,nm+k) + c(i,j+1)*dif(1,3)/dif(i,k)/2.0
! !
!             b(nm+i,k)     = -c(nm+i, j )*dif(1,3)/dif(i,k)/2.0
!             d(nm+i,k)     = -c(nm+i,j+1)*dif(1,3)/dif(i,k)/2.0
!             b(nm+i,i)     =  b(nm+i,i) + c(nm+k,j )*dif(1,3)/dif(i,k)/2.0
!             d(nm+i,i)     =  d(nm+i,i) + c(nm+k,j+1)*dif(1,3)/dif(i,k)/2.0
!           endif
!         enddo
!       enddo
!
!     else                                              ! boundary conditions at bulk
!       ! do i = 3,nm+1
!       !   g(nm+i-2)       = c(nm+i,j) - c(nm+i,j)       ! set concentration or potential c(8:11) - x_Cl, x_Na, x_O2, Phi
!       !   b(nm+i-2,nm+i)  = 1.0
!       ! enddo
!       g(3)   = c(3,j) - c(3,j)
!       b(3,3) = 1.0
!       ! print*, g(1:4)
!       ! print*, c(1,j), c(2,j), c(3,j)*ctot*1e3, c(4,j)*ctot*1e3
!     endif
    !
    ! if (j == 2) then
    !     print*, 'h', h
    !     print*, 'A'
    !     do i = 1, N
    !       write(*, matrix_fmt) A(i,:)
    !     end do
    !
    !     print*,
    !     print*, 'B'
    !     do i = 1, N
    !       write(*, matrix_fmt) B(i,:)
    !     end do
    !
    !     print*,
    !     print*, 'D'
    !     do i = 1, N
    !       write(*, matrix_fmt) D(i,1:N)
    !     end do
    !
    !     print*,
    !     print*, 'G'
    !     do i = 1, N
    !       write(*, matrix_fmt) G(i)
    !     end do
    !
    !     stop
    ! end if
!
    call band(j)

  end do

    nerr=0
    do j = 1,nj
      do i = 1,n
        if(dabs(c(i,j)).gt.1.d-10*dabs(cold (i,j)).and.dabs(c(i,j)).gt.1.d-18) nerr = nerr+1
      enddo
    enddo

    c = cold + c

! ----------------------------------------------------------------------------------------------------------------------------------
    print *, nerr,'fluxes ',(c(i,1),i=1,nm)
    ! if ((nerr.gt.0).and.(iteration.lt.2100)) go to 5
    ! print *, ' z    OH-    Cl-    Na+    O2'
    ! print 999, (h*dble(j-1), ctot*1000.*c(7,j), ctot*1000.*c(8,j), ctot*1000.*c(9,j), ctot*1000.*c(10,j), j=1,nj)
    print *, 'solution number ', iteration

    header_fmt = '(1(A5,2X), 1(A12,2X),   20(A15,2X)) '
    data_fmt   = '(1(I5,2X), 1(F12.5,2X), 20(ES15.8,2X))'


! ***** WRITE TO SCREEN **** !
    write(*, header_fmt) 'Iter', 'Position', 'C_O2', 'C_H2O', 'Flux_O2', 'Flux_H2O'
    write(*, header_fmt) '#',    'η',       'mol/L', 'mol/L', 'Flux_O2', 'Flux_H2O'
    do j = 1, NJ                                    !
      ! c_H2O   = ctot*1000.*c(6,j)
      ! c_OH    = ctot*1000.*c(7,j)
      ! c_Cl    = ctot*1000.*c(8,j)
      ! c_Na    = ctot*1000.*c(9,j)
      ! c_O2    = ctot*1000.*c(10,j)
      ! Phi     = c(11,j) * Rigc*Temp/Fconst
      ! Flux_H2O = c(1,j)
      ! Flux_OH  = c(2,j)
      ! Flux_Cl  = c(3,j)
      ! Flux_Na  = c(4,j)
      ! Flux_O2  = c(5,j)

      ! c_H2O   = ctot*1000.*c(3,j)
      ! c_O2    = ctot*1000.*c(4,j)
      ! Flux_H2O= c(1,j)
      Flux_O2 = c(1,j)
      c_O2    = c(2,j)
                                                    !
      write(*, data_fmt) iteration, h*dble(j-1), c_O2, Flux_O2!, c_OH, c_Na, c_Cl, Phi, c_H2O, Flux_H2O, &
                    ! & Flux_OH, Flux_Cl, Flux_Na, Flux_O2
    end do



! ***** WRITE TO FILE ****** !
    if (iteration == 1) then
      open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

      ! write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Iter', 'Position', 'C_O2', 'C_H2O', 'Flux_O2', 'Flux_H2O'
      write(56, header_fmt) '#',    'η',       'mol/L', 'mol/L', 'Flux_O2', 'Flux_H2O'
    end if

    do j = 1, NJ                                    !
      ! c_H2O   = ctot*1000.*c(6,j)
      ! c_OH    = ctot*1000.*c(7,j)
      ! c_Cl    = ctot*1000.*c(8,j)
      ! c_Na    = ctot*1000.*c(9,j)
      ! c_O2    = ctot*1000.*c(10,j)
      ! Phi     = c(11,j) * Rigc*Temp/Fconst
      ! Flux_H2O = c(1,j)
      ! Flux_OH  = c(2,j)
      ! Flux_Cl  = c(3,j)
      ! Flux_Na  = c(4,j)
      ! Flux_O2  = c(5,j)

      Flux_H2O= c(1,j)
      Flux_O2 = c(2,j)
      c_H2O   = ctot*1000.*c(3,j)
      c_O2    = ctot*1000.*c(4,j)
                                                    !
      write(56, data_fmt) iteration, h*dble(j-1), c_O2, c_H2O, Flux_O2, Flux_H2O
    end do

    if ((nerr.gt.0).and.(iteration.lt.10000)) go to 5
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
