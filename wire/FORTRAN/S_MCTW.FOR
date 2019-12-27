        winapp    
C      PROGRAM    ***   MCTW.FOR   ***
C
C  To calculate the response of a vibrating wire in the
C  mixing chamber of a dilution refrigerator
C
C     V       : molar volume in cm3 per mole
C     MSTAR   : effective mass RATIO
C     VISC(T) : FUNCTION  to give viscosity de Waele form
C     CONC(T) : FUNCTION  to give concentration
C     VMOL    : Molar volume  He3 units!!  cm^3
C     K       : Stokes function K
C     K1      : Stokes function k'
C     RHOWIR  : wire density
C     MALPHA  : MAGIC ALPHA
C     T     : TEMPERATURE in mK
C
C    REQUIRES SUBROUTINES
C       JENSEN  ... mfp corrections to stokes
C       STOKES  ... k and k'
C       DALPHA  ... gets magic alpha
C         WIRE  ... gets wire properties
      REAL*8 MSTAR, K, K1, F0, PI, HBAR, T
      REAL*8  DELTA, MFP, EOVL, GAMMA, BETA, ZETA, DF1, DF2
      REAL*8 DIA, RAD, RHOWIR, RHO3W, MALPHA
      REAL*8 CONC, VISC, VMOL, RHO
      INTEGER ITEMP, ISTART, ISTOP, ISTEP
      CHARACTER*1 TITLE(80)
      CONC(T)=0.066+0.5056*T**2-0.2488*T**3+18.22*T**4-74.22*T**5
      VISC(T)= 0.0000000277/T**2 + 0.00000034/T
C
C care: molar vol and density in units of cm^3 and g/cm^3
C
      VMOL(T)=27.58*(1. + 0.286*CONC(T))/CONC(T)
      RHO(T)=3.016*MSTAR/VMOL(T)
      PI=3.1415927
      MSTAR=2.46
      WRITE(*,*)'Program MCTW'
      WRITE(*,*)'Data output into file MCTW.OUT in current directory'
      WRITE(*,*)'Enter character string to identify output file; 80 max'
      READ(*,98)TITLE
   98 FORMAT(80A1)
      HBAR= 1.054D-34
      WRITE(*,*)'enter resonant FREQUENCY'
      READ(*,*) F0
      OPEN(10, FILE= 'MCTW.OUT')
      WRITE(10,99) TITLE
   99 FORMAT(1X,80A1)
C
C  get wire parameters
C
      WRITE(10,97)
   97 FORMAT(1x,'Program to calculate the response of a vibrating wire
     +in saturated'/1x,'dilute phase in the mixing chamber of a dilution
     + fridge'/1x,' This program assumes the de Waele form for the visco
     +sity'/1x)
      CALL WIRE(DIA,RHOWIR)
C
C  get magic alpha
C
      CALL DALPHA( MALPHA )
      WRITE(*,96) MALPHA
      WRITE(10,96) MALPHA
   96 FORMAT(1x,' Magic alpha = ',F5.2)
      RAD=DIA/2.E3
      WRITE(10,*)'Switching is 10x mean free path and fudge factor for f
     +inite width is 1.14'
      WRITE(*,105) F0
      WRITE(10,105) F0
  105 FORMAT(/1X,'RESONANT FREQUENCY = ', F6.1,' Hz')
      WRITE(10,*)
      WRITE(10,*)'TEMP[mK]    WIDTH      SHIFT '
      WRITE(*,*)'TEMP[mK]    WIDTH      SHIFT '
C
C  set up loop initial parameters
C
      ISTART=500
      ISTOP=200
      ISTEP=10
   10 CONTINUE
      DO 20 ITEMP =ISTART, ISTOP, -ISTEP
        T=ITEMP*0.0001
        DELTA=DSQRT( VISC(T)/(2.*PI*RHO(T)*1000*F0))
C
C   from Fermi liquid theory, the ratio of the viscosity to the
C   mean free path is given by :
C        eta =1/5 n.pf.l
C   numerically, given the definition of pf= hbar(3pi^2n)^4/3
C
        EOVL= 0.2*(6.023E23/(VMOL(T)*1.0E-6))**(4./3.)
     +   *(3.*PI*PI)**(1./3.)*HBAR
        MFP=VISC(T)/EOVL
        ZETA=0.579*MFP
        GAMMA=RAD/DELTA
C
C   Switching mean free path * 10
C
        BETA=ZETA/RAD*(1. + MALPHA*10.*MFP/RAD)/(1. + 10.*MFP/RAD)
        CALL JENSEN( BETA, GAMMA, K, K1)
        RHO3W=RHO(T)/RHOWIR
        DF2= F0 * RHO3W *K1*(1. - 1.14*RHO3W*K)
        DF1= F0 * RHO3W * K * 0.5 * (1. - 0.75*RHO3W*K)
        T=T*1000
        IF(ITEMP.EQ.100)THEN
          WRITE(10,*)
          WRITE(10,*)'TEMP[mK]    WIDTH      SHIFT '
        ENDIF
        WRITE(*,110) T,DF2,DF1
        WRITE(10,110) T,DF2,DF1
  110   FORMAT(2X, F5.2,2(2X,F9.4))
   20 CONTINUE
      IF (ISTEP.EQ.10) THEN
        ISTEP=5
        ISTART=195
        ISTOP=60
      ELSEIF (ISTEP.EQ.5) THEN
        ISTEP=1
        ISTART=59
        ISTOP=20
      ENDIF
      IF (ISTEP.NE.1.OR.T.GT.ISTOP/10)GOTO 10
      CLOSE(10)
      WRITE(*,*)'Data output into file MCTW.OUT in current directory'
      END


      SUBROUTINE JENSEN( BETA, GAMMA, KJ, KJ1)
      REAL*8 BETA, GAMMA, KJ, KJ1, K, K1, DENOM, G1
C
C   To calculate the mean free path corrections to the Stokes functions
C
C   calls subroutine STOKES
C
C   BETA    : Jensen beta parameter
C   GAMMA   : ratio of wire radius to viscous penetration depth
C   K, K1   : Stokes functions
C   KJ, KJ1 : modified functions
C
      CALL STOKES( GAMMA, K, K1 )
      G1=GAMMA**2/4.
      DENOM = (1. + G1 * BETA *K1)**2  +  (G1*BETA*(K-1.))**2
      KJ= 1. + (K-1.)/DENOM
      KJ1=(K1 +G1*BETA*((K-1.)**2 + K1**2))/DENOM
      RETURN
      END


      SUBROUTINE STOKES(G,K,K1)
C
C  evaluates the STOKES functions K and K'
C
C  uses the methods outlined in
C
C  STOKES--- Mathematical and physical papers Vol III.
C
C  for gamma>=3 use equations 113
C            <3               103-105
C
C   G    : gamma
C   K    : Stokes function K
C   K1   : Stokes function K'
C

      DOUBLE PRECISION G, K, K1
      DOUBLE PRECISION M0, M1, AL, E0, E1, G1
      DOUBLE PRECISION PI, A, B, C, D, E, G2, G4, G5, G6, G8, G10, G12
      PI=3.1415927
      IF (G. GE. 3.0) THEN
         G2=G*G
         G4=G2*G2
         G5=G4*G
         A=2.828427 / G
         B=.3535534 / G /G2
         C=.552427 / G5
         D=2.9637719 / G5 / G2
         E=20.4632383 / G5 / G4
         K = 1. + A + B - 0.5 / G4 + C - D + 12.8875857 /G4/G4 -E
         K1 = A + 2./G2 - B + C - 1.625/G4/G2 + D - E
       ELSE
         G1=G/2.0
         G2=G1*G1
         G4=G2*G2
         G6=G4*G2
         G8=G4*G4
         G10=G6*G4
         G12=G6*G6
         M0=G2 - G6/12. + G10/2880.
         M1=G2 - G6/36. + G10/14400.
         E0=G4/2. - G8/144. +G12/86400.
         E1=G4/4. - G8/576. +G12/518400.
         AL=.5772158 + DLOG(G1)
         A=-(AL*M0) + PI/4.*E0 - 0.5*M1 + (G2 - G6/6.545454 +G10
     +    /1261.313)
         B=PI/4.*M0 + AL*E0 - 0.5*(1.-E1) -
     +    (G4*0.75 - G8/69.12 + G12/35265.31)
         C=-(PI/4.*M1) + AL*(1.-E1) + (G4/2.666666 - G8/276.48 +
     +    G12/211591.84)
         D=-(AL*M1) - PI/4.*(1.-E1) + (G2 - G6/19.636363 +G10/6306.57)
         K=1. +2.*(A*C + B*D)/(G2*(C*C + D*D))
         K1= 2.*(B*C - A*D)/(G2*(C*C + D*D))
       ENDIF
      RETURN
      END

      SUBROUTINE DALPHA( MALPHA)
C
C  Gets magic alpha for dilute phase
C
      REAL*8 MALPHA
      CHARACTER ANSWER
      MALPHA=2.2
   10 WRITE(*,*)'Dilute phase MAGIC ALPHA assumed to be 2.2 Is this corr
     +ect? [y/n]'
      READ(*,100) ANSWER
  100 FORMAT(A1)
      IF (ANSWER.NE. 'Y'.AND. ANSWER .NE.'N'.AND. ANSWER.NE. 'y'
     + .AND. ANSWER .NE.'n')GOTO 10
      IF (ANSWER .EQ. 'Y'.OR. ANSWER .EQ. 'y')GOTO 20
        WRITE(*,*)'Enter new value'
        READ(*,*) MALPHA
   20 CONTINUE
      RETURN
      END




      SUBROUTINE WIRE(DIA, RHOWIR)
C
C gets the properties of the vibrating wire
C
      REAL*8 DIA, RHOWIR
      INTEGER N
      CHARACTER*10 METAL
   20   WRITE(*,*)'Enter material: 1 = TANTALUM; 2 = NbTi; 3 = other'
        READ(*,*) N
        IF(N.EQ.1)THEN
          METAL='TANTALUM'
          RHOWIR=16.7
        ELSEIF(N.EQ.2) THEN
          METAL='NbTi'
          RHOWIR=6.05
        ELSEIF(N.EQ.3) THEN
          WRITE(*,*)'INPUT THE NAME OF THE MATERIAL'
          READ(*,101)METAL
  101     FORMAT(A10)
          WRITE(*,*)'INPUT THE DENSITY in g/cm^3'
          READ(*,*)RHOWIR
        ELSEIF(N.NE.1 .AND. N.NE.2 .AND. N.NE.3)THEN
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'
          GOTO 20
        ENDIF
C
C   get wire diameter
C
       WRITE(*,*) 'ENTER wire DIAMETER in mm: 1= 0.124 mm, 2= OTHER'
       READ (*,*) N
       IF (N.NE.1 .AND. N.NE.2)THEN
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'
       ELSEIF(N.EQ.1)THEN
          DIA=0.124
       ELSEIF(N.EQ.2)THEN
          WRITE(*,*)'input diameter in mm'
          READ(*,*) DIA
       ENDIF
      WRITE(*,120)METAL ,RHOWIR, DIA
      WRITE(10,120)METAL ,RHOWIR, DIA
  120 FORMAT(1X,'Wire properties: material is ', A10,
     +' of density ', F8.4,' g/cm^3'
     +/1X,' and diameter ', F7.5,' mm'/1x)
      RETURN
      END




