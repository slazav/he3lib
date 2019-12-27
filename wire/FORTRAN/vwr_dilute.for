      winapp     
     
                     
C     PROGRAM    ***   VWR_Dilute.FOR   ***
C
C  To calculate the response of a vibrating wire in the
C  mixing chamber of a dilution refrigerator
C
C     V       : molar volume in cm3 per mole
C     MSTAR   : effective mass RATIO
C     VISC(T) : FUNCTION  to give viscosity de Waele form
C     CONC(T) : FUNCTION  to give concentration
C     VMOL    : Molar volume  He3 units!!  cm^3
C     K       : Stokes function k
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
C     
      REAL*8 MSTAR, K, K1, F0, PI, HBAR, T
      REAL*8  DELTA, MFP, EOVL, GAMMA, BETA, ZETA, DF1, DF2
      REAL*8 DIA, RAD, RHOWIR, RHO3W, MALPHA
      REAL*8 CONC, VISC, VMOL, RHO
      INTEGER ITEMP, ISTART, ISTOP, ISTEP
      CHARACTER*1 TITLE(80)

      CONC(T)=0.066+0.5056*T**2-0.2488*T**3+18.22*T**4-74.22*T**5
      VISC(T)= 0.000000027/T**2 + 0.00000034/T
C
C care: molar vol and density in units of cm^3 and g/cm^3
C
      VMOL(T)=27.58*(1. + 0.286*CONC(T))/CONC(T)
      RHO(T)=3.016*MSTAR/VMOL(T)
      PI=3.141592653589793D0
      MSTAR=2.46
      WRITE(*,*)'Program VWR_Dilute_Width'    
      WRITE(*,97)             
      
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
   97 FORMAT(//1x,'Program to calculate the response of a vibrating',
     +' wire in saturated'/1x,'dilute phase in the mixing chamber of',
     +' a dilution fridge',
     + //1x,'This program assumes the de Waele form for the viscosity'
     +/1x)
      CALL WIRE(DIA,RHOWIR)
C
C  get magic alpha
C
      CALL DALPHA( MALPHA )
      WRITE(*,96) MALPHA
      WRITE(10,96) MALPHA
   96 FORMAT(1x,' Magic alpha = ',F5.2)
      RAD=DIA/2.E3  
      WRITE(*,*)'Switching is 10x mean free path and fudge factor for f
     +inite width is 1.14'      
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
