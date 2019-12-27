      winapp
C      PROGRAM    ***   VWR_PURE.FOR   ***
C
C  To calculate the response of a vibrating wire in normal
C  and superfluid helium-3.
C
C     P       : PRESSURE in bars
C     TC      : P dependent transition temperature
C     TAB     : A-B transition temp
C     V       : molar volume in cm3 per mole
C     MSTAR   : effective mass RATIO for s/f state
C     F1      : Fermi liquid correction
C     CGAMMA  : coefficient of the specific heat C/nRT
C     A       : viscosity coefficient in eta=(AT^2+B)^-1
C     B       :    "           "          "       "
C     VISC(T) : FUNCTION  as above
C     K       : Stokes function K
C     K1      : Stokes function k'
C     RHO3    : helium-3 density
C     RHO     : normal fluid density below Tc
C     RHOWIR  : wire density
C     MALPHA  : MAGIC ALPHA
C     T       : TEMPERATURE IN mK
C     TR      : reduced TEMP below Tc
C
C    REQUIRES SUBROUTINES
C       YOSHID  ... Yoshida functions a la CHH
C       JENSEN  ... mfp corrections to stokes
C       STOKES  ... k and k'
C       VISCOE  ... Hall's viscosity
C       HE3PRE  ... pressure dependent props of He3
C       REDVIS  ... reduced viscosity from CHH
C        ALPHA  ... gets magic alpha
C         WIRE  ... gets wire properties
C
      REAL*8  P, A, B, MSTAR, V, TC, CGAMMA, K, K1, F0, PI, HBAR, T
      REAL*8 F1, DELTA, MFP, EOVL, ETATC, GAMMA, BETA, ZETA, DF1, DF2
      REAL*8 DIA, RAD, RHOWIR, RHO3, RHO3W, MALPHA,  TR, RHO
      REAL*8 MALPH, ETA, REDVIS, STEP, Y0, Y5, Y6, EC, TAB
      REAL*8 VISC
      LOGICAL APHASE
      CHARACTER*1 TITLE(80)
      INTEGER*4 ITEMP, I
      VISC(T)= 1./( A*T*T + B)
      PI=3.141592653589793D0 
      OPEN(10, FILE= 'HE3TW.OUT')

      WRITE(*,*)'Program VWR_PURE.FOR' 
      WRITE(*,97)
      WRITE(10,97)
   97 FORMAT(//1x,'Program to calculate the response of a vibrating',
     +' wire in normal and superfluid He-3'/1x
     +'Uses Carless, Hall and Hook form for viscosity corrections.',
     +/1x, 'Temperatures on Greywall scale.'/1x)
      WRITE(*,*)'Data output into file HE3TW.OUT in current directory'
      WRITE(*,*)
      WRITE(*,*)'Enter character string to identify output file; 80 max'
      READ(*,98) TITLE
   98 FORMAT(80A1)
      HBAR= 1.054D-34
      WRITE(*,*)'enter pressure in bar'
      READ(*,*) P
      WRITE(*,*)'enter resonant FREQUENCY'
      READ(*,*) F0
      WRITE(10,99) TITLE
   99 FORMAT(1X,80A1)

C
C  get wire parameters
C
      CALL WIRE(DIA,RHOWIR)
      WRITE(*,105) F0
      WRITE(10,105) F0
  105 FORMAT(1X,' RESONANT FREQUENCY = ', F6.1,' Hz')
C
C  get magic alpha
C
      CALL ALPHA( MALPHA )
C
C   get pressure dependent parameters
C
      CALL HE3PRE(P, TC, TAB, V, MSTAR, CGAMMA)
      F1=3.*(MSTAR-1.)
      CALL VISCOE(P, A, B)
      RHO3= 3.016/V
      ETATC=VISC(TC)
      WRITE(*,100) P, TC, TAB, ETATC, RHO3, V , F1, MALPHA
      WRITE(10,100) P, TC, TAB, ETATC, RHO3, V , F1, MALPHA
  100 FORMAT(/1X, ' PRESSURE IS ',F6.3,' BAR',/1X,
     + ' Tc = ',F7.3,' mK   TAB =',F7.3,' mK   viscosity at Tc = ',F9.5,
     + ' Pa s',/1x,
     + ' density of He-3 = ',F8.4,' g/cc        molar volume = ',  F7.3,
     + ' cm^3',/1x,
     + '              F1 = ', F7.3, '               magic alpha = ',
     +  F7.3,/1x)
C
C
C   from Fermi liquid theory, the ratio of the viscosity to the
C   mean free path is given by :
C        eta =1/5 n.pf.l
C   numerically, given the definition of pf= hbar(3pi^2n)^4/3
C
      EOVL= 0.2*(6.023E29/V)**(4./3.)*(3.*PI*PI)**(1./3.)*HBAR
C
C normal state calculation
C

      RHO3W=RHO3/RHOWIR
      RAD=DIA/2.E3
      WRITE(10,*)'FUDGE FACTOR FOR FINITE WIDTH CORRECTION IS 1.14' 
      WRITE(*,*)'FUDGE FACTOR FOR FINITE WIDTH CORRECTION IS 1.14'
      WRITE(10,*)'TEMP[mK]    WIDTH      SHIFT '
      WRITE(*,*)'TEMP[mK]    WIDTH      SHIFT '
      DO 10 ITEMP =200,0,-1
        T=ITEMP
        IF ( T.LT.TC ) T=TC
        DELTA=DSQRT( VISC(T)/(2.*PI*RHO3*1000*F0))
        MFP=VISC(T)/EOVL
        ZETA=0.579*MFP
        GAMMA=RAD/DELTA
C
C  Carless Hall and Hook form for beta
C
        BETA=ZETA/RAD*(1. + MALPHA*MFP/RAD)/(1. + MFP/RAD)
        CALL JENSEN( BETA, GAMMA, K, K1)
        DF2= F0 * RHO3W *K1*(1. - 1.14*RHO3W*K)
        DF1= F0 * RHO3W * K * 0.5 * (1. - 0.75*RHO3W*K)
        WRITE(*,110) T,DF2,DF1
        WRITE(10,110) T,DF2,DF1
  110   FORMAT(2X, F5.2,2(2X,F9.4))
        IF (T.LE.TC) GOTO 20
   10 CONTINUE
   20 CONTINUE
C
C   the superfluid state
C
      WRITE(10,120)
  120 FORMAT(/////1X,   ' THE SUPERFLUID STATE',/1X)
      WRITE(10,*)'TEMP[mK]     WIDTH        SHIFT      ETA/ETA(Tc) '
      WRITE(*,*)'TEMP[mK]     WIDTH        SHIFT      ETA/ETA(Tc) '
      TR=1.0
      STEP=0.02
      APHASE=.TRUE.
      DO 30 I=0,44
C
C   get Yoshida functions
C
        TR=TR-STEP
        T=TR*TC
        IF( APHASE .AND.T.LE.TAB)THEN
          APHASE=.FALSE.
          WRITE(*,*)'Now in B-phase - previously A-phase'
          WRITE(10,*)'Now in B-phase - previously A-phase'
        ENDIF
        CALL YOSHID(TR, Y0, Y5, Y6)
        RHO =RHO3*(1. + F1/3.)*Y0/(1. + F1*Y0/3.)
        RHO3W=RHO/RHOWIR
        EC=REDVIS(TR)
        ETA=ETATC*EC
        ZETA=0.5*Y5*ETA/EOVL
        MFP=ETA/(EOVL*Y6)
        MALPH=1.156*MALPHA/(Y5*Y6)
        DELTA=DSQRT( ETA/(2.*PI*RHO*1000*F0))
        GAMMA=RAD/DELTA
        BETA=ZETA/RAD*(1. + MALPH*MFP/RAD)/(1. + MFP/RAD)
        CALL JENSEN( BETA, GAMMA, K, K1)
        DF2= F0 * RHO3W *K1*(1. - 1.14*RHO3W*K)
        DF1= F0 * RHO3W * K * 0.5 * (1. - 0.75*RHO3W*K)
        WRITE(10,130) T,DF2,DF1, EC
        WRITE(*,130) T,DF2,DF1, EC
  130   FORMAT(2X, F6.3,2(2X,F11.6),2X,F9.4)
   30 CONTINUE
      CLOSE(10)
      WRITE(*,*)'Data output into file HE3TW.OUT in current directory'
      END









