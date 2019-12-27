 %   PROGRAM    ***   HE3TW.FOR   ***\

%  To calculate the response of a vibrating wire in normal\
%  and superfluid helium-3.\
%
%     P       : PRESSURE in bars\
%     TC      : P dependent transition temperature\
%     TAB     : A-B transition temp\
%     V       : molar volume in cm3 per mole\
%     MSTAR   : effective mass RATIO for s/f state\
%     F1      : Fermi liquid correction\
%     CGAMMA  : coefficient of the specific heat C/nRT\
%     A       : viscosity coefficient in eta=(AT^2+B)^-1\
%     B       :    "           "          "       "\
%     VISC(T) : FUNCTION  as above\
%     K       : Stokes function K\
%     K1      : Stokes function k'\
%     RHO3    : helium-3 density\
%     RHO     : normal fluid density below Tc\
%     RHOWIR  : wire density\
%     MALPHA  : MAGIC ALPHA\
%     T       : TEMPERATURE IN mK\
%     TR      : reduced TEMP below Tc\
 
%   REQUIRES SUBROUTINES\
%      YOSHID  ... Yoshida functions a la CHH\
%       JENSEN  ... mfp corrections to stokes\
%       STOKES  ... k and k'\
%       VISCOE  ... Hall's viscosity\
%       HE3PRE  ... pressure dependent props of He3\
%       REDVIS  ... reduced viscosity from CHH\
%        ALPHA  ... gets magic alpha\
%         WIRE  ... gets wire properties\
 
 
      REAL*8  P, A, B, MSTAR, V, TC, CGAMMA, K, K1, F0, PI, HBAR, T\
      REAL*8 F1, DELTA, MFP, EOVL, ETATC, GAMMA, BETA, ZETA, DF1, DF2\
      REAL*8 DIA, RAD, RHOWIR, RHO3, RHO3W, MALPHA, TR, RHO\
      REAL*8 MALPH, ETA, REDVIS, STEP, Y0, Y5, Y6, EC, TAB\
      REAL*8 VISC\
      LOGICAL APHASE\
      CHARACTER*1 TITLE(80)\
      INTEGER*4 ITEMP, I\
      VISC(T)= 1./( A*T*T + B)\
      PI=3.1415926\
      WRITE(*,*)'Program HE3TW'\
      WRITE(*,*)'Data output into file HE3TW.OUT in current directory'\
      WRITE(*,*)'Enter character string to identify output file; 80 max'\
      READ(*,98) TITLE\
   98 FORMAT(80A1)\
      HBAR= 1.054D-34\
      WRITE(*,*)'enter pressure in bar'\
      READ(*,*) P\
      WRITE(*,*)'enter resonant FREQUENCY'\
      READ(*,*) F0\
      OPEN(10, FILE= 'HE3TW.OUT')\
      WRITE(10,99) TITLE\
   99 FORMAT(1X,80A1)\
\
C\
C  get wire parameters\
C\
      CALL WIRE(DIA,RHOWIR)\
      WRITE(*,105) F0\
      WRITE(10,105) F0\
  105 FORMAT(1X,' RESONANT FREQUENCY = ', F6.1,' Hz')\
C\
C  get magic alpha\
C\
      CALL ALPHA( MALPHA )\
C\
C   get pressure dependent parameters\
C\
      CALL HE3PRE(P, TC, TAB, V, MSTAR, CGAMMA)\
      F1=3.*(MSTAR-1.)\
      CALL VISCOE(P, A, B)\
      RHO3= 3.016/V\
      ETATC=VISC(TC)\
      WRITE(*,100) P, TC, TAB, ETATC, RHO3, V , F1, MALPHA\
      WRITE(10,100) P, TC, TAB, ETATC, RHO3, V , F1, MALPHA\
  100 FORMAT(/1X, ' PRESSURE IS ',F6.3,' BAR',/1X,\
     + ' Tc = ',F7.3,' mK   TAB =',F7.3,' mK   viscosity at Tc = ',F9.5,\
     + ' Pa s',/1x,\
     + ' density of He-3 = ',F8.4,' g/cc        molar volume = ',  F7.3,\
     + ' cm^3',/1x,\
     + '              F1 = ', F7.3, '               magic alpha = ',\
     +  F7.3,/1x)\
C\
C\
C   from Fermi liquid theory, the ratio of the viscosity to the\
C   mean free path is given by :\
C        eta =1/5 n.pf.l\
C   numerically, given the definition of pf= hbar(3pi^2n)^4/3\
C\
      EOVL= 0.2*(6.023E29/V)**(4./3.)*(3.*PI*PI)**(1./3.)*HBAR\
C\
C normal state calculation\
C\
\
      RHO3W=RHO3/RHOWIR\
      RAD=DIA/2.E3\
      WRITE(10,*)'FUDGE FACTOR FOR FINITE WIDTH CORRECTION IS 1.14'\
      WRITE(10,*)'TEMP[mK]    WIDTH      SHIFT '\
      WRITE(*,*)'TEMP[mK]    WIDTH      SHIFT '\
      DO 10 ITEMP =50,0,-1\
        T=ITEMP\
        IF ( T.LT.TC ) T=TC\
        DELTA=DSQRT( VISC(T)/(2.*PI*RHO3*1000*F0))\
        MFP=VISC(T)/EOVL\
        ZETA=0.579*MFP\
        GAMMA=RAD/DELTA\
C\
C  Carless Hall and Hook form for beta\
C\
        BETA=ZETA/RAD*(1. + MALPHA*MFP/RAD)/(1. + MFP/RAD)\
        CALL JENSEN( BETA, GAMMA, K, K1)\
        DF2= F0 * RHO3W *K1*(1. - 1.14*RHO3W*K)\
        DF1= F0 * RHO3W * K * 0.5 * (1. - 0.75*RHO3W*K)\
        WRITE(*,110) T,DF2,DF1\
        WRITE(10,110) T,DF2,DF1\
  110   FORMAT(2X, F5.2,2(2X,F9.4))\
        IF (T.LE.TC) GOTO 20\
   10 CONTINUE\
   20 CONTINUE\
C\
C   the superfluid state\
C\
      WRITE(10,120)\
  120 FORMAT(/////1X,   ' THE SUPERFLUID STATE',/1X)\
      WRITE(10,*)'TEMP[mK]     WIDTH        SHIFT      ETA/ETA(Tc) '\
      WRITE(*,*)'TEMP[mK]     WIDTH        SHIFT      ETA/ETA(Tc) '\
      TR=1.0\
      STEP=0.02\
      APHASE=.TRUE.\
      DO 30 I=0,44\
C\
C   get Yoshida functions\
C\
        TR=TR-STEP\
        T=TR*TC\
        IF( APHASE .AND.T.LE.TAB)THEN\
          APHASE=.FALSE.\
          WRITE(*,*)'Now in B-phase - previously A-phase'\
          WRITE(10,*)'Now in B-phase - previously A-phase'\
        ENDIF\
        CALL YOSHID(TR, Y0, Y5, Y6)\
        RHO =RHO3*(1. + F1/3.)*Y0/(1. + F1*Y0/3.)\
        RHO3W=RHO/RHOWIR\
        EC=REDVIS(TR)\
        ETA=ETATC*EC\
        ZETA=0.5*Y5*ETA/EOVL\
        MFP=ETA/(EOVL*Y6)\
        MALPH=1.156*MALPHA/(Y5*Y6)\
        DELTA=DSQRT( ETA/(2.*PI*RHO*1000*F0))\
        GAMMA=RAD/DELTA\
        BETA=ZETA/RAD*(1. + MALPH*MFP/RAD)/(1. + MFP/RAD)\
        CALL JENSEN( BETA, GAMMA, K, K1)\
        DF2= F0 * RHO3W *K1*(1. - 1.14*RHO3W*K)\
        DF1= F0 * RHO3W * K * 0.5 * (1. - 0.75*RHO3W*K)\
        WRITE(10,130) T,DF2,DF1, EC\
        WRITE(*,130) T,DF2,DF1, EC\
  130   FORMAT(2X, F6.3,2(2X,F11.6),2X,F9.4)\
   30 CONTINUE\
      CLOSE(10)\
      WRITE(*,*)'Data output into file HE3TW.OUT in current directory'\
      END\
\
\
\
      SUBROUTINE STOKES(G,K,K1)\
C\
C  evaluates the STOKES functions K and K'\
C\
C  uses the methods outlined in\
C\
C  STOKES--- Mathematical and physical papers Vol III.\
C\
C  for gamma>=3 use equations 113\
C            <3               103-105\
C\
C   G    : gamma\
C   K    : Stokes function K\
C   K1   : Stokes function K'\
C\
C *** requires the following NAG subroutines\
C\
C    XO1AAF - returns the valur of PI\
C\
      DOUBLE PRECISION G, K, K1\
      DOUBLE PRECISION M0, M1, AL, E0, E1, G1\
      DOUBLE PRECISION PI, A, B, C, D, E, G2, G4, G5, G6, G8, G10, G12\
      PI=3.1415926\
      IF (G. GE. 3.0) THEN\
         G2=G*G\
         G4=G2*G2\
         G5=G4*G\
         A=2.828427 / G\
         B=.3535534 / G /G2\
         C=.552427 / G5\
         D=2.9637719 / G5 / G2\
         E=20.4632383 / G5 / G4\
         K = 1. + A + B - 0.5 / G4 + C - D + 12.8875857 /G4/G4 -E\
         K1 = A + 2./G2 - B + C - 1.625/G4/G2 + D - E\
       ELSE\
         G1=G/2.0\
         G2=G1*G1\
         G4=G2*G2\
         G6=G4*G2\
         G8=G4*G4\
         G10=G6*G4\
         G12=G6*G6\
         M0=G2 - G6/12. + G10/2880.\
         M1=G2 - G6/36. + G10/14400.\
         E0=G4/2. - G8/144. +G12/86400.\
         E1=G4/4. - G8/576. +G12/518400.\
         AL=.5772158 + DLOG(G1)\
         A=-(AL*M0) + PI/4.*E0 - 0.5*M1 + (G2 - G6/6.545454 +G10\
     +    /1261.313)\
         B=PI/4.*M0 + AL*E0 - 0.5*(1.-E1) -\
     +    (G4*0.75 - G8/69.12 + G12/35265.31)\
         C=-(PI/4.*M1) + AL*(1.-E1) + (G4/2.666666 - G8/276.48 +\
     +    G12/211591.84)\
         D=-(AL*M1) - PI/4.*(1.-E1) + (G2 - G6/19.636363 +G10/6306.57)\
         K=1. +2.*(A*C + B*D)/(G2*(C*C + D*D))\
         K1= 2.*(B*C - A*D)/(G2*(C*C + D*D))\
       ENDIF\
      RETURN\
      END\
\
      SUBROUTINE JENSEN( BETA, GAMMA, KJ, KJ1)\
      REAL*8 BETA, GAMMA, KJ, KJ1, K, K1, DENOM, G1\
C\
C   To calculate the mean free path corrections to the Stokes functions\
C\
C   calls subroutine STOKES\
C\
C   BETA    : Jensen beta parameter\
C   GAMMA   : ratio of wire radius to viscous penetration depth\
C   K, K1   : Stokes functions\
C   KJ, KJ1 : modified functions\
C\
      CALL STOKES( GAMMA, K, K1 )\
      G1=GAMMA**2/4.\
      DENOM = (1. + G1 * BETA *K1)**2  +  (G1*BETA*(K-1.))**2\
      KJ= 1. + (K-1.)/DENOM\
      KJ1=(K1 +G1*BETA*((K-1.)**2 + K1**2))/DENOM\
      RETURN\
      END\
\
C\
C    Yoshida functions as given in Carless et al JLTP 50,605,1983\
C\
      SUBROUTINE YOSHID(T, Y0, Y5, Y6)\
C\
C  Y0    : Y0\
C  Y5    : combination CHH eq A6\
C  Y6    :                    A7\
C  T     : reduced temp T/Tc\
C  TS    : SCALED reduced temp\
C\
      REAL*8 T, TS, Y0, Y5, Y6\
      TS=T*(0.9074 - 0.0075*T - 0.0216*T*T + 0.1396*T**3 - 0.0611*T**4)\
      IF(T.LE.0.94)THEN\
          Y0=DEXP(-1.76388/TS)*(3.454 - 0.88*TS + 4.625*TS*TS -\
     +      1.367*TS**3)/DSQRT(TS)\
        ELSE\
          Y0=1.985*TS - 0.985\
      ENDIF\
      IF( T .LE. 0.8)THEN\
          Y5=DEXP(1.76388/TS)*(0.10177 + 1.1958*TS - 1.425*TS*TS +\
     +      0.392*TS**3)/DSQRT(TS)\
        ELSE\
          Y5=DEXP(1.76388/TS)*(0.19847 +0.335*DSQRT(1.- TS))/DSQRT(TS)\
      ENDIF\
      IF( T .LE. 0.9)THEN\
          Y6=DEXP(-1.76388/TS)*(2.402 + 0.4467*TS - 2.117*TS*TS +\
     +      4.1*TS**3)\
        ELSE\
          Y6=1. - DSQRT(1.- TS)*(-4.517 + 13.275*TS -7.5*TS*TS)\
      ENDIF\
      RETURN\
      END\
\
\
\
      SUBROUTINE VISCOE(P,A,B)\
c\
c  subroutine to calculate the viscosity coefficients\
c\
c      eta= 1/( AT^2 + B)\
c\
c  ( with eta in Pa s, T in mK.)\
c  using the data given in Carless, Hall and Hook; JLTP 50,583,83\
c  table 1 , page 593 as smoothed by AMG.\
c  The data, converted to a Greywall T scale\
c  by multiplication of T by 1.12.\
c\
      REAL*8 P, A, B\
      IF (P.LE.1.28)THEN\
       A=0.38-.007*(1.28-P)/1.18\
       B=0.06-0.02*(1.28-P)/1.18\
      ENDIF\
      IF (P.GT. 1.28 .AND. P.LE. 4.65) THEN\
       A=0.424-.044*(4.65-P)/3.37\
       B=0.19-0.13*(4.65-P)/3.37\
      ENDIF\
      IF (P.GT. 4.65 .AND. P.LE. 9.89) THEN\
       A=0.495-.071*(9.89-P)/5.24\
       B=0.43-0.24*(9.89-P)/5.24\
      ENDIF\
      IF (P.GT. 9.89 .AND. P.LE.19.89) THEN\
       A=0.603-.108*(19.89-P)/10\
       B=0.94-0.56*(19.89-P)/10\
      ENDIF\
      IF (P.GT.19.89)THEN\
       A=0.603+0.107*(P-19.89)/9.45\
       B=0.94+0.56*(P-19.89)/9.45\
      ENDIF\
      A=A*10*1.12*1.12\
      B=B*10\
      RETURN\
      END\
\
      SUBROUTINE HE3PRE( P, TC, TAB, V, MSTAR, GAMMA )\
c  subroutine to give  various physical parameters for\
c  Helium-3 which depend on the pressure\
c\
c  no external subroutines required\
c\
c  P        : Pressure in BARS\
c  TC       : Transition temperature in mK\
c  TAB      : A-B transition temperature\
c  V        : molar volume in cm^3 per mole\
c  MSTAR    : effective mass ratio\
c  GAMMA    : C/nRT coefficient of the specific heat\
c\
      REAL*8 P, PINT, TC, V, MSTAR, GAMMA, PFERMI, TAB, PPCP\
      REAL*8 A(6), T(6) ,G(6)\
      INTEGER*4 I\
c  coefficients for the molar volume\
      A(1)=36.837231\
      A(2)=-0.11803474D01\
      A(3)=0.83421417D-1\
      A(4)=-0.38859562D-2\
      A(5)=0.94759780D-4\
      A(6)=-0.91253577D-6\
c  coefficient for Tc\
      T(1)=0.92938375\
      T(2)=0.13867188\
      T(3)=-0.69302185D-2\
      T(4)=0.25685169D-3\
      T(5)=-0.57248644D-5\
      T(6)=0.53010918D-7\
c  coefficients for gamma\
      G(1)=0.27840464D1\
      G(2)=0.69575243D-1\
      G(3)=-0.14738303D-2\
      G(4)=0.46153498D-4\
      G(5)=-0.53785385D-6\
      G(6)=0.0D1\
c\
c  set up initial constants for the summations\
c\
      TC = 0.0D0\
      V = 0.0D0\
      GAMMA = 0.0D0\
      PINT=1.0D0\
c  the loop\
      DO 10 I=1,6\
      V = V + A(I)*PINT\
      TC = TC + T(I)*PINT\
      GAMMA = GAMMA +G(I)*PINT\
      PINT = PINT*P\
   10 CONTINUE\
      PFERMI = 2.7551D-19/V**(.3333333333)\
      MSTAR=3.06413D-18*GAMMA/( V *PFERMI )\
c  A-B transition temp calculation\
      PPCP=21.22\
      TAB=0.0D0\
      IF (P.LE.PPCP)THEN\
        TAB=TC\
      ELSE\
        A(1)=-0.10322623D-1\
        A(2)=-0.53633181D-2\
        A(3)=0.83437032D-3\
        A(4)=-0.61709783D-4\
        A(5)=0.17038992D-5\
        DO 20 I=1,5\
         TAB=TAB+A(I)*(P-PPCP)**I\
   20   CONTINUE\
        TAB=TAB+2.273D0\
      ENDIF\
      RETURN\
      END\
\
      REAL*8 FUNCTION REDVIS(TR)\
C\
C  a plot of CHH  reduced viscosity gives the following reasonable\
C  form.   Basically, LOG( redvis) viz sqrt(1-T)^.5 reasonably slow\
C\
\
      REAL*8 TR, SLOPE, N0, D0\
      IF(TR.LE.0.6)THEN\
       REDVIS=0.11\
       GOTO 10\
      ELSEIF(TR.LE.0.7)THEN\
       SLOPE=-0.8562\
       D0=0.4\
       N0=-0.9586\
      ELSEIF(TR.LE.0.8)THEN\
       SLOPE=-0.6183\
       D0=0.3\
       N0=-0.8861\
      ELSEIF(TR.LE.0.9)THEN\
       SLOPE=-1.4172\
       D0=0.2\
       N0=-0.8239\
      ELSEIF(TR.LE.0.95)THEN\
       SLOPE=-1.7352\
       D0=0.1\
       N0=-0.6383\
      ELSEIF(TR.LE.0.975)THEN\
       SLOPE=-1.6177\
       D0=0.05\
       N0=-0.4776\
      ELSEIF(TR.LE.1.)THEN\
       SLOPE=-2.3503\
       D0=0.025\
       N0=-0.3716\
      ENDIF\
      REDVIS=10**(N0+SLOPE*(DSQRT(1.-TR)-DSQRT(D0)))\
   10 CONTINUE\
      END\
\
\
      SUBROUTINE ALPHA( MALPHA)\
C\
C  Gets magic alpha\
C\
      REAL*8 MALPHA\
      CHARACTER ANSWER\
      MALPHA=1.9\
   10 WRITE(*,*)'MAGIC ALPHA assumed to be 1.9 Is this correct? [y/n]'\
      READ(*,100) ANSWER\
  100 FORMAT(A1)\
      IF (ANSWER.NE. 'Y'.AND. ANSWER .NE.'N'.AND. ANSWER.NE. 'y'\
     + .AND. ANSWER .NE.'n')GOTO 10\
      IF (ANSWER .EQ. 'Y'.OR. ANSWER .EQ. 'y')GOTO 20\
        WRITE(*,*)'Enter new value'\
        READ(*,*) MALPHA\
   20 CONTINUE\
      RETURN\
      END\
\
\
\
\
      SUBROUTINE WIRE(DIA, RHOWIR)\
C\
C gets the properties of the vibrating wire\
C\
      REAL*8 DIA, RHOWIR\
      INTEGER N\
      CHARACTER*10 METAL\
   20   WRITE(*,*)'Enter material: 1 = TANTALUM; 2 = NbTi; 3 = other'\
        READ(*,*) N\
        IF(N.EQ.1)THEN\
          METAL='TANTALUM'\
          RHOWIR=16.7\
        ELSEIF(N.EQ.2) THEN\
          METAL='NbTi'\
          RHOWIR=6.05\
        ELSEIF(N.EQ.3) THEN\
          WRITE(*,*)'INPUT THE NAME OF THE MATERIAL'\
          READ(*,101)METAL\
  101     FORMAT(A10)\
          WRITE(*,*)'INPUT THE DENSITY in g/cm^3'\
          READ(*,*)RHOWIR\
        ELSEIF(N.NE.1 .AND. N.NE.2 .AND. N.NE.3)THEN\
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'\
          GOTO 20\
        ENDIF\
C\
C   get wire diameter\
C\
       WRITE(*,*) 'ENTER wire DIAMETER in mm: 1= 0.124 mm, 2= OTHER'\
       READ (*,*) N\
       IF (N.NE.1 .AND. N.NE.2)THEN\
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'\
       ELSEIF(N.EQ.1)THEN\
          DIA=0.124\
       ELSEIF(N.EQ.2)THEN\
          WRITE(*,*)'input diameter in mm'\
          READ(*,*) DIA\
       ENDIF\
      WRITE(*,120)METAL ,RHOWIR, DIA\
      WRITE(10,120)METAL ,RHOWIR, DIA\
  120 FORMAT(1X,'Wire properties: material is ', A10,\
     +' of density ', F8.4,' g/cm^3'\
     +/1X,' and diameter ', F7.5,' mm'/1x)\
      RETURN\
      END\
\
\
\
* -----------------------------------------------------------------\
 \
\
\

\par }
 