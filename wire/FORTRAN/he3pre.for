
      SUBROUTINE HE3PRE( P, TC, TAB, V, MSTAR, GAMMA )
c  subroutine to give  various physical parameters for
c  Helium-3 which depend on the pressure
c
c  no external subroutines required
c
c  P        : Pressure in BARS
c  TC       : Transition temperature in mK
c  TAB      : A-B transition temperature
c  V        : molar volume in cm^3 per mole
c  MSTAR    : effective mass ratio
c  GAMMA    : C/nRT coefficient of the specific heat
c
      REAL*8 P, PINT, TC, V, MSTAR, GAMMA, PFERMI, TAB, PPCP
      REAL*8 A(6), T(6) ,G(6)
      INTEGER*4 I
c  coefficients for the molar volume
      A(1)=36.837231
      A(2)=-0.11803474D01
      A(3)=0.83421417D-1
      A(4)=-0.38859562D-2
      A(5)=0.94759780D-4
      A(6)=-0.91253577D-6
c  coefficient for Tc
      T(1)=0.92938375
      T(2)=0.13867188
      T(3)=-0.69302185D-2
      T(4)=0.25685169D-3
      T(5)=-0.57248644D-5
      T(6)=0.53010918D-7
c  coefficients for gamma
      G(1)=0.27840464D1
      G(2)=0.69575243D-1
      G(3)=-0.14738303D-2
      G(4)=0.46153498D-4
      G(5)=-0.53785385D-6
      G(6)=0.0D1
c
c  set up initial constants for the summations
c
      TC = 0.0D0
      V = 0.0D0
      GAMMA = 0.0D0
      PINT=1.0D0
c  the loop
      DO 10 I=1,6
      V = V + A(I)*PINT
      TC = TC + T(I)*PINT
      GAMMA = GAMMA +G(I)*PINT
      PINT = PINT*P
   10 CONTINUE
      PFERMI = 2.7551D-19/V**(.3333333333)
      MSTAR=3.06413D-18*GAMMA/( V *PFERMI )
c  A-B transition temp calculation
      PPCP=21.22
      TAB=0.0D0
      IF (P.LE.PPCP)THEN
        TAB=TC
      ELSE
        A(1)=-0.10322623D-1
        A(2)=-0.53633181D-2
        A(3)=0.83437032D-3
        A(4)=-0.61709783D-4
        A(5)=0.17038992D-5
        DO 20 I=1,5
         TAB=TAB+A(I)*(P-PPCP)**I
   20   CONTINUE
        TAB=TAB+2.273D0
      ENDIF
      RETURN
      END