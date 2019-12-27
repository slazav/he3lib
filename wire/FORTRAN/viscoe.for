
      SUBROUTINE VISCOE(P,A,B)
c
c  subroutine to calculate the viscosity coefficients
c
c      eta= 1/( AT^2 + B)
c
c  ( with eta in Pa s, T in mK.)
c  using the data given in Carless, Hall and Hook; JLTP 50,583,83
c  table 1 , page 593 as smoothed by AMG.
c  The data, converted to a Greywall T scale
c  by multiplication of T by 1.12.
c
      REAL*8 P, A, B
      IF (P.LE.1.28)THEN
       A=0.38-.007*(1.28-P)/1.18
       B=0.06-0.02*(1.28-P)/1.18
      ENDIF
      IF (P.GT. 1.28 .AND. P.LE. 4.65) THEN
       A=0.424-.044*(4.65-P)/3.37
       B=0.19-0.13*(4.65-P)/3.37
      ENDIF
      IF (P.GT. 4.65 .AND. P.LE. 9.89) THEN
       A=0.495-.071*(9.89-P)/5.24
       B=0.43-0.24*(9.89-P)/5.24
      ENDIF
      IF (P.GT. 9.89 .AND. P.LE.19.89) THEN
       A=0.603-.108*(19.89-P)/10
       B=0.94-0.56*(19.89-P)/10
      ENDIF
      IF (P.GT.19.89)THEN
       A=0.603+0.107*(P-19.89)/9.45
       B=0.94+0.56*(P-19.89)/9.45
      ENDIF
      A=A*10*1.12*1.12
      B=B*10
      RETURN
      END