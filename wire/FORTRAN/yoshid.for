C      
C    Yoshida functions as given in Carless et al JLTP 50,605,1983
C
      SUBROUTINE YOSHID(T, Y0, Y5, Y6)
C
C  Y0    : Y0
C  Y5    : combination CHH eq A6
C  Y6    :                    A7
C  T     : reduced temp T/Tc
C  TS    : SCALED reduced temp
C
      REAL*8 T, TS, Y0, Y5, Y6
      TS=T*(0.9074 - 0.0075*T - 0.0216*T*T + 0.1396*T**3 - 0.0611*T**4)
      IF(T.LE.0.94)THEN
          Y0=DEXP(-1.76388/TS)*(3.454 - 0.88*TS + 4.625*TS*TS -
     +      1.367*TS**3)/DSQRT(TS)
        ELSE
          Y0=1.985*TS - 0.985
      ENDIF
      IF( T .LE. 0.8)THEN
          Y5=DEXP(1.76388/TS)*(0.10177 + 1.1958*TS - 1.425*TS*TS +
     +      0.392*TS**3)/DSQRT(TS)
        ELSE
          Y5=DEXP(1.76388/TS)*(0.19847 +0.335*DSQRT(1.- TS))/DSQRT(TS)
      ENDIF
      IF( T .LE. 0.9)THEN
          Y6=DEXP(-1.76388/TS)*(2.402 + 0.4467*TS - 2.117*TS*TS +
     +      4.1*TS**3)
        ELSE
          Y6=1. - DSQRT(1.- TS)*(-4.517 + 13.275*TS -7.5*TS*TS)
      ENDIF
      RETURN
      END

