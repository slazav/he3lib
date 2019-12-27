
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
      PI=3.141592653589793D0
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