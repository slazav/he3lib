
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