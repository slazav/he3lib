      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EPS    = 5.D-1
      OLDEPS = 1.0D0
10    EPS1   = 1.D0 + EPS
      IF (EPS1 .EQ. 1.D0) GO TO 20
      OLDEPS = EPS
      EPS    = EPS * 5.D-1
      GO TO 10
20    WRITE (*,*) ' D1MACH(3) = ', OLDEPS
      STOP
      END
