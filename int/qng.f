       INTEGER LIMIT, LENW
       PARAMETER (LIMIT = 2000)
       PARAMETER (LENW = 4*LIMIT)
       REAL*8 A,ABSERR,B,F,EPSABS,EPSREL,RESULT,WORK(LENW)
       INTEGER IER,NEVAL,LAST,IWORK(LIMIT)
       EXTERNAL F
       A = 0.0D0
       B = 1.0D0
       EPSABS = 0.0D0
       EPSREL = 1.0D-3
       CALL DQAGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER,
     .      LIMIT,LENW,LAST,IWORK,WORK)

       write(6,1) result
 1     format(' result=',1pe15.7)
       STOP
       END

       REAL*8 FUNCTION F(X)
       REAL*8 X
       F = DEXP(X)/(X*X+0.1D+01)
       RETURN
       END
