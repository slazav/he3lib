        winapp      
*       plot CHH viscosity
        real*8 eta, t
        integer i
        
        open (333, file="DIBfitvisc.txt")
        
        do 222 i=1,100
        t = i*.01
        eta=redvis(t)
        write(*,444) t,eta
        write(333,444) t,eta 
              
 444    format(2F10.5)     




 222    Continue
        end



      REAL*8 FUNCTION REDVIS(TR)
C
C  a plot of CHH  reduced viscosity gives the following reasonable
C  form.   Basically, LOG( redvis) viz sqrt(1-T)^.5 reasonably slow
C
  
      REAL*8 TR, SLOPE, N0, D0
      IF(TR.LE.0.6)THEN
       REDVIS=0.11
       GOTO 10
      ELSEIF(TR.LE.0.7)THEN
       SLOPE=-0.8562
       D0=0.4
       N0=-0.9586
      ELSEIF(TR.LE.0.8)THEN
       SLOPE=-0.6183
       D0=0.3
       N0=-0.8861
      ELSEIF(TR.LE.0.9)THEN
       SLOPE=-1.4172
       D0=0.2
       N0=-0.8239
      ELSEIF(TR.LE.0.95)THEN
       SLOPE=-1.7352
       D0=0.1
       N0=-0.6383
      ELSEIF(TR.LE.0.975)THEN
       SLOPE=-1.6177
       D0=0.05
       N0=-0.4776
      ELSEIF(TR.LE.1.)THEN
       SLOPE=-2.3503
       D0=0.025
       N0=-0.3716
      ENDIF
      REDVIS=10**(N0+SLOPE*(DSQRT(1.-TR)-DSQRT(D0)))
   10 CONTINUE
      END
