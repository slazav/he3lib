
      SUBROUTINE DALPHA( MALPHA)
C
C  Gets magic alpha for dilute phase
C
      REAL*8 MALPHA
      CHARACTER ANSWER
      MALPHA=2.2
   10 WRITE(*,*)'Dilute phase MAGIC ALPHA assumed to be 2.2 Is this corr
     +ect? [y/n]'
      READ(*,100) ANSWER
  100 FORMAT(A1)
      IF (ANSWER.NE. 'Y'.AND. ANSWER .NE.'N'.AND. ANSWER.NE. 'y'
     + .AND. ANSWER .NE.'n')GOTO 10
      IF (ANSWER .EQ. 'Y'.OR. ANSWER .EQ. 'y')GOTO 20
        WRITE(*,*)'Enter new value'
        READ(*,*) MALPHA
   20 CONTINUE
      RETURN
      END