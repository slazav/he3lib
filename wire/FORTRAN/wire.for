
      SUBROUTINE WIRE(DIA, RHOWIR)
C
C gets the properties of the vibrating wire
C
      REAL*8 DIA, RHOWIR
      INTEGER N
      CHARACTER*10 METAL
   20   WRITE(*,*)'Enter material: 1 = TANTALUM; 2 = NbTi; 3 = other'
        READ(*,*) N
        IF(N.EQ.1)THEN
          METAL='TANTALUM'
          RHOWIR=16.7
        ELSEIF(N.EQ.2) THEN
          METAL='NbTi'
          RHOWIR=6.05
        ELSEIF(N.EQ.3) THEN
          WRITE(*,*)'INPUT THE NAME OF THE MATERIAL'
          READ(*,101)METAL
  101     FORMAT(A10)
          WRITE(*,*)'INPUT THE DENSITY in g/cm^3'
          READ(*,*)RHOWIR
        ELSEIF(N.NE.1 .AND. N.NE.2 .AND. N.NE.3)THEN
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'
          GOTO 20
        ENDIF
C
C   get wire diameter
C
       WRITE(*,*) 'ENTER wire DIAMETER in mm: 1= 0.124 mm, 2= OTHER'
       READ (*,*) N
       IF (N.NE.1 .AND. N.NE.2)THEN
          WRITE(*,*)'********* ERROR --- TRY AGAIN *********'
       ELSEIF(N.EQ.1)THEN
          DIA=0.124
       ELSEIF(N.EQ.2)THEN
          WRITE(*,*)'input diameter in mm'
          READ(*,*) DIA
       ENDIF
      WRITE(*,120)METAL ,RHOWIR, DIA
      WRITE(10,120)METAL ,RHOWIR, DIA
  120 FORMAT(1X,'Wire properties: material is ', A10,
     +' of density ', F8.4,' g/cm^3'
     +/1X,' and diameter ', F7.5,' mm'/1x)
      RETURN
      END
