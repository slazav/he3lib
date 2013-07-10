*UPTODATE M01AGETEXT
      SUBROUTINE M01AGE(A, NR, NC, IC, T, TT, IFAIL)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     SORTS THE ROWS OF AN INTEGER MATRIX INTO ASCENDING ORDER OF
C     AN INDEX
C     COLUMN.
      INTEGER P01AAE, NR, NC, IC, M, I, J, K, IJ, N, L, IFAIL,
     * A(NR,NC), T(NC), TT(NC), IU(22), IL(22)
      CHARACTER SRNAME*8
      DATA SRNAME /'M01AGE'/
      IF (NR) 820, 820, 20
   20 IF (NC) 800, 800, 40
   40 IF (IC) 780, 780, 60
   60 IF (IC-NC) 80, 80, 760
   80 M = 1
      I = 1
      J = NR
  100 IF (I-J) 120, 500, 500
  120 K = I
      IJ = (J+I)/2
      DO 140 N=1,NC
         T(N) = A(IJ,N)
  140 CONTINUE
      IF (A(I,IC)-T(IC)) 200, 200, 160
  160 DO 180 N=1,NC
         A(IJ,N) = A(I,N)
         A(I,N) = T(N)
         T(N) = A(IJ,N)
  180 CONTINUE
  200 L = J
      IF (A(J,IC)-T(IC)) 220, 340, 340
  220 DO 240 N=1,NC
         A(IJ,N) = A(J,N)
         A(J,N) = T(N)
         T(N) = A(IJ,N)
  240 CONTINUE
      IF (A(I,IC)-T(IC)) 340, 340, 260
  260 DO 280 N=1,NC
         A(IJ,N) = A(I,N)
         A(I,N) = T(N)
         T(N) = A(IJ,N)
  280 CONTINUE
      GO TO 340
  300 DO 320 N=1,NC
         A(L,N) = A(K,N)
         A(K,N) = TT(N)
  320 CONTINUE
  340 L = L - 1
      IF (A(L,IC)-T(IC)) 360, 360, 340
  360 DO 380 N=1,NC
         TT(N) = A(L,N)
  380 CONTINUE
  400 K = K + 1
      IF (A(K,IC)-T(IC)) 400, 420, 420
  420 IF (K-L) 300, 300, 440
  440 IF (L-I-J+K) 480, 480, 460
  460 IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 540
  480 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 540
  500 M = M - 1
      IF (M) 520, 840, 520
  520 I = IL(M)
      J = IU(M)
  540 IF (J-I-11) 560, 120, 120
  560 IF (I-1) 580, 100, 580
  580 I = I - 1
  600 I = I + 1
      IF (I-J) 620, 500, 620
  620 DO 640 N=1,NC
         T(N) = A(I+1,N)
  640 CONTINUE
      IF (A(I,IC)-T(IC)) 600, 600, 660
  660 K = I
  680 DO 700 N=1,NC
         A(K+1,N) = A(K,N)
  700 CONTINUE
      K = K - 1
      IF (T(IC)-A(K,IC)) 680, 720, 720
  720 DO 740 N=1,NC
         A(K+1,N) = T(N)
  740 CONTINUE
      GO TO 600
  760 IFAIL = P01AAE(IFAIL,4,SRNAME)
      GO TO 860
  780 IFAIL = P01AAE(IFAIL,3,SRNAME)
      GO TO 860
  800 IFAIL = P01AAE(IFAIL,2,SRNAME)
      GO TO 860
  820 IFAIL = P01AAE(IFAIL,1,SRNAME)
      GO TO 860
  840 IFAIL = 0
  860 RETURN
      END
** END OF M01AGETEXT
