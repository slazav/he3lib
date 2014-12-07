!!! ROTA specific functions

! Nuclear stage heat capacity [J/K] vs T[K] and I[A]
      function rota_c_ns(t, i)
        implicit none
        include 'he3.fh'
        real*8 t, i, k, kmag
        k = 9.66D-5     ! J (K/T)^2 - in heat capacity
        kmag = 0.113D0  ! T/A - magnet constant
        rota_c_ns = k * (kmag * i/t)**2
      end

! Calibration of fork N, T/Tc, vs width (Hz) and P (bar)
      function rota_fork_cal(w, p, n)
        implicit none
        include 'he3.fh'
        real*8 w,p,n,a, ttc,ttc1,e
        integer cnt
        a=NaN
        if (nint(n).eq.1) a = 11700D0 ! fork K, 30.4.2010, 29 bar
        if (nint(n).eq.2) a = 17543D0 ! fork E, 30.4.2010, 29 bar
        a = a * (he3_pf(p)/he3_pf(29D0))**4 ! a = pf^4 alpha
        ttc = 0D0
        e = 1D0
        cnt=100
        do while (e.gt.1D-6.and.cnt.gt.0)
          ttc1 = he3_gap(ttc, p)/dlog(a/w)
          e = dabs(ttc-ttc1)
          ttc = ttc1
          cnt = cnt-1
         enddo
        rota_fork_cal=ttc
      end

! Q value of the nmrA spectrometer vs frequency (measured)
      function rota_nmra_q(f0)
        implicit none
        include 'he3.fh'
        real*8 f0
        rota_nmra_q = 31.796D-6*f0 + 110.03979D0
      end

! frequencies of nmrA spectrometer,kHz for n=1..8 (use real*8 n!)
      function rota_nmra_f(n)
        implicit none
        include 'he3.fh'
        real*8 n
        integer ni, f(8)
        ni = int(n);
        f = (/553, 588, 623, 674, 710, 742, 789, 833/)
        if (ni.ge.1.and.ni.le.8) then
          rota_nmra_f = dble(f(ni))
        else
          rota_nmra_f = 0D0
        endif
      end

! Bz field profile of the A spectrometer
      function rota_bza(I, Imin, r, z)
        implicit none
        include 'he3.fh'
        real*8 I,Imin,r,z
        real*8 Bmin, Imin0

        ! Field of the pinch coil at 1A
        Bmin = rota_hmina_n *(loop_bz(rota_hmina_r, r, z)
     .                      - loop_bz(rota_hmina_r, 0D0, 0D0))
     .       + rota_hmina

        ! Effective pinch coil current due to distortion of the
        ! NMR field by the pinch coil
        Imin0 = I*rota_hmina_i0

        rota_bza =
     .     I*rota_nmra       ! field of the NMR solenoid
     .   - (Imin+Imin0)*Bmin ! field of the pinch coil
      end
