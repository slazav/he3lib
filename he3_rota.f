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
        Imin0 = I*rota_hmina_i0i

        rota_bza =
     .     I*rota_nmra       ! field of the NMR solenoid
     .   - (Imin+Imin0)*Bmin ! field of the pinch coil
      end

      function rota_bza_xyz(I, Imin, x,y,z)
      end

! normal phase spectrum
      function rota_nspeca(f0, I, Imin)
        implicit none
        include 'he3.fh'
        real*8 f0,I,Imin

        ! we need a grid with an isotropic spacing,
        ! x=-R..R, y=0..R, z=0..2R
        integer N, nx,ny,nz, ix,iy,iz
        parameter (N=20, nx=2*N-1, ny=N, nz=3*N-1)
        real*8 dx, Hres, F(nx,ny,nz), Fx, Fy, Fz, Fg
        real*8 x(nx), y(ny), z(nz), r
        real*8 e, dd, g, r1
        real*8 nmra_r, nmra_d
        data nmra_r /0.4D0/, nmra_d /0.66D0/

        dx = rota_rcell/dble(N-1)
        Hres = const_2pi*f0/he3_gyro ! resonance field

        ! fill arrays
        do ix=1,nx
          x(ix) = rota_rcell * (2D0*dble(ix-1)/dble(nx-1) - 1D0)
          do iy=1,ny
            y(iy) = rota_rcell * dble(iy-1)/dble(ny-1)
            r = dsqrt(x(ix)**2 + y(iy)**2)
            do iz=1,nz
              z(iz) = rota_rcell * 3D0*dble(iz-1)/dble(nz-1)
              F(ix,iy,iz) = rota_bza(I,Imin, r,z(iz))
     .                       + x(ix)*0.2D0 - Hres
            enddo
          enddo
        enddo

        ! calculate integral
        rota_nspeca=0D0
        do ix=1,nx
          do iy=1,ny
            do iz=1,nz
              if (ix.lt.nx) then
                Fx = (F(ix+1,iy,iz) - F(ix,iy,iz))/dx
              else
                Fx = (F(ix,iy,iz) - F(ix-1,iy,iz))/dx
              endif
              if (iy.lt.ny) then
                Fy = (F(ix,iy+1,iz) - F(ix,iy,iz))/dx
              else
                Fy = (F(ix,iy,iz) - F(ix,iy-1,iz))/dx
              endif
              if (iz.lt.nz) then
                Fz = (F(ix,iy,iz+1) - F(ix,iy,iz))/dx
              else
                Fz = (F(ix,iy,iz) - F(ix,iy,iz-1))/dx
              endif
              Fg = dsqrt(Fx**2 + Fy**2 + Fz**2)
              e=dx*(dabs(Fx)+dabs(Fy)+dabs(Fz))/Fg
              dd = (1D0-dabs(F(ix,iy,iz)/e))/e
              if (dd.lt.0D0) dd=0D0
              if (x(ix)**2+y(iy)**2.le.rota_rcell) then
                ! sensitivity of the NRM pick-up coil
                r1 = dsqrt(y(iy)**2 + z(iz)**2)
                g = loop_bz(nmra_r, r1, x(ix)-nmra_d)
     .             +loop_bz(nmra_r, r1, x(ix)+nmra_d)
              else
                g = 0D0
              endif
              rota_nspeca = rota_nspeca + g*dd*dx**3
            enddo
          enddo
        enddo

      end
