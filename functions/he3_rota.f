!HH> ROTA-specific functions

!> Nuclear stage heat capacity [J/K] vs T[K] and I[A]
      function rota_c_ns(t, i) !F>
        implicit none
        include 'he3.fh'
        real*8 t, i, k, kmag
        k = 9.66D-5     ! J (K/T)^2 - in heat capacity
        kmag = 0.113D0  ! T/A - magnet constant
        rota_c_ns = k * (kmag * i/t)**2
      end

!> Calibration of fork N, T/Tc, vs width (Hz) and P (bar)
      function rota_fork_cal(w, p, n) !F>
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

!> Q value of the nmrA spectrometer vs frequency (measured)
      function rota_nmra_q(f0) !F>
        implicit none
        include 'he3.fh'
        real*8 f0
        rota_nmra_q = 31.796D-6*f0 + 110.03979D0
      end

!> ROTA: frequencies of nmrA spectrometer,kHz for n=1..8 (use real*8 n!)
      function rota_nmra_f(n) !F>
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

!> ROTA: Bz field profile of the A spectrometer
      function rota_bza(I, Imin, r, z) !F>
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

!> ROTA: normal phase spectrum
      function rota_nspeca(f0, I, Imin) !F>
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Q-balls in the zero temperature limit

!> ROTA: Derivative of the textural angle beta_N in the center of the cell (rota-specific, measured), [rad/cm]
      function rota_qball_dbetan(p, f0) !F>
        implicit none
        include 'he3.fh'
        real*8 p, f0, A,B,C, cOb
        A =  9.144220D-6*p**4 - 1.202768D-3*p**3
     .      +4.139982D-2*p**2 - 6.613861D-2*p - 4.830108D-1
        B = -9.742400D-6*p**4 + 1.570559D-3*p**3
     .      -5.013987D-2*p**2 - 3.998610D-1*p - 8.127855D-1
        C = -1.165609D-5*p**4 + 6.445247D-4*p**3
     .      -2.218588D-2*p**2 + 7.691508D-1*p + 5.337443D+0
        cOb = A*(f0/1D6)**2 + B*(f0/1D6) + C
        rota_qball_dbetan = 1D9*cOb/
     .     (const_2pi*he3_nu_b(0D0,p)*he3_cperp(0D0,p))
      end

!> ROTA: nu_z (1/2 of distance between visible axial levels) (no rotation, rota-specific, measured), [Hz]
      function rota_qball_fz0(p,f0,imin) !F>
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin + rota_hmina_i0f*f0
        rota_qball_fz0 = he3_cpar(0D0,p)/const_2pi
     .           * sqrt( 8D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

!> ROTA: nu_r (1/2 of distance between visible radial levels) (no rotation, rota-specific, measured), [Hz]
      function rota_qball_fr0(p,f0,imin) !F>
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin + rota_hmina_i0f*f0
        rota_qball_fr0 = he3_cperp(0D0,p)/const_2pi
     .    * sqrt(2D0*(rota_qball_dbetan(p,f0)*he3_nu_b(0D0,p)/f0)**2
     .             - 4D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

!> ROTA: z size of the magnon condensate (no rotation, rota-specific, measured), [cm]
      function rota_qball_az0(P,f0,imin) !F>
        implicit none
        include 'he3.fh'
        real*8 P,f0,imin
        rota_qball_az0 = he3_cpar(0D0,P)/const_2pi
     .           * sqrt(2D0/rota_qball_fz0(P,f0,imin)/f0)
      end

!> ROTA: r size of the magnon condensate (no rotation, rota-specific, measured), [cm]
      function rota_qball_ar0(P,f0,imin) !F>
        implicit none
        include 'he3.fh'
        real*8 P,f0,imin
        rota_qball_ar0 = he3_cperp(0D0,P)/const_2pi
     .           * sqrt(2D0/rota_qball_fr0(P,f0,imin)/f0)
      end

!> ROTA: tau_RD for the magnon condensate with given radial and axial frequencies (rota-specific, measured) [s]
      function rota_qball_trd(P, f0,  fr, fz) !F>
        implicit none
        include 'he3.fh'
        real*8 P,f0, fr,fz,ar,az, krd,chi,oB,H

        chi = he3_chi_b(0D0,P) * he3_chi_n(P)
        oB  = const_2pi * he3_nu_b(0D0,P)
        H   = const_2pi * f0 / he3_gyro

        az = he3_cpar(0D0,p)/const_2pi * sqrt(2D0/fz/f0)
        ar = he3_cperp(0D0,p)/const_2pi * sqrt(2D0/fr/f0)
        krd = 8D0*const_pi**1.5D0*chi*ar**2*az*H*rota_nmra_q(f0)
        rota_qball_trd = rota_rrda / krd
      end

!> ROTA: tau_RD for the magnon condensate (no rotation, rota-specific, measured) [s]
      function rota_qball_trd0(P,f0,imin) !F>
        implicit none
        include 'he3.fh'
        real*8 P,f0, imin, fr,fz
        fr = rota_qball_fr0(P, f0, imin)
        fz = rota_qball_fz0(P, f0, imin)
        rota_qball_trd0 = rota_qball_trd(P, f0,  fr, fz)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
