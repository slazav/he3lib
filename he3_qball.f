!!! Q-balls in the zero temperature limit

! Leggett frequency [Hz] (measurements, same as he3_nu_b(0,p) )
      function qball_nu_b(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_nu_b = (214.5D0 + 3.730D0*p
     .             - 646D0/(p+7.17D0))*1D3
      end

! cpar [cm/s] (measured)
      function qball_cpar(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_cpar = (40400D0/(p+21.0D0)+565D0)
     .         * sqrt(1.652D0/rota_hmina_mr)
      end

! cper [cm/s] (measured)
      function qball_cper(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_cper = (36320D0/(p+21.6D0)+488D0)
     .         * sqrt(1.652D0/rota_hmina_mr)
      end

! Derivative of the textural angle beta_N in the center
!  of the cell (rota-specific, measured), [rad/cm]
      function qball_dbetan(p, f0)
        implicit none
        include 'he3.fh'
        real*8 p, f0, A,B,C, cOb
        A = - 7.238979D-4*p**3 + 3.407278D-2*p**2
     .      - 3.239890D-2*p - 5.758475D-1
        B = + 1.089347D-3*p**3 - 4.355957D-2*p**2
     .      - 4.251184D-1*p - 6.850477D-1
        C = - 1.281463D-4*p**3 - 5.982353D-3*p**2
     .      + 6.663006D-1*p + 5.294395D0
        cOb = A*(f0/1D6)**2 + B*(f0/1D6) + C
        qball_dbetan = 1D9*cOb/
     .     (const_2pi*qball_nu_b(p)*qball_cper(p))
      end

! nu_z (1/2 of distance between visible axial levels) (no rotation, rota-specific, measured), [Hz]
      function qball_fz0(p,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin - rota_hmina_i0*f0
        qball_fz0 = qball_cpar(p)/const_2pi
     .           * sqrt( 8D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

! nu_r (1/2 of distance between visible radial levels) (no rotation, rota-specific, measured), [Hz]
      function qball_fr0(p,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin - rota_hmina_i0*f0
        qball_fr0 = qball_cper(p)/const_2pi
     .      * sqrt( 2D0*(qball_dbetan(p,f0)*qball_nu_b(p)/f0)**2
     .             - 4D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

! z size of the magnon condensate (no rotation, rota-specific, measured), [cm]
      function qball_az0(P,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 P,f0,imin
        qball_az0 = qball_cpar(P)/const_2pi
     .           * sqrt(2D0/qball_fz0(P,f0,imin)/f0)
      end

! r size of the magnon condensate (no rotation, rota-specific, measured), [cm]
      function qball_ar0(P,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 P,f0,imin
        qball_ar0 = qball_cper(P)/const_2pi
     .           * sqrt(2D0/qball_fr0(P,f0,imin)/f0)
      end

!  tau_RD for the magnon condensate with given radial and axial frequencies (rota-specific, measured) [s]
      function qball_trd(P, f0,  fr, fz)
        implicit none
        include 'he3.fh'
        real*8 P,f0, fr,fz,ar,az, krd,chi,oB,H

        chi = he3_chi_b(0D0,P) * he3_chi_n(P)
        oB  = const_2pi * he3_nu_b(0D0,P)
        H   = const_2pi * f0 / he3_gyro

        az = qball_cpar(p)/const_2pi * sqrt(2D0/fz/f0)
        ar = qball_cper(p)/const_2pi * sqrt(2D0/fr/f0)
        krd = 8D0*const_pi**1.5D0*chi*ar**2*az*H*rota_nmra_q(f0)
        qball_trd = rota_rrda / krd
      end

!  tau_RD for the magnon condensate (no rotation, rota-specific, measured) [s]
      function qball_trd0(P,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 P,f0, imin, fr,fz
        fr = qball_fr0(P, f0, imin)
        fz = qball_fz0(P, f0, imin)
        qball_trd0 = qball_trd(P, f0,  fr, fz)
      end
