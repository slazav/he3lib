!!! Q-balls in the zero temperature limit

! Leggett frequency [Hz] (measurements)
      function qball_nu_b(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_nu_b = (304.1D0 + 2.077D0*p
     .              - 2693D0/(p+13.84D0))*1D3
      end

! cpar [cm/s] (measured)
      function qball_cpar(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_cpar = (38666D0/(20.544D0+P) + 593.7D0)
     .         * sqrt(1.652D0/rota_hmina_mr)
      end

! cper [cm/s] (measured)
      function qball_cper(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_cper = (34849D0/(21.138D0+P) + 512.7D0)
     .         * sqrt(1.652D0/rota_hmina_mr)
      end

! lambda_g1 (measured)
      function qball_lg1(p)
        implicit none
        include 'he3.fh'
        real*8 p, chi
        chi = he3_chi_b(0D0, p)*he3_chi_n(p)
        qball_lg1 = chi/(he3_gyro**2) *
     .     (qball_cpar(p)**2 - qball_cper(p)**2)
      end
! lambda_g2 (measured)
      function qball_lg2(p)
        implicit none
        include 'he3.fh'
        real*8 p, chi
        chi = he3_chi_b(0D0, p)*he3_chi_n(p)
        qball_lg2 = chi/(he3_gyro**2) *
     .     (2D0*qball_cper(p)**2 - qball_cpar(p)**2)/4D0
      end

! Derivative of the textural angle beta_N in the center
!  of the cell (rota-specific, measured), [rad/cm]
      function qball_dbetan(p, f0)
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
        qball_dbetan = 1D9*cOb/
     .     (const_2pi*qball_nu_b(p)*qball_cper(p))
      end

! nu_z (1/2 of distance between visible axial levels) (no rotation, rota-specific, measured), [Hz]
      function qball_fz0(p,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin + rota_hmina_i0f*f0
        qball_fz0 = qball_cpar(p)/const_2pi
     .           * sqrt( 8D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

! nu_r (1/2 of distance between visible radial levels) (no rotation, rota-specific, measured), [Hz]
      function qball_fr0(p,f0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,f0,imin,w0,imin1
        w0=const_2pi*f0
        imin1 = imin + rota_hmina_i0f*f0
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

!  measured textural parameter a
      function qball_text_a(P)
        implicit none
        include 'he3.fh'
        real*8 P
        qball_text_a = 1D-13 * (
     .      3.748354D-5*P**3 - 1.252762D-3*P**2
     .    + 7.177360D-2*P + 7.563687D-1)
      end

!  measured textural parameter d
      function qball_text_d(P)
        implicit none
        include 'he3.fh'
        real*8 P
        qball_text_d = 1D-13*(
     .    - 7.801035D-5*P**3 + 2.234899D-3*P**2
     .    + 2.046272D-2*P + 7.348175D-1)
      end

!  measured textural parameter lambda_SG
      function qball_text_lsg(P)
        implicit none
        include 'he3.fh'
        real*8 P
        qball_text_lsg =
     .    - 1.655711D-3*P**3 + 7.259754D-2*P**2
     .    + 7.128994D-2*P +1.170009D+1
      end
