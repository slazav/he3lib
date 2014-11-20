!!! Q-balls in the zero temperature limit

! Leggett frequency [Hz] (measurements, same as he3_nu_b(0,p) )
      function qball_nu_b(p)
        implicit none
        include 'he3.fh'
        real*8 p
        qball_nu_b = (214.5D0 + 3.730D0*p
     .             - 646D0/(p+7.17D0))*1000
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
      function qball_dbetan(p, nu0)
        implicit none
        include 'he3.fh'
        real*8 p, nu0, A,B,C, cOb
        A = - 4.322539D-4*p**3 + 2.555245D-2*p**2
     .      - 1.338065D-3*p - 7.022825D-1
        B =   6.848111D-4*p**3 - 3.144707D-2*p**2
     .      - 4.791693D-1*p - 4.604034D-1
        C =   5.573567D-6*p**3 - 1.008726D-2*p**2
     .      + 6.884537D-1*p + 5.199754D0
        cOb = A*(nu0/1D6)**2 + B*(nu0/1D6) + C
        qball_dbetan = 1D9*cOb/
     .     (const_2pi*qball_nu_b(p)*qball_cper(p))
      end

! nu_z (1/2 of distance between visible axial levels) (rota-specific, measured), [Hz]
      function qball_fz(p,nu0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,nu0,imin,w0,imin1
        w0=const_2pi*nu0
        imin1 = imin - rota_hmina_i0*nu0
        qball_fz = qball_cpar(p)/const_2pi
     .           * sqrt( 8D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

! nu_r (1/2 of distance between visible radial levels) (rota-specific, measured), [Hz]
      function qball_fr(p,nu0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,nu0,imin,w0,imin1
        w0=const_2pi*nu0
        imin1 = imin - rota_hmina_i0*nu0
        qball_fr = qball_cper(p)/const_2pi
     .      * sqrt( 2D0*(qball_dbetan(p,nu0)*qball_nu_b(p)/nu0)**2
     .             - 4D0*he3_gyro*rota_hmina_mr*imin1/w0)
      end

! z size of the magnon condensate (rota-specific, measured), [cm]
      function qball_az(p,nu0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,nu0,imin
        qball_az = qball_cpar(p)/const_2pi
     .           * sqrt(2D0/qball_fz(p,nu0,imin)/nu0)
      end

! r size of the magnon condensate (rota-specific, measured), [cm]
      function qball_ar(p,nu0,imin)
        implicit none
        include 'he3.fh'
        real*8 p,nu0,imin
        qball_ar = qball_cper(p)/const_2pi
     .           * sqrt(2D0/qball_fr(p,nu0,imin)/nu0)
      end



