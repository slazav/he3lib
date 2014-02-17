!!! He3 other textural parameters

! Textural parameter a, erg/cm^3 1/G^2
! See Thuneberg-2001 f.25 and f.6
      function he3_text_a(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        real*8 gd,gap,f0a,z3,z5,y0


        gd  = he3_gd(p)
        gap = he3_trivgap(ttc,p)
        f0a = he3_f0a(p)
        z3  = he3_z3(ttc,gap)
        z5  = he3_z5(ttc,gap)
        y0  = 1-z3

        he3_text_a = NaN
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          he3_text_a = 2.5D0*gd
     .     * (0.5D0*const_hbar*he3_gyro/(1D0+f0a*(2D0+y0)/3D0))**2
     .     * (5D0-3D0*z5/z3)
        elseif (ttc.eq.1D0) then
          he3_text_a = 2.5D0*gd
     .     * 5D0 * (0.5D0*const_hbar*he3_gyro/(1D0+f0a))**2
        endif
      end

! Textural parameter lambda_{DV}, erg/cm^3 1/(cm/s)^2
! See Thuneberg-2001 f.26 and f.7
      function he3_text_ldv(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        real*8 gd,gap,f1s,z3,z5,y0
        gd  = he3_gd(p)
        gap = he3_trivgap(ttc,p)
        f1s = he3_f1s(p)
        z3  = he3_z3(ttc,gap)
        z5  = he3_z5(ttc,gap)
        y0  = 1-z3

        he3_text_ldv=NaN
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          he3_text_ldv = 5D0*gd
     .     * ( he3_meff(p) *he3_vf(p) / (1D0 + f1s*y0/3D0))**2
     .     * (1D0-1.5D0*z5/z3)
        elseif (ttc.eq.1D0) then
          he3_text_ldv = 5D0*gd
     .     * ( he3_meff(p) *he3_vf(p) / (1D0 + f1s/3D0))**2
        endif
      end

! Textural parameter lambda_{HV}, g/(cm^3 Gauss^2)
! See Thuneberg-2001 f.27 and f.8
      function he3_text_lhv(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        real*8 gap, y0, f1s,f0a, z3,z5,z7, gape
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc,gap, 0D0)
        f1s = he3_f1s(p)
        f0a = he3_f0a(p)
        z3  = he3_z3(ttc,gap)
        z5  = he3_z5(ttc,gap)
        z7  = he3_z7(ttc,gap)
        gape = gap*const_kb*1D-3*he3_tc(p)

        he3_text_lhv=NaN
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          he3_text_lhv = he3_rho(p) / gape**2
     .     * (1D0 + f1s/3D0) / (1D0 + f1s*Y0/3D0)**2
     .     * (0.5D0*const_hbar*he3_gyro /
     .            (1D0 + f0a*(2D0 + Y0)/3D0))**2
     .     * (z3 - 0.9D0*z5 + 0.9D0 * z5**2/z3 - 1.5D0*z7)
        elseif (ttc.eq.1D0) then
          he3_text_lhv = he3_rho(p)
     .     / (1D0 + f1s/3D0)
     .     * (0.5D0*const_hbar*he3_gyro / (1D0 + f0a))**2
     .     / (const_kb*1D-3*he3_tc(p))**2 * 2D0/3D0/const_pi
           ! note: z3/gap^2 = 2/3/pi at ttc=1
        endif
        ! negative values can appear because of inaccurate
        ! z3,z5,z7 at low temperatures
        if (he3_text_lhv.lt.0D0) he3_text_lhv=0D0
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Surface energy coefficient d, erg/(cm^2 Gauss^2)
! Came from ROTA texture library.
! Some G-L extrapolation is used
      function he3_text_d(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        real*8 gap,y0,f0a,n0,xi,d0

        he3_text_d = NaN
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          gap = he3_trivgap(ttc,p)
          y0  = he3_yosida(ttc,gap,0D0)
          f0a = he3_f0a(p)
          n0  = he3_2n0(p)/2D0
          xi  = he3_xigl(ttc,p)
          d0  = 2.2D0 + p*0.5D0/34.39D0
          he3_text_d = n0*(const_hbar*he3_gyro)**2
     .     * xi*d0*(1-y0)
     .     / (4D0*(1D0+f0a)*(3D0+f0a*(2D0+y0)))
        elseif (ttc.eq.1D0) then
          he3_text_d = 0D0
        endif
      end

! Vortex energy coefficient \lambda_{LH}
! Came from ROTA texture library. (difference: 5/2a)
! See Thuneberg-2001 f.30 and Kopu-2007 f.5
! Some G-L extrapolation is used
      function he3_text_llh(ttc,p, omega)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,omega
        real*8 h2m, rc, ri

        h2m = const_hbar/2D0/he3_amass
        rc = he3_xigl(ttc,p)
        ri = dsqrt(h2m/omega)
        he3_text_llh = h2m * omega * (dlog(ri/rc)-0.75D0)
     .    * he3_text_lhv(ttc,p)
      end

! lambda/omega value used in texture library
      function he3_text_lo(ttc,p, omega)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,omega
        he3_text_lo = 2.5D0
     .   * he3_text_llh(ttc,p,omega) / he3_text_a(ttc,p) / omega
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Magnetic length, cm
! see Thuneberg-2001, p.662
      function he3_text_xih(ttc, p, h)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,h
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          he3_text_xih =
     .      dsqrt(65D0*he3_text_lg2(ttc, p)
     .            /8D0/he3_text_a(ttc, p)/h**2)
        elseif (ttc.eq.1D0) then
            he3_text_xih = 0D0
        endif
      end

! Dipole length, cm
! see Thuneberg-2001, p.662
      function he3_text_xid(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        he3_text_xid =
     .    dsqrt(13D0/12D0 * he3_text_lg2(ttc, p)/he3_ld(ttc, p))
      end

! Dipole velocity vd in cm/s
      function he3_text_vd(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        he3_text_vd =
     .    dsqrt(0.4D0*he3_text_a(ttc,p)/he3_text_lhv(ttc,p))
      end


