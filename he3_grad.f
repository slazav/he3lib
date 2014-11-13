!!! He3 gradient energy and spin wave velocity

! K1=K2=K3 in a simple approximation, see VW 7.23m
      function he3_grad_K0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, Y0, gape
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        gape = gap*const_kb*1D-3*he3_tc(p)
        he3_grad_K0 = 1/(5D0*gape**2)
     .     * (const_hbar/2D0/he3_amass)**2
     .     * (1D0-Y0) * he3_rho(p)
      end

! c parameter, see VW 7.25
      function he3_grad_c(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p, gap, y0,f1a,f1s
        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0D0)
        f1a = he3_f1a(p)
        f1s = he3_f1s(p)
        he3_grad_c = -(1D0-y0)*he3_rho(p)/10D0 *
     .    (3D0 + f1a)/(3D0 + f1s) /
     .    (1D0 + f1a*(2D0+3D0*y0)/15D0)
      end

! delta parameter, see VW 7.25
      function he3_grad_delta(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, Y0, f1a
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        f1a = he3_f1a(p)
        he3_grad_delta = f1a*(1D0-Y0)/(3D0+f1a*Y0)
      end

! K1=K2 with Fermi-liquid corrections
      function he3_grad_K12(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap
        gap = he3_trivgap(ttc, p)*const_kb*1D-3*he3_tc(p)
        he3_grad_K12 = - 2D0/gap**2
     .    * (const_hbar/2D0/he3_amass)**2
     .    * he3_grad_c(ttc, p)
      end

! K3 with Fermi-liquid corrections
      function he3_grad_K3(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap
        gap = he3_trivgap(ttc, p)*const_kb*1D-3*he3_tc(p)
        he3_grad_K3 = - 2D0/gap**2
     .    * (const_hbar/2D0/he3_amass)**2
     .    * (1D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
      end

! K = 2K1+K2+K3
      function he3_grad_K(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap
        gap = he3_trivgap(ttc, p)*const_kb*1D-3*he3_tc(p)
        he3_grad_K = - 2D0/gap**2
     .    * (const_hbar/2D0/he3_amass)**2
     .    * (4D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
      end

! K' = K2+K3
      function he3_grad_Kp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap
        gap = he3_trivgap(ttc, p)*const_kb*1D-3*he3_tc(p)
        he3_grad_Kp = - 2D0/gap**2
     .    * (const_hbar/2D0/he3_amass)**2
     .    * (2D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
      end

! lambda_G1 = Delta^2/2 (K2+K3)
      function he3_grad_lg1(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        he3_grad_lg1 =
     .    - (const_hbar/2D0/he3_amass)**2
     .    * (2D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
      end
! lambda_G2 = Delta^2/2 K1
      function he3_grad_lg2(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        he3_grad_lg2 =
     .    - (const_hbar/2D0/he3_amass)**2
     .    * he3_grad_c(ttc, p)
      end

! lambda_SG^b = Delta^2 K2
      function he3_grad_lsgb(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        he3_grad_lsgb =
     .    - 2D0*(const_hbar/2D0/he3_amass)**2
     .    * he3_grad_c(ttc, p)
      end

! Fomin's spin wave velocities (Fomin-1980 f.51)
      function he3_cpar(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, chi, ret
        chi = he3_chi_b(ttc, p) * he3_chi_n(p)
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          ret = - 2D0*he3_gyro**2/chi
     .      * (const_hbar/2D0/he3_amass)**2
     .      * (4D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
          ret = dsqrt(ret)
        elseif (ttc.eq.1D0) then
          ret = 0D0
        else
          ret = NaN
        endif
        he3_cpar = ret
      end
      function he3_cperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, chi, ret
        chi = he3_chi_b(ttc, p) * he3_chi_n(p)
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          ret = - 2D0*he3_gyro**2/chi
     .      * (const_hbar/2D0/he3_amass)**2
     .      * (3D0 + he3_grad_delta(ttc, p)/2D0) * he3_grad_c(ttc, p)
          ret = dsqrt(ret)
        elseif (ttc.eq.1D0) then
          ret = 0D0
        else
          ret = NaN
        endif
        he3_cperp = ret
      end

! Leggett's spin wave velocities (Leggett-1975 XII.B)
      function he3_clpar(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, chi, ret
        chi = he3_chi_b(ttc, p) * he3_chi_n(p)
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          ret = - 2D0*he3_gyro**2/chi
     .      * (const_hbar/2D0/he3_amass)**2
     .      * 2D0 * he3_grad_c(ttc, p)
          ret = dsqrt(ret)
        elseif (ttc.eq.1D0) then
          ret = 0D0
        else
          ret = NaN
        endif
        he3_clpar = ret
      end
      function he3_clperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, chi, ret
        chi = he3_chi_b(ttc, p) * he3_chi_n(p)
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          ret = - 2D0*he3_gyro**2/chi
     .      * (const_hbar/2D0/he3_amass)**2
     .      * (4D0 + he3_grad_delta(ttc, p)) * he3_grad_c(ttc, p)
          ret = dsqrt(ret)
        elseif (ttc.eq.1D0) then
          ret = 0D0
        else
          ret = NaN
        endif
        he3_clperp = ret
      end


