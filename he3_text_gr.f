!!! He3 gradient energy and spin wave velocity

! Textural parameter lambda_{G2}, erg/cm
! See Thuneberg-2001 f.28 and f.10
      function he3_text_lg2(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p, gap, y0,f1a

        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0D0)
        f1a = he3_f1a(p)
        he3_text_lg2 = const_hbar**2 * he3_rho(p)
     .  / (40D0*he3_meff(p) * he3_amass)
     .  * (1D0 + f1a/3D0)*(1D0-y0)/(1D0+f1a*(2D0+3D0*y0)/15D0)
      end

! Textural parameter lambda_{G1}, erg/cm
! See Thuneberg-2001 f.28 and f.10
      function he3_text_lg1(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        he3_text_lg1 = he3_text_lg2(ttc,p)
     .    * (2D0 + he3_text_delta(ttc,p))
      end

! Textural parameter delta
      function he3_text_delta(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, Y0, f1a
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        f1a = he3_f1a(p)
        he3_text_delta = f1a*(1D0-Y0)/(3D0+f1a*Y0)
      end

! Perpendicular spin wave velocity, cm/s
! See doc_tech/egrad.pdf
      function he3_text_cperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, lg2, chi, d
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          lg2 = he3_text_lg2(ttc, p)
          d   = he3_text_delta(ttc, p)
          chi = he3_chi_b(ttc, p) * he3_chi_n(p)
          he3_text_cperp = he3_gyro * sqrt((6D0+d)*lg2/chi)
        elseif (ttc.eq.1D0) then
          he3_text_cperp = 0D0
        endif
      end

! Parallel spin wave velocity, cm/s
! See doc_tech/egrad.pdf
      function he3_text_cpar(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, lg2, chi, d
        if (ttc.ge.0D0.and.ttc.lt.1D0) then
          lg2 = he3_text_lg2(ttc, p)
          d   = he3_text_delta(ttc, p)
          chi = he3_chi_b(ttc, p) * he3_chi_n(p)
          he3_text_cpar = he3_gyro * sqrt((8D0+2D0*d)*lg2/chi)
        elseif (ttc.eq.1D0) then
          he3_text_cpar = 0D0
        endif
      end

!! Bending stiffness coefficient c, erg/cm -- not used!
!! See Hakonen-1989 f.28
!      function he3_text_c(ttc, p)
!        implicit none
!        include 'he3.fh'
!        real*8 ttc, p
!        he3_text_c = 65D0/8D0 * he3_text_lg2(ttc, p)
!      end
