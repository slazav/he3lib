!!! He3 textural parameters

! Textural parameter a in 10^(-14) erg/(Gauss^2 cm^3)
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

        he3_text_a = 2.5D0*gd*(
     .   0.5*const_hbar*he3_gyro/(1D0+f0a*(2D0+y0)/3D0))**2
     .   * (5D0-3D0*z5/z3)

!        if (dabs(1-ttc)<1D-8) then
!          he3_text_a = he3_text_a * 5D0
!        else
!          he3_text_a = he3_text_a * (5D0-3D0*z5/z3)
!        end if
      end

! Textural parameter delta
      function he3_text_delta(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, Y0, f1a
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        f1a = he3_f1a(p)
        he3_text_delta = f1a*(1-Y0)/(3+f1a*Y0)
      end

! Textural parameter lhv in 10^(-8) g/(cm^3 Gauss^2)
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

        gape = gap*const_kb*he3_tc(p)

        he3_text_lhv = he3_rho(p) / gape**2
     .   * (1D0 + f1s/3D0) / (1D0 + f1s*Y0/3D0)**2
     .   * (0.5D0*const_hbar*he3_gyro /
     .          (1D0 + f0a*(2D0 + Y0)/3D0))**2
     .   * (z3 - 0.9D0*z5 + 0.9D0 * z5**2/z3 - 1.5D0*z7)
      end

! Textural parameter lg2 in erg/cm
      function he3_text_lg2(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p, gap, y0,f1a

        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0)
        f1a = he3_f1a(p)
        he3_text_lg2 = const_hbar**2 * const_na
     .  * (1D0 + f1a/3D0)*(1D0-y0)/(1D0+f1a*(2D0+3D0*y0)/15D0)
     .  / (40D0*he3_meff(p) * he3_vm(p))
      end


! Dipole velocity vd in cm/s
      function he3_text_vd(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        he3_text_vd =
     .    dsqrt(0.4D0*he3_text_a(ttc,p)/he3_text_lhv(ttc,p))
      end

! Dipole length, cm
! see Thuneberg-2001, p.662
      function he3_text_ksid(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        he3_text_ksid =
     .    dsqrt(he3_text_lg2(ttc, p)/he3_ld(ttc, p))
      end

