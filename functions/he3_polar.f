!HH> He3 Polar phase

!> susceptibility component along d vector chi_par / chi_n
!> chi_perp/chi_0 = 1
!> see Leggett-1975, VIIID f.7.53 and f.7.54
      function he3p_chi_par(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap,f0a,Y0
        f0a = he3_f0a(p)
        gap = he3_gap(ttc,p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        he3p_chi_par =
     .    (1D0 + f0a) * Y0 / (1D0 + f0a * Y0)
      end

!> dipolar length perpendicular to the l vector, 4g_d/K1
      function he3p_xid_perp(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,K
        K = he3_grad_k12(ttc,p)
        he3p_xid_perp = dsqrt(K/he3_gd(p))/2D0
      end
!> dipolar length parallel to the l vector, 4g_d/(K1+K2+K3)
      function he3p_xid_par(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,K
        K = 2D0*he3_grad_k12(ttc,p) + he3_grad_k3(ttc,p)
        he3p_xid_par = dsqrt(K/he3_gd(p))/2D0
      end

!> magnetic length perpendicular to the l vector (4g_d/K1)
      function he3p_xih_perp(ttc, p, h) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,h, gap,dchi,K
        gap = he3_gap(ttc,p)*const_kb*1D-3*he3_tc(p)
        dchi = he3_chi_n(p) * (1D0-he3p_chi_par(ttc,p))
        K = he3_grad_k12(ttc,p)
        he3p_xih_perp = (gap/h) * dsqrt(K/dchi)
      end

!> magnetic length parallel to the l vector (4g_d/(K1+K2+K3))
      function he3p_xih_par(ttc, p, h) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,h, gap,dchi,K
        gap = he3_gap(ttc,p)*const_kb*1D-3*he3_tc(p)
        dchi = he3_chi_n(p) * (1D0-he3p_chi_par(ttc,p))
        K = 2D0*he3_grad_k12(ttc,p)+he3_grad_k3(ttc,p)
        he3p_xih_par = (gap/h) * dsqrt(K/dchi)
      end
