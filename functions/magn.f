!HH> Magnetic models

!H> Curie-Weiss magnet

!> See <a href="https://arxiv.org/pdf/1301.2141.pdf">[Kochmansky]</a>.

!> Magnetization of Curie-Weiss magnet,  M/mu vs T/Tc and muB/kTc
!> Solving equation m = \tanh((m+btc)/ttc) by Newton method.
!> <br>Works for positive and negative btc.
      function magn_cw_m(ttc, btc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, btc
        real*8 s, y, dy, F,Fp
        s=1D0
        if (btc.lt.0D0) s=-1D0
        y=1D0
        dy=1D0
        do while (dy.GT.1D-10)
          F = s*y - dtanh((s*y+btc)/ttc);
          Fp = s - (1D0 - dtanh((s*y+btc)/ttc)**2)*s/ttc
          dy = F/Fp
          y=y-dy
        enddo
        magn_cw_m = s*y
      end

!> Magnetic susceptibility of Curie-Weiss magnet,  chi*kTc/mu^2 vs T/Tc and muB/kTc
      function magn_cw_chi(ttc,btc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, btc
        real*8 y
        y = magn_cw_m(ttc, btc)
        magn_cw_chi = (1D0-y**2)/(ttc-1D0+y**2)
      end

!> Entropy of Curie-Weiss magnet, S/R vs T/Tc and muB/kTc
      function magn_cw_s(ttc,btc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, btc
        real*8 y
        y = magn_cw_m(ttc, btc)
        magn_cw_s = dlog(2D0*dexp(-y**2/2D0/ttc)*dcosh((y+btc)/ttc))
     .    - (y/2D0 + btc)*y/ttc
      end

!> Heat capacity of Curie-Weiss magnet,  C/R vs T/Tc and muB/kTc
        function magn_cw_c(ttc,btc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, btc
        real*8 y
        y = magn_cw_m(ttc, btc)
        magn_cw_c = (1D0-y**2)*(y+btc)**2 / (ttc-1+y**2) / ttc
      end

!><p> <img src="img/magn.cw.png">
