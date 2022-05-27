!HH> Magnetic models

!H> Curie-Weiss magnet with S=1/2

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

        ! note: (dlog(dcosh(x))' = dtanh(x)
        magn_cw_s = dlog(2D0) + dlog(dcosh((y+btc)/ttc))
     .            - y*(y + btc)/ttc

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

!><p> <img src="img/magn_cw.png">



!H> Paramegnetetic material with internal field

!> See Pobell book f9.15

!> Use muB/kT = gyro*hbar*sqrt(Bext^2 + Bint^2) / kT

!> Magnetization of paramagnetic material, M/Ms vs muB/kT and spin
      function magn_par_m(bt, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 bt, spin
        real*8 x,y
        x = bt/2D0
        y = (2D0*spin + 1D0)*x
        magn_par_m = (2D0*spin+1)/(2D0*spin)/dtanh(y)
     .             - 1D0/(2D0*spin)/dtanh(x)
      end

!> Entropy of paramagnetic material, S/R vs muB/kT and spin
      function magn_par_s(bt, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 bt, spin
        real*8 x,y
        x = bt/2D0
        y = (2D0*spin + 1D0)*x
!        C = (x/sinh(x)**2 - (y/sinh(y))**2
        magn_par_s = x/dtanh(x) - y/dtanh(y) + dlog(sinh(y)/sinh(x))
      end

!> Heat capacity of paramagnetic material, C/R vs muB/kT and spin
      function magn_par_c(bt, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 bt, spin
        real*8 x,y
        x = bt/2D0
        y = (2D0*spin + 1D0)*x
        magn_par_c = (x/dsinh(x))**2 - (y/dsinh(y))**2
      end

!> <p>Example for copper nuclei.  Internal field $B_i$ = 0.36e-3 T, gyromagnetic ratio $\gamma$ = 71.118e6 rad/s/T,
!> spin $J$ = 3/2. First argument for calling functions is $\gamma \hbar \sqrt{B^2 + B_i^2} / k_B T$.

!><p> <img src="img/magn_par.png">
