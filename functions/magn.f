!HH> Magnetic models

!> <p>TODO: it's better to re-chack magnetization units. It's checked that magnetization
!> is consistent with susceptibility and entropy is consistent with heat capacity and
!> demagnetization cooling effect.


!H> Curie-Weiss magnet with spin 1/2

!> See <a href="https://arxiv.org/pdf/1301.2141.pdf">[Kochmansky]</a>.


!> y-function, dimensionless magnetization of S=1/2 Curie-Weiss magnet,  M/mu vs T/Tc and muB/kTc
!> Solving equation m = \tanh((m+btc)/ttc) by Newton method.
!> Works for positive and negative field.
      function magn_cw_y(ttc, btc) !F>
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
        magn_cw_y = s*y
      end

!> Molar magnetization of S=1/2 Curie-Weiss magnet,  M[J/T/mole] vs T[K], B[T], Tc[K], gyro[rad/s/T]
      function magn_cw_m(T, B, Tc, gyro) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Tc, gyro
        real*8 ttc, btc, mu,y

        mu = 0.5*gyro*const_hbar
        ttc = T/Tc
        btc = mu*B / (const_kb*Tc)
        y=magn_cw_y(ttc, btc)

        magn_cw_m = y * 0.5*gyro*const_hbar * const_na
      end

!> Molar magnetic susceptibility of S=1/2 Curie-Weiss magnet, chi [J/T^2/mole] vs T[K], B[T], Tc[K], gyro[rad/s/T]
      function magn_cw_chi(T, B, Tc, gyro) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Tc, gyro
        real*8 ttc, btc, mu, y

        mu = 0.5*gyro*const_hbar
        ttc = T/Tc
        btc = mu*B / (const_kb*Tc)
        y = magn_cw_y(ttc, btc)

        magn_cw_chi = (1D0-y**2)/(ttc-1D0+y**2)
        magn_cw_chi = magn_cw_chi * mu**2/(const_kb*Tc) * const_na
      end

!> Entropy of S=1/2 Curie-Weiss magnet, S/R vs T[K], B[T], Tc[K], gyro[rad/s/T]
      function magn_cw_s(T, B, Tc, gyro) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Tc, gyro
        real*8 ttc, btc, mu, y

        mu = 0.5*gyro*const_hbar
        ttc = T/Tc
        btc = mu*B / (const_kb*Tc)
        y = magn_cw_y(ttc, btc)

        ! note: (dlog(dcosh(x))' = dtanh(x)
        magn_cw_s = dlog(2D0) + dlog(dcosh((y+btc)/ttc))
     .            - y*(y + btc)/ttc

      end

!> Heat capacity of S=1/2 Curie-Weiss magnet,  C/R vs T[K], B[T], Tc[K], gyro[rad/s/T]
      function magn_cw_c(T, B, Tc, gyro) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Tc, gyro
        real*8 ttc, btc, mu, y

        mu = 0.5*gyro*const_hbar
        ttc = T/Tc
        btc = mu*B / (const_kb*Tc)
        y = magn_cw_y(ttc, btc)

        magn_cw_c = (1D0-y**2)*(y+btc)**2 / (ttc-1+y**2) / ttc
      end


! y =  tanh((y+btc)/ttc) => dy = (1-y^2)/(ttc-1+y^2) [ dbtc - (y+btc)/ttc*dttc]
! D = (dS/dB)/(dS/dT) = (dy/dB)/(dy/dT)   [because S=S(y), dS=S'dy]
! dy/dB = dy/dbtc * mu/kTc = (1-y^2)/(ttc-1+y^2) *  mu/kTc
! dy/dT = dy/dttc * 1/Tc  = - (1-y^2)/(ttc-1+y^2) * (y+btc)/ttc *1/Tc
! D = - mu/k ttc/(y+btc)

!> Cooling effect of demagnetization D[K/T] vs T[K], B[T], Tc[K], gyro[rad/s/T]
!> In the demagnetization process $dQ = T\,dS = T(dS/dT)\,dT + T(dS/dB)\,dB$.
!> <br>Then $dT = dQ/C + D\,dB$, where $D = (dS/dB)/(dS/dT)$
      function magn_cw_d(T, B, Tc, gyro) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Tc, gyro
        real*8 ttc, btc, mu, y

        mu = 0.5*gyro*const_hbar
        ttc = T/Tc
        btc = mu*B / (const_kb*Tc)
        y = magn_cw_y(ttc, btc)

        magn_cw_d = - ttc/(y+btc) * mu/const_kB
      end

!><p>Example for Curie-Weiss material with Curie temperature $T_c=0.5$ mK and
!>gyromagnetic ratio $\gamma=203.789\cdot10^6$ rad/s/T:

!><p> <img src="img/magn_cw.png">



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Paramegnetetic material with internal field

!> See Pobell book f9.15

!> Molar magnetization of paramagnetic material, M[J/T/mole] vs T[K], B[T], Bi[T], gyro[rad/s/T], spin[half-int]
      function magn_par_m(T, B, Bi, gyro, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Bi, gyro, spin
        real*8 x,y
        x = gyro*const_hbar*dsqrt(B**2+Bi**2)/(const_kb*T) / 2D0
        y = (2D0*spin + 1D0)*x
        magn_par_m = (2D0*spin+1)/(2D0*spin)/dtanh(y)
     .             - 1D0/(2D0*spin)/dtanh(x)
        magn_par_m = magn_par_m * spin*gyro*const_hbar*const_na
      end

!> Molar magnetic susceptibility of paramagnetic material, chi[J/T^2/mole] vs T[K], B[T], Bi[T], gyro[rad/s/T], spin[half-int]
      function magn_par_chi(T, B, Bi, gyro, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Bi, gyro, spin
        real*8 x,y
        x = gyro*const_hbar*dsqrt(B**2+Bi**2)/(const_kb*T) / 2D0
        y = (2D0*spin + 1D0)*x
        magn_par_m = (2D0*spin+1)/(2D0*spin)/dtanh(y)
     .             - 1D0/(2D0*spin)/dtanh(x)

        magn_par_chi = spin*gyro*const_hbar*const_na
     .               * gyro*const_hbar/(2D0*const_kb*T)
     .               * B/dsqrt(B**2+Bi**2)
     .               * ( - (2D0*spin+1)**2/(2D0*spin)/dsinh(y)**2
     .                   + 1D0/(2D0*spin)/dsinh(x)**2 )
      end

!> Entropy of paramagnetic material, S/R vs T[K], B[T], Bi[T], gyro[rad/s/T], spin[half-int]
      function magn_par_s(T, B, Bi, gyro, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Bi, gyro, spin
        real*8 x,y
        x = gyro*const_hbar*dsqrt(B**2+Bi**2)/(const_kb*T) / 2D0
        y = (2D0*spin + 1D0)*x
        magn_par_s = x/dtanh(x) - y/dtanh(y) + dlog(sinh(y)/sinh(x))
      end

!> Heat capacity of paramagnetic material, C/R vs T[K], B[T], Bi[T], gyro[rad/s/T], spin[half-int]
      function magn_par_c(T, B, Bi, gyro, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Bi, gyro, spin
        real*8 x,y
        x = gyro*const_hbar*dsqrt(B**2+Bi**2)/(const_kb*T) / 2D0
        y = (2D0*spin + 1D0)*x
        magn_par_c = (x/dsinh(x))**2 - (y/dsinh(y))**2
      end

! We want to find D = (dS/dB)/(dS/dT). S depends only on A = T/sqrt(B*B + Bint*Bint)
! Then D = (dA/dB) / (dA/dT) = - T*B/(B*B + Bint*Bint)

!> Cooling effect of demagnetization D[K/T] vs T[K], B[T], Bi[T], gyro[rad/s/T], spin[half-int]
!> In the demagnetization process $dQ = T\,dS = T(dS/dT)\,dT + T(dS/dB)\,dB$.
!> <br>Then $dT = dQ/C + D\,dB$, where $D = (dS/dB)/(dS/dT)$
      function magn_par_d(T, B, Bi, gyro, spin) !F>
        implicit none
        include 'he3.fh'
        real*8 T, B, Bi, gyro, spin
        magn_par_d = - T*B/(B**2 + Bi**2)
      end

!> <p>Example for copper nuclei.  Internal field $B_i = 0.36\cdot 10^{-3}$ T,
!> gyromagnetic ratio $\gamma = 71.118\cdot10^6$ rad/s/T, spin $J$ = 3/2:

!><p> <img src="img/magn_par.png">


