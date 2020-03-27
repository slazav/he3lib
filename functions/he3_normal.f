!HH> Normal 3He liquid parameters beyond zero-temperature limit

! TODO -- Vm, Cv, Cp, ... from ??

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Heat capacity, Cv/R vs T [K], Vm [cm^3/mol]
!> Original formula from Greywall-1983.
!> Note that Cv = Cp up to terms (T/Tf)^3.
      function he3_cv_n(t, v) !F>
        implicit none
        include 'he3.fh'
        real*8 t,v,a,b,c,d
        real*8 s1,s2,s3
        integer i,j
        dimension a(5,4), b(4,3), c(3,3), d(3)
        a(1,1) = -2.9190414D0
        a(1,2) =  5.2893401D+2
        a(1,3) = -1.8869641D+4
        a(1,4) =  2.6031315D+5
        a(2,1) =  0D0
        a(2,2) =  0D0
        a(2,3) =  0D0
        a(2,4) =  0D0
        a(3,1) = -2.4752597D+3
        a(3,2) =  1.8377260D+5
        a(3,3) = -3.4946553D+6
        a(3,4) =  0D0
        a(4,1) =  3.8887481D+4
        a(4,2) = -2.8649769D+6
        a(4,3) =  5.2526785D+7
        a(4,4) =  0D0
        a(5,1) = -1.7505655D+5
        a(5,2) =  1.2809001D+7
        a(5,3) = -2.3037701D+8
        a(5,4) =  0D0

        b(1,1) = -6.5521193D-2
        b(1,2) =  1.3502371D-2
        b(1,3) =  0D0
        b(2,1) =  4.1359033D-2
        b(2,2) =  3.8233755D-4
        b(2,3) = -5.3468396D-5
        b(3,1) =  5.7976786D-3
        b(3,2) = -6.5611532D-4
        b(3,3) =  1.2689707D-5
        b(4,1) = -3.8374623D-4
        b(4,2) =  3.2072581D-5
        b(4,3) = -5.3038906D-7

        c(1,1) = -2.5482958D+1
        c(1,2) =  1.6416936D+0
        c(1,3) = -1.5110378D-2
        c(2,1) =  3.7882751D+1
        c(2,2) = -2.8769188D+0
        c(2,3) =  3.5751181D-2
        c(3,1) =  2.4412956D+1
        c(3,2) = -2.4244083D+0
        c(3,3) =  6.7775905D-2

        d(1) = -7.1613436D+0
        d(2) =  6.0525139D-1
        d(3) = -7.1295855D-3

        if (t.lt.0.1D0) then
          s1=0D0
          do i=1,5
            do j=0,3
              s1 = s1 + a(i,j+1) * v**(-j) * t**i
            enddo
          enddo
          he3_cv_n = s1
          return
        endif

        if (t.ge.0.1D0.and.t.lt.2.5D0) then
          s1=0D0
          do i=0,3
            do j=0,2
              s1 = s1 + b(i+1,j+1) * v**j * t**(-i)
            enddo
          enddo
          s2=0D0
          do i=1,3
            do j=0,2
              s2 = s2 + c(i,j+1) * v**j * t**(-i)
            enddo
          enddo
          s3=0D0
          do j=0,2
            s3 = s3 + d(j+1) * v**j
          enddo
          he3_cv_n = s1 + dexp(-s3/t) * s2
          return
        endif
        he3_cv_n = NaN
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Heat conductivity, K [erg/s cm K] vs T [K] and Vm [cm^3/mol]
! Original formula from Greywall-1984 paper
      function he3_tcond_n_greywall(t, vm)
        implicit none
        real*8 he3_tcond_n_greywall
        include 'he3.fh'
        real*8 t,vm

        if (t.lt.0.05D0) then
          he3_tcond_n_greywall = 1D0/t/(
     .      1D0/(- 4.1884746D1
     .           + 1.9262839D0 * vm)
     .      + t/(- 1.8546379D0
     .           + 2.3695190D-1 *vm
     .           - 6.8284756D-3 *vm**2)
     .   + t**2/(+ 4.3617792D-1
     .           - 4.2101673D-2 *vm
     .           + 1.0050221D-3 *vm**2)
     .   + t**3/(- 9.4328296D-2
     .           + 8.9196267D-3 *vm
     .           - 2.0903165D-4 *vm**2))
        elseif (t.le.1.3D0) then
          he3_tcond_n_greywall =
     .      + 2.5498997D0  /t**2
     .      - 1.1861905D-1 /t**2 * vm
     .      + 1.7187787D-3 /t**2 * vm**2
     .      - 1.4861472D2  /t
     .      + 7.2176329D0  /t * vm
     .      - 7.5439157D-2 /t * vm**2
     .      + 1.0311239D3
     .      - 4.1084636D1  *vm
     .      + 6.8188534D-1 *vm**2
     .      - 3.3746517D3  *t
     .      + 2.2612612D2  *t *vm
     .      - 3.4207801D0  *t *vm**2
     .      + 2.5913792D3  *t**2
     .      - 1.4574998D2  *t**2 *vm
     .      + 2.1389643D0  *t**2 *vm**2
        else
          he3_tcond_n_greywall = NaN
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> He3-n thermal conductivity, K [erg/s cm K] vs T [K] and P [par].
!>
!> Dyugaev-1985. Measurements from Greywall-1984 (7mK-1K) and
!> Kerrisk,Keller-1969 (1.5K-Tcr) are used to obtain some semi-theoretical
!> model for thermal conductivity and viscosity (see below).
!>
!> <img src="img/1984_greywall_tcond.png">
!>
      function he3_tcond_n(t, p) !F>
        implicit none
        include 'he3.fh'
        real*8 t,p
        real*8 tk, k0, dk0

        tk = 0.121830D0 - p*1.459840D-3 + 1.321430D0/(p+10.23327D0)
        k0 = 77.37037D0 - p*0.489282D0 + 623.451D0/(p+16.32097D0)
        dk0 = 3.103759D0 + p*0.0249499D0
     .      - p**2*1.82331D-3 + p**3*4.035088D-5
!        dki = 2.797180D0 - p*0.0107957D0
!     .      + p**2*3.24248D-3 - p**3*1.219298D-4

        if (t.lt.he3_tcr) then
          he3_tcond_n = k0* (tk/t + t/tk + dk0)
        else
          he3_tcond_n = NaN
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Emery factor vs T [K] and P [par].
!>
!> <p>In normal He3 transport properties (viscosity, spin diffusion, thermal conductivity?)
!> are suppressed just above Tc because of some fluctuation effects (Emery-1978).
!> This is clearly seen at viscosity measurements (Parpia-1978, Carless-1983, Nakagawa-1996).
!> The factor has form 1 - $G(1 - \theta/\alpha*\atan(\altha/\theta))$, where $\theta = \sqrt{T/T_c - 1}$.
!>
!> <p>Values of $G$ (pressure independent) and $\alpha$ were obtained by fitting
!> viscosity data from Carless-1983 and Nakagawa-1996 in assumption that at high
!> temperatures they should follow Dyugaev-1985 model.
!>
!> <p>Note that in Carless-1983 temperature scale Alvesalo-80 is used. To convert
!> it to Greywall-86 scale one should multiply temperature by $k=0.893$.
      function he3_emery_factor(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p, A,G,th

        if (ttc.lt.1D0) then
          he3_emery_factor = NaN
        else
          ! Fitting done in data/1983_carless_visc/process_data.m
          G = 2.318858D-1
          A = -4.067317D-2 + 4.841729D1/(6.0755 + p)
          th = dsqrt(ttc-1D0)
          he3_emery_factor  = 1D0 - G*(1D0 - th/A*datan(A/th));
        endif
      end


!> He3-n viscosity, eta [poise] vs T [K] and Vm [cm^3/mol].
!> Pure Dyugaev-1985 model without Emery effect. See function <tt>he3_visc_n</tt> below.
!>
      function he3_visc_n0(t, p) !F>
        implicit none
        include 'he3.fh'
        real*8 t,p
        real*8 te, e0

        te = 0.064795D0 + p*2.949606D-5 + 5.351244D0/(p+16.87556D0)
        e0 = 22.3125D0  + p*0.375D0 - 36.09375D0/(p+7.5D0)

        if (t.lt.he3_tcr) then
          he3_visc_n0 = 1D-6*e0* ((te/t)**2 + 1.41D0*te/t + 1D0)
        else
          he3_visc_n0 = NaN
        endif
      end

!> He3-n viscosity, eta [poise] vs T [K] and Vm [cm^3/mol].
!>
!> Model from Dyugaev-1985, it uses thermal conductivity experimental data to get viscosity.
!> At low temperature Emery effect, reduction of viscosity close to $T_c$ due to fluctuation
!> effects, is taken into account.
!> Very good agreement with Betts-1963,1965 at high temperatures, and with Carless-1983, Nakagawa-1996
!> at low tempeatures (temperature scale correction for Carless-1983 is needed).
!>
!> <p><img src="img/1965_betts_visc_fig1.png">
!> <p><img src="img/1983_carless_visc_fig4.png">
!>
!> <p> There is also a complete viscosity model in Huang-2012, but for me it does
!> not look as good as this one.
!>
      function he3_visc_n(t, p) !F>
        implicit none
        include 'he3.fh'
        real*8 t,p
        he3_visc_n = he3_visc_n0(t,p)
     .             * he3_emery_factor(t*1D3/he3_tc(p),p)
      end
