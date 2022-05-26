!HH> Mathematics

! Integration of a real*8 function
! from xmin to xmax using imax points
! Gauss-2pt quadrature
      function math_dint(func, xmin, xmax, nx)
        implicit none
        include 'he3_math.fh'
        real*8 func, xmin, xmax
        real*8 dx, xp, xm
        integer i, nx

        dx=(xmax-xmin)/dble(nx)
        math_dint = 0D0
        do i=1,nx
          xp = xmin + dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = xmin + dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          math_dint = math_dint
     .       + (func(xp) + func(xm)) * dx/2D0
        enddo
      end

! Integration of complex*16 function
! from xmin to xmax using imax points
! Gauss-2pt quadrature
      function math_cint(func, xmin, xmax, nx)
        implicit none
        include 'he3_math.fh'
        complex*16 func
        real*8 xmin, xmax
        real*8 dx, xp, xm
        integer i, nx

        dx=(xmax-xmin)/dble(nx)
        math_cint = (0D0, 0D0)
        do i=1,nx
          xp = xmin + dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = xmin + dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          math_cint = math_cint
     .       + (func(xp) + func(xm))
     .           * dcmplx(dx/2D0, 0D0)
        enddo
      end

! 2D integration of real*8 function
! from xmin to xmax using imax points
! Gauss-2pt quadrature
      function math_dint2d(func,
     .         xmin, xmax, nx, ymin, ymax, ny)
        implicit none
        include 'he3_math.fh'
        real*8 func
        real*8 xmin, xmax, ymin, ymax
        real*8 s1, s2, dx, xp, xm, dy, yp, ym
        integer ix, iy, nx, ny
        s1 = 0.5D0 - 0.5D0/dsqrt(3D0)
        s2 = 0.5D0 + 0.5D0/dsqrt(3D0)
        dx=(xmax-xmin)/dble(nx)
        dy=(ymax-ymin)/dble(ny)
        math_dint2d = 0D0
        do ix=1,nx
          do iy=1,ny
            xp = xmin + dx * (dble(ix) - s1)
            xm = xmin + dx * (dble(ix) - s2)
            yp = ymin + dy * (dble(iy) - s1)
            ym = ymin + dy * (dble(iy) - s2)
            math_dint2d = math_dint2d
     .         + (func(xm, ym) + func(xp, ym)
     .          + func(xm, yp) + func(xp, yp))
     .             * dx*dy/4D0
          enddo
        enddo
      end

! 2D integration of complex*16 function
! from xmin to xmax using imax points
! Gauss-2pt quadrature
      function math_cint2d(func,
     .         xmin, xmax, nx, ymin, ymax, ny)
        implicit none
        include 'he3_math.fh'
        complex*16 func
        real*8 xmin, xmax, ymin, ymax
        real*8 s1, s2, dx, xp, xm, dy, yp, ym
        integer ix, iy, nx, ny

        s1 = 0.5D0 - 0.5D0/dsqrt(3D0)
        s2 = 0.5D0 + 0.5D0/dsqrt(3D0)
        dx=(xmax-xmin)/dble(nx)
        dy=(ymax-ymin)/dble(ny)
        math_cint2d = (0D0, 0D0)
        do ix=1,nx
          do iy=1,ny
            xp = xmin + dx * (dble(ix) - s1)
            xm = xmin + dx * (dble(ix) - s2)
            yp = ymin + dy * (dble(iy) - s1)
            ym = ymin + dy * (dble(iy) - s2)
            math_cint2d = math_cint2d
     .         + (func(xm, ym) + func(xp, ym)
     .          + func(xm, yp) + func(xp, yp))
     .             * dcmplx(dx*dy/4D0, 0D0)
          enddo
        enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Coordinates and weights for Gauss-7pt+Kronrod-15pt quadrature:
      block data math_intgk_block
        common /math_intgk/ crd,wg,wk
        real*8 crd(15), wg(7), wk(15)
        data
     .    crd /-0.991455371120813D0, -0.949107912342759D0,
     .         -0.864864423359769D0, -0.741531185599394D0,
     .         -0.586087235467691D0, -0.405845151377397D0,
     .         -0.207784955007898D0,  0.000000000000000D0,
     .          0.207784955007898D0,  0.405845151377397D0,
     .          0.586087235467691D0,  0.741531185599394D0,
     .          0.864864423359769D0,  0.949107912342759D0,
     .          0.991455371120813D0/,
     .    wg /0.129484966168870D0, 0.279705391489277D0,
     .        0.381830050505119D0, 0.417959183673469D0,
     .        0.381830050505119D0, 0.279705391489277D0,
     .        0.129484966168870D0/,
     .    wk /0.022935322010529D0, 0.063092092629979D0,
     .        0.104790010322250D0, 0.140653259715525D0,
     .        0.169004726639267D0, 0.190350578064785D0,
     .        0.204432940075298D0, 0.209482141084728D0,
     .        0.204432940075298D0, 0.190350578064785D0,
     .        0.169004726639267D0, 0.140653259715525D0,
     .        0.104790010322250D0, 0.063092092629979D0,
     .        0.022935322010529D0/
      end

! 1D integration of real*8 function from xmin to xmax using imax points
! Gauss-7pt+Kronrad-15pt quadrature with error estimation
      function math_dint_gk(func,
     .         xmin, xmax, nx, aerr)
        implicit none
        include 'he3_math.fh'

        common /math_intgk/ crd,wg,wk
        real*8 crd(15), wg(7), wk(15)
        real*8 func
        real*8 xmin, xmax, aerr
        real*8 dx, x0, f, intk, intg
        integer ix, ixq, nx
        external func
        dx=(xmax-xmin)/dble(nx)
        intk=0D0
        intg=0D0
        do ix=1,nx
          x0 = xmin + dx*(dble(ix)-0.5D0)
          do ixq=1,15
            f=func(x0 + 0.5D0*dx*crd(ixq))
            intk=intk + wk(ixq)*f
            if (mod(ixq,2).eq.0) then
              intg=intg + wg(ixq/2)*f
            endif
          enddo
        enddo
        intk=intk*dx/2D0
        intg=intg*dx/2D0
        aerr=(200D0*dabs(intk-intg))**1.5D0
        math_dint_gk = intk
      end

! 1D adaptive integration of real*8 function
! Gauss-7pt+Kronrad-15pt quadrature
! one of aerr_lim and rerr_lim should be positive
      function math_dint_gka(func,
     .    xmin, xmax, aerr_lim, rerr_lim)
        ! wrapper for the recoursive adaptive integration
        implicit none
        include 'he3_math.fh'
        real*8 xmin, xmax, res
        real*8 aerr_lim, rerr_lim, aerr,rerr
        real*8 func
        external func, math_dint_gka_int

        ! calculate first approximation of the integral
        math_dint_gka = math_dint_gk(func, xmin, xmax, 10, aerr)
        rerr=dabs(aerr/math_dint_gka)

        ! finish if this accuracy is enough
        if (aerr.lt.aerr_lim) return
        if (rerr.lt.rerr_lim) return

        ! calculate desired absolute error limit
        ! (we want recoursion with only absolute limit because
        !  of zero regions in our functions)
        aerr = aerr_lim
        rerr = dabs(rerr_lim*math_dint_gka)
        if (rerr.gt.aerr) aerr = rerr

        ! run recoursion
        math_dint_gka=0D0
        call math_dint_gka_int(math_dint_gka_int, func,
     .    xmin, xmax, aerr, math_dint_gka, 0)
      end

!     internal subroutine for the integration
      subroutine math_dint_gka_int(myself, func,
     .         xmin, xmax, aerr_lim, res, depth)
        implicit none
        include 'he3_math.fh'

        common /math_intgk/ crd,wg,wk
        real*8 crd(15), wg(7), wk(15)

        real*8 func
        real*8 xmin, xmax, res
        real*8 dx,x0, f, intk, intg
        real*8 aerr,rerr, aerr_lim
        real*8 NaN
        integer ixq, depth
        external func

        intk=0D0
        intg=0D0
        NaN=0D0/0D0
        x0 = (xmin + xmax)/2D0
        if (x0.eq.xmin.or.x0.eq.xmax) return ! not enough accuracy dx<<x
        dx=(xmax-xmin)/2D0
        do ixq=1,15
          f=func(x0 + dx*crd(ixq))
          intk=intk + wk(ixq)*f
          if (mod(ixq,2).eq.0) then
            intg=intg + wg(ixq/2)*f
          endif
        enddo
        intk=intk*dx
        intg=intg*dx
        aerr=(200D0*dabs(intk-intg))**1.5D0
        ! write(*,*) '>', xmax-xmin, intk, aerr, aerr_lim

        if (aerr_lim.le.0D0.or.aerr.le.aerr_lim) then
          res = res + intk
          return
        endif
        if (depth.ge.50) then
          write(*,*) 'math_dint_gka warning: ',
     .               'max depth is reached: ', depth
          return
        endif
        call myself(myself, func,
     .       xmin, x0, aerr_lim/2D0, res, depth+1)
        call myself(myself, func,
     .       x0, xmax, aerr_lim/2D0, res, depth+1)
        return
      end


! 2D integration of real*8 function from xmin to xmax using imax points
! Gauss-7pt+Kronrad-15pt quadrature with error estimation
      function math_dint2d_gk(func,
     .         xmin, xmax, nx, ymin, ymax, ny, aerr)
        implicit none
        include 'he3_math.fh'

        common /math_intgk/ crd,wg,wk
        real*8 crd(15), wg(7), wk(15)
        real*8 func
        real*8 xmin, xmax, ymin, ymax, aerr
        real*8 dx,dy, x0,y0, f, intk, intg
        integer ix, iy, ixq, iyq, nx, ny
        external func
        dx=(xmax-xmin)/dble(nx)
        dy=(ymax-ymin)/dble(ny)
        intk=0D0
        intg=0D0
        do ix=1,nx
          do iy=1,ny
            x0 = xmin + dx*(dble(ix)-0.5D0)
            y0 = ymin + dy*(dble(iy)-0.5D0)
            do ixq=1,15
              do iyq=1,15
                f=func(x0 + 0.5D0*dx*crd(ixq),
     .                 y0 + 0.5D0*dy*crd(iyq))
                intk=intk + wk(ixq)*wk(iyq)*f
                if (mod(ixq,2).eq.0.and.mod(iyq,2).eq.0) then
                  intg=intg + wg(ixq/2)*wg(iyq/2)*f
                endif
              enddo
            enddo
          enddo
        enddo
        intk=intk*dx*dy/4D0
        intg=intg*dx*dy/4D0
        aerr=(200D0*dabs(intk-intg))**1.5D0
        math_dint2d_gk = intk
      end

! 2D adaptive integration of real*8 function
! Gauss-7pt+Kronrad-15pt quadrature
      subroutine math_dint2d_gka(myself, func,
     .         xmin, xmax, ymin, ymax, aerr_lim, rerr_lim, res)
        implicit none
        include 'he3_math.fh'

        common /math_intgk/ crd,wg,wk
        real*8 crd(15), wg(7), wk(15)

        real*8 func
        real*8 xmin, xmax, ymin, ymax, res
        real*8 dx,dy, x0,y0, f, intk, intg
        real*8 aerr,rerr, aerr_lim, rerr_lim
        integer ixq, iyq
        external func

        intk=0D0
        intg=0D0
        x0 = (xmin + xmax)/2D0
        y0 = (ymin + ymax)/2D0
        dx=(xmax-xmin)/2D0
        dy=(ymax-ymin)/2D0
        do ixq=1,15
          do iyq=1,15
            f=func(x0 + dx*crd(ixq),
     .             y0 + dy*crd(iyq))
            intk=intk + wk(ixq)*wk(iyq)*f
            if (mod(ixq,2).eq.0.and.mod(iyq,2).eq.0) then
              intg=intg + wg(ixq/2)*wg(iyq/2)*f
            endif
          enddo
        enddo
        intk=intk*dx*dy
        intg=intg*dx*dy
        aerr=(200D0*dabs(intk-intg))**1.5D0
        rerr=aerr/intk
        if (aerr_lim.gt.0D0.and.dabs(aerr).le.aerr_lim) then
          res = res + intk
!        write(*,*) xmax-xmin, ymax-ymin, xmin, ymin, intk, rerr
          return
        endif
        if (rerr_lim.gt.0D0.and.dabs(rerr).le.rerr_lim) then
          res = res + intk
!        write(*,*) xmax-xmin, ymax-ymin, intk, rerr
          return
        endif
        call myself(myself, func,
     .       xmin, x0, ymin, y0,
     .       aerr_lim, rerr_lim, res)
        call myself(myself, func,
     .       x0, xmax, ymin, y0,
     .       aerr_lim, rerr_lim, res)
        call myself(myself, func,
     .       xmin, x0, y0, ymax,
     .       aerr_lim, rerr_lim, res)
        call myself(myself, func,
     .       x0, xmax, y0, ymax,
     .       aerr_lim, rerr_lim, res)
        return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate int( f1 \delta(f2) )
! Discretization of Dirac delta functions in level set methods
! Journal of Computational Physics, Volume 207, Issue 1, 20 July 2005, Pages 28-51
! Bjo:rn Engquist, Anna-Karin Tornberg, Richard Tsai

      function math_dint2d_delta(func1, func2,
     .         xmin, xmax, nx, ymin, ymax, ny)
        implicit none
        include 'he3_math.fh'
        integer nx,ny,ix,iy
        real*8  xmin, xmax, ymin, ymax
        real*8  dx,dy, x0,y0, f1,f2,f2x,f2y,e,de
        real*8  func1, func2
        external func1, func2

        dx=(xmax-xmin)/dble(nx)
        dy=(ymax-ymin)/dble(ny)
        math_dint2d_delta=0D0
        do ix=1,nx
          do iy=1,ny
            x0 = xmin + dx*(dble(ix)-0.5D0)
            y0 = ymin + dy*(dble(iy)-0.5D0)
            f1  = func1(x0,y0)
            f2  = func2(x0,y0)
            f2x = (func2(x0+dx/2D0,y0)-func2(x0-dx/2D0,y0))/dx
            f2y = (func2(x0,y0+dy/2D0)-func2(x0,y0-dy/2D0))/dy
            ! in the case of f with |nabla|_2 = 1 (distance)
            ! e0=dx=dy is good. For large gradient of f we want smaller e0
            ! e0 -> e0/|nabla|_2
            ! e =  |nabla|_1/|nabla|_2 * e0
            e   = (abs(f2x)+abs(f2y)) * (dx+dy)/2D0
            de = (1-abs(f2/e))/e ! hat approximation of delta
            if (de.gt.0D0) then
              math_dint2d_delta = math_dint2d_delta
     .                          + f1*de*dx*dy
            endif
          enddo
        enddo
      end

      function math_dint3d_delta(func1, func2,
     .         xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz)
        implicit none
        include 'he3_math.fh'
        integer nx,ny,nz,ix,iy,iz
        real*8  xmin, xmax, ymin, ymax, zmin, zmax
        real*8  dx,dy,dz, x0,y0,z0, f1,f2,f2x,f2y,f2z,e,de
        real*8  func1, func2
        external func1, func2

        dx=(xmax-xmin)/dble(nx)
        dy=(ymax-ymin)/dble(ny)
        dz=(zmax-zmin)/dble(nz)
        math_dint3d_delta=0D0
        do ix=1,nx
          do iy=1,ny
            do iz=1,nz
              x0 = xmin + dx*(dble(ix)-0.5D0)
              y0 = ymin + dy*(dble(iy)-0.5D0)
              z0 = zmin + dz*(dble(iz)-0.5D0)
              f1  = func1(x0,y0,z0)
              f2  = func2(x0,y0,z0)
              f2x = (func2(x0+dx/2D0,y0,z0)-func2(x0-dx/2D0,y0,z0))/dx
              f2y = (func2(x0,y0+dy/2D0,z0)-func2(x0,y0-dy/2D0,z0))/dy
              f2z = (func2(x0,y0,z0+dz/2D0)-func2(x0,y0,z0-dz/2D0))/dz
              e   = (abs(f2x)+abs(f2y)+abs(f2z)) * (dx+dy+dz)/3D0
              de = (1-abs(f2/e))/e ! hat approximation of delta
              if (de.gt.0D0) then
                math_dint3d_delta = math_dint3d_delta
     .                            + f1*de*dx*dy*dz
              endif
            enddo
          enddo
        enddo
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> complete elliptic integral E(x)
!> from IMSL library
      function math_ele(x) !F>
        implicit none
        real*8 x, A(10),B(10), suma, sumb, e, eta, math_ele
        integer i

        data A /0.1494662175718132D-03, 0.2468503330460722D-02,
     .          0.8638442173604073D-02, 0.1077063503986645D-01,
     .          0.7820404060959553D-02, 0.7595093422559432D-02,
     .          0.115695957452954D-01,  0.2183181167613048D-01,
     .          0.5680519456755915D-01, 0.4431471805608895D00/
        data B /0.3185919565550157D-04, 0.9898332846225384D-03,
     .          0.6432146586438302D-02, 0.1680402334636338D-01,
     .          0.2614501470031388D-01, 0.3347894366576162D-01,
     .          0.4271789054738309D-01, 0.5859366125553149D-01,
     .          0.9374999972120313D-01, 0.2499999999999017D00/
        if (x.lt.0D0.or.x.gt.1D0) then
          math_ele = 0D0/0D0 !NaN
          return
        endif
        e = 1D-12 ! precision
        eta = 1D0 - x
        if (eta .GE. e) then
           suma = 0D0
           sumb = 0D0
           do I=1,10
             suma = (suma+A(I))*eta
             sumb = (sumb+B(I))*eta
           enddo
           math_ele = suma - dlog(eta)*sumb
           math_ele = math_ele + 1D0 + e
        else
           math_ele = 1D0
        endif
      end

!> complete elliptic integral K(x)
!> from IMSL library
      function math_elk(x) !F>
        implicit none
        real*8 x, A(11),B(11), suma, sumb, e, eta, math_elk
        integer i
        data A /0.1393087857006646D-03, 0.2296634898396958D-02,
     .          0.8003003980649985D-02, 0.9848929322176892D-02,
     .          0.6847909282624505D-02, 0.6179627446053317D-02,
     .          0.8789801874555064D-02, 0.1493801353268716D-01,
     .          0.3088514627130518D-01, 0.9657359028085625D-01,
     .          0.138629436111989D01/
        data B /0.2970028096655561D-04, 0.9215546349632497D-03,
     .          0.5973904299155429D-02, 0.155309416319772D-01,
     .          0.2393191332311079D-01, 0.3012484901289892D-01,
     .          0.373777397586236D-01, 0.4882804190686239D-01,
     .          0.7031249973903835D-01, 0.124999999999908D00, 0.5D00/
        if (x.lt.0D0.or.x.gt.1D0) then
          math_elk = 0D0/0D0 !NaN
          return
        endif
        e = 1D-12 ! precision
        eta = 1D0 - x
        if (eta .GE. e) then
          suma = A(1)
          sumb = B(1)
          do I=2,11
            suma = suma*eta + A(I)
            sumb = sumb*eta + B(I)
          enddo
          math_elk = suma - dlog(eta)*sumb
        else
          math_elk = A(11) - dlog(eta)*B(11)
        endif
      end

!> magnetic field Bz of a current loop
      function loop_bz(Rl, r, z) !F>
        implicit none
        include 'he3.fh'
        real*8 Rl,r,z
        real*8 k,ele,elk,pre

        k = 4D0*Rl*r/((Rl+r)**2+z**2)
        ele = math_ele(k)
        elk = math_elk(k)

        pre = const_mu0/(2D0*const_pi)/dsqrt((Rl+r)**2+z**2)
        loop_bz = pre *
     .    (ele*(Rl**2-r**2-z**2)/((Rl-r)**2+z**2) + elk)
      end

!> magnetic field Bz of a current loop
      function loop_br(Rl, r, z) !F>
        implicit none
        include 'he3.fh'
        real*8 Rl,r,z
        real*8 k,ele,elk,pre

        k = 4D0*Rl*r/((Rl+r)**2+z**2)
        ele = math_ele(k)
        elk = math_elk(k)

        pre = const_mu0/(2D0*const_pi)/dsqrt((Rl+r)**2+z**2)
        loop_br = pre * z/r *
     .    (ele*(Rl**2+r**2+z**2)/((Rl-r)**2+z**2) - elk)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve a cubic equation A3 x^3 + A2 x^2 + A1*x + A0 = 0
!
! Calculation of Densities from Cubic Equations of State: Revisited
! Ulrich K. Deiters*, Ricardo Macias-Salina
! Ind. Eng. Chem. Res., 2014, 53 (6), pp 2529�2536

      subroutine solve_cubic(A3,A2,A1,A0, x1,x2,x3)
        implicit none
        real*8 A3,A2,A1,A0, x1,x2,x3
        real*8 lm, b0,b1,b2, t, xi,yi, y,yp,ypp,dx, c0,c1, D

        if (A3.eq.0D0) then
          call solve_quadr(A2,A1,A0, x1,x2)
          x3=1D0/0D0 ! NaN
          return
        endif

        ! normalize and scaling
        b0 = A0/A3
        b1 = A1/A3
        b2 = A2/A3

        lm = dabs(b0)**(1D0/3D0)
        if (lm.lt.sqrt(dabs(b1))) lm=sqrt(dabs(b1))
        if (lm.lt.dabs(b2))       lm=dabs(b2)
        b0 = b0/lm**3
        b1 = b1/lm**2
        b2 = b2/lm

        if (b0.eq.0D0) then
          x1 = 0D0
          call solve_quadr(1D0,b2,b1, x2,x3)
          if (x1.gt.x2) then ! swap x1 and x2
            t=x2
            x2=x1
            x1=t
          endif
          if (x2.gt.x3) then ! swap x2 and x3
            t=x3
            x3=x2
            x2=t
          endif
          goto 30
        endif

        ! inflection point
        xi = -b2/3D0
        yi = xi**3 + b2*xi**2 + b1*xi + b0
        if (yi.eq.0D0) then
          x2 = xi
          c1 = x2+b2
          c0 = c1*x2 + b1
          call solve_quadr(1D0,c1,c0, x1,x3)
          goto 30
        endif

        D = b2**2 - 3D0*b1
        if (D.eq.0D0) then
          x1 = xi - yi**(1D0/3D0)
          x2 = 1D0/0D0 ! NaN
          x3 = 1D0/0D0 ! NaN
          goto 30
        endif

        x1 = xi
        if (D.gt.0D0) x1 = xi - yi/abs(yi) *2D0/3D0*dsqrt(D)


        ! LOOP
20      y  = x1**3 + b2*x1**2 + b1*x1 + b0
        yp = 3D0*x1**2 + 2D0*b2*x1 + b1
        ypp = 6D0*x1 + 2D0*b2
        dx = y*yp/(yp**2-0.5D0*y*ypp)
        x1 = x1 - dx
        if (dabs(dx/x1).gt.1D-10) goto 20

        if (D.gt.0D0) then
          c1 = x1 + b2
          c0 = c1*x1 + b1
          call solve_quadr(1D0,c1,c0, x2,x3)
          if (x1.gt.x2) then ! swap x1 and x2
            t=x2
            x2=x1
            x1=t
          endif
          if (x2.gt.x3) then ! swap x2 and x3
            t=x3
            x3=x2
            x2=t
          endif
          goto 30
        endif
        x2 = 1D0/0D0 ! NaN
        x3 = 1D0/0D0 ! NaN

30      x1=x1*lm
        x2=x2*lm
        x3=x3*lm

        return
      end

! solve quadratic equation A2 x^2 + A1*x + A0 = 0
      subroutine solve_quadr(A2,A1,A0, x1,x2)
        implicit none
        real*8 A2,A1,A0, x1,x2
        real*8 D
        D = A1**2 - 4D0*A2*A0
        if (D.lt.0D0) then
          x1 = 1D0/0D0 ! NaN
          x2 = 1D0/0D0 ! NaN
          return
        endif
        x1 = (-A1-dsqrt(D))/(2D0*A2)
        x2 = (-A1+dsqrt(D))/(2D0*A2)
        return
      end


!>  Evaluates complex Stokes function K + 1i*K'
!>  Uses the methods outlined in STOKES - Mathematical and physical papers Vol III.
!>  For G>=3 use eq.113, for G<3 eq.103-105.
!>  Code is taken from Lancaster ULT wire calibration program.
      function math_stokes(g)  !FC>
        implicit none
        include 'he3.fh'
        real*8 g, k, k1
        real*8 M0, M1, AL, E0, E1, G1
        real*8 A, B, C, D, E, G2, G4, G5, G6, G8, G10, G12
        if (G. GE. 3.0) then
          G2=G*G
          G4=G2*G2
          G5=G4*G
          A=2.828427D0  / G
          B=.3535534D0  / G /G2
          C=.552427D0   / G5
          D=2.9637719D0 / G5 / G2
          E=20.4632383D0 / G5 / G4
          K = 1D0 + A + B - 0.5D0 / G4
     +        + C - D + 12.8875857D0/G4/G4 - E
          K1 = A + 2D0/G2 - B + C - 1.625D0/G4/G2 + D - E
        else
          G1=G/2D0
          G2=G1*G1
          G4=G2*G2
          G6=G4*G2
          G8=G4*G4
          G10=G6*G4
          G12=G6*G6
          M0 = G2 - G6/12D0 + G10/2880D0
          M1 = G2 - G6/36D0 + G10/14400D0
          E0 = G4/2. - G8/144. +G12/86400D0
          E1 = G4/4. - G8/576. +G12/518400D0
          AL = 0.5772158D0 + dlog(G1)
          A = -(AL*M0) + const_pi/4D0*E0 - 0.5D0*M1
     +        + (G2 - G6/6.545454D0 + G10/1261.313D0)
          B = const_pi/4D0*M0 + AL*E0 - 0.5D0*(1D0-E1)
     +        - (G4*0.75D0 - G8/69.12D0 + G12/35265.31D0)
          C = - (const_pi/4D0*M1) + AL*(1D0-E1)
     +        + (G4/2.666666D0 - G8/276.48D0 + G12/211591.84D0)
          D = -(AL*M1) - const_pi/4D0*(1D0-E1)
     +        + (G2 - G6/19.636363D0 +G10/6306.57D0)
          K  = 1D0 + 2D0*(A*C + B*D)/(G2*(C*C + D*D))
          K1 = 2D0*(B*C - A*D)/(G2*(C*C + D*D))
        endif
        math_stokes = dcmplx(K,K1)
      end

!> Stokes K function, real component of math_stokes(g).
      function math_stokes_k(g) !F>
        implicit none
        include 'he3.fh'
        real*8 g
        math_stokes_k = real(math_stokes(g))
      end

!> Stokes K' function, imag component of math_stokes(g).
      function math_stokes_kp(g) !F>
        implicit none
        include 'he3.fh'
        real*8 g
        math_stokes_kp = imag(math_stokes(g))
      end
