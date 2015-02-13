

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
! complete elliptic integral E(x), from IMSL library
      function math_ele(x)
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
! complete elliptic integral K(x), from IMSL library
      function math_elk(x)
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

! magnetic field Bz of a current loop
      function loop_bz(Rl, r, z)
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

! magnetic field Bz of a current loop
      function loop_br(Rl, r, z)
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
