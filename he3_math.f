

! Integration of a real*8 function
! from xmin to xmax using imax points
! Gauss-2pt quadrature
      function math_dint(func, xmin, xmax, nx)
        implicit none
        include 'he3.fh'
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
        include 'he3.fh'
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
        include 'he3.fh'
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
        include 'he3.fh'
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

!! Adaptive integration of a real*8 function
!! using slatec library
!      function math_dint_slatec(func, xmin, xmax, epsabs, epsrel)
!        implicit none
!        include 'he3.fh'
!        real*8 func,xmin,xmax,epsabs,epsrel
!        external F
!        integer limit, lenw
!        parameter (limit = 1000)
!        parameter (lenw = 4*LIMIT)
!        real*8 abserr,res,work(lenw)
!        integer ier,neval,last,iwork(limit)
!        call dqags(func,xmin,xmax,epsabs,epsrel,res,
!     .      abserr,neval,ier, limit,lenw,last,iwork,work)
!        math_dint_slatec = res
!      end

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


! 2D integration of real*8 function from xmin to xmax using imax points
! Gauss-7pt+Kronrad-15pt quadrature with error estimation
      function math_dint2d_gk(func,
     .         xmin, xmax, nx, ymin, ymax, ny, aerr)
        implicit none
        include 'he3.fh'

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
        include 'he3.fh'

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
