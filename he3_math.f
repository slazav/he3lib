
! Intergation of a real*8 function
! from xmin to xmax using imax points Gaussian quadrature
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

      function math_dint_slatec(func, xmin, xmax, epsabs, epsrel)
        implicit none
        include 'he3.fh'
        real*8 func,xmin,xmax,epsabs,epsrel
        external F
        integer limit, lenw
        parameter (limit = 1000)
        parameter (lenw = 4*LIMIT)
        real*8 abserr,res,work(lenw)
        integer ier,neval,last,iwork(limit)
        call dqags(func,xmin,xmax,epsabs,epsrel,res,
     .      abserr,neval,ier, limit,lenw,last,iwork,work)
        math_dint_slatec = res
      end

! Intergation of complex*16 function
! from xmin to xmax using imax points Gaussian quadrature
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

! 2D intergation of real*8 function
! from xmin to xmax using imax points Gaussian quadrature
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

! 2D intergation of complex*16 function
! from xmin to xmax using imax points Gaussian quadrature
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



