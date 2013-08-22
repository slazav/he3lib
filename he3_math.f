
! Intergation of a real*8 function
! from xmin to xmax using imax points Gaussian quadrature
      function math_dint(func, xmin, xmax, maxi, args)
        implicit none
        include 'he3.fh'
        real*8 func, xmin, xmax, args(*)
        real*8 dx, xp, xm
        integer i, maxi

        dx=(xmax-xmin)/dble(maxi)
        math_dint = 0D0
        do i=1,maxi
          xp = xmin + dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = xmin + dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          math_dint = math_dint
     .       + (func(xp, args) + func(xm, args)) * dx/2D0
        enddo
      end

! Intergation of complex*16 function
! from xmin to xmax using imax points Gaussian quadrature
      function math_cint(func, xmin, xmax, maxi, args)
        implicit none
        include 'he3.fh'
        complex*16 func
        real*8 xmin, xmax, args(*)
        real*8 dx, xp, xm
        integer i, maxi

        dx=(xmax-xmin)/dble(maxi)
        math_cint = 0D0
        do i=1,maxi
          xp = xmin + dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = xmin + dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          math_cint = math_cint
     .       + (func(xp, args) + func(xm, args))
     .           * dcmplx(dx/2D0, 0D0)
        enddo
      end
