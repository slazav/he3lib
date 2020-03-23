!     Test program:
!     calculate area of a sphere by integrating
!       delta(distance to the surface)
      program test2d
        implicit none
        real*8 func1,func2, s0
        external func1,func2
        integer i
        include "../../he3_math.fh"

        do i=10,120,2
          s0 = math_dint3d_delta(func1, func2,
     .             -1D0, 1D0, i, -1D0, 1D0, i, -1D0, 1D0, i)
          write (*,*) i, s0
        enddo
      end

      function func1(x,y,z)
        implicit none
        real*8 x,y,z, func1
        func1 = 1D0
      end

      function func2(x,y,z)
        implicit none
        real*8 x,y,z, func2
        real*8 x0,y0,z0,r0,sc
        x0=0.11D0
        y0=0.02D0
        y0=0.11D0
        r0=0.7D0
        sc=4D0
        func2 = (r0 - dsqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)) * sc
      end
