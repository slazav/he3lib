      program test_dint
        implicit none
        real*8 func, k, s0,s1
        external func
        integer i
        common /func_par/ k
        include "../he3_math.fh"
        do i=1,1000
          k = dble(i)/10D0
          s0 = 2D0 * datan(k)/k
          s1 = math_dint(func, -1D0, 1D0, 50)
          write (*,*) k, (s0-s1)/s0
        enddo
      end

!     test function
!      - f(x) = 1/(1+(kx)^2)
!      - int_{-1}^{1} f(x)dx = 2/k atan(k)
!
      function func(x)
        implicit none
        real*8 x, k, func
        common /func_par/ k
        func = 1D0/(1D0+(k*x)**2)
      end
