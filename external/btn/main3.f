c********************************************************************
c BTN: Sample problem 3 (complex usage of BTN, CHKDER)
c main program
c last changed: 10/10/90
c********************************************************************
c
      program          main
      implicit         double precision (a-h,o-z)
      parameter        (nn = 100, ndmx = 10,
     *                 lw = (3+7*ndmx)*(nn+ndmx))
      double precision x(nn), g(nn), w(lw)
      external         sfun
c
c Read data from file
c
      open (7, file = 'main3.d', status = 'old')
      read (7,800) n
c
c Set default values of parameters
c
      call chkpar (n, kmax, msglvl, iunit, istart, iend, istep)
      call btnpar (nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *    initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *    ireset, indstp)
c
c Read in new values of selected parameters
c
      read (7,800) ichk
      read (7,800) istart
      read (7,800) iend
      read (7,800) istep
      read (7,800) msglvl
      read (7,800) iunit
      read (7,800) kmax
c
      read (7,810) rnktol
      read (7,810) accrcy
c
c check derivative values (if desired)
c
      call xstart (x, n)
      if (ichk .ne. 0) then
          call chkder (n, x, f, g, w, lw, sfun, iflag,
     *                 kmax, errmax, imax,
     *                 msglvl, iunit, istart, iend, istep)
          write (iunit,815) errmax, imax
          stop
      end if
c
c solve optimization problem
c
      call btn (n, x, f, g, w, lw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, sfun)
c
c write results to file
c
      write (iunit,820)
      write (iunit,830) iflag
      write (iunit,840) f
      write (iunit,850)
      write (iunit,860) (x(i), i = 1,n)
      write (iunit,870)
      write (iunit,860) (g(i), i = 1,n)
c
      stop
800   format (10x, i5)
810   format (10x, d16.4)
815   format (' CHKDER: Max. error in gradient = ', 1pd12.4,
     *        ' at index ', i5)
820   format (' Results of BTN', /)
830   format (' Error code = ', i4)
840   format (' Function value = ', 1pd24.16)
850   format (/, ' x = ')
860   format (3(1pd16.8))
870   format (/, ' g = ')
      end
c
c----------------------------------------------------------------------
c initial guess
c----------------------------------------------------------------------
c
      subroutine xstart (x, n)
      double precision x(n)
c
c specify initial values of the variables
c
      do 10 i = 1,n
          x(i) = 0.d0
10    continue
c
      return
      end
c
c----------------------------------------------------------------------
c function evaluation
c----------------------------------------------------------------------
c
      subroutine sfun (n ,x, f, g)
c
c evaluate the nonlinear function and its gradient
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*)
C$    PRIVATE
c
c compute function value and gradient
c
      f = 0.d0
      do 30 i = 1,n
          b = i
          if (i .gt. n/2) b = -b
          f = f + 5.d-1*i*x(i)*x(i) - b*x(i)
          g(i) = i*x(i) - b
30    continue
c
      return
      end
