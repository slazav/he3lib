c********************************************************************
c BTN: Sample problem 4 (using BTN to solve a constrained problem)
c
c A logarithmic barrier method is used to minimize a nonlinear
c function subject to simple bounds on the variables ("box" constraints)
c
c main program
c last changed: 05/20/91
c********************************************************************
c
      program          main
      implicit         double precision (a-h,o-z)
      parameter        (n = 100, ndmx = 10,
     *                 lw = (3+7*ndmx)*(n+ndmx))
      double precision x(n), g(n), w(lw), lb(n), ub(n), mu
      common /bnd/     lb, ub
      common /bar/     mu
      external         sfun, maxstp
c
c Initialize nonlinear parameters x and barrier parameter mu
c
      call xstart (x, n)
      mu    = 1.d0
c
c set up parameters for derivative checker
c
      call chkpar (n, kmax, msglvl, iunit, istart, iend, istep)
c
c check derivative values at initial point
c
      call chkder (n, x, f, g, w, lw, sfun, iflag,
     *             kmax, errmax, imax,
     *             msglvl, iunit, istart, iend, istep)
      if (iflag .ne. 0) stop
c
c Set customizing parameters for BTN (the user-supplied routine
c maxstp will be used to bound the step length in the line search)
c
      call btnpar (nodemx, kmax, maxit, msglvl, iprec, nlmax, initv,
     *          tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *          ireset, indstp)
      indstp = 1
c
c solve optimization problem
c
      scale = 10.d0
10    write (*,860) mu
      call btn (n, x, f, g, w, lw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, maxstp)
c
c Print results
c
      if (iflag .eq. 999) then
          write (*,*) ' Fatal error in BTN'
          stop
      end if
c
c update barrier parameter
c
      mu = mu / scale
      if (mu .gt. 5.d-6) go to 10
      write (iunit,810)
      if (iflag .ne.   0) then
          write (    *,*) ' Error code = ', iflag
          write (iunit,*) ' Error code = ', iflag
      end if
      ig   = idamax (n, g, 1)
      write (iunit,820) f
      write (iunit,830) abs(g(ig))
      write (iunit,840)
      write (iunit,850) (x(i), i = 1,n)
      stop
810   format (//, ' Results of optimization problem', /)
820   format (' Final function value = ', 1pd24.16)
830   format (' Norm of gradient     = ', 1pd12.4)
840   format (' Parameters')
850   format (4(1pd16.8))
860   format (//,' mu =', 1pd12.4)
      end
c
c----------------------------------------------------------------------
c initial guess
c----------------------------------------------------------------------
c
      subroutine xstart (x, n)
      double precision x(n)
      double precision lb(100), ub(100), mu
      common /bnd/     lb, ub
      common /bar/     mu
c
c specify initial values of variables and their bounds
c
      do 10 i = 1,n
          x(i)  = 0.5d0
          lb(i) = 0.d0
          ub(i) = 1.d0
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
c evaluate nonlinear function and gradient
c
      implicit         double precision (a-h,o-z)
      double precision x(n), f, g(n)
      double precision lb(100), ub(100), mu
      common /bnd/     lb, ub
      common /bar/     mu
C$    PRIVATE
c
      eps = 1.d-20
      f = 0.d0
      do 30 i = 1,n
          b = 1.d0
          if (i .gt. n/2) b = -b
          d1 = x(i) - lb(i)
          d2 = ub(i) - x(i)
          if (d1 .lt. eps) write (*,800) d1
          if (d2 .lt. eps) write (*,800) d2
          f = f + b * x(i)
     *       - mu * (log (d1)) - mu * (log (d2))
          g(i) =  b - mu / d1 + mu / d2
30    continue
c
      return
800   format (' Warning: negative slack value ', d16.8, /,
     *        ' Parameter mu may be too small')
      end
c
c
      subroutine maxstp (n, x, d, stepmx)
c
c computes maximum stepsize moving from the (feasible) point x
c in direction d, under box constraints
c
c Parameters:
c n       --> integer, number of variables
c x       --> double precision, size n, current and new vector
c             of parameters
c d       --> double precision, size n, search direction vector
c stepmx <--  double precision, maximum step allowed
c
      implicit         double precision (a-h, o-z)
      double precision d(*), x(*), lb(100), ub(100)
      common /bnd/     lb, ub
c
c  set up
c
      alpha   = 1.d8
      ep = 1.d-8
      i1 = 1
      i2 = n
      do 10 i = i1,i2
          if (d(i) .gt. ep) then
              t = (ub(i) - x(i)) / d(i)
           else if (d(i) .lt. -ep) then
              t = (x(i) - lb(i)) /(-d(i))
           endif
           if (t .lt. alpha) alpha = t
10    continue
      stepmx = alpha*0.99999d0
c
      return
      end
