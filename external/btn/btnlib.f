c********************************************************************
c BTN: parallel software for unconstrained optimization
c last changed: 10/02/91
c For: Sequent Balance (parallel) and scalar computers
c********************************************************************
c
      subroutine btnez (n, x, f, g, w, lw, sfun, iflag)
c
c parallel block truncated-Newton routine for unconstrained optimization
c
c Parameters
c n        -> integer, number of variables
c x       <-> double precision, size n, nonlinear variables (initial guess
c                 and solution)
c f       <-  double precision, nonlinear function value (output)
c g       <-  double precision, size n, gradient (output)
c w        -  double precision, size lw, work space; currently must be of
c                 length at least 3n + 3k + 4k^2 + 7(n*k);
c                 where k = # of processors being used
c lw       -> integer, length of w
c sfun     -  subroutine to evaluation nonlinear function:
c                 subroutine sfun (n, x, f, g)
c             Routine sfun must be declared EXTERNAL in the calling program.
c             The parameters n, x, f, and g are as above; x and n must
c             not be modified by routine sfun; sfun must return the
c             nonlinear function value in f and the gradient vector in g.
c             To run in parallel on the Sequent, this routine (and any
c             routine that it calls) must have the command
c             C$    PRIVATE
c             as its first executable statement.
c iflag   <-  integer, error code:
c                 0 => normal return: successful termination
c                 1 => linesearch failed to find lower point
c                 2 => search direction points uphill
c                 3 => function appears to be unbounded below
c                 4 => more than maxfun evaluations of f(x) in linesearch
c                 5 => not converged after nlmax iterations
c                 6 => n <= 0
c               900 => insufficient work space provided
c               999 => disaster (see file fort.iunit for message)
c             See the user manual for more information about the meaning
c             of these error codes.
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      external         sfun
c
c set parameters for BTN
c
      call btnpar (nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *    initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *    ireset, indstp)
c
c solve optimization problem
c
      call btn (n, x, f, g, w, lw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, sfun)
c
      return
      end
c
c
      subroutine btnpar (nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *    initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *    ireset, indstp)
      implicit double precision (a-h,o-z)
C$    integer  cpus_online
c
c Sets parameters for BTN
c
c Parameters:
c nodemx  <-  integer, maximum number of processors used
c             default:   nodemx = [# of processors]/2 on Sequent
c kmax    <-  integer, block size
c             default:   kmax = nodemx
c maxit   <-  integer, maximum number of inner iterations allowed.
c             default:   maxit = 20
c msglvl  <-  integer,  amount of printing desired,
c             (<-1: none
c               -1: warning messages only
c                0: one line per outer iteration)
c             default:   msglvl = 0
c iprec   <-  integer, type of preconditioner
c             (0 = none, 1 = BFGS, 2 = approx. BFGS)
c             default:      iprec = 1.
c nlmax   <-  integer, maximum number of outer iterations allowed
c             default:   nlmax = 500
c initv   <-  integer, type of initialization
c             (1 = random,2 = rhs + random, 3 = limited memory)
c             default:   initv = 3
c tolq    <-  double precision, tolerance for quadratic-based test
c             default:   tolq =  5.d-1
c iunit   <-  integer, output unit # for printing from this processor
c rnktol  <-  double precision, tolerance for rank test in inner iteration
c             default: rnktol = 1.d-9]
c maxfun  <-  integer, maximum number of function evaluations allowed in
c             the linesearch.
c             default:  maxfun = 2*nlmax
c accrcy  <-  double precision, user-chosen convergence tolerance (stop the
c             algorithm if the infinity-norm of the gradient is  <=
c             accrcy*(1.d0 + |f(x)|).
c             default:  accrcy = 0.d0
c stepmx  <-  double precision, maximum step allowed in line search
c             default:   stepmx = 1.d3
c eta     <-  double precision, accuracy parameter for line search (must be
c             strictly between 0 and 1, and it is recommended that
c             it be chosen between .1 and .9)
c             default:  eta = 0.2
c ireset  <-  integer, reset preconditioner after ireset iterations.
c             If ireset <= 0, then no resetting is done.
c             default:  ireset = 0 (no resetting)
c indstp  <-  integer, indicates if the user wishes to control the
c             maximum step size in the line search with a user provided
c             subroutine MAXSTP.  If indstp=0, BTN automatically determines
c             the maximum step.  If indstp<>0, then subroutine MAXSTP must
c             be provided.
c             default:  indstp = 0 (automatic setting of step size)
c
      nodemx  = 5
C$    nodemx  = cpus_online()/2
      kmax    = nodemx
      maxit   = 20
      msglvl  = 0
      iprec   = 1
      nlmax   = 500
      initv   = 3
      tolq    = 5.0d-1
      iunit   = 10
      rnktol  = 1.0d-9
      maxfun  = 2*nlmax
      accrcy  = 0.d0
      stepmx  = 1.d3
      eta     = 0.2d0
      ireset  = 0
      indstp  = 0
c
      return
      end
c
c
      subroutine btn (n, x, f, g, w, lw, sfun, iflag,
     *            nodemx, kmax, maxit, msglvl, iprec, nlmax,
     *            initv, tolq, iunit, rnktol, maxfun, accrcy,
     *            stepmx, eta, ireset, indstp, maxstp)
c
c parallel block truncated-Newton routine
c
c Parameters
c n        -> integer, number of variables
c x        -> double precision, size n, nonlinear variables (initial
c                  guess and solution)
c f       <-  double precision, nonlinear function value (output)
c g       <-  double precision, size n, gradient (output)
c w        -  double precision, size lw, work space; currently must be of
c                 length at least 3n + 3kmax + 4kmax^2 + 7(n*kmax);
c lw       -> integer, length of w
c sfun     -  subroutine to evaluation nonlinear function:
c                 subroutine sfun (n, x, f, g)
c             Routine sfun must be declared EXTERNAL in the calling program.
c             The parameters n, x, f, and g are as above; x and n must
c             not be modified by routine sfun; sfun must return the
c             nonlinear function value in f and the gradient vector in g.
c             To run in parallel on the Sequent, this routine (and any
c             routine that it calls) must have the command
c             C$    PRIVATE
c             as its first executable statement.
c iflag   <-  integer, error code:
c                 0 => ok
c                 1 => linesearch failed to find lower point
c                 2 => search direction points uphill
c                 3 => function appears to be unbounded below
c                 4 => more than maxfun evaluations of f(x)
c                 5 => not converged after nlmax iterations
c                 6 => n <= 0
c               900 => insufficient work space provided
c               999 => disaster (see file fort.iunit for message)
c            See the user manual for more information about the meaning
c            of these error codes.
c nodemx   -> integer, number of processors used
c kmax     -> integer, block size [usually equal to nodemx]
c maxit    -> integer, maximum number of inner iterations allowed
c                 [typical value: 20]
c msglvl   -> integer, amount of printing desired (<0 = none;
c                 0 = 1 line per outer iteration)
c iprec    -> integer, type of preconditioner (0 = none, 1 = BFGS,
c                 2 = approx. BFGS) [it is strongly suggested that
c                 iprec = 1 be used]
c nlmax    -> integer, maximum number of outer iterations allowed
c                 [typical value: 500]
c initv    -> integer, type of initialization (1 = random,
c                 (2 = rhs + random, 3 = limited memory) [it is
c                 suggested that initv = 3 be used]
c tolq     -> double precision, tolerance for quadratic-based test
c                 [typical value: 5.d-1]
c iunit    -> integer, output unit # for printing from this processor
c rnktol   -> double precision, tolerance for rank test in inner iteration
c                 [typical value: 1.d-9]
c maxfun   -> integer, maximum number of function evaluations allowed in
c                 the linesearch [typical value: 2*nlmax]
c accrcy   -> double precision, user-chosen convergence tolerance (stop
c                 the algorithm if the infinity-norm of the gradient is
c                 <= accrcy*(1.d0 + |f(x)|).  Under normal circumstances,
c                 the user should set accrcy=0.0 and let the stringent
c                 convergence test built-in to BTN be used.  If the function
c                 f(x) or its gradient cannot be obtained accurately, then
c                 the user may wish to override this test by setting
c                 accrcy = 1.d-6 (say), or perhaps some slightly larger value.
c                 [typical value: 0.d0]
c stepmx   -> double precision, maximum step allowed in line search
c                 [typical value: 1.d3]
c eta      -> double precision, accuracy parameter for line search (must be
c                 strictly between 0 and 1, and it is recommended that
c                 it be chosen between .1 and .9) [typical value: 0.2]
c ireset   -> integer, reset preconditioner after ireset iterations; if
c                 no resetting is desired, set ireset=0; normally it
c                 will not be necessary to reset the preconditioner based
c                 on the iteration count, but it can occasionally improve
c                 performance on certain specially structured problems; see
c                 the user manual for further guidance.
c                 [recommended value: 0]
c indstp   -> integer, indicates if the user wishes to select the maximum
c                 step in the line search manually at every iteration.
c                 See sample main program 4 for an example.
c                 If indstp=0, then BTN sets the maximum step automatically.
c                 If indstp<>0, then the user-supplied routine MAXSTP is
c                 called (see below).
c                 [typical value: 0]
c maxstp   -  (optional) subroutine to set the maximum step for the line
c             search.  It must have the following calling sequence:
c                 subroutine maxstp (n, x, d, stepmx)
c             where n and x are as for subroutine sfun, d is the current
c             search direction (the new vector of parameters will be of
c             the form x+a*d for some scalar a), and stepmx is the largest
c             allowable step (computed by the user within subroutine
c             maxstp).  The values of n, x, and d must not be modified by
c             the user, and stepmx must be defined by the user.  Subroutine
c             maxstp must be declared EXTERNAL in the calling program.
c             NOTE: most applications will not require this feature, and
c             should set indstp=0 when calling BTN (see above); when
c             indstp=0, BTN will not attempt to call maxstp, and no subroutine
c             need be provided.  In this case, the call to BTN can have the
c             form:  call btn (..., indstp, sfun), where sfun is the name
c             of the function evaluation routine described above.  See the
c             sample main programs for further examples.
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      integer          flgbcg, cubesz
C$    integer*4        m_get_numprocs
      logical          msg
      external         sfun
c
c check for inappropriate parameter settings
c
      if (maxit .lt. 0) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: maxit < 0'
              write (*,*) '       maxit = ', maxit
              write (*,*) '       Resetting maxit = 20'
          end if
          maxit = 20
      end if
      if (iprec .lt. 0 .or. iprec .gt. 2) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: iprec out of range'
              write (*,*) '       iprec = ', iprec
              write (*,*) '       Resetting iprec = 1'
          end if
          iprec = 1
      end if
      if (initv .lt. 1 .or. initv .gt. 3) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: initv out of range'
              write (*,*) '       initv = ', initv
              write (*,*) '       Resetting initv = 3'
          end if
          initv = 3
      end if
      if (tolq .le. 0.d0 .or. tolq .ge. 1.d0) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: tolq out of range'
              write (*,*) '       tolq = ', tolq
              write (*,*) '       Resetting tolq = 0.5'
          end if
          tolq = 5.d-1
      end if
      if (rnktol .lt. 0.d0) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: rnktol out of range'
              write (*,*) '       rnktol = ', rnktol
              write (*,*) '       Resetting rnktol = 1.d-9'
          end if
          rnktol = 1.d-9
      end if
      if (eta .le. 0.d0 .or. eta .ge. 1.d0) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: eta out of range'
              write (*,*) '       eta = ', eta
              write (*,*) '       Resetting eta = 0.2'
          end if
          eta = 2.d-1
      end if
c
c set up
c
      if (n .le. 0) then
          iflag = 6
          go to 500
      end if
      cubesz = kmax
C$    cubesz = m_get_numprocs()
c
c compare n with blocksize
c
      if (n .lt. kmax) then
          if (msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: blocksize too small'
              write (*,*) '       kmax = ', kmax
              write (*,*) '       n    = ', n
              write (*,*) '       Setting kmax = n'
          end if
          kmax = n
      end if
      maxit1   = maxit
      kxk      = kmax * kmax
      k        = kmax
      nk       = n * k
      nf       = 1
c
c subdivide the work array
c
      lp      = 1
      lpre    = lp     + n
      ivinit  = lpre   + n
      iv0     = ivinit + nk
      iv1     = iv0    + nk
      ir      = iv1    + nk
      iar     = ir     + nk
      ip      = iar    + nk
      ialpha  = ip     + nk
      ibeta   = ialpha + kxk
      il1     = ibeta  + kxk
      il2     = il1    + kxk
      id1     = il2    + kxk
      id2     = id1    + kmax
      is      = id2    + kmax
      iwk     = is     + kmax
      lprenw  = iwk    + nk
      lwtest  = lprenw + n - 1
      if (lw .lt. lwtest) then
          call disast ('BTN', 'INSUFFICIENT STORAGE', iunit, iflag)
          write (iunit,*) ' # of locations required = ', lwtest
          iflag = 900
          go to 500
      end if
c
      msg     = msglvl .ge.  0
      ncg     = 0
      iter    = 0
      iflag   = 0
      flgbcg  = 0
      intflg  = 0
      itsave  = 0
c
c initial printouts
c
      call sfun (n, x, f, g)
      fsave = abs(f)
      gnrm = pdnrmi (n, g)
      if (msg) then
          write (*,800)
          write (*,810) iter, nf, ncg, f, gnrm
      end if
c
c set tolerances for convergence tests
c
      call settol (tol0,tol1,tol2,tol3,tol4)
c
c convergence test at initial point
c
      ftest = 1.d0 + dabs(f)
      if (gnrm .lt. tol0*ftest) go to 500
c
c initialize preconditioner to the identity
c
      if (iprec .gt. 0) call pdfill (n, 1.d0, w(lprenw))
c
c  main loop
c
      do 10 iter = 1,nlmax
c
c  set the preconditioner
c
          itsave = itsave + 1
          if (iprec .eq. 0) then
              call pdfill (n, 1.d0, w(lpre))
          else
              call pdcopy (n, w(lprenw), w(lpre))
          end if
c
c  get search direction
c
          call precg (w(iv1), n, n, g, kmax, iter,
     *        initv, flgbcg)
          call bcg (n, w(lp), g, gnrm, k, kmax, flgbcg,
     *        maxit1, x, nbl, w(lpre), w(lprenw),
     *        tolq, n, kmax, w(iv0), w(iv1), w(ir),
     *        w(iar), w(ip), w(ialpha), w(ibeta), w(il1),
     *        w(il2), w(id1), w(id2), w(is), w(iwk),
     *        nodemx, iprec, sfun, iunit, rnktol)
          k = kmax
          if (iflag .eq. 999) go to 500
          call pdneg (n, w(lp))
          dg = pddot (n, w(lp), g)
          ncg = ncg + nbl
          if (initv .eq. 3 .and.
     *        flgbcg .gt. 0 .and. intflg .eq. flgbcg)
     *        initv = 2
          call postcg (w(iv1), n, n, k, iter, initv,
     *        w(ivinit), n, w(lp), g, gnrm, flgbcg, cubesz)
          if (flgbcg .ge. 0) intflg = flgbcg
          if (indstp .ne. 0) call maxstp (n, x, w(lp), stepmx)
c
c parallel line search
c
          oldf = f
          call psrch (iflag, n, x, f, g, w(lp), w(ip), w(ir),
     *               w(iar), sfun, nf, eta, stepmx, kmax,
     *               cubesz, dg, w(iv0), w(iwk), alpha)
          if (iflag .lt. 0) then
              maxit1 = 1 + 3/kmax
          else
              maxit1 = maxit
          end if
          if (iflag .lt. 0) iflag = 0
          gnrm = pdnrmi (n, g)
          if (msg) write (*,810) iter, nf, ncg, f, gnrm
          if (iflag .ne. 0) go to 500
c
c test for convergence
c
          diff  = abs(oldf - f)
          ftest = 1.d0 + dabs(f)
          xnrm = pdnrmi (n, x)
          xtest = 1.d0 + xnrm
          pnrm = pdnrmi (n, w(lp))
c
          if ((alpha*pnrm  .lt.   tol1*xtest
     *         .and. diff  .lt.   tol2*ftest
     *         .and. gnrm  .lt.   tol3*ftest)
     *         .or.  gnrm  .lt.   tol4*ftest
     *         .or.  gnrm  .lt. accrcy*ftest) go to 500
          if (nf .gt. maxfun) then
              iflag = 4
              go to 500
          end if
          ftest = abs(f)
          if (ftest .le. fsave*1.d-5 .or. itsave .eq. ireset) then
              if (iprec .gt. 0) call pdfill (n, 1.d0, w(lprenw))
              fsave = ftest
              itsave = 0
          end if
10    continue
      iflag = 5
c
500   if (msg) write (*,820) iflag
      return
800   format (/, 3x, 'it', 3x, 'nf', 2x, 'ncg', 11x, 'f', 14x, '|g|')
810   format (' ', i4, 1x, i4, 1x, i4, 2x, 1pd16.8, 2x, 1pd12.4)
820   format (' BTN terminating with error code = ', i4)
      end
c
c
      subroutine settol (tol0,tol1,tol2,tol3,tol4)
      implicit double precision (a-h,o-z)
c
c set tolerances for convergence tests
c
      eps    = d1mach(3)
      if (1.d0 + 2.d0*eps .eq. 1.d0) then
          write (*,*) ' ### BTN: ERROR ###'
          write (*,*) '     Machine epsilon too small'
          write (*,*) '     Value = ', eps
          write (*,*) '     Modify routine D1MACH'
          write (*,*) '     Terminating execution'
          stop
      end if
      rteps  = sqrt(eps)
      rtol   = 1.d1*rteps
      rtolsq = rtol*rtol
      peps   = eps**0.6666d0
c
      tol0 = 1.d-2  * rteps
      tol1 = rtol   + rteps
      tol2 = rtolsq + eps
      tol3 = sqrt(peps)
      tol4 = 1.d-2*rtol
c
      return
      end
c
c
      subroutine bcg (n, x, b, bnrm, k, kmax, iflag, maxit,
     *                xnl, ncg, PC, PCnew,
     *                tolq, nn, kk, V0, V1, R, AR, U,
     *                alpha, beta, L1, L2, D1, D2, s,
     *                wk, nodemx, iprec, sfun, iunit, rnktol)
c
c Parallel block conjugate-gradient iteration (to compute a
c search direction for a truncated-Newton optimization algorithm)
c
c PARAMETERS
c n          -> integer, dimension of problem
c x         <-  double precision, size n, search direction
c b          -> double precision, size n, right-hand side vector
c bnrm       -> double precision, infinity-norm of right-hand side
c k         <-> integer, current block size
c kmax       -> integer, maximum block size
c iflag     <-  integer, flag:
c                       0 => okay
c                      -1 => steepest-descent direction
c                       i => linear depend. in col. i at iteration 1
c                     999 => disaster (message to unit iunit)
c maxit      -> integer, maximum number of inner iterations allowed
c                     (if maxit<0, then nonlinearity test is also used)
c xnl        -> double precision, size n, vector of nonlinear parameters
c ncg       <-  integer, counter for number of inner iterations
c PC        <-> double precision, size n, current and new preconditioner
c PCnew      -  double precision, size n, work array (new PC)
c tolq       -> double precision, tolerance for quadratic test
c nn         -> integer, leading dimension of n*k matrices
c kk         -> integer, leading dimension of k*k matrices
c V0         -  double precision, size n*k, work array (old Lanczos
c               vectors)
c V1         -> double precision, size n*k, work array (current Lanczos
c               vectors)
c R          -  double precision, size n*k, work array (preconditioned V1)
c AR         -  double precision, size n*k, work array (Hessian times R)
c U          -  double precision, size n*k, work array (inner search
c               directions)
c alpha      -  double precision, size k*k, work array (diagonal block
c               of V'AV)
c beta       -  double precision, size k*k, work array (off-diag block
c               of V'AV)
c L1         -  double precision, size k*k, work array (Cholesky factor
c               of V'AV)
c L2         -  double precision, size k*k, work array (Cholesky factor
c               of V'AV)
c D1         -  double precision, size   k, work array (Cholesky factor
c               of V'AV)
c D2         -  double precision, size   k, work array (Cholesky factor
c               of V'AV)
c s          -  double precision, size   k, work array (inner step sizes)
c wk         -  double precision, size n*k, work array
c nodemx     -> integer, maximum number of processors available
c iprec      -> integer, choice of preconditioner:
c                       0 => no preconditioner
c                       1 => LMQN diagonal PC
c                       2 => approximate LMQN diagonal PC
c sfun       -> external, subroutine sfun (n,x,f,g), evaluate f(x)
c               To run in parallel on the Sequent, this routine (and any
c               routine that it calls) must have the command
c               C$    PRIVATE
c               as its first executable statement.
c iunit      -> integer, output unit # for printing
c rnktol     -> double precision, tolerance for rank test in inner iteration
c
      implicit         double precision (a-h,o-z)
      double precision b(n), xnl(n), x(n), PC(n), PCnew(n),
     *                 V0(nn,kk), V1(nn,kk), R(nn,kk), AR(nn,kk),
     *                 U(nn,kk), alpha(kk,kk), beta(kk,kk),
     *                 L1(kk,kk), L2(kk,kk), D1(kk), D2(kk), s(kk),
     *                 wk(nn,kk)
      integer          kmax, cubesz
      logical          indef
      external         sfun
c
c set up
c
      cubesz = kmax
      kold   = k
      iflag  = 0
      info   = 0
      iend   = 0
      indef  = .false.
      xAx    = 0.d0
      nblock = n/k
      if (nblock*k .lt. n) nblock = nblock + 1
c
c**********************************************************************
c Initialization
c**********************************************************************
c
      call pdfill (n, 0.d0, x)
      call dfill (k, 0.d0, s, 1)
      if (bnrm .eq. 0.d0) then
          call disast ('BCG', 'bnrm = 0', iunit, iflag)
          go to 500
      end if
      call pdcopy (n, b, V1)
      call pdscal (n, 1.d0/bnrm, V1)
c
c**********************************************************************
c Main loop
c**********************************************************************
c
      maxnbl = nblock
      if (maxit  .eq.     0) maxit  = 1
      if (maxnbl .gt. maxit) maxnbl = maxit
c
      do 100 nbl = 1, maxnbl
          ncg  = nbl
          kold = k
c
c**********************************************************************
c Get QR of V1
c**********************************************************************
c
          call msolve (V1, nn, n, k, R, nn, PC)
          call pgetch (V1, nn, R, nn, beta, kk, n, kold,
     *            k, iflag, rnktol, iunit)
          if (iflag .eq. 999) return
          if (nbl .eq. 1) then
              if (k .eq. 0) then
                  call disast ('BCG', 'K=0 AT ITER=1', iunit, iflag)
                  iflag = -1
                  call pdcopy  (n, b, x)
                  return
              else
                  if (k. lt. kold) iflag = k+1
                  kold = k
              endif
          else
              if (k .eq. 0) go to 500
          endif
c
c Compute V1 x (Beta')^(-1)
c
          call vltnv1 (V1, nn, n, k, beta, kk, iunit, iflag)
          if (iflag .eq. 999) return
          call msolve (V1, nn, n, k, R, nn, PC)
          if (nbl .eq. 1) then
              s(1) = pddot (n, b, R(1,1))
              s(1) = s(1) / bnrm
          end if
c
c**********************************************************************
c Compute AR, alpha
c**********************************************************************
c
          call getar (R, nn, k, AR, nn, cubesz, wk, nn, n, xnl,
     *                b, sfun)
          call AtrBA (R, nn, k, AR, nn, alpha, kk, n)
c
c update diagonal preconditioner (if desired)
c
          if (iprec .eq. 1) then
              call frmpc (PCnew, alpha, R, AR, nn, n, kk, k, wk)
          endif
          if (iprec .eq. 2) then
              call frmpc1 (PCnew, alpha, R, AR, nn, n, kk, k, wk)
          endif
c
c**********************************************************************
c Cholesky for L1
c set D1 = D2 (kold * kold), then solve for L1 (kold * k)
c L2 * D1 * L1 = beta, where L2, D1 are of dimension kold
c**********************************************************************
c
          if (nbl .gt. 1) then
              call dcopy (kold, D2, 1, D1, 1)
              do 29 i = 1,k
                  if (L2(i,i) .eq. 0.d0) then
                      call disast ('LSOL', 'l(i,i) = 0',
     *                             iunit, iflag)
                      return
                  end if
29            continue
C$DOACROSS SHARE(L2,L1,beta,D1,kk,kold,iunit,iflag)
              do 30 i = 1,k
                  call lsol (L2(i,i), kk, kold-i+1, beta(i,i),
     *                  L1(i,i), iunit, iflag)
                  call dvdiv (kold-i+1,L1(i,i),1,D1(i),1,L1(i,i),1)
30            continue
          end if
c
c**********************************************************************
c Cholesky for alpha
c Cholesky factorization of block tridiagonal matrix.
c In the following L2 will be the factor of the diagonal block.
c All processors compute the Cholesky simultaneously.
c**********************************************************************
c
          call matcpy (alpha, kk, L2, kk, k, k)
          if (nbl .gt. 1) then
C$DOACROSS SHARE(L2,L1,D1,kold), LOCAL(j,l)
              do 50 i = 1,k
              do 50 j = 1,i
                  do 40 l = i,kold
                      L2(i,j) = L2(i,j) - L1(l,i) * D1(l) * L1(l,j)
40                continue
                  if (j .lt. i) L2(j,i) = L2 (i,j)
50            continue
          end if
          call dpofa2 (L2, kk, k, info, iunit, iflag, rnktol)
          if (iflag .eq. 999) return
          if (info .ne. 0) then
              indef = .true.
              k = abs(info) - 1
              if (nbl .eq. 1 .and. k .eq. 0) then
                  call pdcopy (n, b, x)
                  iflag = -1
                  return
              endif
              if (nbl .gt. 1 .and. k .eq. 0) go to 500
          endif
c
c get L2 and D2
c
          call getld2 (L2, kk, k, D2, iunit, iflag)
          if (iflag .eq. 999) return
c
c**********************************************************************
c Search direction U
c**********************************************************************
c compute  V1 - U L1
c     V1 (n * k), U (n * kold), L1 (kold * k) lower triangular
c     V1 - U L1 will be stored in R
c     U*L1 is stored in U (=U)
c
          if (nbl .gt. 1) then
              call timsl1 (U, nn, n, kold, L1, kk, k)
              call aminsb (R, nn, U, nn, R, nn, n, k)
          end if
c
c solve for U:
c    U * L2' = V1
c    V1 is n * k
c    L2 is k * k lower triangular
c
          call matcpy (R, nn, U, nn, n, k)
          call vltnv2 (U, nn, n, k, L2, kk, iunit, iflag)
          if (iflag .eq. 999) return
c
c**********************************************************************
c Step size
c**********************************************************************
c
          if (nbl .eq. 1) then
              call lsol  (L2, kk, k, s, s, iunit, iflag)
              if (iflag .eq. 999) return
              call dvdiv (k, s, 1, D2, 1, s, 1)
          else
              call dvmul  (kold, D1, 1, s, 1, s, 1)
              call premlt (s, nn, kold, 1, L1, kk, k, wk, n)
              call lsol   (L2, kk, k, wk, s, iunit, iflag)
              if (iflag .eq. 999) return
              call dvdiv  (k, s, 1, D2, 1, s, 1)
              call dscal  (k, -1.d0, s, 1)
          end if
c
c**********************************************************************
c Update x
c**********************************************************************
c
          call maxpy (x, n, n, 1, 1.d0, U, nn, k, s, kk)
c
c**********************************************************************
c compute termination criterion using quadratic test
c**********************************************************************
c
          if (indef) go to 500
          if (nbl .eq. maxnbl) go to 500
          if (nbl .gt. 1 .and. k .eq. 0) go to 500
c
c quadratic test
c
          call dvmul  (k, D2, 1, s, 1, wk, 1)
          xAx  = xAx + ddot(k, s, 1, wk, 1)
          bx = pddot (n, b, x)
          qnew = xAx * 0.5d0 - bx
          if (nbl .gt. 1) then
              diff = qnew - qold
              if (qnew .eq. 0.d0) then
                  call disast ('BCG','Q(P) = 0',iunit, iflag)
                  go to 500
              end if
              if ((nbl * diff / qnew) .lt. tolq) iend = 1
          end if
          qold = qnew
c
          if (iend .eq. 1) go to 500
c
c**********************************************************************
c Update for next iteration:  form new V
c**********************************************************************
c AR:= AR - V0 * beta; V0 (n * kold), beta (kold * k)
c
          if (nbl .gt. 1) then
              call timsl1 (V0, nn, n, kold, beta, kk, k)
              call aminsb (AR, nn, V0, nn, AR, nn, n, k)
          end if
c
c V1:= AR - V1 * alpha; update V0
c
          call matcpy (V1, nn,  V0, nn, n, k)
          call atimsb (V0, nn, n, k, alpha, kk, k, V1, nn)
          call aminsb (AR, nn, V1, nn, V1, nn, n, k)
c
100   continue
c
500   call pdscal (n, bnrm, x)
      return
      end
c
c
      subroutine pgetch (V, ldv, Av, ldav, Beta, ldb, n, k,
     *             kr, iflag, rnktol, iunit)
c
c Cholesky factorization of V'A V, where
c
c PARAMETERS
c V         -> double precision, size n*k, input matrix
c ldv       -> integer, leading dimension of V
c Av        -> double precision, size n*k, input matrix A*V
c ldav      -> integer, leading dimension of AV
c Beta     <-  double precision, size k*k, Cholesky factor
c ldb       -> integer, leading dimension of Beta
c n         -> integer, size of arrays
c k         -> integer, block size (current)
c kr       <-  integer, rank of Beta (new block size)
c iflag    <-  integer, flag (=999 in case of disaster)
c rnktol    -> double precision, tolerance for rank test in inner
c              iteration
c iunit     -> integer, unit # for disaster messages
c
      implicit    double precision (a-h,o-z)
      dimension   V(ldv,*), Av(ldav,*), Beta(ldb,*)
c
c form V'AV
c
      call AtrBA (V, ldv, k, Av, ldav, Beta, ldb, n)
c
c find Cholesky factorization of V'AV (simultaneously on all processors)
c and determine new block size (i.e., rank of Beta)
c
      call dpofa2 (Beta, ldb, k, info, iunit, iflag, rnktol)
      if (iflag .eq. 999) return
      if (info .ne. 0) then
          kr = abs(info) - 1
      else if (info .eq. 0) then
          kr = k
      endif
c
      return
      end
c
c
      subroutine AtrBA (A, lda, k, BA, ldb, C, ldc, n)
c
c computes C = A'* BA, where BA is B*A with B symmetric.
c A (n x k), BA (n * k) and C (k * k)
c The result is explicitly symmetrized, and both halves are computed.
c
c PARAMETERS
c A         -> double precision, size n*k, input matrix
c lda       -> integer, leading dimension of A
c k         -> integer, block size
c BA        -> double precision, size n*k, input matrix B*A
c ldb       -> integer, leading dimension of BA
c C        <-  double precision, size k*k, result matrix
c ldc       -> integer, leading dimension of C
c n         -> integer, size of arrays
c
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), BA(ldb,*), C(ldc,*)
c
c form piece of A'BA
c
C$DOACROSS SHARE(C,A,BA,k,n), LOCAL(j)
      do 10 i = 1, k
          do 5 j = 1, k
              C(i,j) = ddot(n, A(1,i), 1, BA(1,j), 1)
5         continue
10    continue
c
c symmetrize result
c
C$DOACROSS SHARE(C), LOCAL(j)
      do 20 i = 1,k
          do 15 j = 1,i-1
              C(j,i) = (C(i,j)+C(j,i))*5.d-1
              C(i,j) = C(j,i)
15        continue
20    continue
c
      return
      end
c
c
      subroutine getar (R, ldr, k, AR, ldar, kmax, wk, ldwk,
     *                  n, xnl, g, sfun)
c
c getar: forms A*R, where R is a block matrix, by finite differencing
c
c PARAMETERS
c
c R             -> double precision, size n*k, original block matrix
c ldr           -> integer, leading dimension of R
c k             -> integer, block size
c AR           <-  double precision, size n*k, result: A*R (A = Hessian)
c ldar          -> integer, leading dimension of AR
c kmax          -> integer, # of active processors
c wk            -  double precision, size n*k, work array
c ldwk          -> integer, leading dimension of wk
c n             -> integer, # of variables
c xnl           -> double precision, current estimate of nonlinear
c                  variables
c g             -> double precision, current gradient
c sfun          -  subroutine to evaluation nonlinear function
c                  To run in parallel on the Sequent, this routine (and
c                  any routine that it calls) must have the command
c                  C$    PRIVATE
c                  as its first executable statement.
c
      implicit         double precision (a-h,o-z)
      double precision R(ldr,*), AR(ldar,*), wk(ldwk,*), g(*), xnl(*)
      integer          kmax
      external         sfun
c
C$DOACROSS SHARE(R,AR,xnl,g,wk,n,k,sfun)
      do 20 j = 1,kmax
          call atimes (R(1,j), AR(1,j), n, xnl, g, wk(1,j),
     *                 sfun, k)
20    continue
c
      return
      end
c
c
      subroutine atimes (v, av, n, x, g, wk, sfun, k)
      implicit         double precision  (a-h,o-z)
      double precision v(*), av(*), x(*), g(*), wk(*)
C$    private
c
c compute matrix/vector product via finite-differencing
c
c PARAMETERS
c v       -> double precision, size n, input for matrix/vector product
c av      <- double precision, size n, result of matrix/vector product
c n       -> integer, dimension of problem
c x       -> double precision, size n, current nonlinear parameter vector
c g       -> double precision, size n, current nonlinear gradient
c wk      -  double precision, size n, work array
c sfun    -> subroutine sfun (n,x,f,g) to evaluate nonlinear function
c            To run in parallel on the Sequent, this routine (and any
c            routine that it calls) must have the command
c            C$    PRIVATE
c            as its first executable statement.
c k       -> integer, the current block size
c
c REQUIRES
c d1mach -  function to specify machine epsilon [d1mach(3)]
c BLAS   -  basic linear algebra subroutines
c
      eps = d1mach(3)
      h   = 1.d1*sqrt(eps)
      rh  = 1.d0 / h
      call dcopy (n, x, 1, wk, 1)
      call daxpy (n, h, v, 1, wk, 1)
      call sfun  (n, wk, fw, av)
      call daxpy (n, -1.d0, g, 1, av, 1)
      call dscal (n, rh, av, 1)
c
      return
      end
c
c
      subroutine frmpc  (D, alpha, R, AR, nn, n, kk, k, wk)
c
c Update the diagonal preconditioner, based on BFGS formula
c
c PARAMETERS
c D        <-> double precision, size n, diagonal preconditioner
c alpha     -> double precision, size k*k, matrix from routine BCG
c R         -> double precision, size n*k, matrix from routine BCG
c AR        -> double precision, size n*k, matrix from routine BCG
c nn        -> integer, leading dimension of R and AR
c n         -> integer, size of arrays
c kk        -> integer, leading dimension of alpha
c k         -> integer, current block size
c wk        -  double precision, size 1, work variable
c
      implicit         double precision (a-h,o-z)
      double precision D(*), R(nn,*), AR(nn,*), alpha(kk,*),
     *                 wk(*)
c
c update D for each vector in the block
c
      do 30 ind = 1,k
c
c form v'Dv
c
          vdv = 0.d0
C$DOACROSS SHARE(R,D,ind), REDUCTION(vdv)
          do 10 i = 1,n
              vdv = vdv + R(i,ind)*D(i)*R(i,ind)
10        continue
c
c update D
c
          vgv = alpha(ind,ind)
          if (abs(vdv) .le. 1.d-12) go to 30
          if (abs(vgv) .le. 1.d-12) go to 30
C$DOACROSS SHARE(R,AR,D,vdv,vgv,ind), LOCAL(t1,t2)
          do 20 i = 1,n
              t1 = D(i)*R(i,ind)
              t2 = AR(i,ind)
              D(i) = D(i) - t1*t1/vdv + t2*t2/vgv
              if (D(i) .le. 1.d-6) D(i) = 1.d0
20        continue
c
30    continue
c
      return
      end
c
c
      subroutine frmpc1 (D, alpha, R, AR, nn, n, kk, k, wk)
c
c Update the diagonal preconditioner, BFGS (approximate) formula
c
c PARAMETERS
c D        <-> double precision, size n, diagonal preconditioner
c alpha     -> double precision, size k*k, matrix from routine BCG
c R         -> double precision, size n*k, matrix from routine BCG
c AR        -> double precision, size n*k, matrix from routine BCG
c nn        -> integer, leading dimension of R and AR
c n         -> integer, size of arrays
c kk        -> integer, leading dimension of alpha
c k         -> integer, current block size
c wk        -  double precision, size k, work array
c
      implicit         double precision (a-h,o-z)
      double precision D(*), R(nn,*), AR(nn,*), alpha(kk,*),
     *                 wk(*)
c
c form v'Dv for this processor, and then produce global result
c
C$DOACROSS SHARE(r,D,wk,n), LOCAL(vdv,i)
      do 20 ind = 1,k
          vdv = 0.d0
          do 10 i = 1,n
              vdv = vdv + r(i,ind)*D(i)*r(i,ind)
10        continue
          wk(ind) = vdv
20    continue
c
c update D
c
C$DOACROSS SHARE(R,AR,alpha,wk,D,k), LOCAL(vgv,vdv,t1,t2,ind)
      do 40 i = 1,n
          do 30  ind = 1,k
              vgv = alpha(ind,ind)
              vdv = wk(ind)
              if (abs(vdv) .le. 1.d-12) go to 30
              if (abs(vgv) .le. 1.d-12) go to 30
              t1 = D(i)*R(i,ind)
              t2 = AR(i,ind)
              D(i) = D(i) - t1*t1/vdv + t2*t2/vgv
              if (D(i) .le. 1.d-6) D(i) = 1.d0
30        continue
40    continue
c
      return
      end
c
c
      subroutine precg (V, ldv, n, b, kmax, iter, initv, iflag)
c
c initialize the matrix V
c      initv = 1  col 1 = b, all others random
c      initv = 2  col 1 = b, col 2 = previous direc, all others random
c      initv = 3  col 1 = b, alternate previous direc and previous grad
c
c PARAMETERS
c V      <-> double precision, size n*kmax, initial matrix for block
c            CG iteration
c ldv     -> integer, leading dimension of V
c n       -> integer, size of arrays
c b       -> double precision, size n, right-hand side vector
c kmax    -> integer, # of active processors
c iter    -> integer, outer iteration number
c initv   -> integer, type of initialization desired
c iflag   -> integer, indicator flag (tests if initv=3 is failing)
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*), b(*)
c
c Store b as first column of V
c
      call pdcopy (n, b, V)
c
c INITV=1,2,3:  Dynamic (and random) initialization
c
      krand = 2
      if (initv .eq. 2) then
          if (iter .gt. 1) krand = 3
      else if (initv .eq. 3) then
          if (iflag .eq. -1) then
              krand = 3
              go to 20
          endif
          krand = 2*iter
          if (krand .gt. kmax) go to 40
      endif
20    do 30 j = krand,kmax
          call rancol (V, ldv, j, n)
30    continue
c
40    return
      end
c
c
      subroutine postcg (V, ldv, n, k, iter, initv,
     *                  Vinit, ldvi, x, g, gnrm, iflag, kmax)
c
c initialize the matrix V
c    initv = 2  col 1 = b, col 2 = previous direc(x), all others random
c    initv = 3  col 1 = b, alternate previous direcs(x) and grads(g)
c If directions are linearly dependent (at the first inner iteration)
c the offending column is replaced by a random column.  If this happens
c for the same column two outer iterations in a row, subroutine btn
c switches to initialization strategy 3.
c
c V        <-  double precision, size n*k, new initialization matrix
c ldv       -> integer, leading dimension of V
c n         -> integer, size of problem
c k         -> integer, block size
c iter      -> integer, outer iteration number
c initv     -> integer, choice of initialization scheme
c Vinit    <-> double precision, size n*k, storage of initialization
c              vectors
c ldvi      -> integer, leading dimension of Vinit
c x         -> double precision, size n, search-direction vector
c g         -> double precision, size n, current gradient
c gnrm      -> double precision, infinity-norm of g
c iflag     -> integer, flag from BCG (indicates linear dependence)
c kmax      -> integer, number of active processors
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*), Vinit(ldvi,*), x(*), g(*)
c
      if (k .eq. 1) return
      if (initv .eq. 1) return
c
c store search direction in column 2
c
      xnrm = pdnrmi (n, x)
      call pdcopy (n, x, V(1,2))
      call pdscal (n, 1.d0/xnrm, V(1,2))
      if (initv .eq.  2) return
      if (iflag .eq. -1) return
      if (k     .eq.  2) return
c
c store current gradient in column 3
c
      call pdcopy (n, g, V(1,3))
      call pdscal (n, 1.d0/gnrm, V(1,3))
      if (k .eq. 3) return
      if (iflag .gt. 0)
     *     call rancol (Vinit, ldvi, iflag, n)
c
c update remaining columns
c
      do 10 j = k, 4, -1
          call pdcopy (n, Vinit(1,j-2), Vinit(1,j))
          call pdcopy (n, Vinit(1,j  ), V(1,j))
10    continue
c
      call pdcopy (n, V(1,2), Vinit(1,2))
      call pdcopy (n, V(1,3), Vinit(1,3))
c
      return
      end
c
c
      subroutine psrch (iflag, n, x, f, g, d, x1, g1, gopt,
     *            sfun, nf, eta, stepmx, kmax, cubesz,
     *            dg, alf, f1, alpha)
c
c Parallel line search
c      Using Armijo convergence test based on function values only
c      See Luenberger (2nd edition), p. 212
c
c Parameters:
c iflag  <--  integer, error code:
c                 0 => okay
c                -1 => okay, but alpha <> 1
c                 1 => no acceptable point found
c                 2 => d is a direction of ascent
c                 3 => function may be unbounded below
c n       --> integer, number of variables
c x      <--> double precision, size n, current and new vector of
c             parameters
c f      <--> double precision, current and new function value
c g      <--> double precision, size n, current and new gradient vector
c d       --> double precision, size n, search direction vector
c x1      --  double precision, size n*k, work array to store temporary x
c g1      --  double precision, size n*k, work array to store temporary g
c gopt    --  double precision, size   n, work array to store temporary g
c             (optimal)
c sfun    --> subroutine to evaluate f(x) (call sfun (n,x,f,g))
c             To run in parallel on the Sequent, this routine (and any
c             routine that it calls) must have the command
c             C$    PRIVATE
c             as its first executable statement.
c nf     <--> integer, total number of function/gradient evaluations
c eta     --> double precision, parameter for line search
c stepmx  --> double precision, maximum step allowed
c kmax    --> integer, number of processors available to evaluate f(x)
c cubesz  --> integer, number of processors
c alf     --  double precision, length cubesz, array for step lengths
c f1      --  double precision, length cubesz, array for function values
c alpha  <--  double precision, step length
c
      implicit         double precision (a-h, o-z)
      double precision d(*), g(*), g1(n,*), gopt(*),
     *                 x(*), x1(n,*), alf(*), f1(*)
      integer          kmax, cubesz
      logical          aup, adown
      external         sfun
c
c  set up
c
      aup     = .true.
      adown   = .true.
      alfopt  = 0.d0
      fopt    = f
      itmax   = 2 + 30/kmax
c
      iflag   = 0
      if (dg .ge. 0.d0) then
          iflag = 2
          return
      endif
      icount = 0
c
c  get maximum step and initialize alpha
c
      if (kmax .eq. 1) then
          alf(1) = 1.d0
          if (stepmx .lt. alf(1)) alf(1) = stepmx
          amax = alf(1)
          amin = amax
      else
          amax = kmax
          if (amax .gt. stepmx) amax = stepmx
          if (amax .gt. 1.d0) then
              amin = 1.d0 / dfloat(kmax)
              if (kmax .eq. 2) amin = 1.d0
          else
              amin = amax / dfloat(kmax)
          end if
          call setalf (amax, amin, alf, kmax, 1)
      endif
c
c  test function values at initial points
c
20    nf = nf + 1
      fmin = f
      imin = 0
      icount = icount + 1
C$DOACROSS SHARE(n,alf,d,x,x1,f1,g1,f,eta,dg,fmin,alpha,
C$&        ftest,imin,sfun),
C$&        LOCAL(f2)
      do 700 id = 1,kmax
          call dsvtvp (n, alf(id), d, 1, x, 1, x1(1,id), 1)
          call feval (n, x1(1,id), f1(id), g1(1,id), sfun)
          ftest = f + eta * alf(id) * dg
c
c determine minimum function value using the subroutine gopf
c
          f2   = f1(id)
          if (f2 .gt. ftest) f2 = fmin + 1.d0
C$        call m_lock
          if (f2 .lt. fmin) then
              fmin  = f2
              alpha = alf(id)
              imin  = id
          endif
C$        call m_unlock
700    continue
c
c Test for failure at this iteration and success at last iteration.
c In this case, use the previous optimal step and exit.
c
      if (fmin .ge. fopt .and. alfopt .ne. 0.d0) then
          imin  = kmax
          alpha = alfopt
          fmin  = fopt
          call pdcopy (n, gopt, g1(1,imin))
          go to 30
      end if
c
c Test for too many iterations
c
      if (icount .ge. itmax .and. imin .eq. 0) then
          iflag = 1
          return
      endif
      if (icount .ge. itmax .and. imin .eq. kmax) then
          iflag = 3
      endif
c
c Test to see if step must be reduced
c
      if (imin .eq. 0 .and. adown) then
          aup  = .false.
          amax = 0.9d0 * amin
          amin = amax / dfloat (kmax)
          if (kmax .eq. 1) then
              amax = amin / 2.d0
              amin = amax
          end if
          call setalf (amax, amin, alf, kmax, 2)
          go to 20
      endif
c
c Test to see if step must be increased
c
      if (kmax .gt. 1 .and. imin .eq. kmax
     *          .and. amax .lt. stepmx .and. aup) then
          adown = .false.
          if (fmin .lt. fopt) then
              alfopt = alpha
              fopt   = fmin
              call pdcopy (n, g1(1,imin), gopt)
          end if
          amin = amax * 1.1d0
          amax = amin * kmax
          if (amax .gt. stepmx) amax = stepmx
          call setalf (amax, amin, alf, kmax, 1)
          go to 20
      endif
c
c form the new x, and retrieve g
c
30    call pdsvtv (n, alpha, d, x, x)
c     call pdsvtvp (n, alpha, d, x, x) [NOTE CHANGE FROM INTEL NAME]
      if (alpha .ne. 1.d0) iflag = -1
      f = fmin
      call pdcopy (n, g1(1,imin), g)
c
      return
      end
c
c
      subroutine setalf (amax, amin, alf, cubesz, ind)
c
c determine a set of steplengths alpha for the line search
c
c Parameters
c amax    --> maximum value of alpha
c amin   <--  minimum value of alpha
c alf     --  double precision, length cubesz, array for step lengths
c cubesz  --> number of processors
c ind     --> indicator: =1 (normal) =2 (shrinking step)
c
      implicit         double precision (a-h, o-z)
      double precision alf(*)
      integer          cubesz
c
c Set up
c
      np1 = cubesz - 1
      if (amin .ge. amax) amin = amax / dfloat(cubesz)
c
c loop
c
      do 200 id = 1,cubesz
         id1 = id - 1
         if (ind .eq. 2) go to 100
c
         if (cubesz .gt. 2) then
             if (amin .lt. 1.d0 .and. amax .gt. 1.d0) then
                 imid = cubesz/2
                 if (id1 .lt. imid) then
                     alf(id) = amin + id1 * (1.d0-amin)/dfloat(imid)
                 else if (id1 .gt. imid) then
                     alf(id) = 1.d0 +
     *                  (id1-imid) * (amax-1.d0)/dfloat(np1-imid)
                 else if (id1 .eq. imid) then
                     alf(id) = 1.d0
                 end if
             else
                 alf(id) = amin + id1 * (amax-amin)/dfloat(np1)
             end if
         else
             if (id .eq. 1) alf(id) = amin
             if (id .eq. 2) alf(id) = amax
         end if
         go to 200
c
c Shrinking of step (more rapid reduction of step to zero)
c
100      if (cubesz .ge. 2) then
             alf(id) = amax / 2.d0**(cubesz-id1-1)
             amin  = amax / 2.d0**(cubesz-1)
         else
             alf(id) = amax
             amin  = amax
         end if
200   continue
c
      return
      end
c
c
      subroutine aminsb (A, lda, B, ldb, C, ldc, m, n)
      double precision A(lda,*), B(ldb,*), C(ldc,*)
c
c computes C = A - B, for m*n matrices A, B, and C
c [C can overwrite either A or B]
c
C$DOACROSS SHARE(A,B,C,m)
      do 10 j = 1, n
          call dvsub (m, A(1,j), 1, B(1,j), 1, C(1,j), 1)
10    continue
c
      return
      end
c
c
      subroutine atimsb (A, lda, n, m, B, ldb, k, C, ldc)
      implicit         double precision (a-h, o-z)
      double precision A(lda,*), B(ldb,*), C(ldc,*)
c
c forms C = A x B,  for A(n*m), B(m*k), and C(n*k)
c
C$DOACROSS SHARE(A,B,C,m,n,lda), LOCAL(i)
      do 10 j = 1,k
      do 10 i = 1,n
          C(i,j) = ddot (m,A(i,1),lda,B(1,j),1)
10    continue
c
      return
      end
c
c
      subroutine disast (routne, msg, iunit, iflag)
c
c print fatal error message on fort.iunit, then return
c
c PARAMETERS
c routne -> character*(*), name of routine where error was detected
c msg       -> character*(*), error message
c iunit       -> integer, unit for output
c iflag <-  integer, flag=999 upon return from disaster
c
      character*(*) routne, msg
c
      write (iunit,800) routne, msg
      iflag = 999
c
      return
800   format (' ********************',         /,
     *        ' ERROR, ERROR, ERROR',          /,
     *        ' Fatal error in routine ', a10, /,
     *        ' ', a40,                        /,
     *        ' Terminating',                  /,
     *        ' ********************')
      end
c
c
      subroutine dpofa2 (A, lda, n, info, iunit, iflag, rnktol)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*),  t, s, tol, mnorm
      integer          info, n
c
c Cholesky factor of a double precision positive definite matrix.
c this is a modification of the LINPACK routine DPOFA, designed to
c produce the lower triangular matrix a column at a time.
c
c info:   integer = 0 for normal return
c                 = k signals an error condition - the leading minor
c                     of order k is not positive definite.
c
      tol = rnktol * mnorm (A, lda, n, n)
c
c  get LDL decomposition from a Cholesky factorization
c
      do 30 i = 1,n
          s = ddot (i-1, A(i,1), lda, A(i,1), lda)
          s = A(i,i) - s
c
          info = i
          if (s .lt. -tol) return
          if (s .le.  tol) then
              info = - info
              return
          endif
          A(i,i) = sqrt(s)
c
          do 20 k = i+1,n
              t = a(k,i)  - ddot (i-1, A(i,1), lda, A(k,1), lda)
              if (A(i,i) .eq. 0.d0) then
                  call disast ('DPOFA2', 'A(i,i) = 0', iunit, iflag)
                  return
              end if
              t = t / A(i,i)
              A(k,i)  = t
20        continue
30    continue
      info = 0
c
      return
      end
c
c
      subroutine getld2 (L, ldl, k, D, iunit, iflag)
      implicit         double precision (a-h, o-z)
      double precision L(ldl,*), D(*)
c
c  get LDL decomposition from a Cholesky factorization
c
      do 10 i = 1,k
          D(i) = L(i,i)
          if (D(i) .eq. 0.d0) then
              call disast ('GETLD2', 'D(i) = 0', iunit, iflag)
              return
          end if
          do 10 j = 1,k
              if (j .gt. i) L(i,j) = 0.d0
              if (j .le. i) L(i,j) = L(i,j) / D(j)
10    continue
      call dvmul (k, D, 1, D, 1, D, 1)
c
      return
      end
c
c
      subroutine lsol (l, ldl, k, y, x, iunit, iflag)
      double precision l(ldl,*), y(*), x(*), ddot
C$    private
c
c solve Lx = y where L is lower triangular
c the vectors x and y may be the same (overwrite solution on rhs)
c
      do 10 i = 1,k
          x(i) = y(i) - ddot (i-1, l(i,1), ldl, x, 1)
          if (l(i,i) .eq. 0.d0) then
              call disast ('LSOL', 'l(i,i) = 0', iunit, iflag)
              return
          end if
          x(i) = x(i) / l(i,i)
10    continue
c
      return
      end
c
c
      subroutine matcpy (A, lda, B, ldb, m, n)
      double precision A(lda,*), B(ldb,*)
c
c copies matrix A  (m x n) onto matrix B
c
C$DOACROSS SHARE(A,B,m)
      do 10 j = 1,n
          call dcopy (m, A(1,j), 1, B(1,j), 1)
10    continue
c
      return
      end
c
c
      subroutine maxpy (y, ldy, n, k, alpha, a, lda, m, x, ldx)
      implicit         double precision (a-h,o-z)
      double precision y(ldy,*), a(lda,*), x(ldx,*)
c
c form y = y + alpha * Ax
c      y      n x k
c      A      n x m
c      x      m x k
c      alpha  scalar
c
C$DOACROSS SHARE(A,x,y,alpha,lda,m,n), LOCAL(j)
      do 20 i = 1,k
          do 10 j = 1,n
              y(j,i) = y(j,i) +
     *                 alpha * ddot(m,A(j,1),lda,x(1,i),1)
10        continue
20    continue
c
      return
      end
c
c
      double precision function mnorm (A, lda, m, n)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), colnrm, zero
      data             zero /0.d0/
c
c compute the 1-norm of the m x n matrix A
c
      mnorm = zero
      do 10 j = 1,n
          colnrm = dasum (m, A(1,j), 1)
          if (mnorm .lt. colnrm) mnorm = colnrm
10    continue
c
      return
      end
c
c
      subroutine msolve (A, lda, n, k, Am, ldam, Rm)
      double precision A(lda,*), Am(ldam,*), Rm(*)
c
c  preconditions the n*k matrix A, using the vector Rm.
c  resulting matrix is in Am [A and Am can be the same]
c
C$DOACROSS SHARE(Am,A,Rm,n)
      do 10 j = 1,k
          call dvdiv (n, A(1,j), 1, Rm, 1, Am(1,j), 1)
10    continue
c
      return
      end
c
c
      subroutine premlt (A, lda, k1, n,  L, ldl, k2, B,ldb)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), L(ldl,*), B(ldb,*)
c
c forms B =  L' x A
c for A(k1*n) and L(k1*k2) lower triangular (with k1 .ge. k2)
c
C$DOACROSS SHARE(A,B,L,k1,k2), LOCAL(i)
      do 10 j = 1,n
      do 10 i = 1,k2
          B(i,j) = ddot (k1-i+1,A(i,j),1,L(i,i),1)
10    continue
c
      return
      end
c
c
      subroutine rancol (V, ldv, j, m)
c
c  generate a random vector for the j-th column of the matrix V
c
c PARAMETERS
c V      <-  double precision, size m*j, initialization matrix
c            (see PREBCG)
c ldv     -> integer, leading dimension of V
c j       -> integer, column to be randomized
c m       -> integer, length of row
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*)
      integer          click, seed, rand
c
c setup
c
      click = 52343
      icol  = j
      seed  = mod(2*(icol)*click+1,65536)
      rand  = seed
c
c generate random column
c
      do 20 i = 1, m
          rand = mod(3125*rand,65536)
          V(i,j) = (rand - 32768.0d0)/16384.0d0
20    continue
c
      return
      end
c
c
      subroutine timsl1 (A, lda, n, k1, L, ldl, k2)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), L(ldl,*)
c
c performs A = A x L
c for A(n*k1) and L(k1*k2) lower triangular (with k1 .ge. k2)
c
      do 40 i = 1,k2
          call dscal (n, L(i,i), A(1,i), 1)
          do 30 k = i+1,k1
              call daxpy (n, L(k,i), A(1,k), 1, A(1,i), 1)
30        continue
40    continue
c
      return
      end
c
c
      subroutine vltnv1 (V, ldv, m, k, L, ldl, iunit, iflag)
      implicit         double precision (a-h,o-z)
      double precision L(ldl,*), V(ldv,*)
c
c Solves P L' = V, for L(k*k) lower triangular, P(m*k) and V(m*k)
c V is overwritten with P
c
      do 40 i = 1,k
          do 20 j = 1,i-1
              call daxpy (m, -L(i,j), V(1,j), 1, V(1,i), 1)
20        continue
          rl = 1.d0 / L(i,i)
          call dscal (m, rl, V(1,i), 1)
40    continue
c
      return
      end
c
c
      subroutine vltnv2 (V, ldv, m, k, L, ldl, iunit, iflag)
      implicit         double precision (a-h,o-z)
      double precision L(ldl,*), V(ldv,*)
c
c Solves P L' = V, for L(k*k) unit lower triangular, P(m*k) and V(m*k)
c V is overwritten with P
c
      do 40 i = 1,k
          do 20 j = 1,i-1
              call daxpy (m, -L(i,j), V(1,j), 1, V(1,i), 1)
20        continue
40    continue
c
      return
      end
c
c
      subroutine feval (n, x, f, g, sfun)
c
c compute function evaluation within line search
c
c PARAMETERS
c n       -> integer, dimension of problem
c x       -> double precision, size n, nonlinear parameter vector
c f      <-  double precision, function value at x
c g      <-  double precision, size n, gradient at x
c sfun    -> subroutine sfun (n,x,f,g) to evaluate nonlinear function
c            To run in parallel on the Sequent, this routine (and any
c            routine that it calls) must have the command
c            C$    PRIVATE
c            as its first executable statement.
c
      implicit         double precision  (a-h,o-z)
      double precision x(*), g(*), f
C$    private
c
c perform function evaluation
c
      call sfun  (n, x, f, g)
c
      return
      end
