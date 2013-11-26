c********************************************************************
c BTN: CHKDER, parallel derivative checker for nonlinear functions
c last changed: 0927/91
c for Sequent Balance (parallel) and scalar computers
c******************************************************************
c
      subroutine chkez (n, x, f, g, w, lw, sfun, iflag, errmax,
     *                  imax)
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      external         sfun
c
c Check derivatives of objective function via finite-differencing.
c Easy to use version
c
c PARAMETERS
c n       -> integer, dimension of problem
c x       -> double precision, size n, point where derivatives are checked
c f      <-  double precision, value of function at the point x (from sfun)
c g      <-  double precision, size n, gradient at the point x (from sfun)
c w       -  double precision, size n*(2k), work array
c lw      -> integer, declared length of array w [not used]
c sfun    -> subroutine sfun (n,x,f,g) to evaluate f(x), g(x);
c            see subroutine BTNEZ for details. Routine SFUN must
c            be declared EXTERNAL in the calling program.
c            To run in parallel on the Sequent, this routine (and any
c            routine that it calls) must have the command
c            C$    PRIVATE
c            as its first executable statement.
c iflag  <-  integer, error code:   0 => gradient values appear accurate
c                                 100 => gradient error > .01
c errmax <-  double precision, size of largest error in the gradient
c imax   <-  integer, index where errmax was obtained
c
c set default values of customizing parameters
c
      call chkpar (n, kmax, msglvl, iunit, istart, iend,
     *      istep)
c
c check selected gradient values
c
      call chkder (n, x, f, g, w, lw, sfun, iflag,
     *             kmax, errmax, imax,
     *             msglvl, iunit, istart, iend, istep)
      return
      end
c
c
      subroutine chkpar (n, kmax, msglvl, iunit, istart, iend,
     *      istep)
C$    integer*4 cpus_online
c
c  sets parameters for  derivative checker
c
c n         -> integer, number of variables
c kmax     <-  integer, block size = # of processors used
c              default: kmax = [# of processors]/2 on Sequent
c msglvl   <-  integer, amount of printing desired,
c              (<-1 =>  none
c              ( -1 =>  warning messages only
c                 0 =>  one line with maximum error
c                 1 =>  individual components are printed. Processor i
c                       will print to iunit + i
c              default: msglvl = 0
c iunit    <-  integer, see discussion for msglvl
c              default: iunit = 10
c istart   <-  integer, index of the first partial derivative
c              to be checked
c              default: istart = 1
c iend     <-  integer, index of the last partial derivative to be
c              checked
c              default: iend = istep*kmax
c istep    <-  integer, gradient is tested every istep components.
c              default: istep = n/kmax
c              Warning: if n is large or if gradient values are
c                       expensive, checking the gradient values can
c                       be very time consuming.  The default values
c                       result in each processor checking exactly
c                       one gradient component.
c
      nodemx  = 5
C$    nodemx  = cpus_online()/2
      kmax    = nodemx
      if (kmax .gt. n) kmax = n
      msglvl  = 0
      iunit   = 10
      istart  = 1
      istep   = n/kmax
      iend    = istart + istep*kmax - 1
      if (iend .gt. n) iend = n
c
      return
      end
c
c
      subroutine chkder (n, x, f, g, w, lw, sfun, iflag,
     *                   kmax, errmax, imax,
     *                   msglvl, iu, istart, iend, istep)
c
c Check derivatives of objective function via finite-differencing.
c
c PARAMETERS
c n       -> integer, dimension of problem
c x       -> double precision, size n, point where derivatives are checked
c f      <-  double precision, value of function at the point x (from sfun)
c g      <-  double precision, size n, gradient at the point x (from sfun)
c w       -  double precision, size n*(2kmax), work array
c lw      -> integer, length of work array w [not used]
c sfun    -> subroutine sfun (n,x,f,g) to evaluate f(x), g(x);
c            see subroutine BTNEZ for details. Routine SFUN must
c            be declared EXTERNAL in the calling program.
c            To run in parallel on the Sequent, this routine (and any
c            routine that it calls) must have the command
c            C$    PRIVATE
c            as its first executable statement.
c iflag  <-  integer, error code:   0 => gradient values appear accurate
c                                 100 => gradient error > .01
c kmax    -> integer, block size [normally equal to number of processors,
c            unless more than one processor is being used to compute each
c            function value; see subroutine BTN for further details]
c errmax <-  double precision, size of largest error in the gradient
c imax   <-  integer, index where errmax was obtained
c msglvl  -> integer, specify amount of printing desired:
c                  <-1 => none
c                   -1 => warning messages only
c                    0 => one line with maximum error
c                    1 => output to unit (iu)
c iu      -> integer, see discussion of parameter msglvl
c istart  -> integer, test gradient starting at component istart
c iend    -> integer, test gradient starting at component iend
c istep   -> integer, test gradient at every istep component
c                    e.g., if istep=2 then test every 2nd component
c
c REQUIRES
c d1mach  -  function to specify machine epsilon [d1mach(3)]
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(n,*)
C$    integer*4        m_get_myid
      logical          msg0, msg1
c
c Check value of machine epsilon
c
      eps = d1mach(3)
      if (1.d0 + 2.d0*eps .eq. 1.d0) then
          write (*,*) ' ### CHKDER: ERROR ###'
          write (*,*) '     Machine epsilon too small'
          write (*,*) '     Value = ', eps
          write (*,*) '     Modify routine D1MACH'
          write (*,*) '     Terminating execution'
          stop
      end if
c
c Setup: obtain processor information
c
      iflag  = 0
      h      = sqrt(eps)
      imax   = 0
      msg0   = msglvl .ge. 0
      msg1   = msglvl .ge. 1
      if (msg1) write (iu,800)
      errmax = 0.d0
c
c evaluate function and gradient at x
c
      call sfun (n, x, f, g)
c
c compute finite-difference estimate of gradient
c
C$DOACROSS SHARE(w,x,g,f,h,errmax,imax,istart,istep,
C$&              n,sfun,kmax),
C$&        LOCAL(fx,gi,erri,i1,i2)
      do 10 i = istart, iend, istep
          i1 = (i-istart)/istep + 1
C$        i1 = m_get_myid() + 1
          i2 = i1 + kmax
          call dcopy (n, x, 1, w(1,i2), 1)
          w(i,i2) = w(i,i2) + h
          call sfun (n, w(1,i2), fx, w(1,i1))
          gi = (fx - f) / h
          erri = abs(gi - g(i)) / (1.d0 + abs(g(i)))
c         if (msg1) write (iu,810) i, g(i), gi, erri
C$        call m_lock
          if (erri .gt. errmax) then
              errmax = erri
              imax   = i
          end if
C$        call m_unlock
10    continue
c
      if (errmax .gt. 1.d-2) iflag = 100
      if (msg0) write ( *,820) errmax, imax
      if (msg1) write (iu,820) errmax, imax
c
      return
800   format (' Testing Derivatives', //,
     *        '   i', 9x, 'g(i)', 14x, 'gi', 13x, 'error')
810   format (' ', i3, 2x, d16.8, 2x, d16.8, 2x, d12.4)
820   format (/, ' CHKDER:  Max. error in gradient = ', d12.4, /,
     *           '          observed at component = ', i6)
      end
