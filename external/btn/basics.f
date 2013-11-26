c********************************************************************
c BTN: BASICS, low level routines (including Linpack BLAS)
c last changed: 11/14/90
c for Sequent Balance (parallel) and scalar computers
c******************************************************************
c
      subroutine dsvtvp (n, a, y, ly, x, lx, x1, lx1)
      double precision y(*), x(*), x1(*), a
C$    private
c
c Form x1 = x + a*y (with strides ly, lx, and lx1, as in Linpack)
c
      iy = 1
      ix = 1
      ix1 = 1
      do 10 i = 1,n
          x1(ix1) = x(ix) + a*y(iy)
          iy = iy + ly
          ix = ix + lx
          ix1 = ix1 + lx1
10    continue
c
      return
      end
c
c
      subroutine dfill (n, a, x, lx)
      double precision x(*), a
c
c Set x(i) = a (with stride lx, as in Linpack)
c
      ix = 1
      do 10 i = 1,n
          x(ix) = a
          ix = ix + lx
10    continue
c
      return
      end
c
c
      subroutine dneg (n, x, lx)
      double precision x(*)
c
c Set x(i) = -x(i) (with stride lx, as in Linpack)
c
      ix = 1
      do 10 i = 1,n
          x(ix) = -x(ix)
          ix = ix + lx
10    continue
c
      return
      end
c
c
      subroutine dvsub (n, x, lx, y, ly, z, lz)
      double precision x(*), y(*), z(*)
C$    private
c
c Form z(i) = x(i) - y(i) (with strides, as in Linpack)
c
      ix = 1
      iy = 1
      iz = 1
      do 10 i = 1,n
          z(iz) = x(ix) - y(iy)
          ix = ix + lx
          iy = iy + ly
          iz = iz + lz
10    continue
c
      return
      end
c
c
      subroutine dvdiv (n, x, lx, y, ly, z, lz)
      double precision x(*), y(*), z(*)
C$    private
c
c Form z(i) = x(i) / y(i) (with strides, as in Linpack)
c
      ix = 1
      iy = 1
      iz = 1
      do 10 i = 1,n
          z(iz) = x(ix) / y(iy)
          ix = ix + lx
          iy = iy + ly
          iz = iz + lz
10    continue
c
      return
      end
c
c
      subroutine dvmul (n, x, lx, y, ly, z, lz)
      double precision x(*), y(*), z(*)
c
c Form z(i) = x(i) * y(i) (with strides, as in Linpack)
c
      ix = 1
      iy = 1
      iz = 1
      do 10 i = 1,n
          z(iz) = x(ix) * y(iy)
          ix = ix + lx
          iy = iy + ly
          iz = iz + lz
10    continue
c
      return
      end
c
c
      subroutine pdsvtv (n, a, y, x, x1)
c     subroutine pdsvtvp (n, a, y, x, x1) [NOTE CHANGE FROM INTEL NAME]
      double precision y(*), x(*), x1(*), a
c
c Form x1 = x + a*y
c
C$DOACROSS SHARE(x1, x, y, a)
      do 10 i = 1,n
          x1(i) = x(i) + a*y(i)
10    continue
c
      return
      end
c
c
      subroutine pdfill (n, a, x)
      double precision x(*), a
c
c Set x(i) = a (with stride = 1)
c
C$DOACROSS SHARE(x, a)
      do 10 i = 1,n
          x(i) = a
10    continue
c
      return
      end
c
c
      subroutine pdneg (n, x)
      double precision x(*)
c
c Set x(i) = -x(i) (with stride = 1)
c
C$DOACROSS SHARE(x)
      do 10 i = 1,n
          x(i) = -x(i)
10    continue
c
      return
      end
c
c
      double precision function pddot (n, x, y)
c
c Parallel inner product routine
c Forms the inner product of its part of x and y in parallel
c
c PARAMETERS
c n         -> integer, size of arrays
c x         -> double precision, size n, input array 1
c y         -> double precision, size n, input array 2
c
      implicit         double precision (a-h,o-z)
      double precision x(*), y(*)
c
      t = 0.d0
C$DOACROSS SHARE(x,y), REDUCTION(t)
      do 10 i = 1,n
          t = t + x(i)*y(i)
10    continue
c
      pddot = t
c
      return
      end
c
c
      double precision function pdnrmi (n, x)
c
c Parallel infinity-norm routine
c
c PARAMETERS
c n         -> integer, size of arrays on each processor
c x         -> double precision, size n, input array
c
      implicit         double precision (a-h,o-z)
      double precision x(*), t
c
      t = 0.d0
C$DOACROSS SHARE(x), LOCAL(absxi), REDUCTION(t)
      do 10 i = 1,n
          absxi = abs(x(i))
          t = max(t,absxi)
10    continue
c
      pdnrmi = t
c
      return
      end
c
c
      subroutine pdcopy (n, x, y)
c
c Parallel vector copy routine
c Copies: y(i) = x(i)
c
c PARAMETERS
c n         -> integer, size of arrays
c x         -> double precision, size n, input array 1
c y         -> double precision, size n, input array 2
c
      implicit         double precision (a-h,o-z)
      double precision x(*), y(*)
c
C$DOACROSS SHARE(x,y)
      do 10 i = 1,n
          y(i) = x(i)
10    continue
c
      return
      end
c
c
      subroutine pdscal (n, a, x)
c
c Parallel vector scaling routine
c Scales: x(i) = a*x(i)
c
c PARAMETERS
c n         -> integer, size of arrays
c a         -> double precision, scaling factor
c x         -> double precision, size n, input/output array
c
      implicit         double precision (a-h,o-z)
      double precision x(*), a
c
C$DOACROSS SHARE(x,a)
      do 10 i = 1,n
          x(i) = a*x(i)
10    continue
c
      return
      end
c
c
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
C$    private
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
c
c***********************************************************************
c
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C$    private
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
c
c***********************************************************************
c
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
C$    private
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
c
c***********************************************************************
c
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C$    private
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
c
c***********************************************************************
c
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dmax
      integer i,incx,ix,n
C$    private
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
c
c***********************************************************************
c
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
C$    private
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
c
c***********************************************************************
c
      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(*), cutlo, cuthi, hitest, sum,
     *                   xmax, zero, one
      data   zero, one /0.0d0, 1.0d0/
C$    private
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
c
c***********************************************************************
c
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
C$    private
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
