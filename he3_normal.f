! Normal 3He liquid parameters
! beyond zero-temperature limit

! TODO -- Vm, Cv, Cp, ... from ??

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     heat capacity, Cv/R [dimensionless]
!     arguments: T [K], Vm [cm^3/mol]
!     (Cv = Cp up to terms (T/Tf)^3)
!     Greywall-1983
      function He3_cv_n(t, v)
        implicit none
        include 'he3.fh'
        real*8 t,v,a,b,c,d
        real*8 s1,s2,s3
        integer i,j
        dimension a(5,4), b(4,3), c(3,3), d(3)
        a(1,1) = -2.9190414D0
        a(1,2) =  5.2893401D+2
        a(1,3) = -1.8869641D+4
        a(1,4) =  2.6031315D+5
        a(2,1) =  0D0
        a(2,2) =  0D0
        a(2,3) =  0D0
        a(2,4) =  0D0
        a(3,1) = -2.4752597D+3
        a(3,2) =  1.8377260D+5
        a(3,3) = -3.4946553D+6
        a(3,4) =  0D0
        a(4,1) =  3.8887481D+4
        a(4,2) = -2.8649769D+6
        a(4,3) =  5.2526785D+7
        a(4,4) =  0D0
        a(5,1) = -1.7505655D+5
        a(5,2) =  1.2809001D+7
        a(5,3) = -2.3037701D+8
        a(5,4) =  0D0

        b(1,1) = -6.5521193D-2
        b(1,2) =  1.3502371D-2
        b(1,3) =  0D0
        b(2,1) =  4.1359033D-2
        b(2,2) =  3.8233755D-4
        b(2,3) = -5.3468396D-5
        b(3,1) =  5.7976786D-3
        b(3,2) = -6.5611532D-4
        b(3,3) =  1.2689707D-5
        b(4,1) = -3.8374623D-4
        b(4,2) =  3.2072581D-5
        b(4,3) = -5.3038906D-7

        c(1,1) = -2.5482958D+1
        c(1,2) =  1.6416936D+0
        c(1,3) = -1.5110378D-2
        c(2,1) =  3.7882751D+1
        c(2,2) = -2.8769188D+0
        c(2,3) =  3.5751181D-2
        c(3,1) =  2.4412956D+1
        c(3,2) = -2.4244083D+0
        c(3,3) =  6.7775905D-2

        d(1) = -7.1613436D+0
        d(2) =  6.0525139D-1
        d(3) = -7.1295855D-3

        if (t.lt.0.1D0) then
          s1=0D0
          do i=1,5
            do j=0,3
              s1 = s1 + a(i,j+1) * v**(-j) * t**i
            enddo
          enddo
          he3_cv_n = s1
          return
        endif

        if (t.ge.0.1D0.and.t.lt.2.5D0) then
          s1=0D0
          do i=0,3
            do j=0,2
              s1 = s1 + b(i+1,j+1) * v**j * t**(-i)
            enddo
          enddo
          s2=0D0
          do i=1,3
            do j=0,2
              s2 = s2 + c(i,j+1) * v**j * t**(-i)
            enddo
          enddo
          s3=0D0
          do j=0,2
            s3 = s3 + d(j+1) * v**j
          enddo
          he3_cv_n = s1 + dexp(-s3/t) * s2
          return
        endif
        he3_cv_n = NaN
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
