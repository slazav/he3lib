! Vapor pressure [bars] vs T [K]
! Sherman, R.H.; Sydoriak, S.G.; Roberts, T.R.
! The 1962 He3 scale of temperatures
! http://archive.org/details/jresv68An6p579

      function He3_Pvap(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.0.2D0.and.T.le.3.324D0) then
          He3_Pvap =
     .     - 2.49174D0 / T
     .     + 4.80386D0
     .     - 0.286001D0 * T
     .     + 0.198608D0 * T**2
     .     - 0.0502237D0 * T**3
     .     + 0.00505486D0 * T**4
          He3_Pvap = dexp(He3_Pvap) * 1.333224D-3 ! torr -> bar
        else
          He3_Pvap = -1D0
        endif
        return
      end


