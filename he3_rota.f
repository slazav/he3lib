!!! ROTA specific functions

! Nuclear stage heat capacity [J/K] vs T[K] and I[A]
      function rota_c_ns(t, i)
        implicit none
        include 'he3.fh'
        real*8 t, i, k, kmag
        k = 9.66D-5     ! J (K/T)^2 - in heat capacity
        kmag = 0.113D0  ! T/A - magnet constant
        rota_c_ns = k * (kmag * i/t)**2
      end
