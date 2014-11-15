      program main
        implicit none
        include '../he3.fh'

        real*8 ttc, p, nuB

        data   ttc /0.2D0/,    ! Temperature, T/Tc
     .         p   /10D0/     ! Pressure, bar

        nuB = he3_nu_b(ttc,p)
        write(*,*) 'nu_b  = ', nuB/1000, 'kHz'
        write(*,*) 'm_3   = ', he3_amass, 'g'

      end
