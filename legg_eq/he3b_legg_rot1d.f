! Leggett equations in rotating frame. 1D calculations along rotation axis.
! Spin currents, spin diffusion, leggett-takagi, extra Mz relaxation.
! input:
!   U(7)   - Mx My Mz Nx Ny Nz Theta
!   Ux(7)  - dU/dx
!   Uxx(7) - d2U/dx2
! output
!   Fv(7)  - dU/dt
! use BLK_LEGG_EQ common block (see he3b_legg_rot1d.hf) to set parameters
! Origin: Dmitriev

      subroutine he3b_legg_rot1d(U,UX,UXX,FV)
        implicit none
        include '../he3.fh'
        include 'he3b_legg_rot1d.fh'

        real*8 U,UX,UXX,FV
        dimension U(7),UX(7),UXX(7),FV(7)
        real*8 Wx, Wz, Wzh, dW
        real*8 UN,UNx,UNy,UNz,UMxm,UMym,UMzm
        real*8 DD45,ST,CT,CTM,CT1,CTG,UT,AUT,AF,DAF,FTN,DFTN,B
        real*8 UJX,UJY,UJZ,DJX,DJY,DJZ ! spinn current J and dJ/dz

C       calculate freq
        Wx = he3_gyro*Hx
        Wz = he3_gyro*Hz
        Wzh = 0.5D0*Wz
        dW = Wz-W0

C       fix n vector length
        UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/UN
        UNy = U(5)/UN
        UNz = U(6)/UN

        UMxm = U(1)-Wx/Wz
        UMym = U(2)
        UMzm = U(3)-1.0D0

        DD45=UNx*UX(5)-UX(4)*UNy       ! Nx Ny` - Nx` Ny
        ST=dsin(U(7))
        CT=dcos(U(7))
        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM      ! ctg(T/2) = sin(T)/(1-cos(T))
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0

        AUT = UT*Flegg**2/Wz*4.0D0*PI**2
        AF  = -Cpar**2/Wz
        DAF = -2D0*Cpar*dCpar/Wz

        FTN=CTM*DD45-ST*UX(6)-UX(7)*UNz
        DFTN=CTM*(UNx*UXX(5)-UXX(4)*UNy)-ST*UXX(6)-UXX(7)*UNz-
     *   CT1*UX(7)*UX(6)+ST*UX(7)*DD45   !!! dFTN/dz

C       components of spin current, Ji
        UJX=2.0D0*(UX(7)*UNx+ST*UX(4)+CTM*(UNy*UX(6)-UX(5)*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*FTN
        UJY=2.0D0*(UX(7)*UNy+ST*UX(5)-CTM*(UNx*UX(6)-UX(4)*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*FTN
        UJZ=2.0D0*(UX(7)*UNz+ST*UX(6)+CTM*(UNx*UX(5)-UX(4)*UNy))+
     *   (CTM*UNz**2+CT)*FTN

C       dJi/dz
        DJX=AF*(2.0D0*(UXX(7)*UNx+CT1*UX(7)*UX(4)+ST*UXX(4)+
     *   ST*UX(7)*(UNy*UX(6)-UX(5)*UNz)+CTM*(UNy*UXX(6)-UXX(5)*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*DFTN+(ST*UX(7)*UNx*UNz+
     *   CTM*(UX(4)*UNz+UNx*UX(6))+UX(5)*ST+UNy*CT*UX(7))*FTN)-
     *   Diff*UXX(1)+DAF*UJX
        DJY=AF*(2.0D0*(UXX(7)*UNy+CT1*UX(7)*UX(5)+ST*UXX(5)-
     *   ST*UX(7)*(UNx*UX(6)-UX(4)*UNz)-CTM*(UNx*UXX(6)-UXX(4)*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*DFTN+(ST*UX(7)*UNy*UNz+
     *   CTM*(UX(5)*UNz+UNy*UX(6))-UX(4)*ST-UNx*CT*UX(7))*FTN)-   !!!!!!!!!
     *   Diff*UXX(2)+DAF*UJY
        DJZ=AF*(2.0D0*(UXX(7)*UNz+CT1*UX(7)*UX(6)+ST*UXX(6)+
     *   ST*UX(7)*DD45+CTM*(UNx*UXX(5)-UXX(4)*UNy))+
     *   (CTM*UNz**2+CT)*DFTN+(ST*UX(7)*UNz**2+
     *   CTM*2.0D0*UNz*UX(6)-ST*UX(7))*FTN)-Diff*UXX(3)+DAF*UJZ

        B = UNx*UMxm + UMym*UNy + UMzm*UNz

C       Leggett equations
        FV(1)=   dW*U(2)           + AUT*UNx - DJX
        FV(2)= - dW*U(1) + Wx*U(3) + AUT*UNy - DJY
        FV(3)=           - Wx*U(2) + AUT*UNz - DJZ - UMzm*T1
        FV(4)= - W0*UNy - Wzh*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)=   W0*UNx - Wzh*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=          - Wzh*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)= Wz*B + UT/Tf

      end
