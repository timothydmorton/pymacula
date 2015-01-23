! =====================================================================================
! MACULA version 1.0
!
PROGRAM maculacall

! COMPILING: g95 -O3 -c macula.f90 && g95 -O3 -o macula maculacall.f90 macula.o
! EXECUTING: ./macula

use maculamod

implicit none

 ! === INPUTS ===
 INTEGER, PARAMETER :: ndata = 1E3	  ! # of data points
 INTEGER, PARAMETER :: Nspot = 5	  ! # of star spots
 INTEGER, PARAMETER :: mmax = 1		  ! # of quarters/instruments
 LOGICAL, PARAMETER :: derivatives=.TRUE. ! Calculate partial derivatives?
 LOGICAL, PARAMETER :: temporal=.TRUE.    ! Calculate temporal derivatives?
 LOGICAL, PARAMETER :: TdeltaV=.TRUE.     ! Calculate transit depth variations?
 INTEGER :: i, k, l
 REAL(8), DIMENSION(ndata) :: t
 REAL(8), DIMENSION(12) :: Theta_star
 REAL(8), DIMENSION(8,Nspot) :: Theta_spot
 REAL(8), DIMENSION(2,mmax) :: Theta_inst
 REAL(8), DIMENSION(mmax) :: Tstart, Tend
 REAL(8), DIMENSION(ndata) :: Fmod, deltaratio
 REAL(8), DIMENSION(ndata,12) :: dFmod_star
 REAL(8), DIMENSION(ndata,8,Nspot) :: dFmod_spot
 REAL(8), DIMENSION(ndata,2,mmax) :: dFmod_inst
 REAL(8), DIMENSION(ndata) :: dFmoddt
 REAL(8), PARAMETER :: degs = 0.017453292519943295D0
 REAL(8), PARAMETER :: Pi = 3.141592653589793D0

 ! Create some dummy times
 DO i=1,ndata
  t(i) = 55000.0D0 + 0.1D0*i
 END DO

 ! Define start and end time of this dummy data set
 Tstart(1) = t(1) - 0.05D0
 Tend(1) = t(ndata) + 0.05D0

 ! Star's parameters
 call random_flat(12,Theta_star)
 Theta_star(1)  = Theta_star(1)*Pi-0.5D0*Pi	! inclination
 Theta_star(2)  = 0.23333D0*(Tend(1)-Tstart(1))*Theta_star(2)+0.1D0*(Tend(1)-Tstart(1))	! P_EQ
 Theta_star(3)  = 0.40D0*Theta_star(3) - 0.20D0	! kappa2
 Theta_star(4)  = 0.40D0*Theta_star(4) - 0.20D0	! kappa4
 Theta_star(5)  = 0.3999D0		!c1
 Theta_star(6)  = 0.4269D0		!c2
 Theta_star(7)  =-0.0227D0		!c3
 Theta_star(8)  =-0.0839D0		!c4
 Theta_star(9)  = 0.3999D0		!d1
 Theta_star(10) = 0.4269D0		!d2
 Theta_star(11) =-0.0227D0		!d3
 Theta_star(12) =-0.0839D0		!d4

 ! Spot parameters
 DO k=1,Nspot
   call random_flat(8,Theta_spot(:,k))
   Theta_spot(1,k) = Theta_spot(1,k)*2.0D0*Pi - Pi			! longitude
   Theta_spot(2,k) = Theta_spot(2,k)*Pi - 0.5D0*Pi			! latitude
   Theta_spot(3,k) = Theta_spot(3,k)*10.0D0*degs			! alpha_max
   Theta_spot(4,k) = Theta_spot(4,k)*0.5D0				! fspot
   Theta_spot(5,k) = Theta_spot(5,k)*(Tend(mmax)-Tstart(1))+Tstart(1)	! tmax
   Theta_spot(6,k) = (Tend(mmax)-Tstart(1))*Theta_spot(6,k)		! lifetime
   Theta_spot(7,k) = (Tend(mmax)-Tstart(1))*Theta_spot(7,k)		! ingress
   Theta_spot(8,k) = (Tend(mmax)-Tstart(1))*Theta_spot(8,k)		! egress
 END DO

 ! Instrumental parameters
 Theta_inst(1,1)  = 1.00D0	! U(1)
 Theta_inst(2,1)  = 1.00D0	! B(1)

 call macula(t,ndata,Nspot,mmax,derivatives,temporal,TdeltaV,&	! Controls
             Theta_star,Theta_spot,Theta_inst,&		! Star, Spot & Instrument parameters
             Tstart,Tend,&				! Times of start/end for each data set
             Fmod,dFmod_star,dFmod_spot,dFmod_inst,dFmoddt,deltaratio) ! Outputs

 ! Output
 open(11,FILE='output.dat',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,ndata
   write(11,*) t(i),Fmod(i),deltaratio(i),dFmoddt(i)!,dFmod_star,dFmod_spot,dFmod_inst
 END DO
 close(11)

CONTAINS

! =======================================================
 SUBROUTINE random_flat(n,r)
 INTEGER i, n
 DOUBLE PRECISION seeda, r(n)

 call random_seed()
 DO i=1,n
   call random_number(seeda)
   r(i)=seeda
 END DO

 END SUBROUTINE random_flat
! =======================================================

END PROGRAM maculacall
