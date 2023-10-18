MODULE PLASMA

!module that contains subroutines that calculates
! density of the solar corona plasma
! all calculations are based on Parkes's model 1958

IMPLICIT NONE

CONTAINS

SUBROUTINE CORONA(dens,sw,x)

USE CONSTANT

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER:: V_1AU = 6.25e+7
! wind velocity at the distance 1AU =625km/s
DOUBLE PRECISION, PARAMETER:: T_e_model=1.0e+6
! the electron temperature of the corona
! only for model !!!!!!!!
DOUBLE PRECISION, PARAMETER:: N_1AU = 6.59D0
! electron number density at 1AU =6.59cm^(-3)
! the constant defined experimentally
! Mann, G et al A&A, 348, 614-620 (1999)

DOUBLE PRECISION, INTENT(in)::x
DOUBLE PRECISION, INTENT(OUT)::dens, sw
DOUBLE PRECISION:: vj, vk, v0, v_old
!service variables
DOUBLE PRECISION:: C
!constant that defined by experimental data
DOUBLE PRECISION :: x_c, v_T0;
! critical distance and velocity

!C = AU*AU*N_1AU*V_1AU
! constant from the paper 
c = 6.3e+34

v_T0 = sqrt(T_e_model*k_b/(mu*m_p))
! critical velocity

x_c = G*M_s/(2.0*v_T0*v_T0) 
! critical distance

IF (x > x_c) then
Vj = 5e+9
Vk = v_T0
ELSE
vk = 1e-12
vj = v_T0
END IF
! velocity of the solar wind at the distance of 1AU 

  v0=1e+10

  DO WHILE (abs(PARKER(x,v0))> 1e-10) 
    v_old=v0
    v0 = vj - PARKER(x,vj)*(vj-vk)/(PARKER(x,vj)-PARKER(x,vk))
    IF (Parker(x,v0)*Parker(x,vk)< 0.) THEN
      vj=v0
    ELSE
      vk=v0
    END IF
  END DO

sw =v0
!final solar wind velocity

Dens= C/(x*x*v0) 
! final density

CONTAINS

DOUBLE PRECISION FUNCTION PARKER(R,v_wind)

! Function of a density and solar wind velocity based
! on Parker's model of corona
 USE CONSTANT
 
IMPLICIT NONE 
 
  DOUBLE PRECISION, INTENT(IN)::R,v_wind;
  
  Parker = (v_wind/v_T0)**2 - log(v_wind*v_wind/(v_T0*v_T0))-4.0*log(R/x_c)- 4.0*x_c/R + 3.0

  END FUNCTION PARKER
  
END SUBROUTINE CORONA

!********************************************************

INTEGER FUNCTION X_NUMBER(x_min,x_max)

USE CONSTANT
USE PARAMS

IMPLICIT NONE 

DOUBLE PRECISION, INTENT(IN)::x_min,X_max;

INTEGER :: xstep;
DOUBLE PRECISION:: x_local,d_x,dxmin;
DOUBLE PRECISION::omega_local,density_local,wind_local;

  dxmin =d/15.d0  
!the minimum step in x
   
  x_local = xmin; 
  xstep =1;

  DO WHILE (x_local < X_max)
    
  CALL CORONA(density_local, wind_local,x_local+R_s+x_0)
  Omega_local = sqrt (4.0*Pi*e*e*Density_local/m)
  d_x = d*40.*1e+6/(Omega_local/(2*pi))

  IF (d_x .LE. dxmin) d_x = dxmin 

  xstep =xstep + 1
  x_local = x_local+d_x

 END DO    

 X_Number=xstep
 
!the number of x points
WRITE(output,*) 'Mesh dimensions calculated            -  OK'
END FUNCTION 

!***************************************

SUBROUTINE LOOP_DEN(x_min,x_max,Density_x,Omega_x,R_x)

  USE CONSTANT
  USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) ::x_min,x_max;
DOUBLE PRECISION, DIMENSION(:),INTENT(OUT)::Density_x,Omega_x,R_x;
DOUBLE PRECISION:: x_local,dx,lgrad,lcons;

dx =d/10.d0  
x_local = x_min;

lgrad=(dlog(1.d8)-dlog(1d13))/(dlog(5d9)-dlog(2d8))!-2.86135!
lcons=dlog(1d8)-lgrad*dlog(5d9)
! powerlaw of increasing density from x0 to xmax
! constant density from xmin to x0
!lgrad is (n(x0)-n(xend))/(x0-xend)

 DO i=1, Nx
	Density_x(i)=dexp(lcons+lgrad*dlog(x_0-x_local))
	Omega_x(i) = sqrt (4.0*Pi*e*e*Density_x(i)/m)
	R_x(i) = x_local
	x_local = x_local+dx
  END DO
  
  write(output,*) "Density calculated                    -  OK"

END SUBROUTINE LOOP_DEN

SUBROUTINE NEWLOOP_DEN(x_min,x_max,Density_x,Omega_x,R_x)

  USE CONSTANT
  USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) ::x_min,x_max;
DOUBLE PRECISION, DIMENSION(:),INTENT(OUT)::Density_x,Omega_x,R_x;
DOUBLE PRECISION:: x_local,dx;

dx =d/10.d0  
x_local = x_min;

! combination of aschwanden et al. 2002 for a RHESSI loop and Kontar et al. 2008 for Val C

 DO i=1, Nx
	!Density_x(i)=1.16d17*exp((x_local-x_0)/1.43d7)+1d13*((x_0-x_local)/1.0d8)**(-2.5)
	
	Density_x(i)=1.16d17*exp((x_local-x_0)/5.43d7)+2d9*(1.+0.1*sin(x_local/2e6))
	
	Omega_x(i) = sqrt (4.0*Pi*e*e*Density_x(i)/m)
	R_x(i) = x_local
	x_local = x_local+dx
  END DO
  
  write(output,*) "NEW density calculated                    -  OK"

END SUBROUTINE NEWLOOP_DEN

SUBROUTINE CORONA_PARAMETRS(x_min,x_max,Density_x,Omega_x,Wind_x,R_x)

  USE CONSTANT
  USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) ::x_min,x_max;
DOUBLE PRECISION, DIMENSION(:),INTENT(OUT)::Density_x,Omega_x,Wind_x,R_x;
DOUBLE PRECISION:: x_local,d_x,dxmin;

  dxmin =d/15.d0  
  !the minimum step in x
   
  x_local = x_min;

  DO i=1, Nx
 
  CALL CORONA(density_x(i), wind_x(i),x_local+R_s+x_0)

  Omega_x(i) = sqrt (4.0*Pi*e*e*Density_x(i)/m)

  d_x = d*40.*1e+6/(Omega_x(i)/(2*pi))
 !d_x decreases with frequency
  d_x =d/15.d0
 !fixed d_x
  R_x(i) = x_local

  IF (d_x .LE. dxmin) d_x = dxmin

  x_local = x_local+d_x

  END DO
  ! smooth profile is ready
  
!  density_x=density_x*(1.+(R_s/(R_x+R_s))**2*.1*sin(R_x/(Pi*d)))  
!  Omega_x = sqrt (4.0*Pi*e*e*Density_x/m)
  
  !variations are added 
  
  write(output,*) "Density calculated                    -  OK"

END SUBROUTINE CORONA_PARAMETRS

!-------------------------------------------------------

SUBROUTINE KNUMBER(vx,kx)
! subroutine to construct k-mesh

  USE CONSTANT
  USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:),INTENT(OUT)::vx ,kx;
DOUBLE PRECISION ::dk_add;

v_T = sqrt(k_b*T_e/m);
!electron thermal velocity

v_Ti= sqrt(k_b*T_i/m_p);
!ion thermal velocity

Omega_ce=e*B/(m*v_c);
! electron cyclotron frequency

vmin= V_T;
! minimal velocity

dv =(vmax1 - vmin)/(nv1-1)
dk_add=2*v_t/(vmax1*(2*nv-2*nv1+1))
!dv2=(vmax  - vmax1)/(nv-nv1-1)
! velocity step in the resonance region

DO I=1, Nv1
vx(i)= vmin+dv*(i-1)
kx(i)= v_t/vx(i)
vx(2*nv-i+1)= -vx(i)
kx(2*nv-i+1)= -kx(i)
END DO

DO I=Nv1+1,2*Nv-Nv1
kx(i)=kx(i-1)-dk_add
vx(i)=v_t/kx(i)
END DO


OPEN(UNIT=20, FILE=kv_mesh, STATUS='REPLACE',ACTION ='WRITE')
DO I=1, 2*Nv
WRITE(20,'(F8.5,A,ES14.4)')kx(i),' ',vx(i)/vbeam
END DO
CLOSE(20)
 
END SUBROUTINE KNUMBER

END MODULE PLASMA
