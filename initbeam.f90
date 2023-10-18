MODULE INITBEAM
 

   USE CONSTANT

   USE PARAMS
 
   IMPLICIT NONE

   CONTAINS

 DOUBLE PRECISION FUNCTION F0(v_local, x_local)


   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: v_local, x_local;  
   DOUBLE PRECISION :: a, del_pow;
   	    
	del_pow=6.0
	a=Nbeam*(1.0-del_pow)/(vbeam**(1.0-del_pow)-V_T**(1.0-del_pow))

if ((v_local <= 5.*V_T).AND.(v_local >=  v_T)) &
F0=a*exp(-x_local*x_local/(d*d))*(5.*v_T)**(-1.0*del_pow)

if ((v_local <= vbeam).AND.(v_local >= 5.*v_T)) &

F0=a*exp(-x_local*x_local/(d*d))*(v_local)**(-1.0*del_pow)
if (v_local >  vbeam) F0 = 1e-4*a*exp(-x_local*x_local/(d*d))*(v_local)**(-1.0*del_pow)
		  
END FUNCTION F0

!**********************************************************

SUBROUTINE initial_values(N_b,e_b,R_x,Density_x,Omega_x,Tqv_x,a2_x,gamma_ei,gamma_nov,b2_x)

USE CONSTANT
USE PARAMS

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:),INTENT(OUT)::Tqv_x,a2_x,gamma_ei,gamma_nov,b2_x;
DOUBLE PRECISION, DIMENSION(:),INTENT(IN)::R_x,Density_x,Omega_x;
DOUBLE PRECISION, INTENT(OUT)::E_b,N_b;
DOUBLE PRECISION :: nvlots,dvlots;
!beam density and energy density

  E_b=0.
  N_b=0.
  x=0.0D0
  
 nvlots=1e+6
dvlots=(vmax1 - vmin)/(nvlots-1.0)
  
! DO  i=2, Nv1;
!   v = vmin+(i)*dv  
!   N_b= N_b+(F0(v,x)+F0(v-dv,x))*dv/2.d0
!   E_b = E_b + (F0(v,x)*v*v+F0(v-dv,x)*(v-dv)*(v-dv))*dv/4.d0
DO  i=2, nvlots;
  v = vmin+(i)*dvlots  
  N_b= N_b+(F0(v,x)+F0(v-dvlots,x))*dvlots/2.d0
  E_b = E_b + (F0(v,x)*v*v+F0(v-dvlots,x)*(v-dvlots)*(v-dvlots))*dvlots/4.d0

  END DO
  ! calculation of electron beam density and 
  ! energy energy of electron beam
 ! write(*,*) n_b, e_b

! Calculation of the Landau Damping cutoff threshold
two_rootpi = 2*sqrt(pi)
landau_cutoff = 4*v_t
landau_min = exp(-landau_cutoff**2/(V_T*V_T))/(abs(landau_cutoff/V_T)**3)


!---------------- calculation of quasilinear time ---------

Tqv_x = DENSITY_x/(OMEGA_x*nbeam)

!---------------------------------------------------------  

VALUE1=MINLOC(ABS(R_x))
OMEGA_0 = omega_x(VALUE1(1))
!Minimal value of distance - close to zero
ne=m*OMEGA_0*OMEGA_0/(4.*pi*e*e)

! --- calculation of dimensionless coeficients starts here ----

   k1 = nbeam/vbeam;
   ! electron distribution dimensionl parameter
   k2 = m*nbeam*Vbeam**3/Omega_0;
   ! spectral energy density dimensional parameter
    
              
   b1 = 4.*PI*PI*e*e/m;
   ! First equation coeficient for resonance

    b2_x = b1*Omega_0/(omega_x*Vbeam**4);
    ! Second equation coeficient for resonance

   a1 = Vbeam **3*b1*nbeam/Omega_0;
   ! First equation coeficient
 
   !a2=(4*pi*pi*e*e*k2/m)/m;     
   a2_x = Pi/(Tqv_x*Vbeam);
    ! Second equation coeficient


! Coulomb Loagrithm term and Colomb Damping terms \gamma
! assuimng always T_e > 0.42MK
!lnA=15.89+dlog10(T_e)-0.5*dlog(Density_x)
gamma_ei=4.*PI*e**4*Density_x*(15.89+dlog(T_e)-0.5*dlog(Density_x))/(m*m*V_T*V_T*V_T)
gamma_nov=4.*PI*e**4*Density_x*(15.89+dlog(T_e)-0.5*dlog(Density_x))/(m*m)
 
write(output,*)'Dimension parameters calculated       -  OK'

END SUBROUTINE initial_values

!**********************************************************
SUBROUTINE initial_FW(v_x,F_vx,W_vx,Wt_vx,Wtx_vx,R_x,Omega_x)
USE CONSTANT
USE PARAMS

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:),INTENT(OUT)::F_vx,W_vx,Wt_vx,Wtx_vx;
DOUBLE PRECISION, DIMENSION(:),INTENT(IN)::R_x,v_x,Omega_x;


DO j=1, Nx
  DO i=1, Nv1
    v = vmin+(i-1)*dv  
    F_vx(i,j) = F0(v,R_x(j))/k1
  ! electron distribution function
  END DO
!boundary conditions
    F_vx(1,j) =0.d0
    F_vx(Nv,j)=0.d0
END DO
!-------------------------------
DO j=1, Nx
  DO i=1, 2*Nv
   !v = v_x(i)
   
   W_vx(i,j)  =k_b*T_e*Omega_x(j)**2/(2*pi*pi*v_t**2)/k2;
   Wt_vx(i,j) =k_b*T_e*Omega_x(j)**2/(2*pi*pi*v_t**2)/k2;
   Wtx_vx(i,j)=k_b*T_e*Omega_x(j)**2/(2*pi*pi*v_t**2)/k2;
	! what it was previously: 1e-3*k_b*T_e/((2*pi)**3*k2);
	!k2 is scaling/normalisation
  ! Spectral energy density of the thermal level
  END DO
!boundary conditions
W_vx(1,j) =0.
W_vx(Nv,j) =0.
Wt_vx(1,j) =0.
Wt_vx(Nv,j) =0.
END DO

write(output,*)'Initial values given                  -  OK'
END SUBROUTINE initial_FW

!**********************************************************
END MODULE INITBEAM
