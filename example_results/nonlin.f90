MODULE NONLIN


CONTAINS

SUBROUTINE W_SCATT(W_v,W_t,W_tx,kx,Omega,grad)

 ! Langmuir wave scattering off ions
 ! only induced terms


USE CONSTANT
USE PARAMS

   IMPLICIT NONE

   DOUBLE PRECISION,DIMENSION(:) :: W_v,W_t,W_tx,kx;
   DOUBLE PRECISION, INTENT(IN) :: Omega,grad;
   DOUBLE PRECISION :: Omega1,Omega0, Omega1_t,Omega0_t,alpha,k_m;
   DOUBLE PRECISION :: Omega1_tx,Omega0_tx;
   DOUBLE PRECISION :: ksi, ksi2,ksi3,w_ion,wlt_t,wlt_l;
   DOUBLE PRECISION :: ksi2x,ksi3x,wlt_tx,wlt_lx;

   INTEGER::ik,kk;
   DOUBLE PRECISION :: k0,k01,dk,dk2;

   DOUBLE PRECISION :: ne_local;

!call writeProfiles(W_v,W_v,W_k,99999);
!----------------------------------------------------
ne_local=Omega**2*m/(4*Pi*e*e)

Alpha = sqrt(2*Pi)*Omega**2/(4.*(2.*Pi)**3*ne_local*v_Ti*(1+T_e/T_i)**2)

DO kk=2, 2*Nv-1
   
   w_ion  =0. 
   wlt_t  =0.
   wlt_l  =0.
   wlt_tx  =0.
   wlt_lx  =0.

   k0  = kx(kk)*Omega/v_t

   dk2= abs(kx(kk)-kx(kk-1))*Omega/v_t

   Omega0  =Omega+3.*v_T**2*k0**2/(2.*Omega)
   Omega0_t =Omega+V_c*v_c*k0*k0/(2.*Omega)
   Omega0_tx=Omega+0.5*omega_ce+V_c*v_c*k0*k0/(2.*Omega)

   DO ik=2, 2*Nv

   IF (ik/=kk) THEN

   k01= kx(ik)*Omega/v_T
  
   k_m =abs(k0-k01)

   Omega1  =Omega+3.*v_T*v_T*k01*k01/(2.*Omega)
   Omega1_t =Omega+V_c*v_c*k01*k01/(2.*Omega)
   Omega1_tx=Omega+0.5*omega_ce+V_c*v_c*k01*k01/(2.*Omega)
    

   !ksi  = (3.*V_T*v_T*(k0+k01))/(2.*sqrt(2.)*v_Ti*Omega)
   ksi = (Omega0-Omega1)/(sqrt(2.)*k_m*V_Ti)
   ksi2 = (Omega0_t-Omega1)/(sqrt(2.)*k_m*V_Ti)
   ksi3 = (Omega0-Omega1_t)/(sqrt(2.)*k_m*V_Ti)

   ksi2x = (Omega0_tx-Omega1)/(sqrt(2.)*k_m*V_Ti)
   ksi3x = (Omega0-Omega1_tx)/(sqrt(2.)*k_m*V_Ti)
  
   dk=abs(kx(ik)-kx(ik-1))*Omega/v_T

! L-wave scattering
   w_ion=w_ion+dk*exp(-ksi**2)*(Omega0*W_v(ik)/Omega1-W_v(kk)-&
   (2.*Pi)**3*(Omega0-Omega1)*W_v(kk)*W_v(ik)*k2/(k_b*T_i*Omega1))/k_m;

! o-mode emission

   wlt_t=wlt_t+dk*exp(-ksi2**2)*(Omega0_t*W_v(ik)/Omega1-W_t(kk)-&
   (2.*Pi)**3*(Omega0_t-Omega1)*W_t(kk)*W_v(ik)*k2/(k_b*T_i*Omega1))/k_m;

   wlt_L=wlt_L+dk*exp(-ksi3**2)*(Omega0*W_t(ik)/Omega1_t-W_v(kk)-&
   (2.*Pi)**3*(Omega0-Omega1_t)*W_v(kk)*W_t(ik)*k2/(k_b*T_i*Omega1_t))/k_m;

! x-mode emission   
   wlt_tx=wlt_tx+dk*exp(-ksi2x**2)*(Omega0_tx*W_v(ik)/Omega1-W_tx(kk)-&
   (2.*Pi)**3*(Omega0_tx-Omega1)*W_tx(kk)*W_v(ik)*k2/(k_b*T_i*Omega1))/k_m;
   
   wlt_Lx=wlt_Lx+dk*exp(-ksi3x**2)*(Omega0*W_tx(ik)/Omega1_tx-W_v(kk)-&
   (2.*Pi)**3*(Omega0-Omega1_tx)*W_v(kk)*W_tx(ik)*k2/(k_b*T_i*Omega1_tx))/k_m;

   END IF    

   END DO

W_v(kk) = W_v(kk)+dt2*w_ion*Alpha+dt2*wlt_L*Alpha+dt2*wlt_Lx*Alpha

W_t(kk) = W_t(kk)+dt2*wlt_t*Alpha*POLAR(omega,omega0_t,+1.0d0)

W_tx(kk)= W_tx(kk)+dt2*wlt_tx*Alpha*POLAR(omega,omega0_tx,-1.0d0)

END DO

CONTAINS

DOUBLE PRECISION FUNCTION POLAR(omega_pe,omega_t,sigma)
USE CONSTANT
USE PARAMS
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: omega_pe,Omega_t,sigma
! sigma is the sign of polarisation
! SIGMA =+1 for o-mode, sigma=-1 for x mode

DOUBLE PRECISION:: xx,yy;

xx=omega_pe*omega_pe/(omega_t*omega_t)
yy=omega_ce/omega_t

Polar =DTAN(fi)*DTAN(fi)*(yy+sigma*(1-xx)*DCOS(fi))**2/&
((1-xx)*(2.*(1.+sigma*yy*DCOS(fi))**2-sigma*xx*yy*DCOS(fi)))

END FUNCTION POLAR


      
END SUBROUTINE W_SCATT

END MODULE NONLIN