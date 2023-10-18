MODULE SOLVER
! subroutines solving diff equations

CONTAINS


SUBROUTINE LandauDamping(W_v,vx,Omega_x)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  LandauDamping
!
!PURPOSE:
!
! This subroutine performs the Landau Damping term on
! the spectral energy density array
! CALLS: Constant.f90 and Params.f90
!
! INPUTS:  
!	W(v) 		; spectral energy density this x, t as a function of velocity
!	velocity	; array of velocities
!	Omega_x		; omega_pe for this x
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: 
!	W(v) 		; spectral energy density this x, t as a function of velocity
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	


USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: w_v,vx
DOUBLE PRECISION, INTENT(IN) :: Omega_x
DOUBLE PRECISION Landau

DO I=2*Nv-1,2,-1
	IF (abs(vx(i))<(landau_cutoff)) THEN
! 	Landau =2.*sqrt(PI)*Omega_x*exp(-4.**2)/(abs(4.)**3)
		Landau = two_rootpi*Omega_x*landau_min
	ELSE
		Landau = two_rootpi*Omega_x*exp(-vx(i)**2/(V_T*V_T))/(abs(vx(i)/V_T)**3)
	END IF

	if (Landau*dt >= 1.0) then	
		w_v(i) = 1e-20
	else
		w_v(i) = w_v(i) - w_v(i)*Landau*dt
	endif
END DO

END SUBROUTINE LandauDamping



SUBROUTINE W_EQUATION(W_v,k,Omega_x,d_Omega_x,dx)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  W_Equation
!
!PURPOSE:
!
! This subroutine takes the k-space step for the term
! \frac{\partial \omega}{\partial x}\frac{\partial W}{\partial k}
!
! CALLS:  Constant.f90 and Params.f90
!
! INPUTS:  
!	W(v) 		; array of spectral energy density this (x,t) as a f'n of velocity
!	k		; array of wave numbers
!	Omega_x		; omega_pe for this x
!	d_Omega_x	; difference between omega_pe at x and x_+1
!	dx		; difference between x and x_+1
!
! OPTIONAL INPUTS:
!  
! OUTPUTS:
!	W(v) 		; array of spectral energy density at (x_+1,t) as a f'n of velocity
! 
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	

USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: w_v,k;
DOUBLE PRECISION, INTENT(IN)::Omega_x,d_Omega_x,dx;
DOUBLE PRECISION :: aa_const

aa_const = - 100.*dt*v_t*d_Omega_x/(Omega_x*dx)

DO I=2*Nv-1,2,-1
! shift in K space
	aa = aa_const / (abs(k(i)-k(i-1)))

! 	aa = -dt*V_T*d_omega/(Omega_x*abs(k(i)-k(i-1))*dx)
! 	aa = +dt*V_T/(abs(kx(i)-kx(i-1))*1e+9)
! 	aa = dt*v(i)*v(i)*d_omega / (Omega_x*dx*dv))

	IF(abs(aa)>0.5) write (*,*) 'W_EQUATION:  plasma gradient is too large .... ', aa
	IF (aa < 0.) then
		w1 = aa*(W_v(i)-W_v(i-1))
	ELSE
		w1 = aa*(W_v(i+1)-W_v(i))
	end if
	w_v(i) = w_v(i) + w1

END DO

! shift in K space
END SUBROUTINE W_EQUATION
!----------------------------------------------------------------

SUBROUTINE VanLeer(Fminus2,Fminus1,Fzero,Fplus1,vx,Xminus2,Xminus1,Xzero,Xplus1)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  W_Equation
!
!PURPOSE:
!
! This subroutine is a monotonic transport finite difference method for the equation
!
! \frac{\partial F/{\partial t}+\gamma \frac{\partial F}{\partial X}=0
! 
!  taken from 
! "Towards the ultimate conservative difference scheme. IV. 
!  A new approach to numerical convection" 
! B. van Leer J. Comput. Phys. 23 276- 299(1977)
! 
! In this subroutine this is done for an array of "velocity" (\gamma above) stepping 
! each F_0 forward in time.
!
! CALLS:  Constant.f90 and Params.f90
!
! INPUTS:  
!  	F_-2(v), F_-1(v), F_0(v), F_+1(v),	; function to be stepped forward in time 
!	v					; array of velocities
!	X_-2, X_-1, X_0, X_+1			; X positions corresponding to F indexed positions
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: 
!	F_0(v) 					; Function stpped forward to time t+Dt as f'n of velocity
! 
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!
! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	


USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: Fminus2,Fminus1,Fzero,Fplus1;
DOUBLE PRECISION,DIMENSION(:), INTENT (IN) :: vx;
DOUBLE PRECISION,INTENT(IN) :: Xminus2,Xminus1,Xzero,Xplus1;
DOUBLE PRECISION:: d1,d2,DF_plus,DF_minus,a_minus,a_plus,v_local,a_minusconst,a_plusconst

a_minusconst	= dt / (Xzero-Xminus1)
a_plusconst	= dt / (Xplus1-Xzero)

DO i=2, Nv1-1
! velocity loop
	v_local = vx(i)
	d1 = (Fzero(i)-Fminus1(i)) * (Fplus1(i)-Fzero(i))
	IF ( d1 > 0.) THEN
		DF_plus = d1*(Xplus1-Xminus1) / ((Xplus1-Xzero)*(Fplus1(i)-Fminus1(i)))
	ELSE
		DF_plus = 0.
	END IF
	d2 = (Fminus1(i)-Fminus2(i)) * (Fzero(i)-Fminus1(i))
	IF ( d2 > 0.) THEN
		DF_minus = d2*(Xzero-Xminus2) / ((Xzero-Xminus1)*(Fzero(i)-Fminus2(i)))
	ELSE
		DF_minus = 0.
	END IF
	a_minus = v_local*a_minusconst
	a_plus = v_local*a_plusconst
	Fzero(i) = Fzero(i) -a_minus * (Fzero(i)-Fminus1(i) + (1.-a_plus)*DF_plus/2. - (1-a_minus)*DF_minus/2.)

END DO


END SUBROUTINE VanLeer

SUBROUTINE QUASILINEAR(F_x,W_x,velocity,Weq_coeff)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  Quasilinear
!
!PURPOSE:
!
! This subroutine takes care of the quasilinear terms  in the 1D equations
! of quasi-linear relaxation in an inhomogenous plasma (the right-hand terms)

! CALLS: Constant.f90 and Params.f90
!
! INPUTS:  
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
!	velocity	; array of velocities
!	Weq_koef	; "a2" or "\alpha-2" constant for spectral engergy density at this x
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: 
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	


USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: F_x,W_x;
DOUBLE PRECISION,DIMENSION(:), INTENT (IN) :: velocity;
DOUBLE PRECISION,INTENT (IN):: Weq_coeff;
DOUBLE PRECISION :: dvdv,Fconst,Wconst

dvdv 	= dv*dv
Fconst	= dt*a1
Wconst	= dt*Weq_coeff/dv

DO i=2, Nv1-1
       	! velocity loop
        v = velocity(i)
        m1 = W_x(i) * (F_x(i+1)-F_x(i)) / (v*dvdv)
        if (i==2) then
        	m2 = 0.
        else
        	m2 = W_x(i-1)*(F_x(i)-F_x(i-1))/((v-dv)*dvdv)
        end if
        F_x(i) = F_x(i) + Fconst*(m1-m2)
        W_x(i) = W_x(i) + Wconst*v*v*W_x(i) * (F_x(i+1)-F_x(i))
END DO



END SUBROUTINE QUASILINEAR

SUBROUTINE UPWIND(Fminus1,Fzero,velocity,Xminus1,Xzero)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  Simple first order finite difference operator. Produces a spatial step:
!
! \frac{\partial F/{\partial t}+\gamma \frac{\partial F}{\partial X} 
!
! In this subroutine this is done for an array of "velocity" (\gamma above) stepping 
! each F_0 forward in X.
!
!PURPOSE:
!
! This subroutine is 

! CALLS: Constant.f90 and Params.f90
!
! INPUTS:  
!  	F_-1(v), F_0(v)			; function to be stepped forward in time 
!	v				; array of velocities
!	X_-1, X_0,			; X positions corresponding to F indexed positions
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: F_0(v) 		; Function stepped forward to X t+Dt
! 
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	

USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: Fminus1,Fzero;
DOUBLE PRECISION,DIMENSION(:), INTENT (IN) ::velocity;
DOUBLE PRECISION,INTENT (IN) ::Xminus1,Xzero;

Fzero = Fzero - velocity*(Fzero-Fminus1)*dt/(Xzero-Xminus1)

END SUBROUTINE UPWIND

SUBROUTINE Spontaneous(F_x,W_x,velocity,b2)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  Spontaneous
!
!PURPOSE:
!
! This subroutine takes care of the spontaneous emission from resonant interaction of
! the electorns and waves.

! CALLS: Constant.f90 and Params.f90
!
! INPUTS:  
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
!	velocity	; array of velocities
!	b2		; "b2" or "\beta-2" constant for spectral engergy density at this x
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: 
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	
!  08/2008:	fixed some of the nonsense but correct??????


USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: F_x,W_x;
DOUBLE PRECISION,DIMENSION(:), INTENT (IN) :: velocity;
DOUBLE PRECISION,INTENT (IN):: b2;

DO i=1, Nv1-1
      	! velocity loop
        v = velocity(i)
	F_x(i) = F_x(i) + dt*b1*(F_x(i+1)-F_x(i))/dv
        W_x(i) = W_x(i) + dt*b2*v*v*v*F_x(i)
END DO

END SUBROUTINE Spontaneous

SUBROUTINE CoulombCollisions(F_x,W_x,velocity,gamF,gamW)
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  CoulombCollisions
!
!PURPOSE:
!
! This subroutine takes care of Coloumb collision for both 
! the electorns and waves.

! CALLS: Constant.f90 and Params.f90
!
! INPUTS:  
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
!	velocity	; array of velocities
!	gamF 		; the F equation coeff  (without V) for this x
!	gamW		; the W equation coeff (gamm_ei) for this x
!
! OPTIONAL INPUTS:
!  
! OUTPUTS: 
!	F(v) 		; electron dist at this x, t as a function of velocity
!	W(v) 		; spectral energy density this x, t as a function of velocity
! OPTIONAL OUTPUTS:
!   
!
! KEYWORDS:
!   
!
! COMMON BLOCKS:
!   
!
! SIDE EFFECTS:
!

! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: 	Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  	eduard(at)astro.gla.ac.uk: module structure added
!  06/2008:	Solver.f90 from reorganising Nonlin.f90	
!  08/2008:	fixed all the nonsense IGH


USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(:) :: F_x,W_x;
DOUBLE PRECISION,DIMENSION(:), INTENT (IN) :: velocity;
DOUBLE PRECISION,INTENT (IN):: gamF, gamW;
DOUBLE PRECISION :: vdv, fcolterm, wcolterm;

!W_x=W_x-gamW*dt*W_x

DO i=1, Nv1-1
      	! velocity loop
        v = velocity(i)
	vdv=velocity(i+1)
	! fcolterm should always be >0
	fcolterm=(gamF*dt/dv)*((F_x(i+1)/(vdv*vdv))-(F_x(i)/(v*v)))

	if (fcolterm > F_x(i)) then 
		F_x(i)=1.0E-20
	else
        	F_x(i) = F_x(i) + fcolterm
	endif
	wcolterm=gamW*dt*W_x(i)
	if (wcolterm > W_x(i)) then 
		W_x(i)=1.0E-20
	else
        	W_x(i) = W_x(i) - wcolterm
	endif

END DO

END SUBROUTINE CoulombCollisions

END MODULE SOLVER