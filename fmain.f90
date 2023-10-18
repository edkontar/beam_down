PROGRAM QBEAM_PLASMA
!
! PROJECT:
!   Beam-plasma interaction
!
! NAME:
!
!  QBEAM_PLASMA
!
!PURPOSE:
!
! This program solves the one-dimensional kinetic equations
! of quasilinear relaxation and weak turbulence using explicit
! difference schemes (second order approximation over coordinate
! and velocity, and first order approximation over time)
!
! CALLS:  various subroutines
!
! INPUTS:  see   init.par file
!
! OPTIONAL INPUTS:
!  none
! OUTPUTS:
! F       - electron distribution function
! W     - wave spectral energy density
! The results of the code are returned as data files

! OPTIONAL OUTPUTS:
!   none
!
! KEYWORDS:
!   none
!
! COMMON BLOCKS:
!   none
!
! SIDE EFFECTS:
!
! RESTRICTIONS:

! MODIFICATION HISTORY:
!  1999-2002: Version 1 written by eduardk(at)astro.uio.no 
!  02/2002:  eduard(at)astro.gla.ac.uk: module structure added
!  uploded to github October 2023


  USE CONSTANT

  USE PARAMS  

  USE READER

  USE WRITER  
  
  USE INITBEAM
  
  USE PLASMA

  USE SERVICE
  
  USE NONLIN
  
  USE SOLVER
  !USE 'MPIF'
  
  IMPLICIT NONE

! include 'mpif.h'

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: F, W, W_F,W_Fx;
	! Electron distribution function and
	! Spectral energy density
	
    !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::Fp, Wp,W_Fp,W_Fxp;
    ! arrays for paralel calculations
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: Gamma_ei, Gamma_nov;
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: Rx,DENSITY,OMEGA,Wind,Tqv,a2,b2;
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: Vx,Kx;
    ! Velocity and wave number arrays

  !  INTEGER:: status,psize,my_rank,p,Np,rstat
    !mpi variables
    	  	 
!------------  Main plasma values definition ---------------------
LAST_save=0.
!time_control=0

! reading initial configuration
! in reader.f90 and reads init.par
CALL ReadCONFIG

! calculating the number of x points 
Nx=(xmax-xmin)/(d/10.) 

!----------- ALOCATING MAIN ARRAYS --------------------------
IF (.NOT.ALLOCATED(Rx))      ALLOCATE (Rx(Nx))
IF (.NOT.ALLOCATED(Density)) ALLOCATE (Density(Nx))
IF (.NOT.ALLOCATED(Omega))   ALLOCATE (Omega(Nx))
IF (.NOT.ALLOCATED(Gamma_ei))ALLOCATE (Gamma_ei(Nx))
IF (.NOT.ALLOCATED(Gamma_nov))ALLOCATE (Gamma_nov(Nx))
IF (.NOT.ALLOCATED(Wind))    ALLOCATE (Wind(Nx))
IF (.NOT.ALLOCATED(F))       ALLOCATE (F(2*Nv,Nx))
IF (.NOT.ALLOCATED(W))       ALLOCATE (W(2*Nv,Nx))
IF (.NOT.ALLOCATED(W_F))     ALLOCATE (W_F(2*Nv,Nx))
IF (.NOT.ALLOCATED(W_Fx))    ALLOCATE (W_Fx(2*Nv,Nx))
IF (.NOT.ALLOCATED(Tqv))     ALLOCATE (Tqv(Nx))
IF (.NOT.ALLOCATED(a2))      ALLOCATE (a2(Nx))
IF (.NOT.ALLOCATED(b2))      ALLOCATE (b2(Nx))
IF (.NOT.ALLOCATED(Vx))      ALLOCATE (Vx(2*Nv))
IF (.NOT.ALLOCATED(Kx))      ALLOCATE (Kx(2*Nv))
    
!--------------------------------------------------------------



! In plasma.f90
! calculates denisty_x and omega_x
CALL newloop_den(xmin,xmax,Density,Omega,Rx)

!In plasma.f90
! constructs k-mesh
CALL KNUMBER(Vx,Kx)

  
! In initbeam.f90
! calculates beam density and energy and the coefficients
CALL initial_values(nb,eb,Rx,Density,Omega,Tqv,a2,Gamma_ei,Gamma_nov,b2)
 
! In initbeam.f90
! initial F and W made
CALL initial_FW(vx,F,W,W_F,W_Fx,Rx,Omega)

! ! In writer.f90
! ! write out the initial denisty
CALL write_Density(Rx,density,tqv,Gamma_ei)
! 
! ! In writer.f90
! ! write out initial profile of F and W
 CALL WriteProfiles(Rx,F,W,W_F,W_Fx,0,OMEGA,kx,vx)
!  
! In solver.f90
! calculates energy in the waves and particles
CALL ENERGY(Rx,F,W,W_F,W_Fx,OMEGA,kx,vx,0)

! In writer.f90
! show the parameters (sent to output)
CALL DATA_SHOW(Omega);

! -------------------------------------------------
    call system_clock(count=itime, count_rate=rate, count_max=max)
    time_begin  = Float(itime)/rate
    Time_control= time_save/10.
 
!------------- MAIN LOOP starts here ---------------
 Main_loop : DO WHILE (t .LE. tmax)
! ---------------- TIME STEP ----------------
 
! fixed dt has to be less that tqv, so depend on beam density,
! approx n~10^4 cm^-3 needs dt~10^-5s
! approx n~10^6 cm^-3 needs dt~10^-7s  etc
! dt now read in from init.par
t=t+dt
 
DISPLAY_TIME: IF ((t-last_control) .GE.Time_control) THEN
	! accumulating to the main array
	CALL ENERGY(Rx,F,W,W_F,W_Fx,OMEGA,kx,vx, FLOOR(t/dt))
	CALL CHECK_TIME 
	SAVE_TIME: IF ((t-last_save).GE.time_save) THEN    
		CALL WriteProfiles(Rx,F,W,W_F,W_Fx,FLOOR(t/time_save),OMEGA,kx,vx)
		!saving data to disk
		last_save = time_save*FLOOR(t/time_save);
	END IF SAVE_TIME
	last_control = time_control*FLOOR(t/time_control); 
 END IF DISPLAY_TIME

X_LOOP: DO j=3,Nx-2
	! x-coordinate loop
	! all in solver.f90
  	CALL W_EQUATION(W(:,j),kx,Omega(j),Omega(j+1)-Omega(j),Rx(j+1)-Rx(j))
	CALL VanLeer(F(:,j-2),F(:,j-1),F(:,j),F(:,j+1),vx,Rx(j-2),Rx(j-1),Rx(j),Rx(j+1))
 	CALL VanLeer(W(:,j-2),W(:,j-1),W(:,j),W(:,j+1),3.*v_T*v_T/vx,Rx(j-2),Rx(j-1),Rx(j),Rx(j+1))
 	CALL QUASILINEAR(F(:,j),W(:,j),vx,a2(j))
 	CALL LandauDamping(W(:,j),vx,Omega(j))
	CALL CoulombCollisions(F(:,j),W(:,j),vx,gamma_nov(j),gamma_ei(j))
!	CALL Spontaneous(F(:,j),W(:,j),vx,b2(j))
END DO X_LOOP

! boundary conditions at the bottom of the box
F(:,1)=F(:,2)
F(:,2)=F(:,3)
W(:,2)=W(:,3)
W(:,1)=W(:,2)

! boundary conditions at the top of the box
F(:,Nx)=F(:,Nx-1)   
W(:,Nx)=W(:,Nx-1)

END DO Main_loop
!End of main loop
 
! ----------- main array dealocation --------------------------
DEALLOCATE (Rx)
DEALLOCATE (Density)
DEALLOCATE (Omega)
DEALLOCATE (Wind)
DEALLOCATE (F)
DEALLOCATE (W,W_F,W_Fx)
DEALLOCATE (Tqv)
DEALLOCATE (a2,b2)
DEALLOCATE (Gamma_ei, Gamma_nov)
DEALLOCATE (Vx,Kx)

WRITE(output,*) 'Main arrays deallocated                 -  OK'   
!--------------------------------------------------------------
WRITE(output,*) "CALCULATIONS COMPLETED     -    OK !"
IF (output/=0) close(output);

END PROGRAM
