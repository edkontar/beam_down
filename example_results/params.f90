MODULE PARAMS

! This module contains parameters of calculations
! and variables 


USE CONSTANT

	IMPLICIT NONE

	INTEGER :: i,j;  
	INTEGER :: time_flush,Nx,output;
	! save data Time_flush number of time steps

	DOUBLE PRECISION :: w0,w1, w2, m1, m2,k1,k2,a1,b1;
	! coeficients of kinetic equations

	DOUBLE PRECISION :: v,x, x_0, t, dt, dt2, dv,vmin;
	! intrinsic variables 

	DOUBLE PRECISION :: ne, Eb,Nb;
	! plasma density, quasilinear time, 
	! beam density, and beam energy density

	DOUBLE PRECISION ::aa;        

	DOUBLE PRECISION:: NuPlasma, OMEGA_0;
	!Plasma frequency - defines plasma density
	! Omega = NuPlasma *2*Pi
	
	DOUBLE PRECISION:: Nbeam;
	! Electron beam density

	DOUBLE PRECISION:: Tmax;
	! Time of the calculations

	DOUBLE PRECISION:: Vbeam
	! Velocity of electron beam

	DOUBLE PRECISION:: Xmin, Xmax;
	! define the geometry of calculation area
	! Starting and End of simulation box

	DOUBLE PRECISION:: d;
	! Spatial size of the initial cloud
	
	DOUBLE PRECISION:: Time_save;
	! Time in seconds of each save

	DOUBLE PRECISION:: V_T,t_scat,V_Ti,omega_ce;
	! Thermal electron beam velocity

	DOUBLE PRECISION:: landau_cutoff,landau_min, two_rootpi
;
	! Damping  values

	INTEGER:: allocate_stat, io_stat;

	DOUBLE PRECISION:: time_begin;

	DOUBLE PRECISION:: time_end, last_save,last_control;
	! time diagnostic variables

	DOUBLE PRECISION::  remain_time, time_est, Time_control
	!service variables

	integer :: itime, rate, max, time_step;
	! CPU_time  variables

	INTEGER, DIMENSION(1) :: value1;
	INTEGER, DIMENSION(2) :: Value2;
	INTEGER :: Max_loc,Max_loc_p;
	DOUBLE PRECISION :: max_f, max_f_p;

END MODULE PARAMS