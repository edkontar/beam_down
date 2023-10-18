 Module constant

! This module contains constant of calculations
! There are no inputs to this file and all varialbes
! are parameterised.

	DOUBLE PRECISION, PARAMETER::PI=3.141593D0  

	DOUBLE PRECISION, PARAMETER:: tmin=0;
	! Initial time

	DOUBLE PRECISION, PARAMETER:: vmax1=2.8e+10;
	! Maximal velocity

	DOUBLE PRECISION, PARAMETER:: vmax=3.0e+10;
	! Maximal velocity

	DOUBLE PRECISION, PARAMETER:: e=4.8d-10;
	! Electron charge
	
	DOUBLE PRECISION, PARAMETER::v_c=2.9978e+10;
	!speed of light

	DOUBLE PRECISION, PARAMETER::  m=9.1d-28;
	! electron mass
	
	DOUBLE PRECISION, PARAMETER:: m_p= 1.6726485e-24;
	!proton mass
	
        DOUBLE PRECISION, PARAMETER:: k_b = 1.380662e-16;
	!Boltzman constant
	
	DOUBLE PRECISION, PARAMETER::  G=6.6720e-8;
	! gravitational constant
	
	DOUBLE PRECISION, PARAMETER:: mu=0.6;
	! mean molecular value
	! Priest E.R.(1982) Solar Magnetohydrodynamics. Reidel,Dordrecht

	DOUBLE PRECISION, PARAMETER:: M_s=1.99e+33;
	!Solar mass
	
	DOUBLE PRECISION, PARAMETER:: R_s=6.958e+10 
	!Solar radius
	
	DOUBLE PRECISION, PARAMETER:: T_e =1e+6
	!Electron plasma temperature

	DOUBLE PRECISION, PARAMETER:: T_i =1e+6	
        ! Ion plasma temperature
 
        DOUBLE PRECISION, PARAMETER:: B =0.1d0
	! magnetric field 0.1 G
	
	DOUBLE PRECISION, PARAMETER:: AU = 1.5e+13
	! 1 Astronomical Unit in cm
		
	INTEGER, PARAMETER:: nt=10000;
	! number of time steps
	
	INTEGER, PARAMETER:: nv1= 60
	! Number of V cells

	INTEGER, PARAMETER:: nv = 120
	! Number of V cells

        DOUBLE PRECISION, PARAMETER :: fi=0.174d4
        ! fi =10*PI/180 = 10 degrees

        Character(*),parameter:: vxfwd_sequence = 'fwd'
        Character(*),parameter:: fxv_matrix =  'fxv'
        Character(*),parameter:: wxv_matrix =  'wxv'
        Character(*),parameter:: wtxv_matrix = 'txv'
    	character(*),parameter:: corona_file = 'xdf.dat'
        character(*),parameter:: Tau_qv_file = 'xTau_qv.dat'
    	character(*),parameter:: kv_mesh = 'kx_vx.dat'
	! OUTPUT File names

        Character(*),parameter:: init_params = 'init.par'
        ! Input file name 

	Character(*),parameter:: outfile = 'output.log'
        ! Input file name 

        CHARACTER (8):: xvfwn, st; 

End  Module constant
