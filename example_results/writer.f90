MODULE writer

IMPLICIT NONE

CONTAINS

Subroutine  WriteProfiles (Rx,Fvx, Wvx,WT_vx,WTx_vx, Tstep,Omegap,kx,v_x)
! The procedure is intended to write profiles
! of electron distribution function and  
! spectral energy density of Langmuir waves to the 
! respective disk files

  USE CONSTANT
! General definition of constants
  USE PARAMS
! Load certain global variables
  USE SERVICE

  IMPLICIT NONE

  DOUBLE PRECISION, dimension (:,:),INTENT(IN):: Fvx, Wvx,WT_vx, WTx_vx;
  DOUBLE PRECISION, dimension (:),  INTENT(IN):: Rx, Omegap;
  DOUBLE PRECISION, DIMENSION (:),  INTENT(IN):: Kx,V_x;
  DOUBLE PRECISION, dimension (:),allocatable :: n_beam,E_beam,E_wave,E_wave_o,E_wave_x,E_waveT;

  integer :: WriteStat
! Service variables 
  integer, intent(in)::Tstep
  DOUBLE PRECISION ::Norm_x, Norm_Fvx, Norm_Wvx, Norm_WT_vx,Norm_dens;
  DOUBLE PRECISION ::Norm_Fx, Norm_Fo, Norm_eb, Norm_ew,Norm_wt,Tb_o,Tb_x;

  character (LEN= 32):: vxfwd_flush, fxv_flush, wxv_flush, wTxv_flush
  character (LEN= 5) :: Findex
! Findex - the output file numbers
   
  IF (.NOT.ALLOCATED(N_beam))  ALLOCATE(N_beam(Nx),STAT=allocate_stat)
  IF (.NOT.ALLOCATED(E_beam))  ALLOCATE(E_beam(Nx),STAT=allocate_stat)
  IF (.NOT.ALLOCATED(E_wave))  ALLOCATE(E_wave(Nx),STAT=allocate_stat)
  IF (.NOT.ALLOCATED(E_waveT)) ALLOCATE(E_waveT(Nx),STAT=allocate_stat)
  IF (.NOT.ALLOCATED(E_wave_x)) ALLOCATE(E_wave_x(Nx),STAT=allocate_stat)
  IF (.NOT.ALLOCATED(E_wave_o)) ALLOCATE(E_wave_o(Nx),STAT=allocate_stat)

  IF(allocate_stat /= 0) STOP 'ARRAYS: N_beam, E_beam, etc. ARE NOT ALLOCATED'

  write(Findex,'(i5.5)')  Tstep
  Findex=trim(Findex)

  vxfwd_flush= trim(vxfwd_sequence//Findex//'.dat')  
  ! file with energy distribution 

  fxv_flush= trim(fxv_matrix//Findex//'.dat')  
  ! file with F(v,x)
  
  wxv_flush= trim(wxv_matrix//Findex//'.dat')  
  ! file with W(v,x)

  wTxv_flush= trim(wTxv_matrix//Findex//'.dat')  
  ! file with Wt_F(v,x)
 
 
  OPEN(UNIT=11, FILE = vxfwd_flush, STATUS = 'REPLACE', ACTION='WRITE' ) 
  OPEN(UNIT=12, FILE = fxv_flush, STATUS = 'REPLACE', ACTION='WRITE',RECL=2435456) 
  OPEN(UNIT=13, FILE = wxv_flush, STATUS = 'REPLACE', ACTION='WRITE',RECL=2435456) 
!   OPEN(UNIT=14, FILE = wTxv_flush, STATUS = 'REPLACE', ACTION='WRITE',RECL=2435456) 
 ! open file with density profile
 
  DO j=1, Nx

  Norm_x     = Rx(j) /d;

  N_beam(j)=k1*INTEGRAL(Fvx(:,j),V_x(:),1,Nv1)
  ! Calculation of the local density
  
  E_beam(j)=0.5d0*k1*INTEGRAL(Fvx(:,j)*V_x(:)*V_x(:),V_x(:),1,Nv1)
   ! local energy density of electrons

  E_wave(j) = k2*Omegap(j)*INTEGRAL(Wvx(:,j),kx(:),1,2*Nv)/v_T
  
   ! local energy density of Langmuir waves
  
  
  E_wave_o(j)= k2*Omegap(j)*INTEGRAL(WT_vx(:,j),kx(:),1,2*Nv)/v_T
  E_wave_x(j)= k2*Omegap(j)*INTEGRAL(WTx_vx(:,j),kx(:),1,2*Nv)/v_T
  
  E_waveT(j)= E_wave_x(j)+E_wave_x(j)
   ! local energy of electromagnetic waves

        Norm_Dens  = N_beam(j)/nb;
        Norm_eb    = e_beam(j)/Eb
	Norm_ew    = e_wave(j)/(m*Eb)

	Norm_wT    = e_waveT(j)/(m*Eb)
	Norm_Fx    = e_wave_x(j)/(m*Eb)
	Norm_Fo    = e_wave_o(j)/(m*Eb)
	
        Tb_o = ((2.*PI)**3)*MaxVal(WT_vx(:,j))*k2/(k_b*T_e)	
        Tb_x = ((2.*PI)**3)*MaxVal(WTx_vx(:,j))*k2/(k_b*T_e)

     if (Norm_eb  <10e-38)   Norm_eb=0.0
     if (Norm_ew  <10e-38)   Norm_ew=0.0
     if (Norm_Dens<10e-38) Norm_Dens=0.0

 WRITE (11,'(ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4)',IOSTAT = WriteStat ) Norm_x,' ', &
 Norm_dens,' ',Norm_eb,' ',Norm_ew,' ',Norm_wT,' ',Norm_Fo,' ',Norm_Fx, ' ',Tb_o,' ',Tb_x
 IF (WriteStat /= 0) THEN 
 WRITE (output,*) Norm_x,' ',Norm_dens,' ',Norm_eb,' ',Norm_ew,' ',Norm_wT
 STOP 'Error while writing F, W, and density profile'
 END IF

  DO i=1, NV1

     Norm_Fvx  = Fvx(i,j)!*k1*vbeam/nbeam;

     IF (Norm_Fvx<10e-10) Norm_Fvx=0.0

  WRITE (12,'(ES14.4, A)', ADVANCE="no", IOSTAT = WriteStat ) Norm_Fvx,' '
  IF (WriteStat /= 0) stop 'Error while writing F(x,v)' 
  ! Write in matrix form  - line by line 
  ! Each line - range in velocity
  ! Diagnostics to make sure that disk writes are successful

  END DO
  WRITE (12, '(A1)', ADVANCE ='yes') ' '
  ! new line 

  DO i=1, 2*NV
       Norm_Wvx  = 2.*pi**2*v_t**2*Wvx(i,j)*k2/(k_b*T_e*Omegap(j)**2);!((2.*PI)**3)*Wvx(i,j)*k2/(k_b*T_e)
       ! brightness temperature

! 	   IF (Norm_Wvx<10e-10) Norm_Wvx=0.0      	  
       WRITE (13,'(ES14.4, A)', ADVANCE="NO", IOSTAT = WriteStat ) Norm_Wvx,' '
       IF (WriteStat /= 0) stop 'Error while writing W(x,v)'

  END DO
  WRITE (13, '(A1)', ADVANCE ='yes') ' '

!   DO i=1, 2*NV
! 
        Norm_Wt_vx = (2.*PI)**3*Wt_vx(i,j)*k2/(k_b*T_e)
!        ! brightness temperature
! 
! 	   IF (Norm_WT_vx<1e-10) Norm_Wt_vx=0.0
! 	         	  
!        WRITE (14,'(ES14.4, A)', ADVANCE="NO", IOSTAT = WriteStat ) Norm_Wt_vx,' '
!        IF (WriteStat /= 0) stop 'Error while writing W_F(x,v)' 
!        
!   END DO
!   WRITE (14, '(A1)', ADVANCE ='yes') ' '

  END DO

  CLOSE (11)
  CLOSE (12) 
  CLOSE (13)
!   CLOSE (14)     

! Close files 

   DEALLOCATE(n_beam, E_beam, E_wave, E_wave_o, E_wave_x);

WRITE(output,'(A,F9.3,A)') 'Data saved at T=',t,' sec ...................... OK';
   
END SUBROUTINE WriteProfiles

!----------------------------------------------------------------------------

SUBROUTINE WRITE_DENSITY(Rx,dens_x, tau_qv,Gamma_ei)

USE CONSTANT;

USE PARAMS;

IMPLICIT NONE

 DOUBLE PRECISION, DIMENSION(:),INTENT(IN)::Rx,dens_x, tau_qv,Gamma_ei;
 DOUBLE PRECISION :: xlocal,Omegap,L,Gamma_ei_x;
 INTEGER::WriteStat;

 
 open (UNIT=13, FILE = corona_file, STATUS = 'REPLACE', ACTION='WRITE' ) 
 open (UNIT=19, FILE = tau_qv_file, STATUS = 'REPLACE', ACTION='WRITE' ) 

  L=0.

 do i=1, nx

  xlocal= x_0- Rx(i)
  !height about photo  
  Omegap =sqrt(4.0*pi*e*e*dens_x(i)/m)
  
  gamma_ei_x= gamma_ei(i)
  ! damping rate for Langmuir waves

  L= sqrt(Dens_x(i))*(Rx(i+1)-Rx(i))/(sqrt(dens_x(i+1))-sqrt(dens_x(i)))
      
  write (13,'(ES14.4,A,ES14.4,A,ES14.4,A,ES14.4,A,ES14.4)', IOSTAT=WRITESTAT) xlocal/1e+9,' ',dens_x(i),' ',Omegap/1e+6,' ', &
  gamma_ei_x,' ',L/R_s
  ! Write  - line by line
 
  write (19,*, IOSTAT=WRITESTAT) xlocal/R_S,' ',(xlocal-R_S-x_0)/d,' ',Tau_qv(i)
  If (WriteStat /= 0) THEN 
  write (output,*) 'STEP=',i, xlocal/R_s,' ',xlocal/d,' ',Tau_qv(i)
  stop 'Error while writing quasiliner time' 
  END IF
  ! Write  - line by line
 
 end do 
 
 close (13)
 close (19)
 
 WRITE(output,*)'Plasma parameters saved to disk       -  OK'
 END SUBROUTINE WRITE_DENSITY

!************************************************************
SUBROUTINE DATA_SHOW(Omega_x);

USE CONSTANT
USE PARAMS

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN)::Omega_x

write(output,*) '------------------------------------------'
write(output,'(A,I6,A,I6,A,F6.1)')'Nx =',Nx,'  Nv=',Nv,' Xmax =', Xmax/d
write(output,'(A,F9.3,A,F9.3,A)')'F_max =', Omega_0/(2.*PI*1e9), ' GHz F_min =',Omega_x(Nx)/(2.*PI*1e9),' GHz'
write(output,'(A,F9.3,A)')'Omega_ce /2Pi    =', 1e-6*omega_ce/(2.*PI), ' MHz'
write(output,'(A,F9.2,A)')'Thermal velocity =', V_t/1e+5, ' km/s'
write(output,'(A,F9.2)')'V_beam/v_Te      =', Vbeam/v_T
write(output,'(A,ES14.4,A)')'Beam density     =', nbeam, ' cm^{-3}'
write(output,'(A,ES14.4,A)')'Plasma density   =', ne, ' cm^{-3}'   
write(output,'(A,ES14.4,A)')'Quasilinear time =', ne/(omega_0*nbeam),' sec - INITIAL'
write(output,*) '------------------------------------------'

END SUBROUTINE DATA_SHOW

END MODULE WRITER
