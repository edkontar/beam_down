MODULE SERVICE

IMPLICIT NONE

CONTAINS

SUBROUTINE CHECK_TIME
 ! subroutine that calculates remaining time

  USE PARAMS

  IMPLICIT NONE

  DOUBLE PRECISION::timeacc;
       
    time_step = time_step + 1
     
    CALL SYSTEM_CLOCK(count=itime, count_rate=rate, count_max=max)
  
    time_end = float(itime)/rate

    time_est=time_end - time_begin
         
    if (time_est<0) time_est = time_est + max/rate
    ! condition checks if counter started from zero again 

    timeacc = timeacc + time_est
   	remain_time= (Tmax-T)*timeacc/(T-Tmin)
      ! Find remaining time   

Write (output , '(A,F10.4,A, F5.1,A,F10.1,A)') 'Passed T=',t ,'sec (',t*100/(tmax-tmin),'% completed) for ', time_est,' seconds'				
!   Translate remaining time to hrs.mm.sec and write to screen
 	 time_begin = time_end
	   
END SUBROUTINE CHECK_TIME
!-----------------------------------------------------------
SUBROUTINE WriteHRS(time0)

 USE PARAMS
!Translating remaining time to hrs.mm.sec and write to screen
    
      DOUBLE PRECISION, intent(in):: time0
      DOUBLE PRECISION :: timex

      integer:: hrs,mins,secs
      
      timex=time0
      hrs=Floor(timex/3600)
      timex=timex-hrs*3600

      mins=Floor(timex/60)
      timex=timex-mins*60

      secs=timex
  
 WRITE(output,'(A,I5,A, I2, A, I2, A)')  'Remaining time:', hrs,' hrs ', mins,' mins ',secs,' sec.' 

 END SUBROUTINE WriteHRS

!-----------------------------------------------------------

SUBROUTINE  ENERGY(Rx,Fvx,Wvx,Wt_vx,Wtx_vx,Omegap,kx,v_x,Time)
! The procedure calculates the energy of waves 
!and particles

  USE constant
! General definition of constants
  USE params
! Load certain global variables

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION (:,:),INTENT(IN):: Fvx, Wvx, Wt_vx,Wtx_vx;
  DOUBLE PRECISION, DIMENSION (:),  INTENT(IN):: Rx, Omegap;
  DOUBLE PRECISION, DIMENSION (:),  INTENT(IN):: Kx,V_x;
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: n_beam,E_beam,E_wave,E_Twave;

  DOUBLE PRECISION ::N_electrons,E_electrons,E_waves,E_Twaves,E_total;
  DOUBLE PRECISION, save ::Ne0, Etotal0 
  INTEGER, INTENT(IN):: TIME

  IF (.NOT.ALLOCATED(N_beam))  ALLOCATE(N_beam(Nx))
  IF (.NOT.ALLOCATED(E_beam))  ALLOCATE(E_beam(Nx))
  IF (.NOT.ALLOCATED(E_wave))  ALLOCATE(E_wave(Nx))
  IF (.NOT.ALLOCATED(E_Twave)) ALLOCATE(E_Twave(Nx))

  DO j=2, Nx

  N_beam(j)=k1*INTEGRAL(Fvx(:,j),V_x(:),1,Nv1)
  ! Calculation of the local density
  
  E_beam(j)=0.5d0*k1*INTEGRAL(Fvx(:,j)*V_x(:)*V_x(:),V_x(:),1,Nv1)
   ! local energy density of electrons

  E_wave(j)=k2*Omegap(j)*INTEGRAL(Wvx(:,j),kx(:),1,2*Nv)/(v_T)
   ! local energy density of Langmuir waves
  
!  E_Twave(j)=k2*Omegap(j)*INTEGRAL((Wt_vx(:,j)+Wtx_vx(:,j)),kx(:),1,2*Nv)/v_T
   ! local energy density of electromagnetic waves

  END DO

  N_electrons= INTEGRAL(N_beam(:),Rx(:),1,Nx)/(d*Nbeam)
  !total number of electrons

  E_electrons= INTEGRAL(E_beam(:),Rx(:),1,Nx)/(d*Eb)
  ! energy of electrons

  E_waves = INTEGRAL(E_wave(:),Rx(:),1,Nx)/(d*Eb*m)
  ! energy of Langmuir waves
  
  E_Twaves = INTEGRAL(E_Twave(:),Rx(:),1,Nx)/(d*Eb*m)
  ! energy of electromagnetic waves

  E_total= E_electrons+E_waves+E_Twaves
  ! total energy in the system
  
  IF (Time==0) THEN
  Ne0=N_electrons
  Etotal0 =E_total
  END IF
! write(output,'(ES14.4,ES14.4)') N_electrons, Ne0
  WRITE(output,'(A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F8.5)')'Ne=',N_electrons/Ne0,' Ee=',E_electrons/Etotal0,' EL=',&
  E_waves/Etotal0,' EF=',E_Twaves/Etotal0, ' E=',E_total/Etotal0
  ! writes output to the screen

  DEALLOCATE(N_beam,E_beam,E_wave, E_Twave)

END SUBROUTINE ENERGY

DOUBLE PRECISION FUNCTION INTEGRAL(F_integral,X_integral,n_begin,n_end)
! Function returns integral of 1D array F_integral
! over array X_integral using trapezoid method
! each array has Npoints elements of DOUBLE PRECISION
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:),INTENT(IN):: F_integral,X_integral;
  INTEGER,INTENT(IN)::n_begin,n_end;
  
  !IF (.NOT.ALLOCATED(F_integral)) ALLOCATE(F_integral(n_end))
  !IF (.NOT.ALLOCATED(X_integral)) ALLOCATE(X_integral(n_end))

  DOUBLE PRECISION:: sum;
  INTEGER::i
  
  sum=0.0d0;

  DO i=n_begin+1, n_end
  sum = sum + 0.5d0*(F_integral(i)+F_integral(i-1))*abs(X_integral(i)-X_integral(i-1))
  END DO

  INTEGRAL=sum;

  !DEALLOCATE(F_integral,X_integral)
END FUNCTION INTEGRAL 

END MODULE SERVICE
