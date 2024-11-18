C======================================================================
      SUBROUTINE HOS_WAVE_3D(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,
     1     TRIGSX,TRIGSY,PEX,PEY,TIME,DT,WVN,ETA,EX,EY,VPS,FETA,FVPS,
     1     US,VS,WS,PA,FR2,PERIOD,NCPU)

C     BY DI YANG, 08/2007

C     MAIN SUBROUTINE FOR H.O.S SIMULATION OF 3D WAVE
C     4TH ORDER RUNGE-KUTTA SCHEME IS USED FOR TIME INTEGRATION

      IMPLICIT NONE

      INTEGER IRK
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU,ISMOOTH2

      REAL PEX,PEY,TIME,DT,FR2,PERIOD

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),ETA0(NXMAX,NYMAX/NCPU)
      REAL EX(NXMAX,NYMAX/NCPU),EY(NXMAX,NYMAX/NCPU)
      REAL VPS(NXMAX,*),VPS0(NXMAX,NYMAX/NCPU)
      REAL VPSX(NXMAX,NYMAX/NCPU),VPSY(NXMAX,NYMAX/NCPU)
      REAL FETA(NXMAX,NYMAX/NCPU,4),FVPS(NXMAX,NYMAX/NCPU,4)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)
      REAL PA(NXMAX,*)
      REAL WVN(NXMAX,NYMAX/NCPU,NPW),ZP(NXMAX,NYMAX/NCPU,NPW)
      REAL R(NXMAX,NYMAX/NCPU,NPW)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C-------------------------------------------------
C     SAVE ETA AND VPS AT TIME STEP N
C-------------------------------------------------
      CALL SURF_RK4(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,ETA0,VPS0,
     1     FETA,FVPS,DT,1,NCPU)
C-----END HERE

C------------------------------------------------------------------
C     4TH ORDER RUNGE-KUTTA SCHEME FOR EVOLUTION OF ETA AND VPS
C------------------------------------------------------------------

      CALL RIGH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ETA,EX,EY,VPSX,VPSY,WS,FETA,FVPS,PA,1,FR2,NCPU)

      DO 100 IRK=2,4
         
         CALL SURF_RK4(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,ETA0,VPS0,
     1        FETA,FVPS,DT,IRK,NCPU)

         CALL DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1        PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)

         CALL ZETA(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,ETA,ZP,WVN,NCPU)

         CALL BOUNDVP(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,VPS,R,ZP,WVN,NCPU)

         CALL WSURF(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,WS,R,ZP,WVN,NCPU)

         CALL RIGH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1        PEX,PEY,ETA,EX,EY,VPSX,VPSY,WS,FETA,FVPS,PA,IRK,FR2,NCPU)

 100  CONTINUE
      
      CALL SURF_UPDATE(NXMOD,NYMOD,NXMAX,NYMAX,DT,ETA,VPS,ETA0,VPS0,
     1     FETA,FVPS,NCPU)
      
C-----END HERE

C----------------------------------------------------------------
C     VELOCITY, PD(ETA)/DT, AND PD(VPS)/DT AT TIME STEP N+1
C----------------------------------------------------------------

      CALL SMOOTH1(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,WVN,ETA,VPS,0.95,NCPU)

      CALL DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)

      CALL SMOOTH2(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,ETA,VPS,EX,EY,TIME/PERIOD,ISMOOTH2,NCPU)

      IF(ISMOOTH2.EQ.1) THEN
         CALL DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1        PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)
      ENDIF
      
      CALL ZETA(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,ETA,ZP,WVN,NCPU)
      
      CALL BOUNDVP(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,VPS,R,ZP,WVN,NCPU)
      
      CALL WSURF(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,WS,R,ZP,WVN,NCPU)
      
      CALL UVSURF(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,US,VS,WS,EX,EY,VPSX,VPSY,NCPU)

      RETURN
      END

C=====SUBROUTINE HOS_WAVE_3D END HERE






C======================================================================
      SUBROUTINE SURF_RK4(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,ETA0,VPS0,
     1     FETA,FVPS,DT,NRK4,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE DF/DT BY 4TH ORDER RUNGE-KUTTA SCHEME

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NRK4,NCPU
      
      REAL FAC,DT

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA0(NXMAX,NYMAX/NCPU),ETA(NXMAX,NYMAX/NCPU)
      REAL VPS0(NXMAX,NYMAX/NCPU),VPS(NXMAX,NYMAX/NCPU)
      REAL FETA(NXMAX,NYMAX/NCPU,4),FVPS(NXMAX,NYMAX/NCPU,4)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      
      IF(NRK4.EQ.1) THEN
         
         DO I=1,NXMOD
            DO J=1,NYMAX/NCPU
               ETA0(I,J)=ETA(I,J)
               VPS0(I,J)=VPS(I,J)
            ENDDO
         ENDDO

      ELSE IF(NRK4.EQ.2.OR.NRK4.EQ.3.OR.NRK4.EQ.4) THEN

         FAC=0.5*DT
         IF(NRK4.EQ.4) FAC=DT

         DO I=1,NXMOD
            DO J=1,NYMAX/NCPU
               ETA(I,J)=ETA0(I,J)+FAC*FETA(I,J,NRK4-1)
               VPS(I,J)=VPS0(I,J)+FAC*FVPS(I,J,NRK4-1)
            ENDDO
         ENDDO

      ELSE
         
         PRINT*, "Invalid value for IRK in SURF_RK4 !"
         PRINT*, "IRK=",NRK4
         STOP

      ENDIF

      RETURN
      END

C=====SUBROUTINE SURF_RK4 END HERE






C=======================================================================
      SUBROUTINE SURF_UPDATE(NXMOD,NYMOD,NXMAX,NYMAX,DT,ETA,VPS,ETA0,
     1     VPS0,FETA,FVPS,NCPU)

C     BY DI YANG, 08/2007

C     UPDATE SURFACE DATA AFTER 4TH STEP RUNGE-KUTTA

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU

      REAL DT

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA0(NXMAX,NYMAX/NCPU),ETA(NXMAX,NYMAX/NCPU)
      REAL VPS0(NXMAX,NYMAX/NCPU),VPS(NXMAX,NYMAX/NCPU)
      REAL FETA(NXMAX,NYMAX/NCPU,4),FVPS(NXMAX,NYMAX/NCPU,4)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      DO I=1,NXMOD
         DO J=1,NYMAX/NCPU
            ETA(I,J)=ETA0(I,J)+DT*(FETA(I,J,1)+2.*FETA(I,J,2)
     1           +2.*FETA(I,J,3)+FETA(I,J,4))/6.
            VPS(I,J)=VPS0(I,J)+DT*(FVPS(I,J,1)+2.*FVPS(I,J,2)
     1           +2.*FVPS(I,J,3)+FVPS(I,J,4))/6.
         ENDDO
      ENDDO
      RETURN
      END

C=====SUBROUTINE SURF_UPDATE END HERE







C====================================================================
      SUBROUTINE ZETA(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,ETA,ZP,WVN,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE (ETA^K)/(K!) FOR TAYLOR EXPANSION

      IMPLICIT NONE
      
      INTEGER I,J,K
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU
      
      REAL PEX,PEY
      
      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,NYMAX/NCPU)
      REAL ZP(NXMAX,NYMAX/NCPU,NPW)
      REAL WVN(NXMAX,NYMAX/NCPU,NPW)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            ZP(I,J,1)=ETA(I,J)
         ENDDO
      ENDDO
      
      DO K=2,NPW-1
         DO J=1,NYMAX/NCPU
            DO I=1,NXMOD
               ZP(I,J,K)=ZP(I,J,K-1)*ZP(I,J,1)/FLOAT(K)
            ENDDO
         ENDDO
         CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1        TRIGSY,IFAX,ZP(1,1,K),NCPU)
      ENDDO

      RETURN
      END

C=====SUBROUTINE ZETA END HERE







C=======================================================================
      SUBROUTINE BOUNDVP(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,VPS,R,ZP,WVN,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE VALUE OF PHI AT MEAN SURFACE LEVEL

      IMPLICIT NONE
      
      INTEGER I,J,L,M,K,K1,K2,MODEX,MODEY
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU
      
      REAL PEX,PEY,AX,AY,AN

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL VPS(NXMAX,NYMAX/NCPU)
      REAL R(NXMAX,NYMAX/NCPU,NPW),DR(NXMAX,NYMAX/NCPU)
      REAL ZP(NXMAX,NYMAX/NCPU,NPW),WVN(NXMAX,NYMAX/NCPU,NPW)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DO I=1,NXMOD
         DO J=1,NYMAX/NCPU
            R(I,J,1)=VPS(I,J)
         ENDDO
      ENDDO

C-----------------------------------------
C     CALCULATION IN SPECTRAL SPACE
C-----------------------------------------

      CALL FFTXY_PARA(R(1,1,1),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1     NYMOD,NXMAX,NYMAX,NCPU,-1)

      DO K=2,NPW

         DO I=1,NXMAX
            DO J=1,NYMAX/NCPU
               R(I,J,K)=0.
            ENDDO
         ENDDO
         
         DO K1=1,K-1
            
            K2=K-K1

            DO M=1,NYMAX/NCPU
               DO L=1,NXMOD+2
                  DR(L,M)=R(L,M,K2)*WVN(L,M,K1)
               ENDDO
            ENDDO
            
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     APPLY PSEUDO-SPECTRAL METHOD
C     DO MULTIPLICATION IN PHYSICAL SPACE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            CALL FFTXY_PARA(DR,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,
     1           NXMAX,NYMAX,NCPU,1)

            DO I=1,NXMOD
               DO J=1,NYMAX/NCPU
                  R(I,J,K)=R(I,J,K)-ZP(I,J,K1)*DR(I,J)
               ENDDO
            ENDDO

C~~~~~END HERE

         ENDDO

         CALL FFTXY_PARA(R(1,1,K),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1        NYMOD,NXMAX,NYMAX,NCPU,-1)
         
      ENDDO

      DO K=1,NPW
         CALL FFTXY_PARA(R(1,1,K),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1        NYMOD,NXMAX,NYMAX,NCPU,1)
      ENDDO

C-----END HERE

      RETURN
      END

C=====SUBROUTINE BOUNDVP END HERE







C=======================================================================
      SUBROUTINE WAVENUM(NXMOD,NYMOD,NXMAX,NYMAX,NPW,PEX,PEY,WVN,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE EXPONENTIAL OF WAVENUMBER

      IMPLICIT NONE

      INTEGER L,M,K,MODEX,MODEY
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU

      REAL PEX,PEY,AX,AY,AN

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL WVN(NXMAX,NYMAX/NCPU,NPW)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      
      DO K=1,NPW
         DO M=1,NYMAX/NCPU
            MODEY=(MYID*NYMAX/NCPU+M-1)/2
            AY=PEY*MODEY
            DO L=1,NXMOD+2
               MODEX=(L-1)/2
               AX=PEX*MODEX
               AN=(AX**2+AY**2)**0.5
               WVN(L,M,K)=AN**K
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
      
C=====SUBROUTINE WVNUM END HERE








C======================================================================
      SUBROUTINE WSURF(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,WS,R,ZP,WVN,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE VERTICAL VELOCITY AT SURFACE

      IMPLICIT NONE
      
      INTEGER I,J,L,M,K,K1,K2
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU

      REAL PEX,PEY

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL R(NXMAX,NYMAX/NCPU,NPW),PK(NXMAX,NYMAX/NCPU,NPW)
      REAL WS(NXMAX,NYMAX/NCPU)
      REAL ZP(NXMAX,NYMAX/NCPU,NPW),WVN(NXMAX,NYMAX/NCPU,NPW)
      REAL TMP(NXMAX,NYMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C----------------------------------
C     SUM OF PERTURBATION MODE
C----------------------------------

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            PK(I,J,1)=R(I,J,1)
         ENDDO
      ENDDO

      DO K=2,NPW
         DO J=1,NYMAX/NCPU
            DO I=1,NXMOD
               PK(I,J,K)=PK(I,J,K-1)+R(I,J,K)
            ENDDO
         ENDDO
      ENDDO

C-----END HERE

C---------------------------------------
C     VERTICAL VELOCITY AT SURFACE
C---------------------------------------

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     INITIALIZE VERTICAL VELOCITY
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            WS(I,J)=0.
         ENDDO
      ENDDO

C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     SUM OF HIGH ORDER TERM
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO K=1,NPW-1

         K1=NPW-K
         
         CALL FFTXY_PARA(PK(1,1,K1),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1        NYMOD,NXMAX,NYMAX,NCPU,-1)

         DO M=1,NYMAX/NCPU
            DO L=1,NXMOD+2
               TMP(L,M)=PK(L,M,K1)*WVN(L,M,K+1)
            ENDDO
         ENDDO

         CALL FFTXY_PARA(TMP,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,
     1        NXMAX,NYMAX,NCPU,1)

         DO J=1,NYMAX/NCPU
            DO I=1,NXMOD
               WS(I,J)=WS(I,J)+ZP(I,J,K)*TMP(I,J)
            ENDDO
         ENDDO

      ENDDO
      
C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     ADD LEADING ORDER TERM
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL FFTXY_PARA(PK(1,1,NPW),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1     NYMOD,NXMAX,NYMAX,NCPU,-1)
      
      DO M=1,NYMAX/NCPU
         DO L=1,NXMOD+2
            TMP(L,M)=PK(L,M,NPW)*WVN(L,M,1)
         ENDDO
      ENDDO
      
      CALL FFTXY_PARA(TMP,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,
     1     NXMAX,NYMAX,NCPU,1)
      
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            WS(I,J)=WS(I,J)+TMP(I,J)
         ENDDO
      ENDDO

C~~~~~END HERE

      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,WS,NCPU)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE WSURF END HERE






C==========================================================================
      SUBROUTINE DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)

C     BY DI YANG, 08/2007
      
C     CALCULATE HORIZONTAL DERIVETIVE OF ETA AND VPS

      IMPLICIT NONE
      
      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU

      REAL PEX,PEY
      
      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),EX(NXMAX,*),EY(NXMAX,*)
      REAL VPS(NXMAX,*),VPSX(NXMAX,*),VPSY(NXMAX,*)
      REAL FTMP(NXMAX,NYMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      CALL PDFX_PARA(ETA,EX,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,NXMAX,
     1     NYMAX,NCPU)
      CALL PDFY_PARA(ETA,EY,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX,NYMAX,NCPU)
      CALL PDFX_PARA(VPS,VPSX,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,NXMAX,
     1     NYMAX,NCPU)
      CALL PDFY_PARA(VPS,VPSY,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,
     1     NXMAX,NYMAX,NCPU)

      RETURN
      END

C=====SUBROUTINE DERIVH END HERE






C==========================================================================
      SUBROUTINE UVSURF(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,US,VS,WS,EX,EY,VPSX,VPSY,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE HORIZONTAL VELOCITY AT SURFACE

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      
      REAL PEX,PEY
      
      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR
      
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)
      REAL EX(NXMAX,*),EY(NXMAX,*)
      REAL VPSX(NXMAX,*),VPSY(NXMAX,*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            US(I,J)=VPSX(I,J)-WS(I,J)*EX(I,J)
            VS(I,J)=VPSY(I,J)-WS(I,J)*EY(I,J)
         ENDDO
      ENDDO

      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,US,NCPU)
      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,VS,NCPU)

      RETURN
      END
      
C=====SUBROUTINE UVSURF END HERE





C======================================================================
      SUBROUTINE RIGH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ETA,EX,EY,VPSX,VPSY,WS,FETA,FVPS,PA,NRK4,FR2,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE RIGHT HAND SIDE OF BOUNDARY EVOLUTION EQUATION
C     R1 FOR KBC
C     R2 FOR DBC
C     PA PRESURE OF AIR AT SURFACE

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NRK4,NCPU
      
      REAL PEX,PEY
      
      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL FETA(NXMAX,NYMAX/NCPU,4),FVPS(NXMAX,NYMAX/NCPU,4)
      REAL ETA(NXMAX,*),EX(NXMAX,*),EY(NXMAX,*),WS(NXMAX,*)
      REAL VPSX(NXMAX,*),VPSY(NXMAX,*)
      REAL PA(NXMAX,*)
      REAL T1(NXMAX,NYMAX/NCPU),T2(NXMAX,NYMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      REAL FR2

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C---------------------------------
C     CHECK RUNGE-KUTTA STEP
C---------------------------------

      IF(NRK4.LE.0.OR.NRK4.GT.4) THEN
         PRINT*, "Invalid value for IRK in RIGH !"
         PRINT*, "IRK=",NRK4
         STOP
      ENDIF

C-----END HERE

C--------------------------------
C     CALCULATE 1+EX^2+EY^2
C--------------------------------

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            T1(I,J)=1.+EX(I,J)**2+EY(I,J)**2
         ENDDO
      ENDDO
      
      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,T1,NCPU)

C-----END HERE

C----------------------
C     CALCULATE R1
C----------------------

      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            FETA(I,J,NRK4)=-(VPSX(I,J)*EX(I,J)+VPSY(I,J)*EY(I,J))
     1           +T1(I,J)*WS(I,J)
         ENDDO
      ENDDO

      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,FETA(1,1,NRK4),NCPU)

C-----END HERE

C-------------------------
C     CALCULATE WS^2
C-------------------------
      
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            T2(I,J)=WS(I,J)**2
         ENDDO
      ENDDO

      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,T2,NCPU)

C-----END HERE

C----------------------
C     CALCULATE R2
C----------------------
      
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            FVPS(I,J,NRK4)=-ETA(I,J)/FR2-0.5*(VPSX(I,J)**2+VPSY(I,J)**2) 
     1           +0.5*T1(I,J)*T2(I,J)-PA(I,J)
         ENDDO
      ENDDO

      CALL DEALIASXY_PARA(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
     1     TRIGSY,IFAX,FVPS(1,1,NRK4),NCPU)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE RIGH END HERE
