      PROGRAM HS_TEST

C     BY DI YANG, 08/2007

      IMPLICIT NONE

      INTEGER I,J,K,IT,NTIME
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU,NSWAVEX,NSWAVEY
      INTEGER NWORK,NTRIGS
      INTEGER IOUT,NOUT
      INTEGER IBLOWUP
      INTEGER ISTART
      INTEGER IVEL
      INTEGER NTP
      INTEGER ISWELL
      
      REAL PEX,PEY,TIME,DT,XL,X,TWOPI,FR21,AKAX,AKAY,ENER0,ENER
      REAL P0
      REAL RE, FR2
      REAL LAMBDA_P, USTAR, OMEGA_P

      PARAMETER (NXMAX=260,NYMAX=260,NPW=3)
      PARAMETER (NXMOD=256,NYMOD=256)
      PARAMETER (NWORK=2*NXMAX*NYMAX)
      PARAMETER (NTRIGS=4*MAX(NXMAX,NYMAX))
      PARAMETER (NTIME=20)
      PARAMETER (NOUT=5) ! MAKE OUTPUT EVERY NOUT TIME STEPS
CC      PARAMETER (DT=0.05)
      PARAMETER (NTP=100) ! NUMBER OF TIME STEPS PER PEAK WAVE PERIOD
      PARAMETER (PEX=1.0,PEY=1.0)
CC      PARAMETER (FR21=0.013)
CC      PARAMETER (FR2=0.013)
      PARAMETER (AKAX=0.1,NSWAVEX=25)
C     AKAX: WAVE SLOPE IN X-DIRECTION.  ONLY USED FOR SINUOSOIDAL WAVES
C     NSWAVEX: NUMBER OF PEAK WAVES IN X-DIRECTION.
      PARAMETER (AKAY=0.1,NSWAVEY=6)
C     AKAY: WAVE SLOPE FOR SWELL
C     NSWAVEY: NUMBER OF SWELLS IN X-DIRECTION
      PARAMETER (NCPU=10)
      PARAMETER (P0=0.001)
      PARAMETER (RE=1000)

      PARAMETER (ISTART=0)
C     ISTART=0: START WAVE SIMULATION FROM INITIAL CONDITION
C     ISTART=1: READ IN RESTART FILES FROM fort.160*

      PARAMETER (ISWELL=1)
C     ISWELL: =1, ADD SWELL; =0, NO SWELL

C---------------------------------------------------
C     VARIABLES FOR LEVEL-SET INFLOW CONDITION
C---------------------------------------------------

      PARAMETER (IVEL=1)
C     IF IVEL=0, ONLY PERFORM HOS SIMULATION OF WAVE
C     IF IVEL=1, ALSO CALCULATE VELOCITY BELOW THE WATER SURFACE

      INTEGER NZMAXL,NZL,NZMAXI,NZI
      PARAMETER (NZMAXL=130,NZL=129)
      PARAMETER (NZMAXI=194,NZI=193)
      REAL DEPTH
      REAL DZL(NZMAXI),ZLS(NZMAXI)

C-----END HERE

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)
      
      REAL ETA(NXMAX,NYMAX/NCPU),FETA(NXMAX,NYMAX/NCPU,4)
      REAL EX(NXMAX,NYMAX/NCPU),EY(NXMAX,NYMAX/NCPU)
      REAL VPS(NXMAX,NYMAX/NCPU),FVPS(NXMAX,NYMAX/NCPU,4)
      REAL VPSX(NXMAX,NYMAX/NCPU),VPSY(NXMAX,NYMAX/NCPU)
      REAL US(NXMAX,NYMAX/NCPU),VS(NXMAX,NYMAX/NCPU)
      REAL WS(NXMAX,NYMAX/NCPU),PA(NXMAX,NYMAX/NCPU)
      REAL WVN(NXMAX,NYMAX/NCPU,NPW)
      REAL R(NXMAX,NYMAX/NCPU,NPW),ZP(NXMAX,NYMAX/NCPU,NPW)
      REAL ENALL(NXMAX,NYMAX/NCPU), ENALL0(NXMAX,NYMAX/NCPU)
      REAL ID(NXMAX,NYMAX/NCPU)

      REAL WORK(NWORK)
      REAL TRIGSX(NTRIGS),TRIGSY(NTRIGS),TRIGSZ(NTRIGS)
      INTEGER IFAX(19)

      REAL TMP
      REAL PERIOD, PHAV, SIGMA, BETA, T0
      REAL XP(NXMAX)
      REAL PI
      REAL COEF
      REAL NTMP
      REAL ETMP(NXMAX,NYMAX/NCPU)

      REAL ALPHA, OMEGA0, NU, TA, NN
      REAL F, U10, G, USS
      INTEGER FID

      INTEGER MOVENX, MOVENY
CC
      REAL T1(NXMAX,NYMAX/NCPU),SPETA(NXMAX,NYMAX/NCPU)
      REAL SPX(NXMAX)
      REAL WAVENUMBX,DWK,WAVENUMBY
      INTEGER N,NMAX,M

      INTEGER ICON
      DATA ICON/16/
CC

      CALL MPI_INIT(IERR)     
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      
      PI = ACOS(-1.)

C-------------------------------------------------------------------------
C     INITIALIZE VERTICAL GRIDS FOR CALCULATING VELOCITY FOR LEVEL-SET
C-------------------------------------------------------------------------

      IF(IVEL.EQ.1) THEN

         TWOPI=2.*ACOS(-1.)
         XL=TWOPI/PEX
         
         DEPTH=XL/2.
         DO K=1,NZI
            DZL(K)=DEPTH/FLOAT(NZL-1)
         ENDDO
         DO K=1,NZI
            ZLS(K)=(K-1)*DZL(1)-DEPTH
         ENDDO

      ENDIF

C-----END HERE

C------------------------------------------
C     2D LINEAR WAVE INITIAL CONDITION
C------------------------------------------

      USTAR=0.01343
C     USTAR IS THE DIMENSIONALESS WIND FRICTION VELOCITY.
C     IT IS THE VALUE FROM THE LES.
C     

      U10 = 12.30
C     U10 IS THE DIMENSIONAL WIND VELOCITY AT 10 METER HEIGHT.

      USS=0.45
C     USS IS THE DIMENSIONL WIND FRICTION VELOCITY.

      F = 80000.
C     F IS THE DIMENSIONAL FETCH VALUE.

      G = 9.8
C     G IS THE GRAVITATIONAL ACCELERATION.

      ALPHA = 0.076 * ( U10**2 / F / G )**0.22
C     SPECTRAL COEFFICIENT (PHILLIPS CONSTANT)

      OMEGA_P = 22.*(G**2/U10/F)**(1./3.)
C     OMEGA_P IS THE DIMENSIONAL ANGULAR FREQUENCY AT SPECTRUM PEAK.

      LAMBDA_P = 2*PI*G/OMEGA_P**2
C     LAMBDA_P IS THE DIMENSIONAL WAVELENGTH AT SPECTRUM PEAK.

      FR21 = (USS/USTAR)**2/G/(LAMBDA_P/(2.*PI/PEX/NSWAVEX))
C     FR21 IS THE SQUARE OF THE FROUDE NUMBER.
C     NSWAVEX IS THE NUMBER OF PEAK WAVES IN THE SIMULATION DOMAIN.
C     2*PI/PEX IS THE STREAMWISE SIMULATION DOMAIN SIZE
C     2*PI/PEY IS THE SPANWISE SIMULATION DOMAIN SIZE

      FR2 = FR21

      PRINT*, LAMBDA_P, FR2

      OMEGA0 = ( NSWAVEX * PEX / FR21 )**0.5
C     OMEGA0 IS THE DIMENSIONLESS ANGULAR FREQUENCY AT SPECTRUM PEAK

      PRINT*, 'C/U*=', OMEGA0/(NSWAVEX*PEX)/USTAR

      NU = 3.3
CC      NU=1.0

      WRITE(*,*)'start-code'

C-------------------------------------------------
C     READIN DATA TO CONTINUE THE SIMULATION
C-------------------------------------------------
      
      IF(ISTART.EQ.1) THEN
CC         CALL REREADN_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,PA,ICON,
CC     1        TIME,FR2,NCPU)
         CALL REREAD_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,PA,ICON,
     1        TIME,FR2,NCPU)
         GOTO 500
      ENDIF

      CALL JONSWAP_3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,FR21,ALPHA,OMEGA0,
     1     NU,ETA,VPS,NCPU,WORK,TRIGSX,TRIGSY,IFAX,FR2)

C      stop

CC      DWK = PEX
CC      CALL SPECTK3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,DWK,
CC     1           WORK,TRIGSX,TRIGSY,IFAX,T1,ETA,SPETA)

CC      WRITE(75,*) " VARIABLES=X,Y,SPECT"
CC      WRITE(75,950)G,NXMOD/2-1,NYMOD/2-1
CC 950  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
CC      DO M=1,NYMOD/2-1
CC       DO N=1,NXMOD/2-1
CC        WAVENUMBX=(N-1)*PEX
CC        WAVENUMBY=(M-1)*PEY
CC        WRITE(75,510)WAVENUMBX,WAVENUMBY,SPETA(N,M)
CC       ENDDO
CC      ENDDO
c$$$ 510  FORMAT(3E12.4)
c$$$
c$$$      CALL SPECTKX(NXMOD,NYMOD,NXMAX,NYMAX,PEX,NMAX,DWK,
c$$$     1     WORK,TRIGSX,TRIGSY,IFAX,T1,ETA,SPX)
c$$$      
c$$$      WRITE(76,*) 'VARIABLES=K,SPX'
c$$$      DO M=1,NXMOD/2-1
c$$$         WAVENUMBX=(M-1)*PEX
c$$$         WRITE(76,510) WAVENUMBX,SPX(M)
c$$$      ENDDO
C
      
CC      COEF = 0
C-----END HERE

C------------------------------------------
C     2D LINEAR WAVE INITIAL CONDITION
C------------------------------------------

CC      CALL LINEAR_WAVE_I_3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,AKAX,
CC     1     AKAY,NSWAVEX,NSWAVEY,FR2,ETA,VPS,NCPU)

      IF(ISWELL.EQ.1) THEN

         CALL LINEAR_WAVE_I_2D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,AKAY,
     1        NSWAVEY,FR2,ETA,VPS,NCPU)

      ENDIF
      
C-----END HERE

C---------------------------------------
C     STOKES WAVE INITIAL CONDITION
C---------------------------------------
      
C      READ(11,*) TMP, PERIOD, PHAV
C      READ(11,*) TMP, TMP

C      DT = PERIOD / 100.

C      J = 1
C      DO I = 1, NXMOD
C         READ(11,*) XP(I), ETA(I,J), VPS(I,J)
CC         XP(I) = XP(I) - PI
C      END DO

C      DO J = 2, NYMAX / NCPU
C         DO I = 1, NXMOD
C            ETA(I,J) = ETA(I,1)
C            VPS(I,J) = VPS(I,1)
C         END DO
C      END DO

C      CALL FFTXY_PARA(ETA,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
C     1     NYMAX,NCPU,-1)

C      IF ( MYID .EQ. 0 ) THEN
C         ETA(1,1) = 0
C      END IF

C      CALL FFTXY_PARA(ETA,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
C     1     NYMAX,NCPU,1)

C      CALL SMOOTH3(NXMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,ETA,VPS,NCPU)
      
C--------------------------------
C     EMPTY FIELD
C-------------------------------
C      DO I = 1, NXMOD
C         DO J = 1, NYMAX / NCPU
C           ETA(I,J) = ETA(I,J)/10.
C           VPS(I,J) = VPS(I,J)/10.
C         END DO
C      END DO

CC      SIGMA = ( NSWAVEX * PEX / FR2 )**0.5
CC      BETA = ( RE * SIGMA / 2. )**0.5
CC      PERIOD = 2 * PI / SIGMA
CC      DT = PERIOD / 100.

CC--NONLINEAR TREATMENT
CC      TA = 2 * PERIOD
CC      NN = 2

CC      DO I = 1, NXMOD
CC         XP(I) = 2 * PI / PEX / NXMOD * ( I - 1. )
CC      END DO

C--------------------------------

C----------------------------------
C     NO AIR ABOVE THE SURFACE
C----------------------------------
      
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            PA(I,J)=0.
         ENDDO
      ENDDO
      
C-----END HERE

C--------------------------------
C     INITIALIZE WAVE NUMBER
C--------------------------------

      TIME=0.

 500  CONTINUE
     
      CALL WAVENUM(NXMOD,NYMOD,NXMAX,NYMAX,NPW,PEX,PEY,WVN,NCPU)

C-----END HERE

C---------------------------------------------
C     INITIAL STEP FOR RUNGE-KUTTA SCHEME
C---------------------------------------------
      
      CALL DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)

      CALL ZETA(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,ETA,ZP,WVN,NCPU)      

      CALL BOUNDVP(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,VPS,R,ZP,WVN,NCPU)     

      CALL WSURF(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,WS,R,ZP,WVN,NCPU)     

      CALL UVSURF(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,US,VS,WS,EX,EY,VPSX,VPSY,NCPU)

C      TMP = 0.0
C      DO J=1,NYMOD
C         DO I=1,NXMOD
C            TMP=TMP+ETA(I,J)**2
C         ENDDO
C      ENDDO
C      TMP=TMP/(NXMOD*NYMOD)
C      WRITE(*,*)'TMP2=',TMP,2*(NSWAVEX*PEX)*TMP**0.5
C
C      CALL SMOOTH2(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
C     1     PEX,PEY,ETA,VPS,EX,EY,0.,NCPU)

C      IOUT=NOUT-1

C---------------------------
C     H.O.S SIMULATION
C---------------------------
C      TIME=0.

      SIGMA = ( NSWAVEX * PEX / FR2 )**0.5
C      BETA = ( RE * SIGMA / 2. )**0.5
      PERIOD = 2 * PI / SIGMA
      DT = PERIOD / FLOAT(NTP)
CC      TIME = TMP*PERIOD

      IOUT=0

      DO 100 IT=1,NTIME

         TIME=TIME+DT

         DO J = 1, NYMAX / NCPU
            DO I = 1, NXMOD
               PA(I,J) = 0.
            END DO
         END DO
C
         CALL HOS_WAVE_3D(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,TIME,DT,WVN,ETA,EX,EY,VPS,FETA,FVPS,US,VS,
     1        WS,PA,FR2,PERIOD,NCPU)

         TMP = TIME / PERIOD
CC         TMP = TIME
         CALL OUTSURF(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
     1        FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,IOUT,NOUT,TMP,
     1        NCPU)
         
         IF(IBLOWUP.EQ.1) THEN
            PRINT*, "Code blows up!"
            GOTO 1000
         ENDIF
         
 100  CONTINUE

C------------------------------------------------
C     CALCULATE VELOCITY BELOW WATER SURFACE
C------------------------------------------------
      
      IF(IVEL.EQ.1) THEN

         CALL DERIVH(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1        PEX,PEY,ETA,VPS,EX,EY,VPSX,VPSY,NCPU)
         
         CALL ZETA(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,ETA,ZP,WVN,NCPU)      
         
         CALL BOUNDVP(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,VPS,R,ZP,WVN,NCPU)     
         
         CALL WSURF(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,TRIGSX,
     1        TRIGSY,PEX,PEY,WS,R,ZP,WVN,NCPU)     
         
         CALL UVSURF(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1        PEX,PEY,US,VS,WS,EX,EY,VPSX,VPSY,NCPU)

         CALL VELOCITY_HOS(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,
     1        TRIGSX,TRIGSY,PEX,PEY,R,ETA,WVN,TIME,NZMAXI,NZI,
     1        NZL,DEPTH,DZL,ZLS,US,VS,WS,NSWAVEX,NCPU)

      ENDIF

C-----END HERE


CC      CALL SAVE_A(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
CC     1           FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,IOUT,NOUT,TMP,
CC     1           NCPU)
      CALL SAVE_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,PA,ICON,
     1     TIME,FR2,NCPU)
CC      CALL SAVE_LES_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETA,VPS,PA,ICON,
CC     1     TIME,FR2,NCPU)

 1000 CONTINUE

      CALL MPI_FINALIZE(IERR)
      STOP
      END

C=====PROGRAM HOS_TEST END HERE






C======================================================================
      SUBROUTINE VELOCITY_HOS(NXMOD,NYMOD,NXMAX,NYMAX,NPW,WORK,IFAX,
     1     TRIGSX,TRIGSY,PEX,PEY,R,ETA,WVN,TIME,NZMAXI,NZI,NZL,DEPTH,
     1     DZL,ZLS,US,VS,WS,NSWAVEX,NCPU)

C     BY DI YANG, 05/2008

C     CALCULATE VELOCITY BELOW WAVE SURFACE
C     OUTPUT FLOW FIELD FOR LEVEL SET

      IMPLICIT NONE
      
      INTEGER I,J,L,M,K,K1,K2,II
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NPW,NCPU,NSWAVEX

      REAL PEX,PEY,PMODX,PMODY,TIME

C+++
C     FOR INFLOW CONDITION OF LEVEL SET + IBM

      INTEGER NCREQ,JJ1,JJ2,NYEND

      INTEGER NZMAXI,NZI,KZ,NZS,NZL
C     FROM BOTTOM TO MEAN SURFACE LEVEL: 1~NZL
C     FROM MEAN SURFACE LEVEL: NZL~NZS

      REAL ZLS(*),DZL(*),DEPTH,EMAX,EMAX_ALL,ZG,Z1,Z2,XL,YL,HBAR
      REAL TWOPI,X,Y
C+++     

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL R(NXMAX,NYMAX/NCPU,NPW),PK(NXMAX,NYMAX/NCPU,NPW)
      REAL UL(NXMAX,NYMAX/NCPU,NZMAXI)
      REAL VL(NXMAX,NYMAX/NCPU,NZMAXI)
      REAL WL(NXMAX,NYMAX/NCPU,NZMAXI)
      REAL ZP(NPW),WVN(NXMAX,NYMAX/NCPU,NPW),ETA(NXMAX,NYMAX/NCPU)
      REAL T1(NXMAX,NYMAX/NCPU)
      REAL T2(NXMAX,NYMAX/NCPU)
      REAL T3(NXMAX,NYMAX/NCPU)
      REAL TMP(NXMAX,NYMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      REAL UI(NXMOD,NYMAX/NCPU,NZMAXI)
      REAL VI(NXMOD,NYMAX/NCPU,NZMAXI)
      REAL WI(NXMOD,NYMAX/NCPU,NZMAXI)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)
      REAL PHII(NXMOD,NYMAX/NCPU,NZMAXI)

      REAL AT(NXMAX,NYMAX/NCPU)
      REAL B(NYMAX,NXMAX/NCPU),BT(NYMAX,NXMAX/NCPU)
      REAL BY(NYMAX,NXMAX/NCPU)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NCREQ IS THE # OF CPUS THAT CONTAIN USEFUL ELEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX

      JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
      JJ2=NYMAX/NCPU

      IF(MYID.EQ.NCREQ-1) THEN
         NYEND=JJ1
      ELSE IF(MYID.LT.NCREQ-1) THEN
         NYEND=JJ2
      ELSE
         NYEND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C---------------------------------
C     INITIALIZE PARAMETERS
C---------------------------------

      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      YL=TWOPI/PEY

      EMAX=0.
      
      DO I=1,NXMOD
         DO J=1,NYEND
            IF(ETA(I,J).GT.EMAX) EMAX=ETA(I,J)
         ENDDO
      ENDDO

      CALL MPI_ALLREDUCE(EMAX,EMAX_ALL,1,MPI_DOUBLE_PRECISION,
     1     MPI_MAX,MPI_COMM_WORLD,IERR)

      DO K=1,NZI
         IF(ZLS(K).GT.EMAX_ALL) THEN
            NZS=K-1
            GOTO 50
         ENDIF
      ENDDO
 50   CONTINUE

      IF(NZS.GT.NZI) THEN
         PRINT*, "INVALID NZI"
         PRINT*, "NZS=",NZS
         PRINT*, "NZI=",NZI
         STOP
      ENDIF

C-----END HERE

C----------------------------------
C     SUM OF PERTURBATION MODE
C----------------------------------

C      DO J=1,NYMAX/NCPU
C         DO I=1,NXMOD
C            PK(I,J,1)=R(I,J,1)
C         ENDDO
C      ENDDO

      DO K=1,NPW
         DO J=1,NYMAX/NCPU
            DO I=1,NXMOD
               PK(I,J,K)=R(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,NPW
         CALL FFTXY_PARA(PK(1,1,K),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1        NYMOD,NXMAX,NYMAX,NCPU,-1)
      ENDDO

C-----END HERE

C---------------------------------------
C     VELOCITY BELOW WAVE SURFACE
C---------------------------------------

      DO 100 KZ=1,NZMAXI

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     INITIALIZE VELOCITY
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         DO J=1,NYMAX/NCPU
            DO I=1,NXMOD
               UL(I,J,KZ)=0.
               VL(I,J,KZ)=0.
               WL(I,J,KZ)=0.
            ENDDO
         ENDDO

C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     INITIALIZE TEMPORAL VARIABLES
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         DO K=1,NPW
            ZP(K)=0.
         ENDDO

C~~~~~END HERE
               
C***********************************
C     BELOW MEAN SURFACE LEVEL
C***********************************

         IF (KZ.LE.NZL) THEN

            DO K=1,NPW

               DO I=1,NXMAX
                  DO J=1,NYMAX/NCPU
                     T1(I,J)=0.
                     T2(I,J)=0.
                     T3(I,J)=0.
                  ENDDO
               ENDDO

               DO L=1,NXMOD+1,2
                  PMODX=PEX*(L-1)/2
                  DO M=1,NYMAX/NCPU
                     T1(L,M)=-PMODX*PK(L+1,M,K)*EXP(WVN(L,M,1)*ZLS(KZ))
                     T1(L+1,M)=PMODX*PK(L,M,K)*EXP(WVN(L+1,M,1)*ZLS(KZ))
                  ENDDO
               ENDDO
               
               DO I=1,NXMOD
                  DO J=1,NYMAX/NCPU
                     TMP(I,J)=PK(I,J,K)
                  ENDDO
               ENDDO
               
               CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
            
               DO L=1,NYMOD+1,2
                  PMODY=PEY*(L-1)/2
                  DO M=1,NXMAX/NCPU
                     BY(L,M)=-PMODY*B(L+1,M)
                     BY(L+1,M)=PMODY*B(L,M)
                  ENDDO
               ENDDO
            
               CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
               
               DO L=1,NXMOD
                  DO M=1,NYMAX/NCPU
                     T2(L,M)=T2(L,M)*EXP(WVN(L,M,1)*ZLS(KZ))
CC                     T2(L,M)=PK(L,M,K)*EXP(WVN(L,M,1)*ZLS(KZ))
                  ENDDO
               ENDDO
               
               DO M=1,NYMAX/NCPU
                  DO L=1,NXMOD+2
                     T3(L,M)=PK(L,M,K)*WVN(L,M,1)*EXP(WVN(L,M,1)
     1                    *ZLS(KZ))
                  ENDDO
               ENDDO

               CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1              NYMOD,NXMAX,NYMAX,NCPU,1)
               CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1              NYMOD,NXMAX,NYMAX,NCPU,1)
               CALL FFTXY_PARA(T3,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
     1              NYMOD,NXMAX,NYMAX,NCPU,1)
               
               DO J=1,NYMAX/NCPU
                  DO I=1,NXMOD
                     UL(I,J,KZ)=UL(I,J,KZ)+T1(I,J)
                     VL(I,J,KZ)=VL(I,J,KZ)+T2(I,J)
                     WL(I,J,KZ)=WL(I,J,KZ)+T3(I,J)
                  ENDDO
               ENDDO

            ENDDO


c$$$            DO K=1,NPW-1
c$$$               
c$$$               K1=NPW-K
c$$$               
c$$$               DO L=1,NXMOD+1,2
c$$$                  PMODX=PEX*(L-1)/2
c$$$                  DO M=1,NYMAX/NCPU
c$$$                     T1(L,M)=-PMODX*PK(L+1,M,K1)*WVN(L,M,K)
c$$$     1                    *EXP(WVN(L,M,1)*ZLS(KZ))
c$$$                     T1(L+1,M)=PMODX*PK(L,M,K1)*WVN(L+1,M,K)
c$$$     1                    *EXP(WVN(L+1,M,1)*ZLS(KZ))
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               DO I=1,NXMOD
c$$$                  DO J=1,NYMAX/NCPU
c$$$                     TMP(I,J)=PK(I,J,K1)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
c$$$               
c$$$               DO L=1,NYMOD+1,2
c$$$                  PMODY=PEY*(L-1)/2
c$$$                  DO M=1,NXMAX/NCPU
c$$$                     BY(L,M)=-PMODY*B(L+1,M)
c$$$                     BY(L+1,M)=PMODY*B(L,M)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
c$$$               
c$$$               DO L=1,NXMOD
c$$$                  DO M=1,NYMAX/NCPU
c$$$                     T2(L,M)=T2(L,M)*WVN(L,M,K)*EXP(WVN(L,M,1)*ZLS(KZ))
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               DO M=1,NYMAX/NCPU
c$$$                  DO L=1,NXMOD+2
c$$$                     T3(L,M)=PK(L,M,K1)*WVN(L,M,K+1)*EXP(WVN(L,M,1)
c$$$     1                    *ZLS(KZ))
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               CALL FFTXY_PARA(T3,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               
c$$$               DO J=1,NYMAX/NCPU
c$$$                  DO I=1,NXMOD
c$$$                     UL(I,J,KZ)=UL(I,J,KZ)+ZP(K)*T1(I,J)
c$$$                     VL(I,J,KZ)=VL(I,J,KZ)+ZP(K)*T2(I,J)
c$$$                     WL(I,J,KZ)=WL(I,J,KZ)+ZP(K)*T3(I,J)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$            ENDDO

C     ADD LEADING ORDER TERM

c$$$            DO L=1,NXMOD+1,2
c$$$               PMODX=PEX*(L-1)/2
c$$$               DO M=1,NYMAX/NCPU
c$$$                  T1(L,M)=-PMODX*PK(L+1,M,NPW)*EXP(WVN(L,M,1)*ZLS(KZ))
c$$$                  T1(L+1,M)=PMODX*PK(L,M,NPW)*EXP(WVN(L+1,M,1)*ZLS(KZ))
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            DO I=1,NXMOD
c$$$               DO J=1,NYMAX/NCPU
c$$$                  TMP(I,J)=PK(I,J,NPW)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
c$$$           
c$$$            DO L=1,NYMOD+1,2
c$$$               PMODY=PEY*(L-1)/2
c$$$               DO M=1,NXMAX/NCPU
c$$$                  BY(L,M)=-PMODY*B(L+1,M)
c$$$                  BY(L+1,M)=PMODY*B(L,M)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
c$$$            
c$$$            DO M=1,NYMAX/NCPU
c$$$               DO L=1,NXMOD+2
c$$$                  T2(L,M)=T2(L,M)*EXP(WVN(L,M,1)*ZLS(KZ))
c$$$                  T3(L,M)=PK(L,M,NPW)*WVN(L,M,1)*EXP(WVN(L,M,1)*ZLS(KZ))
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1           NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$            CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1           NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$            CALL FFTXY_PARA(T3,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,
c$$$     1           NXMAX,NYMAX,NCPU,1)
c$$$            
c$$$            DO J=1,NYMAX/NCPU
c$$$               DO I=1,NXMOD
c$$$                  UL(I,J,KZ)=UL(I,J,KZ)+T1(I,J)
c$$$                  VL(I,J,KZ)=VL(I,J,KZ)+T2(I,J)
c$$$                  WL(I,J,KZ)=WL(I,J,KZ)+T3(I,J)
c$$$               ENDDO
c$$$            ENDDO


ccccccccccccccc

c$$$            DO L=1,NXMOD+1,2
c$$$               PMODX=PEX*(L-1)/2
c$$$               DO M=1,NYMAX/NCPU
c$$$                  UL(L,M,KZ)=-PMODX*PK(L+1,M,NPW)*EXP(WVN(L,M,1)
c$$$     1                 *ZLS(KZ))
c$$$                  UL(L+1,M,KZ)=PMODX*PK(L,M,NPW)*EXP(WVN(L+1,M,1)
c$$$     1                 *ZLS(KZ))
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            DO I=1,NXMOD
c$$$               DO J=1,NYMAX/NCPU
c$$$                  TMP(I,J)=PK(I,J,NPW)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
c$$$            
c$$$            DO L=1,NYMOD+1,2
c$$$               PMODY=PEY*(L-1)/2
c$$$               DO M=1,NXMAX/NCPU
c$$$                  BY(L,M)=-PMODY*B(L+1,M)
c$$$                  BY(L+1,M)=PMODY*B(L,M)
c$$$               ENDDO
c$$$            ENDDO
c$$$
c$$$            CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
c$$$            
c$$$            DO L=1,NXMOD
c$$$               DO M=1,NYMAX/NCPU
c$$$                  VL(L,M,KZ)=T2(L,M)*EXP(WVN(L,M,1)*ZLS(KZ))
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            DO M=1,NYMAX/NCPU
c$$$               DO L=1,NXMOD+2
c$$$                  WL(L,M,KZ)=PK(L,M,NPW)*WVN(L,M,1)*EXP(WVN(L,M,1)
c$$$     1                 *ZLS(KZ))
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL FFTXY_PARA(UL(1,1,KZ),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1           NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$            CALL FFTXY_PARA(VL(1,1,KZ),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1           NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$            CALL FFTXY_PARA(WL(1,1,KZ),WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1           NYMOD,NXMAX,NYMAX,NCPU,1)

C*****END HERE

C***********************************************
C     BETWEEN MEAN SURFACE LEVEL AND CREST
C***********************************************

         ELSE IF(KZ.GT.NZL.AND.KZ.LE.NZS) THEN
            
c$$$            ZP(1)=ZLS(KZ)
c$$$            DO K=2,NPW-1
c$$$               ZP(K)=ZP(K-1)*ZP(1)/FLOAT(K)
c$$$            ENDDO

            DO J=1,NYEND
               DO I=1,NXMOD
                  IF(ZLS(KZ).LE.ETA(I,J)) THEN
                     Z1=ETA(I,J)-ZLS(KZ)
                     Z2=ZLS(KZ)-ZLS(NZL)
                     UL(I,J,KZ)=(US(I,J)*Z2+UL(I,J,NZL)*Z1)/(Z1+Z2)
                     VL(I,J,KZ)=(VS(I,J)*Z2+VL(I,J,NZL)*Z1)/(Z1+Z2)
                     WL(I,J,KZ)=(WS(I,J)*Z2+WL(I,J,NZL)*Z1)/(Z1+Z2)
                  ENDIF
               ENDDO
            ENDDO

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     SUM OF HIGH ORDER TERM
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c$$$            DO K=1,NPW-1
c$$$               
c$$$               K1=NPW-K
c$$$               
c$$$               DO L=1,NXMOD+1,2
c$$$                  PMODX=PEX*(L-1)/2
c$$$                  DO M=1,NYMAX/NCPU
c$$$                     T1(L,M)=-PMODX*PK(L+1,M,K1)*WVN(L,M,K)
c$$$                     T1(L+1,M)=PMODX*PK(L,M,K1)*WVN(L+1,M,K)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               DO I=1,NXMOD
c$$$                  DO J=1,NYMAX/NCPU
c$$$                     TMP(I,J)=PK(I,J,K1)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
c$$$               
c$$$               DO L=1,NYMOD+1,2
c$$$                  PMODY=PEY*(L-1)/2
c$$$                  DO M=1,NXMAX/NCPU
c$$$                     BY(L,M)=-PMODY*B(L+1,M)
c$$$                     BY(L+1,M)=PMODY*B(L,M)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
c$$$               
c$$$               DO L=1,NXMOD
c$$$                  DO M=1,NYMAX/NCPU
c$$$                     T2(L,M)=T2(L,M)*WVN(L,M,K)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               DO M=1,NYMAX/NCPU
c$$$                  DO L=1,NXMOD+2
c$$$                     T3(L,M)=PK(L,M,K1)*WVN(L,M,K+1)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$               CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               CALL FFTXY_PARA(T3,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,
c$$$     1              NYMOD,NXMAX,NYMAX,NCPU,1)
c$$$               
c$$$               DO J=1,NYMAX/NCPU
c$$$                  DO I=1,NXMOD
c$$$                     UL(I,J,KZ)=UL(I,J,KZ)+ZP(K)*T1(I,J)
c$$$                     VL(I,J,KZ)=VL(I,J,KZ)+ZP(K)*T2(I,J)
c$$$                     WL(I,J,KZ)=WL(I,J,KZ)+ZP(K)*T3(I,J)
c$$$                  ENDDO
c$$$               ENDDO
c$$$               
c$$$            ENDDO
            
C~~~~~END HERE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     ADD LEADING ORDER TERM
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c$$$            DO L=1,NXMOD+1,2
c$$$               PMODX=PEX*(L-1)/2
c$$$               DO M=1,NYMAX/NCPU
c$$$                  T1(L,M)=-PMODX*PK(L+1,M,NPW)
c$$$                  T1(L+1,M)=PMODX*PK(L,M,NPW)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            DO I=1,NXMOD
c$$$               DO J=1,NYMAX/NCPU
c$$$                  TMP(I,J)=PK(I,J,NPW)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL MPE_2DTRANS(TMP,AT,B,NXMAX,NYMAX/NCPU,NCPU)
c$$$           
c$$$            DO L=1,NYMOD+1,2
c$$$               PMODY=PEY*(L-1)/2
c$$$               DO M=1,NXMAX/NCPU
c$$$                  BY(L,M)=-PMODY*B(L+1,M)
c$$$                  BY(L+1,M)=PMODY*B(L,M)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL MPE_2DTRANS(BY,BT,T2,NYMAX,NXMAX/NCPU,NCPU)
c$$$            
c$$$            DO M=1,NYMAX/NCPU
c$$$               DO L=1,NXMOD+2
c$$$                  T3(L,M)=PK(L,M,NPW)*WVN(L,M,1)
c$$$               ENDDO
c$$$            ENDDO
c$$$            
c$$$            CALL FFTXY_PARA(T3,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,
c$$$     1           NXMAX,NYMAX,NCPU,1)
c$$$            
c$$$            DO J=1,NYMAX/NCPU
c$$$               DO I=1,NXMOD
c$$$                  UL(I,J,KZ)=UL(I,J,KZ)+T1(I,J)
c$$$                  VL(I,J,KZ)=VL(I,J,KZ)+T2(I,J)
c$$$                  WL(I,J,KZ)=WL(I,J,KZ)+T3(I,J)
c$$$               ENDDO
c$$$            ENDDO
            
C~~~~~END HERE

C*****END HERE
               
         ENDIF
         
C         CALL DEALIASXY_PARA2(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
C     1        TRIGSY,IFAX,UL(1,1,KZ),NCPU)
C         CALL DEALIASXY_PARA2(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
C     1        TRIGSY,IFAX,VL(1,1,KZ),NCPU)
C         CALL DEALIASXY_PARA2(NXMOD,NYMOD,1,NXMAX,NYMAX,WORK,TRIGSX,
C     1        TRIGSY,IFAX,WL(1,1,KZ),NCPU)

 100  CONTINUE

C-----END HERE

C---------------------------------------------------------
C     INTERPOLATE VELOCITY TO CARTESIAN COORDINATES
C---------------------------------------------------------

      do i=1,nxmod
         ii=i
         DO J=1,NYEND
            DO K=1,NZI
               PHII(II,J,K)=ETA(I,J)-ZLS(K)
c$$$               IF(ZLS(K).LE.0..AND.ZLS(K).LT.ETA(I,J)) THEN
c$$$                  UI(II,J,K)=UL(I,J,K)
c$$$                  VI(II,J,K)=VL(I,J,K)
c$$$                  WI(II,J,K)=WL(I,J,K)
c$$$               ELSE IF(ZLS(K).GT.0.AND.ZLS(K).LT.ETA(I,J)) THEN
c$$$                  Z1=ZLS(K)
c$$$                  Z2=ETA(I,J)-ZLS(K)
c$$$                  UI(II,J,K)=(UL(I,J,NZL)*Z2+US(I,J)*Z1)/ETA(I,J)
c$$$                  VI(II,J,K)=(VL(I,J,NZL)*Z2+VS(I,J)*Z1)/ETA(I,J)
c$$$                  WI(II,J,K)=(WL(I,J,NZL)*Z2+WS(I,J)*Z1)/ETA(I,J)
               IF(ZLS(K).LT.ETA(I,J)) THEN
                  UI(II,J,K)=UL(I,J,K)
                  VI(II,J,K)=VL(I,J,K)
                  WI(II,J,K)=WL(I,J,K)
               ELSE
                  Z1=ETA(I,J)-ZLS(K)
                  UI(II,J,K)=US(I,J)*EXP(Z1*PEX*NSWAVEX)
                  VI(II,J,K)=VS(I,J)*EXP(Z1*PEX*NSWAVEX)
                  WI(II,J,K)=WS(I,J)*EXP(Z1*PEX*NSWAVEX)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C-----END HERE

C-------------------------------------------------
C     OUTPUT INITIAL CONDITION FOR LEVEL SET
C-------------------------------------------------

      OPEN(UNIT=5000+MYID)
      WRITE(5000+MYID,*) 'VARIABLES=X,Y,Z,U,V,W,PHI'
      WRITE(5000+MYID,*) 'ZONE T="INIT" I=',NXMOD,' J=',NYEND,
     &     ' K=',NZI,' F=POINT'

      DO K=1,NZI
         DO J=1,NYEND
            DO I=1,NXMOD
               X=(I-1)*XL/NXMOD
               y=(j+myid*nymax/ncpu-1)*yl/float(nymod)-yl/2.
               WRITE(5000+MYID,*) X,y,zls(k),UI(I,J,K),
     &              VI(I,J,K),WI(I,J,K),PHII(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      PRINT *,MYID,'IC FOR LEVELSET OUTPUT FINISHED!'

C-----END HERE

      RETURN
      END

C=====SUBROUTINE VELOCITY_HOS END HERE






C--------------------------------------------------------
      SUBROUTINE MOVE_X(NXMAX,NYMAX,NXMOD,NYMOD,F,MOVENUMBER)
CC--BY XIN GUO-----------------------------------------CC
C
C MOVE THE WHOLE DATA TO X DIRECTION WITH MOVENUMBER STEPS
C BY USING THE PERIODIC BOUNDARY CONDITION
C
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@--CC
C
C--INPUT
      INTEGER NXMAX, NYMAX
      INTEGER NXMOD, NYMOD
      INTEGER MOVENUMBER
C--INPUT/OUTPUT
      REAL F(NXMAX,*)
C--INTERNAL
      INTEGER I, J, K
      REAL TMP
C
      DO K = 1, MOVENUMBER
         DO J = 1, NYMOD
            DO I = 1, NXMOD - 1
               TMP = F(I,J)
               F(I,J) = F(I+1,J)
               F(I+1,J) = TMP
            END DO
         END DO
      END DO
C
      END

C--------------------------------------------------------
      SUBROUTINE MOVE_Y(NXMAX,NYMAX,NXMOD,NYMOD,F,MOVENUMBER)
CC--BY XIN GUO-----------------------------------------CC
C
C MOVE THE WHOLE DATA TO X DIRECTION WITH MOVENUMBER STEPS
C BY USING THE PERIODIC BOUNDARY CONDITION
C
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@--CC
C
C--INPUT
      INTEGER NXMAX, NYMAX
      INTEGER NXMOD, NYMOD
      INTEGER MOVENUMBER
C--INPUT/OUTPUT
      REAL F(NXMAX,*)
C--INTERNAL
      INTEGER I, J, K
      REAL TMP
C
      DO K = 1, MOVENUMBER
         DO J = 1, NYMOD - 1
            DO I = 1, NXMOD 
               TMP = F(I,J)
               F(I,J) = F(I,J+1)
               F(I,J+1) = TMP
            END DO
         END DO
      END DO
C
      END

C==========================================================================
      SUBROUTINE OUTSURF(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
     1     FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,IOUT,NOUT,TIME,NCPU)

C     BY DI YANG, 08/2007

      IMPLICIT NONE
      
      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
      INTEGER IOUT,NOUT
      
      REAL PEX,PEY,TWOPI,XL,YL,X,Y,DX,DY,TIME

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL FETA(NXMAX,NYMAX/NCPU,*),FVPS(NXMAX,NYMAX/NCPU,*)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)

      REAL EO(NXMAX,NYMAX),VPO(NXMAX,NYMAX)
      REAL FEO(NXMAX,NYMAX),FVO(NXMAX,NYMAX)
      REAL UO(NXMAX,NYMAX),VO(NXMAX,NYMAX),WO(NXMAX,NYMAX)

      REAL SE(NXMAX), SPF(NXMAX,NYMAX)

      REAL WORK(*), TRIGSX(*), TRIGSY(*)
      INTEGER IFAX(*)

      REAL DWK
      INTEGER NMAX
      REAL T1(NXMAX,NYMAX)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IOUT=IOUT+1
      IF(IOUT.NE.NOUT) GOTO 100
      IOUT=0

C--------------------------------
C     SEND ALL DATA TO CPU1
C--------------------------------

      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,ETA,EO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,VPS,VPO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,FETA(1,1,1),FEO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,FVPS(1,1,1),FVO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,US,UO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,VS,VO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,WS,WO,NCPU)

C-----END HERE

C--------------------------------------
C     LET CPU1 OUTPUT SURFACE DATA
C--------------------------------------

      IF(MYID.NE.0) GOTO 100

      WRITE(*,800) TIME
 800  FORMAT('Output surface data at time=',F11.5)

      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      DX=XL/FLOAT(NXMOD)
      YL=TWOPI/PEY
      DY=YL/FLOAT(NYMOD)

      WRITE(93+MYID*1000,*) " VARIABLES=X,Y,Z,VPS,U,V,W,FETA,FVPS"
CC      WRITE(93+MYID*1000,*) " VARIABLES=X,Y,Z"
      WRITE(93+MYID*1000,900) TIME,NXMOD,NYMOD
      DO J=1,NYMOD
         Y=(J-1)*DY
         DO I=1,NXMOD
            X=(I-1)*DX
            WRITE(93+MYID*1000,901) X,Y,EO(I,J),VPO(I,J),UO(I,J),
     1           VO(I,J),WO(I,J),FEO(I,J),FVO(I,J)
CC            WRITE(93+MYID*1000,901) X,Y,EO(I,J)
         ENDDO
      ENDDO
 900  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 901  FORMAT(25E12.4)

      IF ( MYID .EQ. 0 ) THEN

         DO J=1,NYMAX
            DO I=1,NXMAX
               T1(I,J)=0.
            ENDDO
         ENDDO

         CALL SPECTKX(NXMOD,NYMOD,NXMAX,NYMAX,PEX,NMAX,DWK,WORK,TRIGSX,
     1        TRIGSY,IFAX,T1,EO,SE)

         DO J=1,NYMAX
            DO I=1,NXMAX
               T1(I,J)=0.
            ENDDO
         ENDDO
         
         CALL SPECTK3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,DWK,WORK,TRIGSX,
     1        TRIGSY,IFAX,T1,EO,SPF)

         WRITE(32,976) TIME

         DO I = 1, NXMOD / 2
            X = ( I - 1 ) * PEX
            WRITE(32,*) X, SE(I)
         END DO

         WRITE(33,902) TIME, NXMOD / 2 , NYMOD / 2 
         
         DO J = 1, NYMOD / 2
            Y = ( J - 1 ) * PEY
            DO I = 1, NXMOD / 2
               X = ( I - 1 ) * PEX
               WRITE(33,*) X, Y, SPF(I,J)
            END DO
         END DO
      END IF

 902  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 976  FORMAT(' ZONE T="',F11.5,' "')
 
C-----END HERE

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE OUTSURF END HERE

C==========================================================================
      SUBROUTINE OUTSURFMODE(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,WORK,
     1     TRIGSX,IFAX,NX,FID,ETA,TIME,NCPU)
C     BY XIN GUO
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NYMOD, NXMAX, NYMAX, NCPU
      REAL PEX, PEY
      REAL ETA(NXMAX,*)
      INTEGER NX, FID
      REAL TIME
      REAL WORK(*), TRIGSX(*)
      INTEGER IFAX(*)
CC--INTERNAL
      REAL AMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
C
      CALL FFTX_PARA(ETA,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,-1)
C
      AMP = ( ETA(2*NX-1,1)**2 + ETA(2*NX,1)**2 )**0.5
C
      IF ( MYID .EQ. 0 ) THEN
         WRITE(FID,*) TIME, AMP
      END IF
C
      CALL FFTX_PARA(ETA,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,1)      
C-----END HERE

      RETURN
      END

C=====SUBROUTINE OUTSURFMODE END HERE




C===================================================================
      SUBROUTINE ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,F,FA,NCPU)

C     BY DI YANG, 05/2007

C     COLLECTING 2D PLANE DATA FOR ALL PROCESSORS INTO MYID=0

      IMPLICIT NONE
      
      INTEGER I,J,K,JS
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NUMSEND,NUMRECV
      INTEGER ITAG

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)
      
      REAL F(NXMAX,*),FA(NXMAX,NYMAX)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      
      NUMSEND=NXMAX*NYMAX/NCPU
      NUMRECV=NUMSEND

      IF(NCPU.EQ.1) THEN
         
         DO I=1,NXMOD
            DO J=1,NYMOD
               FA(I,J)=F(I,J)
            ENDDO
         ENDDO
         
      ELSE 

         IF(MYID.NE.0) THEN
         
            ITAG=MYID
            CALL MPI_SEND(F,NUMSEND,MPI_DOUBLE_PRECISION,0,ITAG,
     1           MPI_COMM_WORLD,IERR)

         ELSE IF(MYID.EQ.0) THEN
            
            DO I=1,NXMOD
               DO J=1,NYMAX/NCPU
                  FA(I,J)=F(I,J)
               ENDDO
            ENDDO

            DO K=1,NCPU-1
               JS=K*NYMAX/NCPU+1
               ITAG=K
               CALL MPI_RECV(FA(1,JS),NUMRECV,MPI_DOUBLE_PRECISION,K,
     1              ITAG,MPI_COMM_WORLD,STAT,IERR)
            ENDDO
            
         ENDIF

      ENDIF

      RETURN
      END

C=====SUBROUTINE ALLTOONE END HERE

C===================================================================
      SUBROUTINE ALLTOONE1(NXMOD,NYMOD,NXMAX,NYMAX,F,FA,NCPU)

C     BY XIN GUO, 09/2009

C     COLLECTING 2D PLANE DATA FOR ALL PROCESSORS INTO MYID=0

      IMPLICIT NONE
      
      INTEGER I,J,K,JS
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NUMSEND,NUMRECV
      INTEGER ITAG

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)
      INTEGER STATUS(MPI_STATUS_SIZE)
      
      REAL F(NXMAX,*),FA(NXMAX,NYMAX)
      REAL TMP(NXMAX,NYMAX)

      INTEGER JJ

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DO I = 1, NXMOD
         DO J = 1, NYMOD
            TMP(I,J) = 0
         END DO
      END DO
CC
      DO I = 1, NXMOD
         DO J = 1, NYMAX / NCPU
            JJ = J + MYID * NYMAX / NCPU
            TMP(I,JJ) = F(I,J)
         END DO
      END DO
      
      CALL MPI_ALLREDUCE(F,TMP,NXMAX*NYMAX,MPI_DOUBLE_PRECISION,MPI_SUM,
     1     MPI_COMM_WORLD,STATUS,IERR)
CC
      DO I = 1, NXMAX
         DO J = 1, NYMAX
            FA(I,J) = TMP(I,J)
         END DO
      END DO
CC
      RETURN
      END

C=====SUBROUTINE ALLTOONE END HERE






C===========================================================================
      SUBROUTINE LINEAR_WAVE_I_2D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,AKA,
     1     NSWAVE,FR2,ETA,VPS,NCPU)

C     BY DI YANG, 08/2007
      
C     USE 2D LINEAR WAVE THEORY AS INITIAL CONDITION

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMAX,NYMAX,NXMOD,NYMOD,NCPU,NSWAVE

      REAL PEX,XL,DX,X,FR2
      REAL AKA,AK,AA,OMEG
      REAL TWOPI

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,NYMAX/NCPU),VPS(NXMAX,NYMAX/NCPU)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C--------------------------------
C     INITIALIZE PARAMETERS
C--------------------------------

      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      DX=XL/FLOAT(NXMOD)

      AK=PEX*NSWAVE
      AA=AKA/AK
      OMEG=(AK/FR2)**0.5

C-----END HERE

C----------------------------------------
C     INITIAL CONDITION FOR SURFACE
C----------------------------------------

      DO I=1,NXMOD
         X=(I-1)*DX
         DO J=1,NYMAX/NCPU
            ETA(I,J)=ETA(I,J)+AA*COS(AK*X)
            VPS(I,J)=VPS(I,J)+AA*1./FR2/OMEG*SIN(AK*X)
         ENDDO
      ENDDO

C-----END HERE

      RETURN
      END

C=====SUBROUTINE LINEAR_WAVE_I_2D END HERE






C===========================================================================
      SUBROUTINE LINEAR_WAVE_I_3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,AKAX,
     1     AKAY,NSWAVEX,NSWAVEY,FR2,ETA,VPS,NCPU)

C     BY DI YANG, 08/2007
      
C     USE 3D LINEAR WAVE THEORY AS INITIAL CONDITION

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMAX,NYMAX,NXMOD,NYMOD,NCPU,NSWAVEX,NSWAVEY

      REAL PEX,PEY,XL,DX,X,YL,DY,Y,FR2
      REAL AKAX,AKX,AAX,AKAY,AKY,AAY,OMEGX,OMEGY
      REAL TWOPI

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,NYMAX/NCPU),VPS(NXMAX,NYMAX/NCPU)
      REAL ETR(NYMAX,NXMAX/NCPU),VPST(NYMAX,NXMAX/NCPU)
      REAL AT(NXMAX,NYMAX/NCPU),BT(NYMAX,NXMAX/NCPU)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C--------------------------------
C     INITIALIZE PARAMETERS
C--------------------------------

      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      DX=XL/FLOAT(NXMOD)
      YL=TWOPI/PEY
      DY=YL/FLOAT(NYMOD)

      AKX=PEX*NSWAVEX
      AAX=AKAX/AKX
      OMEGX=(AKX/FR2)**0.5

      AKY=PEY*NSWAVEY
      AAY=AKAY/AKY
      OMEGY=(AKY/FR2)**0.5

C-----END HERE

C----------------------------------------
C     INITIAL CONDITION FOR SURFACE
C----------------------------------------

      DO I=1,NXMOD
         X=(I-1)*DX
         DO J=1,NYMAX/NCPU
            ETA(I,J)=AAX*COS(AKX*X)
            VPS(I,J)=-AAX*1./FR2/OMEGX*SIN(AKX*X)
         ENDDO
      ENDDO

      CALL MPE_2DTRANS(ETA,AT,ETR,NXMAX,NYMAX/NCPU,NCPU)
      CALL MPE_2DTRANS(VPS,AT,VPST,NXMAX,NYMAX/NCPU,NCPU)

      DO I=1,NYMOD
         DO J=1,NXMAX/NCPU
            Y=(I-1)*DY
            ETR(I,J)=ETR(I,J)+AAY*COS(AKY*Y)
            VPST(I,J)=VPST(I,J)-AAY*1./FR2/OMEGY*SIN(AKY*Y)
         ENDDO
      ENDDO

      CALL MPE_2DTRANS(ETR,BT,ETA,NYMAX,NXMAX/NCPU,NCPU)
      CALL MPE_2DTRANS(VPST,BT,VPS,NYMAX,NXMAX/NCPU,NCPU)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE LINEAR_WAVE_I_3D END HERE

C===========================================================================
      SUBROUTINE JONSWAP_2D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,FR2,ALPHA, 
     1     OMEGA0,NU,ETA,VPS,NCPU)
CC--BY XIN GUO---------------------------------------------------------CC
C
C     2D JONSWAP
C
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@--CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL
      INTEGER MYID,NUMPROCS,IERR
CC--INPUT
      INTEGER NXMAX, NYMAX, NXMOD, NYMOD, NCPU
      REAL PEX, FR2, ALPHA, OMEGA0, NU
CC--OUTPUT
      REAL ETA(NXMAX,*), VPS(NXMAX,*)
CC--INTERNAL
      REAL TWOPI, THETA, SW, OMEGA, AA, SIGMA
      REAL XL, DX, X
      INTEGER I, J, K
      INTEGER IDUM
CC--STARTWORKING      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
C--------------------------------
C     INITIALIZE PARAMETERS
C--------------------------------
      TWOPI = 2. * ACOS(-1.)
      XL = TWOPI / PEX
      DX = XL / FLOAT(NXMOD)
    
      DO I = 1, NXMAX
         DO J = 1, NYMAX / NCPU
            ETA(I,J) = 0
            VPS(I,J) = 0
         END DO
      END DO

C-----END HERE

C----------------------------------------
C     INITIAL CONDITION FOR SURFACE
C----------------------------------------

      IDUM = 87953

C      THETA = 0

      DO K = 1, ( NXMOD + 2 ) / 2
C      DO K = 1, 2
         THETA = RAN(IDUM) * TWOPI
         OMEGA = ( K * PEX / FR2 )**0.5
         IF ( OMEGA .LT. OMEGA0 ) THEN
            SIGMA = 0.07
         ELSE
            SIGMA = 0.09
         END IF
         SW = ALPHA / FR2**2 / OMEGA**5 
     1        * EXP(-5./4.*(OMEGA/OMEGA0)**(-4.)) 
     1        * NU**(EXP(-(OMEGA-OMEGA0)**2/2/SIGMA**2/OMEGA0**2))
         SW = SW / FR2 / 2. / OMEGA
         AA = ( 2. * SW * PEX )**0.5
         DO I = 1, NXMOD
            X = ( I - 1 ) * DX
            DO J = 1, NYMAX / NCPU
               ETA(I,J) = ETA(I,J) + AA * COS(K*PEX*X+THETA)
               VPS(I,J) = VPS(I,J) + AA / FR2 / OMEGA 
     1              * SIN(K*PEX*X+THETA)
            END DO
         END DO
      END DO

C-----END HERE

      RETURN
      END

C=====SUBROUTINE LINEAR_WAVE_I_2D END HERE

CC---------------------------------------------------------------------
      SUBROUTINE JONSWAP_3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,FR2,ALPHA, 
     1     OMEGA0,NU,ETA,VPS,NCPU,WORK,TRIGSX,TRIGSY,IFAX,FR22)
CC--BY XIN GUO---------------------------------------------------------CC
C
C     2D JONSWAP
C
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@--CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL
      INTEGER MYID,NUMPROCS,IERR
CC--INPUT
      INTEGER NXMAX, NYMAX, NXMOD, NYMOD, NCPU
      REAL PEX, PEY, FR2, ALPHA, OMEGA0, NU, FR22
CC--OUTPUT
      REAL ETA(NXMAX,*), VPS(NXMAX,*)
      REAL ETA1(NXMAX,NYMAX), VPS1(NXMAX,NYMAX)
CC--FFT
      REAL WORK(*)
      REAL TRIGSX(*), TRIGSY(*)
      INTEGER IFAX(*)
CC--INTERNAL
      REAL TWOPI, THETA, THETA1, SW, OMEGA, AA, SIGMA
      REAL XL, YL, DX, DY, X, Y, K, CK, FR21
      INTEGER I, J, KX, KY, KXR, KYR, KY1, KY2, KYR1, KYR2
      INTEGER KYY1, KYY2
      INTEGER IDUM
      REAL DUM
CC--STARTWORKING      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
C--------------------------------
C     INITIALIZE PARAMETERS
C--------------------------------
      TWOPI = 2. * ACOS(-1.)
      XL = TWOPI / PEX
      YL = TWOPI / PEY
      DX = XL / FLOAT(NXMOD)
      DY = YL / NYMOD
    
      DO I = 1, NXMAX
         DO J = 1, NYMAX / NCPU
            ETA(I,J) = 0
            VPS(I,J) = 0
         END DO
      END DO

      DO I=1,NXMAX
         DO J=1,NYMAX
            ETA1(I,J)=0.
            VPS1(I,J)=0.
         ENDDO
      ENDDO

C-----END HERE

C----------------------------------------
C     INITIAL CONDITION FOR SURFACE
C----------------------------------------

      IDUM = 87953

C      DO I=1,10+MYID
C         DUM=RAN(IDUM)
C      ENDDO
      
C      THETA = 0

CC--TEST
C      WRITE(98,902) TIME, NXMOD / 2 + 1, NYMOD / 2 + 1
C      WRITE(98,901) 0., 0., 0., 0., 0.
CC--@@@@@@@@@@2

      DO KX = 0, NXMOD / 2
         DO KY = 1, NYMAX, 2
CC         DO KY = 1, NYMAX / NCPU, 2
CC
CC            KYR = ( KY + 1 + MYID * NYMAX / NCPU ) / 2 - 1
            KYR = ( KY + 1 ) / 2 - 1
CC
            IF ( KYR .LE. NYMOD / 2 ) THEN
CC
               IF ( KX .NE. 0 .OR. KYR .NE. 0 ) THEN
                  K = ( ( KX * PEX )**2 + ( KYR * PEY )**2 )**0.5
                  CK = KX * PEX / K
                  THETA = RAN(IDUM) * TWOPI
                  THETA1 = RAN(IDUM) * TWOPI
                  OMEGA = ( K / FR2 )**0.5
                  IF ( OMEGA .LT. OMEGA0 ) THEN
                     SIGMA = 0.07
                  ELSE
                     SIGMA = 0.09
                  END IF
                  SW = ALPHA / FR2**2 / OMEGA**5 
     1                 * EXP(-5./4.*(OMEGA/OMEGA0)**(-4.)) 
     1                 * NU**(EXP(-(OMEGA-OMEGA0)**2/2/SIGMA**2
     1                 /OMEGA0**2)) * 4. / TWOPI * CK**2
                  SW = SW / FR2**2 / 2. / OMEGA**3
                  AA = ( 2. * SW * PEX * PEY )**0.5
CC--TEST
C                  WRITE(*,901) KX * PEX, KY * PEY, 
C     1                 SW * FR2**2 * 2 * OMEGA**3, SW, AA
C 902              FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,
C     1                 ' F=POINT')
C 901              FORMAT(5E12.4)
CC--@@@@@@@@@@@@2
                  FR21 = FR22
                  OMEGA = ( K / FR21 )**0.5
C
CC                  KY1 = 1 + MYID * NYMAX / NCPU
CC                  KY2 = NYMAX / NCPU + MYID * NYMAX / NCPU
                  KY1 = 1
                  KY2 = NYMAX
C
                  KYR1 = 2 * KYR + 1
                  KYR2 = 2 * KYR + 2
C
CC                  KYY1 = KYR1 - MYID * NYMAX / NCPU
CC                  KYY2 = KYR2 - MYID * NYMAX / NCPU
                  KYY1 = KYR1
                  KYY2 = KYR2
CC--TEST
                  IF ( KYR1 .EQ. 1 .AND. KYR2 .EQ. 2 ) THEN
                     THETA1 = THETA
                  END IF
CC--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
CC
                  AA = AA / 4.

                  IF(KX.EQ.0.OR.KYR.EQ.0) THEN 
                     AA=AA*SQRT(2.)
                  ENDIF
CC
                  IF ( KYR1 .GE. KY1 .AND. KYR1 .LE. KY2 ) THEN
CC                     ETA(KX*2+1,KYY1) = AA * COS(THETA) + AA 
CC     1                    * COS(THETA1)
CC                     ETA(KX*2+2,KYY1) = AA * SIN(THETA) + AA 
CC     1                    * SIN(THETA1)
CC                     VPS(KX*2+1,KYY1) = AA / FR21 / OMEGA * SIN(THETA) 
CC     1                    + AA / FR21 / OMEGA * SIN(THETA1)
CC                     VPS(KX*2+2,KYY1) = - AA / FR21 / OMEGA * COS(THETA) 
CC     1                    - AA / FR21 / OMEGA * COS(THETA1)
                     ETA1(KX*2+1,KYY1) = AA * COS(THETA) + AA 
     1                    * COS(THETA1)
                     ETA1(KX*2+2,KYY1) = AA * SIN(THETA) + AA 
     1                    * SIN(THETA1)
                     VPS1(KX*2+1,KYY1) = AA / FR21 / OMEGA * SIN(THETA) 
     1                    + AA / FR21 / OMEGA * SIN(THETA1)
                     VPS1(KX*2+2,KYY1) = - AA / FR21 / OMEGA 
     1                    * COS(THETA) - AA / FR21 / OMEGA * COS(THETA1)
                  END IF
CC
                  IF ( KYR2 .GE. KY1 .AND. KYR2 .LE. KY2 ) THEN
CC                     ETA(KX*2+1,KYY2) = AA * SIN(THETA) - AA 
CC     1                    * SIN(THETA1)
CC                     ETA(KX*2+2,KYY2) = - AA * COS(THETA) + AA 
CC     1                    * COS(THETA1)
CC                     VPS(KX*2+1,KYY2) = - AA / FR21 / OMEGA * COS(THETA) 
CC     1                    + AA / FR21 / OMEGA * COS(THETA1)
CC                     VPS(KX*2+2,KYY2) = - AA / FR21 / OMEGA * SIN(THETA) 
CC     1                    + AA / FR21 / OMEGA * SIN(THETA1)
                     ETA1(KX*2+1,KYY2) = AA * SIN(THETA) - AA 
     1                    * SIN(THETA1)
                     ETA1(KX*2+2,KYY2) = - AA * COS(THETA) + AA 
     1                    * COS(THETA1)
                     VPS1(KX*2+1,KYY2) = - AA / FR21 / OMEGA 
     1                    * COS(THETA) + AA / FR21 / OMEGA * COS(THETA1)
                     VPS1(KX*2+2,KYY2) = - AA / FR21 / OMEGA 
     1                    * SIN(THETA) + AA / FR21 / OMEGA * SIN(THETA1)
                  END IF
CC
               END IF
CC
            END IF
         END DO
      END DO

C--------------------------------
C     SEND ALL DATA TO CPU1
C--------------------------------

      CALL ONETOALL(NXMOD,NYMOD,NXMAX,NYMAX,ETA1,ETA,NCPU)
      CALL ONETOALL(NXMOD,NYMOD,NXMAX,NYMAX,VPS1,VPS,NCPU)

C-----END HERE      

      CALL FFTXY_PARA(ETA,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,1)
      CALL FFTXY_PARA(VPS,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,1)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE LINEAR_WAVE_I_3D END HERE
C==========================================================================
      SUBROUTINE OUTSURF_A(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
     1     FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,IOUT,NOUT,TIME,NCPU)

C     BY DI YANG, 08/2007
C     MODIFIED BY HAMID      

      IMPLICIT NONE
      
      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
      INTEGER IOUT,NOUT,ISTART
      
      REAL PEX,PEY,TWOPI,XL,YL,X,Y,DX,DY,TIME

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL FETA(NXMAX,NYMAX/NCPU,*),FVPS(NXMAX,NYMAX/NCPU,*)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)

      REAL EO(NXMAX,NYMAX),VPO(NXMAX,NYMAX)
      REAL FEO(NXMAX,NYMAX),FVO(NXMAX,NYMAX)
      REAL UO(NXMAX,NYMAX),VO(NXMAX,NYMAX),WO(NXMAX,NYMAX)

      REAL SE(NXMAX), SPF(NXMAX,NYMAX)

      REAL WORK(*), TRIGSX(*), TRIGSY(*)
      INTEGER IFAX(*)

      REAL DWK
      INTEGER NMAX

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
C
      IOUT=IOUT+1
      IF(IOUT.NE.NOUT) GOTO 100
      IOUT=0
C
      IF (MYID.EQ.0)THEN
         WRITE(*,800) TIME
 800     FORMAT('Output surface data at time=',F11.5)
      ENDIF
C
      WRITE(93+MYID*1000,900) TIME,NXMOD,NYMOD
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            WRITE(93+MYID*1000,901) ETA(I,J),VPS(I,J),US(I,J),
     1           VS(I,J),WS(I,J),FETA(I,J,1),FVPS(I,J,1)
         ENDDO
      ENDDO
 900  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 901  FORMAT(7E12.4)
 
C-----END HERE

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE OUTSURF END HERE

C==========================================================================
C==========================================================================
      SUBROUTINE REREAD_A(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
     1     FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,ISTART,TIME,NCPU)

C     ADDED BY HAMID      

      IMPLICIT NONE
      
      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
      INTEGER IOUT,NOUT,ISTART
      
      REAL PEX,PEY,TWOPI,XL,YL,X,Y,DX,DY,TIME

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL FETA(NXMAX,NYMAX/NCPU,*),FVPS(NXMAX,NYMAX/NCPU,*)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)

      REAL EO(NXMAX,NYMAX),VPO(NXMAX,NYMAX)
      REAL FEO(NXMAX,NYMAX),FVO(NXMAX,NYMAX)
      REAL UO(NXMAX,NYMAX),VO(NXMAX,NYMAX),WO(NXMAX,NYMAX)

      REAL SE(NXMAX), SPF(NXMAX,NYMAX)

      REAL WORK(*), TRIGSX(*), TRIGSY(*)
      INTEGER IFAX(*)

      REAL DWK
      INTEGER NMAX

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF (ISTART.NE.1) GOTO 100
C
      READ(14+MYID*1000,900) TIME,NXMOD,NYMOD
      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            READ(14+MYID*1000,901) ETA(I,J),VPS(I,J),US(I,J),
     1           VS(I,J),WS(I,J),FETA(I,J,1),FVPS(I,J,1)
         ENDDO
      ENDDO
 900  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 901  FORMAT(7E12.4)
 
C-----END HERE

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE OUTSURF END HERE

C==========================================================================
      SUBROUTINE SAVE_A(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,ETA,VPS,FETA,
     1     FVPS,US,VS,WS,WORK,TRIGSX,TRIGSY,IFAX,IOUT,NOUT,TIME,NCPU)

C     BY DI YANG, 08/2007
C     MODIFIED BY HAMID      

      IMPLICIT NONE
      
      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
      INTEGER IOUT,NOUT,ISTART
      
      REAL PEX,PEY,TWOPI,XL,YL,X,Y,DX,DY,TIME

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL FETA(NXMAX,NYMAX/NCPU,*),FVPS(NXMAX,NYMAX/NCPU,*)
      REAL US(NXMAX,*),VS(NXMAX,*),WS(NXMAX,*)

      REAL EO(NXMAX,NYMAX),VPO(NXMAX,NYMAX)
      REAL FEO(NXMAX,NYMAX),FVO(NXMAX,NYMAX)
      REAL UO(NXMAX,NYMAX),VO(NXMAX,NYMAX),WO(NXMAX,NYMAX)

      REAL SE(NXMAX), SPF(NXMAX,NYMAX)

      REAL WORK(*), TRIGSX(*), TRIGSY(*)
      INTEGER IFAX(*)

      REAL DWK
      INTEGER NMAX

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      NCPU=32
      DO MYID=0,NCPU-1
CC
      WRITE(13+MYID*1000,900) TIME,NXMOD,NYMOD
CC
      WRITE(*,*)MYID,(NYMAX/NCPU)*MYID+1,(NYMAX/NCPU)*(MYID+1)
      DO J=(NYMAX/NCPU)*MYID+1,(NYMAX/NCPU)*(MYID+1)
CC
C      DO J=1,NYMAX/NCPU
         DO I=1,NXMOD
            WRITE(13+MYID*1000,901) ETA(I,J),VPS(I,J),US(I,J),
     1           VS(I,J),WS(I,J),FETA(I,J,1),FVPS(I,J,1)
         ENDDO
      ENDDO
CC
      ENDDO
CC
 900  FORMAT(' ZONE T="',F11.5,'" I=',I4,' J=',I4,' F=POINT')
 901  FORMAT(7E12.4)
 
C-----END HERE

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE OUTSURF END HERE

C==========================================================================







C===========================================================================
      SUBROUTINE SAVE_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETAB,VPS,PA0,ICON,
     1     TIME,FR2,NCPU)

C     SAVE DATA FOR CONTINUING HOS SIMULATION

      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER I,J,ICON
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER MYID,NUMPROCS,IERR

      INTEGER NCREQ,JJ1,JJ2,NYEND

      REAL TIME,FR2

      REAL ETAB(NXMAX,*),VPS(NXMAX,*),PA0(NXMAX,*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NCREQ IS THE # OF CPUS THAT CONTAIN USEFUL ELEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX

      JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
      JJ2=NYMAX/NCPU

      IF(MYID.EQ.NCREQ-1) THEN
         NYEND=JJ1
      ELSE IF(MYID.LT.NCREQ-1) THEN
         NYEND=JJ2
      ELSE
         NYEND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(NYEND.EQ.0) GOTO 100

      ICON=ICON+1

      OPEN(ICON*1000+MYID)

      WRITE(ICON*1000+MYID,*) NXMOD,NYMOD,NCPU,TIME,FR2
      
      DO I=1,NXMOD
         DO J=1,NYEND
            WRITE(ICON*1000+MYID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
         ENDDO
      ENDDO

      CLOSE(ICON*1000+MYID)
      
 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE SAVE_HOS END HERE








C===========================================================================
      SUBROUTINE REREAD_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETAB,VPS,PA0,ICON,
     1     TIME,FR2,NCPU)

C     READIN DATA FOR CONTINUING HOS SIMULATION

      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER I,J,ICON
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NXMOD1,NYMOD1,NCPU1
      INTEGER MYID,NUMPROCS,IERR

      INTEGER NCREQ,JJ1,JJ2,NYEND

      REAL TIME,FR2,FR21

      REAL ETAB(NXMAX,*),VPS(NXMAX,*),PA0(NXMAX,*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NCREQ IS THE # OF CPUS THAT CONTAIN USEFUL ELEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX

      JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
      JJ2=NYMAX/NCPU

      IF(MYID.EQ.NCREQ-1) THEN
         NYEND=JJ1
      ELSE IF(MYID.LT.NCREQ-1) THEN
         NYEND=JJ2
      ELSE
         NYEND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(NYEND.EQ.0) GOTO 100

      OPEN(ICON*1000+MYID)

      READ(ICON*1000+MYID,*) NXMOD1,NYMOD1,NCPU1,TIME,FR21
      
      IF(NXMOD1.NE.NXMOD) THEN
         PRINT*, 'Inconsistent NXMOD:'
         PRINT*, 'NXMOD=',NXMOD
         PRINT*, 'NXMOD1=',NXMOD1
         STOP
      ENDIF
      
      IF(NYMOD1.NE.NYMOD) THEN
         PRINT*, 'Inconsistent NYMOD:'
         PRINT*, 'NYMOD=',NYMOD
         PRINT*, 'NYMOD1=',NYMOD1
         STOP
      ENDIF

      IF(ABS(FR21-FR2).GT.1.E-4) THEN
         PRINT*, 'Inconsistent FR2:'
         PRINT*, 'FR2=',FR2
         PRINT*, 'FR21=',FR21
         STOP
      ENDIF

      DO I=1,NXMOD
         DO J=1,NYEND
            READ(ICON*1000+MYID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
         ENDDO
      ENDDO

      CLOSE(ICON*1000+MYID)
      
 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE REREAD_HOS END HERE








C===========================================================================
      SUBROUTINE SAVE_LES_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETAB,VPS,PA0,
     1     ICON,TIME,FR2,NCPU)

C     SAVE DATA FOR CONTINUING HOS SIMULATION

      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER I,J,ICON
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER MYID,NUMPROCS,IERR

      INTEGER NCREQ,JJ1,JJ2,NYEND

      integer ncpu1,K,id
      parameter(ncpu1=10)

      REAL TIME,FR2

      REAL ETAB(NXMAX,*),VPS(NXMAX,*),PA0(NXMAX,*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NCREQ IS THE # OF CPUS THAT CONTAIN USEFUL ELEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX

      JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
      JJ2=NYMAX/NCPU

      IF(MYID.EQ.NCREQ-1) THEN
         NYEND=JJ1
      ELSE IF(MYID.LT.NCREQ-1) THEN
         NYEND=JJ2
      ELSE
         NYEND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC      IF(NYEND.EQ.0) GOTO 100

      ICON=ICON+1

CC      OPEN(ICON*1000+MYID)

CC      DO I=1,NXMOD
CC         DO J=1,NYEND
CC            WRITE(ICON*1000+MYID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
CC         ENDDO
CC      ENDDO

CC      CLOSE(ICON*1000+MYID)

      DO K=1,NCPU1
         ID=K-1
         OPEN(ICON*1000+ID)

         NCREQ=NCPU1-(NYMAX-NYMOD)*NCPU1/NYMAX
         
         JJ1=NYMAX/NCPU1-MOD(NYMAX-NYMOD,NYMAX/NCPU1)
         JJ2=NYMAX/NCPU1

         IF(ID.EQ.NCREQ-1) THEN
            NYEND=JJ1
         ELSE IF(MYID.LT.NCREQ-1) THEN
            NYEND=JJ2
         ELSE
            NYEND=0
         ENDIF
         
         DO I=1,NXMOD
            DO J=ID*NYMAX/NCPU1+1,ID*NYMAX/NCPU1+NYEND
               WRITE(ICON*1000+ID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
            ENDDO
         ENDDO

         CLOSE(ICON*1000+ID)

      ENDDO

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE SAVE_LES_HOS END HERE







C===========================================================================
      SUBROUTINE REREADN_HOS(NXMOD,NYMOD,NXMAX,NYMAX,ETAB,VPS,PA0,ICON,
     1     TIME,FR2,NCPU)

C     SAVE DATA FOR CONTINUING HOS SIMULATION

      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER I,J,ICON,I1,J1,J0
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU,NCPU1,NXMAX1,NYMAX1
      INTEGER MYID,NUMPROCS,IERR

      INTEGER NCREQ,JJ1,JJ2,NYEND,NYEND2
      INTEGER ID1START,ID1END,NY1START,NY1END,NID,NYSTART1,NYEND1

      INTEGER IDUM
      REAL DUM

      INTEGER NXMOD1,NYMOD1
      REAL TIME,FR2,FR21

      REAL ETAB(NXMAX,*),VPS(NXMAX,*),PA0(NXMAX,*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NCREQ IS THE # OF CPUS THAT CONTAIN USEFUL ELEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      NCREQ=NCPU-(NYMAX-NYMOD)*NCPU/NYMAX

      JJ1=NYMAX/NCPU-MOD(NYMAX-NYMOD,NYMAX/NCPU)
      JJ2=NYMAX/NCPU

      IF(MYID.EQ.NCREQ-1) THEN
         NYEND=JJ1
      ELSE IF(MYID.LT.NCREQ-1) THEN
         NYEND=JJ2
      ELSE
         NYEND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(NYEND.EQ.0) GOTO 100

      OPEN(ICON*1000)
      READ(ICON*1000,*) NXMOD1,NYMOD1,NCPU1,TIME,FR21
      CLOSE(ICON*1000)
      
      IF(NXMOD1.NE.NXMOD) THEN
         PRINT*, 'Inconsistent NXMOD:'
         PRINT*, 'NXMOD=',NXMOD
         PRINT*, 'NXMOD1=',NXMOD1
         STOP
      ENDIF
      
      IF(NYMOD1.NE.NYMOD) THEN
         PRINT*, 'Inconsistent NYMOD:'
         PRINT*, 'NYMOD=',NYMOD
         PRINT*, 'NYMOD1=',NYMOD1
         STOP
      ENDIF

      IF(ABS(FR21-FR2).GT.1.E-4) THEN
         PRINT*, 'Inconsistent FR2:'
         PRINT*, 'FR2=',FR2
         PRINT*, 'FR21=',FR21
         STOP
      ENDIF

      IF(NCPU.NE.NCPU1) GOTO 101

      OPEN(ICON*1000+MYID)
      DO I=1,NXMOD
         DO J=1,NYEND
            READ(ICON*1000+MYID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
         ENDDO
      ENDDO
      CLOSE(ICON*1000+MYID)

 101  CONTINUE

      IF(MOD(NXMOD+2,NCPU1).EQ.0) THEN
         NXMAX1=NXMOD+2
      ELSE
         NXMAX1=NXMOD+2+NCPU1-MOD(NXMOD+2,NCPU1)
      ENDIF
      IF(MOD(NYMOD+2,NCPU1).EQ.0) THEN
         NYMAX1=NYMOD+2
      ELSE
         NYMAX1=NYMOD+2+NCPU1-MOD(NYMOD+2,NCPU1)
      ENDIF

      print*, 'nxmaxs1=',nxmax1
      print*, 'nymaxs1=',nymax1
      print*, 'ncpu1=',ncpu1

      ID1START=NYMAX/NCPU*MYID/(NYMAX1/NCPU1)
      NY1START=NYMAX/NCPU*MYID+1-NYMAX1/NCPU1*ID1START
      ID1END=(NYMAX/NCPU*(MYID+1)-1)/(NYMAX1/NCPU1)
      NY1END=NYMAX/NCPU*(MYID+1)-NYMAX1/NCPU1*ID1END
      J0=0
      DO NID=ID1START,ID1END
         IF(ID1START.EQ.ID1END) THEN
            NYSTART1=NY1START
            NYEND1=NY1END
         ELSE IF(NID.EQ.ID1START) THEN
            NYSTART1=NY1START
            NYEND1=NYMAX1/NCPU1
         ELSEIF(NID.EQ.ID1END) THEN
            NYSTART1=1
            NYEND1=NY1END
         ELSE
            NYSTART1=1
            NYEND1=NYMAX1/NCPU1
         ENDIF

         IF(NID.EQ.NCPU1-1) THEN
            NYEND2=NYMOD-(NYMAX1/NCPU1)*NID
         ELSE
            NYEND2=NYMAX1/NCPU1
         ENDIF

         OPEN(ICON*1000+NID)
         DO I=1,NXMOD
            DO J1=1,NYEND2
               IF(J1.GE.NYSTART1 .AND. J1.LE.NYEND1) THEN
                  J=J0+(J1-NYSTART1)+1
                  READ(ICON*1000+NID,*) ETAB(I,J),VPS(I,J),PA0(I,J)
               ELSE
                  READ(ICON*1000+NID,*) DUM,DUM,DUM
               ENDIF
            ENDDO
         ENDDO
         CLOSE(ICON*1000+NID)
         J0=J0+NYEND1-NYSTART1+1
      ENDDO

 100  CONTINUE

      RETURN
      END

C=====SUBROUTINE REREADN_HOS END HERE
