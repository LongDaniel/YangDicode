C========================================================================
      SUBROUTINE SMOTH1(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,
     1     TRIGSY,ETA,VPS,FAC,NCPU)

C     SIMPLE LOW PASS SPECTRAL FILTER

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NXPAT,NYPAT

      REAL FAC

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR
      
      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL ET(NYMAX,NXMAX/NCPU),VT(NYMAX,NXMAX/NCPU)
      REAL AT(NXMAX,NYMAX/NCPU),BT(NYMAX,NXMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      CALL FFTX_PARA(ETA,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,-1)
      CALL FFTX_PARA(VPS,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,-1)

      NXPAT=NXMOD/2*FAC
      DO I=1,NXMOD+2
         DO J=1,NYMAX/NCPU
            IF(I.GT.NXPAT) THEN
               ETA(I,J)=0.
               VPS(I,J)=0.
            ENDIF
         ENDDO
      ENDDO
      
      CALL MPE_2DTRANS(ETA,AT,ET,NXMAX,NYMAX/NCPU,NCPU)
      CALL MPE_2DTRANS(VPS,AT,VT,NXMAX,NYMAX/NCPU,NCPU)

      CALL FFTX_PARA(ET,WORK,TRIGSY,IFAX,NYMOD,NYMAX,NXMAX,NCPU,-1)
      CALL FFTX_PARA(VT,WORK,TRIGSY,IFAX,NYMOD,NYMAX,NXMAX,NCPU,-1)

      NYPAT=NYMOD/2*FAC
      DO I=1,NYMOD+2
         DO J=1,NXMAX/NCPU
            IF(I.GT.NYPAT) THEN
               ET(I,J)=0.
               VT(I,J)=0.
            ENDIF
         ENDDO
      ENDDO

      CALL FFTX_PARA(ET,WORK,TRIGSY,IFAX,NYMOD,NYMAX,NXMAX,NCPU,1)
      CALL FFTX_PARA(VT,WORK,TRIGSY,IFAX,NYMOD,NYMAX,NXMAX,NCPU,1)
      
      CALL MPE_2DTRANS(ET,BT,ETA,NYMAX,NXMAX/NCPU,NCPU)
      CALL MPE_2DTRANS(VT,BT,VPS,NYMAX,NXMAX/NCPU,NCPU)

      CALL FFTX_PARA(ETA,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,1)
      CALL FFTX_PARA(VPS,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX,NCPU,1)

      RETURN
      END

C=====SUBROUTINE SMOTH1 END HERE






C========================================================================
      SUBROUTINE SMOOTH1(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,WVN,ETA,VPS,RATIO,NCPU)

C     BY DI YANG, 08/2007

C     DAMPING HIGHEST SPECTRAL MODE WITH A SMOOTH FUNCTION

      IMPLICIT NONE

      INTEGER L,M
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU,MODEX,MODEY
      
      REAL PEX,PEY,AX,AY,RATIO
      REAL FACT_FUN,WVF,WVMAX,DELTA,DWVN,FACTOR

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),VPS(NXMAX,*)
      REAL WVN(NXMAX,NYMAX/NCPU,*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C--------------------------------
C     INITIALIZE PARAMETERS
C--------------------------------

      MODEX=((NXMOD+1)/2)/3*2
      MODEY=((NYMOD+1)/2)/3*2
      AX=PEX*MODEX
      AY=PEY*MODEY
      WVMAX=MIN(AX,AY)
      WVF=RATIO*WVMAX
      DELTA=2.*(WVMAX-WVF)

C-----END HERE

C------------------------------------------------
C     DAMPING HIGHEST MODE IN SPECTRAL SPACE
C------------------------------------------------

      CALL FFTXY_PARA(ETA,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)
      CALL FFTXY_PARA(VPS,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)

      DO M=1,NYMAX/NCPU
         DO L=1,NXMOD+2
            IF(WVN(L,M,1).GE.1.E-6) THEN
               DWVN=WVMAX-WVN(L,M,1)
               IF(DWVN.LE.0.) DWVN=0.
               FACTOR=FACT_FUN(DWVN,DELTA)
               ETA(L,M)=ETA(L,M)*FACTOR
               VPS(L,M)=VPS(L,M)*FACTOR
            ENDIF
         ENDDO
      ENDDO

      CALL FFTXY_PARA(ETA,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,1)
      CALL FFTXY_PARA(VPS,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,1)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE SMOOTH1 END HERE





C=======================================================================
      FUNCTION FACT_FUN(T,T0)

      IMPLICIT NONE

      REAL T,T0
      REAL FACT_FUN

      IF(T.GT.T0) THEN
         FACT_FUN=1.
      ELSE
         FACT_FUN=462.*(T/T0)**6-1980.*(T/T0)**7+3465.*(T/T0)**8
     1        -3080.*(T/T0)**9+1386.*(T/T0)**10-252.*(T/T0)**11
      ENDIF

      RETURN
      END

C=====FUNCTION FACT_FUN END HERE







C=========================================================================
      SUBROUTINE SMOOTH2(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,
     1     TRIGSY,PEX,PEY,ETA,VPS,EX,EY,TIME,ISMOOTH2,NCPU)

C     BY DI YANG, 08/2010

C     SMOOTH THE SURFACE IN PHYSICAL SPACE BY 9-POINTS LOCAL AVERAGE

      IMPLICIT NONE

      INTEGER I,J
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU,NSMOOTH,IBREAK
      INTEGER ISMOOTH,ISMOOTH2
     
      REAL PEX,PEY
      REAL ELIMIT,TIME,EXMAX,EYMAX

      DATA ELIMIT/0.7/
 
      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR
      
      REAL ETA(NXMAX,*),EX(NXMAX,*),EY(NXMAX,*),VPS(NXMAX,*)
      REAL EO(NXMAX,NYMAX),EXO(NXMAX,NYMAX),EYO(NXMAX,NYMAX)
      REAL VPO(NXMAX,NYMAX)
      REAL FTMP(NXMAX,NYMAX)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C--------------------------------
C     SEND ALL DATA TO CPU1
C--------------------------------

      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,ETA,EO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,VPS,VPO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,EX,EXO,NCPU)
      CALL ALLTOONE(NXMOD,NYMOD,NXMAX,NYMAX,EY,EYO,NCPU)

C-----END HERE

      ISMOOTH=0
      ISMOOTH2=0
C     ISMOOTH: 0--NO SMOOTH HAS BEEN DONE; 1--SMOOTHED

      IF(MYID.NE.0) GOTO 100

C--------------------------------------
C     INITIALIZE SMOOTHING COUNTER
C--------------------------------------

      NSMOOTH=0

 1    CONTINUE
      IBREAK=0
      NSMOOTH=NSMOOTH+1

C-----END HERE

C---------------------------------------
C     CHECK LOCAL SLOPE ON SURFACE
C     IF TOO STEEP
C     APPLY 9-POINTS LOCAL AVERAGE
C---------------------------------------

      EXMAX=0.
      EYMAX=0.
      DO J=1,NYMOD
         DO I=1,NXMOD
            IF(ABS(EXO(I,J)).GT.EXMAX) EXMAX=ABS(EXO(I,J))
            IF(ABS(EYO(I,J)).GT.EYMAX) EYMAX=ABS(EYO(I,J))
            IF(ABS(EXO(I,J)).GT.ELIMIT.OR.ABS(EYO(I,J)).GT.ELIMIT) THEN
               CALL SMOOTH2_SUB(NXMOD,NYMOD,NXMAX,NYMAX,EO,I,J)
               CALL SMOOTH2_SUB(NXMOD,NYMOD,NXMAX,NYMAX,VPO,I,J)
               IBREAK=1
               ISMOOTH=1
CC               GOTO 50
            ENDIF
         ENDDO
      ENDDO

 50   CONTINUE

      WRITE(22,104) TIME,EXMAX,EYMAX
      WRITE(*,104) TIME,EXMAX,EYMAX
 104  FORMAT(25E12.5)

C-----END HERE

C--------------------------------------------
C     CHECK IF REACH MAX SMOOTHING STEP
C--------------------------------------------

      IF(IBREAK.EQ.1) THEN
         IF(NSMOOTH.GT.50) THEN
            GOTO 100
         ENDIF
         WRITE(22,*) "HOS BREAK"
         WRITE(*,*) "HOS BREAK"
         CALL PDFX(EO,EXO,FTMP,WORK,TRIGSX,IFAX,PEX,NXMOD,NYMOD,NXMAX)
         CALL PDFY(EO,EYO,FTMP,WORK,TRIGSY,IFAX,PEY,NXMOD,NYMOD,NXMAX)
         GOTO 1
      ENDIF

C-----END HERE

 100  CONTINUE

C-----------------------------------------
C     SEND ALL DATA BACK TO EACH CPU
C-----------------------------------------

      CALL ONETOALL(NXMOD,NYMOD,NXMAX,NYMAX,EO,ETA,NCPU)
      CALL ONETOALL(NXMOD,NYMOD,NXMAX,NYMAX,VPO,VPS,NCPU)

      CALL MPI_ALLREDUCE(ISMOOTH,ISMOOTH2,1,MPI_INTEGER,MPI_MAX,
     1     MPI_COMM_WORLD,IERR)

C-----END HERE

      RETURN
      END

C=====SUBROUTINE SMOOTH2 END HERE







C======================================================================
      SUBROUTINE PDFX(F,FX,FTMP,WORK,TRIGS,IFAX,PEX,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FX(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTX(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 L=1,NXMOD+1,2
       PMODX=PEX*(L-1)/2
       DO 20 M=1,NYMOD
        FX(L,M)=-PMODX*FTMP(L+1,M)
        FX(L+1,M)=PMODX*FTMP(L,M)
 20   CONTINUE
      CALL FFTX(FX,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END

C=====SUBROUTINE PDFX END HERE






C=========================================================================
      SUBROUTINE PDFY(F,FY,FTMP,WORK,TRIGS,IFAX,PEY,NXMOD,NYMOD,NXMAX)
C
      REAL F(NXMAX,*),FY(NXMAX,*),FTMP(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
       FTMP(I,J)=F(I,J)
 10   CONTINUE
      CALL FFTY(FTMP,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,-1)
      DO 20 M=1,NYMOD+1,2
       PMODY=PEY*(M-1)/2
       DO 20 L=1,NXMOD
        FY(L,M)=-PMODY*FTMP(L,M+1)
        FY(L,M+1)=PMODY*FTMP(L,M)
 20   CONTINUE
      CALL FFTY(FY,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,1)
      RETURN
      END

C=====SUBROUTINE PDFY END HERE







C=========================================================================
      SUBROUTINE FFTX(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)
C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 20 J=1,NYMOD
      DO 20 I=NXMOD+1,NXMOD+2
 20   F(I,J)=0.
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END

C=====SUBROUTINE FFTX END HERE





C===================================================================
      SUBROUTINE FFTY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)
C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 10 I=1,NXMOD
      DO 10 J=NYMOD+1,NYMOD+2
 10   F(I,J)=0.
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END

C=====SUBROUTINE FFTY END HERE






C=====================================================================
      SUBROUTINE SMOOTH2_SUB(NXMOD,NYMOD,NXMAX,NYMAX,F,I,J)

C     BY DI YANG, 08/2007

C     9-POINTS LOCAL AVERAGE

      IMPLICIT NONE
      
      INTEGER I,J,KS,I1,J1,IP,JP,IM,JM,IDEX,JDEX
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX
      
      DATA KS/32/

      REAL F(NXMAX,*),TMP(NXMAX,NYMAX)

      DO J1=J-KS/2,J+KS/2
         JP=J1
         IF(JP.GT.NYMOD) JP=JP-NYMOD
         IF(JP.LT.1) JP=JP+NYMOD
         DO I1=I-KS/2,I+KS/2
            IP=I1
            IF(IP.GT.NXMOD) IP=IP-NXMOD
            IF(IP.LT.1) IP=IP+NXMOD
            TMP(IP,JP)=0.
            DO JM=-1,1
               DO IM=-1,1
                  IDEX=IP+IM
                  JDEX=JP+JM
                  IF(IDEX.GT.NXMOD) IDEX=IDEX-NXMOD
                  IF(IDEX.LT.1) IDEX=IDEX+NXMOD
                  IF(JDEX.GT.NYMOD) JDEX=JDEX-NYMOD
                  IF(JDEX.LT.1) JDEX=JDEX+NYMOD
                  TMP(IP,JP)=TMP(IP,JP)+F(IDEX,JDEX)
               ENDDO
            ENDDO
            F(IP,JP)=TMP(IP,JP)/9.
         ENDDO
      ENDDO
      
      RETURN
      END

C=====SUBROUTINE SMOOTH2_SUB END HERE





C===================================================================
      SUBROUTINE ONETOALL(NXMOD,NYMOD,NXMAX,NYMAX,FA,F,NCPU)

C     BY DI YANG, MAY 2011
C     DISTRIBUTE 2D PLANE DATA FROM THE FIRST CPU TO ALL CPUS

      IMPLICIT NONE
      
      INTEGER I,J,K,JS
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NCREQ,JJ1,JJ2,NYEND
C      INTEGER NUMSEND,NUMRECV
C      INTEGER ITAG

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)
      
      REAL F(NXMAX,*),F0(NXMAX,NYMAX),FA(NXMAX,NYMAX)

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
      
C-------------------------
C     INITIALIZATION
C-------------------------

      DO J=1,NYMAX/NCPU
         DO I=1,NXMAX
            F(I,J)=0.
         ENDDO
      ENDDO

      DO J=1,NYMAX
         DO I=1,NXMAX
            F0(I,J)=0.
         ENDDO
      ENDDO

      IF(MYID.NE.0) THEN
         DO J=1,NYMAX
            DO I=1,NXMAX
               FA(I,J)=0.
            ENDDO
         ENDDO
      ENDIF

C-----END HERE

C---------------------------------------
C     DISTRIBUTE VALUE TO ALL CPUS
C---------------------------------------

      CALL MPI_ALLREDUCE(FA,F0,NXMAX*NYMAX,MPI_DOUBLE_PRECISION,
     1     MPI_SUM,MPI_COMM_WORLD,IERR)

      DO J=1,NYMAX/NCPU
         JS=MYID*NYMAX/NCPU+J
         DO I=1,NXMAX
            F(I,J)=F0(I,JS)
         ENDDO
      ENDDO
            

C-----END HERE

      RETURN
      END

C=====SUBROUTINE ONETOALL END HERE





C=====================================================================
      SUBROUTINE ONETOALL_OLD(NXMOD,NYMOD,NXMAX,NYMAX,FA,F,NCPU)
      
C     BY DI YANG, 08/2007

C     DISTRABUTING 2D PLANCE DATA FROM MYID=0 TO ALL PROCESSORS

      IMPLICIT NONE
      
      INTEGER I,J,K,JS
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NUMSEND,NUMRECV
      INTEGER ITAG

      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)

      REAL F(NXMAX,NYMAX/NCPU),FA(NXMAX,NYMAX)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NUMSEND=NXMAX*NYMAX/NCPU
      NUMRECV=NUMSEND

      IF(NCPU.EQ.1) THEN
         
         DO I=1,NXMOD
            DO J=1,NYMOD
               F(I,J)=FA(I,J)
            ENDDO
         ENDDO
         
      ELSE

         IF(MYID.EQ.0) THEN
            
            DO I=1,NXMOD
               DO J=1,NYMAX/NCPU
                  F(I,J)=FA(I,J)
               ENDDO
            ENDDO
            
            DO K=1,NCPU-1
               JS=K*NYMAX/NCPU+1
               ITAG=K
               CALL MPI_SEND(FA(1,JS),NUMSEND,MPI_DOUBLE_PRECISION,K,
     1              ITAG,MPI_COMM_WORLD,IERR)
            ENDDO

         ELSE IF(MYID.NE.0) THEN
            
            ITAG=MYID
            CALL MPI_RECV(F,NUMRECV,MPI_DOUBLE_PRECISION,0,ITAG,
     1           MPI_COMM_WORLD,STAT,IERR)
            
         ENDIF
         
      ENDIF

      RETURN
      END

C=====SUBROUTINE ONETOALL END HERE
