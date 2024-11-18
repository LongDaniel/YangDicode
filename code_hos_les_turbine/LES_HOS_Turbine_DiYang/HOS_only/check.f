C=========================================================================
      SUBROUTINE CHECK(NXMOD,NYMOD,NXMAX,NYMAX,WORK,IFAX,TRIGSX,TRIGSY,
     1     PEX,PEY,ITIM,TIME,FR2,ETA,VPS,FETA,ENER0,IBLOWUP,NCPU)

C     BY DI YANG, 08/2007

C     CALCULATE POTENTIAL AND KINETIC ENERGY
C     CALCULATE VOLUME CONSERVATION AND FLUX
C     CHECK IF CODE BLOWS UP OR NOT

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU,ITIM,IBLOWUP
      
      REAL FR2
      REAL PEX,PEY,TIME,ENER,ENER0
      REAL XL,YL,TWOPI
      REAL VOL,FLUX,PE,KE

      INCLUDE "mpif.h"
      
      INTEGER MYID,NUMPROCS,IERR

      REAL ETA(NXMAX,*),FETA(NXMAX,NYMAX/NCPU,*),VPS(NXMAX,*)
      REAL T1(NXMAX,NYMAX/NCPU),T2(NXMAX,NYMAX/NCPU)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C-------------------------------
C     CALCULATE DOMAIN SIZE
C-------------------------------
      
      TWOPI=2.*ACOS(-1.)
      XL=TWOPI/PEX
      YL=TWOPI/PEY

C-----END HERE

C------------------------------
C     VOLUME CONSERVATION
C------------------------------
      
      DO I=1,NXMOD
         DO J=1,NYMAX/NCPU
            T1(I,J)=ETA(I,J)
         ENDDO
      ENDDO

      CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)
      
      VOL=XL*YL*T1(1,1)
      
C-----END HERE

C---------------------
C     VOLUME FLUX
C---------------------

      DO I=1,NXMOD
         DO J=1,NYMAX/NCPU
            T2(I,J)=FETA(I,J,1)
         ENDDO
      ENDDO
      
      CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)
      
      FLUX=XL*YL*T2(1,1)

C-----END HERE

C--------------------------------------
C     POTENTIAL AND KINETIC ENERGY
C--------------------------------------

      DO I=1,NXMOD
         DO J=1,NYMAX/NCPU
            T1(I,J)=ETA(I,J)*ETA(I,J)
            T2(I,J)=VPS(I,J)*FETA(I,J,1)
         ENDDO
      ENDDO

      CALL FFTXY_PARA(T1,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)
      CALL FFTXY_PARA(T2,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
     1     NYMAX,NCPU,-1)

      PE=0.5*XL*YL*T1(1,1)
      KE=0.5*XL*YL*T2(1,1)
      ENER = PE / FR2 + KE
C      IF(ITIM.EQ.1) ENER0=ENER
C      IF ( ITIM .EQ. 0 ) ENER0 = ENER
      IF ( ITIM .EQ. 0 ) ENER0 = 1.

      WRITE(15,900) TIME,VOL,FLUX,PE,KE,ENER,ENER/ENER0
 900  FORMAT(7F11.7)

C-----END HERE

C---------------------------------
C     CHECK IF CODE BLOWS UP
C---------------------------------
      
C      IF(ABS(ENER-ENER0).GE.0.5*ENER0) IBLOWUP=1

C-----END HERE

      RETURN
      END

C=====SUBROUTINE CHECK END HERE
