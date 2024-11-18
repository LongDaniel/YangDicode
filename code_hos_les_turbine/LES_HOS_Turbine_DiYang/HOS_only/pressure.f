C==========================================================================
      SUBROUTINE PRESS1(NXMOD,NXMAX,NYMAX,P0,SIGMA,NSWAVEX,TIME,XP,PA,
     1     NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     STANDING WAVE SHAPE PRESSURE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA
      INTEGER NSWAVEX
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL TMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      DO I = 1, NXMOD
         TMP = P0 * SIN(SIGMA*TIME) * COS(NSWAVEX*XP(I))
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C
      RETURN
      END

C=====SUBROUTINE PRESS1 END HERE

C==========================================================================
      SUBROUTINE PRESS2(NXMOD,NXMAX,NYMAX,P0,SIGMA,NSWAVEX,TIME,XP,PA,
     1     NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     TRAVELLING WAVE SHAPE PRESSURE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA
      INTEGER NSWAVEX
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL TMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      DO I = 1, NXMOD
         TMP = P0 * SIN(NSWAVEX*XP(I)+SIGMA*TIME)
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS2 END HERE

C==========================================================================
      SUBROUTINE PRESS3(NXMOD,NXMAX,NYMAX,PEX,P0,SIGMA,AKAX,NSWAVEX,FR2,
     1     DELTA,TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     DELTA PRESSURE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2, DELTA, PEX
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI, PERIOD
      REAL TMP
      REAL AK, AA
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      PERIOD = 2 * PI / SIGMA
      AK = PEX * NSWAVEX
      AA = AKAX / AK
CC---STANDING WAVE
      DO I = 1, NXMOD 
C         IF ( TIME .LE. PERIOD / 2. 
C     1        .AND. TIME .GE. PERIOD / 2. - 2 * DELTA ) THEN
C            TMP = - AA / FR2 / SIGMA * 1. / 2. / DELTA 
C     1           * ( 1 - COS(PI/DELTA*TIME) ) * SIN(NSWAVEX*XP(I))
C         ELSE
C            TMP = 0
C         END IF
         IF ( TIME .LE. 0 ) THEN
            TMP = AKAX / FR2 / SIGMA * 1. / 2. / DELTA 
     1           * ( 1 - COS(PI/DELTA*TIME) ) * COS(NSWAVEX*XP(I)) 
         ELSE
            TMP = 0
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
CC--TWO STANDING WAVE--TRAVELING WAVE
C      DO I = 1, NXMOD
C         IF ( TIME .LE. 0 ) THEN
C            TMP = AKAX / FR2 / SIGMA * 1. / 2. / DELTA 
C     1           * ( 1 - COS(PI/DELTA*TIME) ) * COS(NSWAVEX*XP(I)) 
C         ELSE IF ( TIME .GE. PERIOD / 4. - DELTA 
C     1           .AND. TIME .LE. PERIOD / 4. ) THEN
C            TMP = - AKAX / FR2 / SIGMA * 1. / 2. / DELTA 
C     1           * ( 1 - COS(PI/DELTA*TIME) ) * SIN(NSWAVEX*XP(I))            
C         ELSE
C            TMP = 0
C         END IF
C         DO J = 1, NYMAX / NCPU
C            PA(I,J) = TMP
C         END DO
C      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS3 END HERE

C==========================================================================
      SUBROUTINE PRESS4(NXMOD,NXMAX,NYMAX,P0,SIGMA,NSWAVEX,TIME,XP,PA,
     1     NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     LONGUET-HIGGINS' PRESSURE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI
      REAL TMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
C
      DO I = 1, NXMOD
         IF ( TIME .LE. PI ) THEN
            TMP = P0 * SIN(TIME) * SIN(NSWAVEX*XP(I)+SIGMA*TIME)
         ELSE
            TMP = 0
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS4 END HERE

C==========================================================================
      SUBROUTINE PRESS5(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,N,
     1     TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     GENERAL STANDING WAVE GENERATOR
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2
      REAL N
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI
      REAL TMP
      REAL T0
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      T0 = PI / N / SIGMA
C
      DO I = 1, NXMOD
         IF ( TIME .LE. T0 ) THEN
            TMP = - N * AKAX / 4. / FR2 
     1           * ( ( N + 2. ) * SIN(SIGMA*((N+1.)*TIME-T0)) 
     1           - ( N - 2. ) * SIN(SIGMA*((N-1.)*TIME+T0)) )
     1           * COS(NSWAVEX*XP(I)) 
         ELSE
            TMP = 0
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS5 END HERE

C==========================================================================
      SUBROUTINE PRESS6(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,N,
     1     TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     GENERAL TREVELING WAVE GENERATOR
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2
      REAL N
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI
      REAL TMP
      REAL T0
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      T0 = PI / N / SIGMA
C
      DO I = 1, NXMOD
         IF ( TIME .LE. T0 ) THEN
            TMP = - N * AKAX / 4. / FR2 
     1           * ( ( N + 2. ) * SIN(SIGMA*((N+1.)*TIME+NSWAVEX*XP(I))) 
     1           - ( N - 2. ) * SIN(SIGMA*((N-1.)*TIME-NSWAVEX*XP(I))) )
         ELSE
            TMP = 0
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS6 END HERE

C==========================================================================
      SUBROUTINE PRESS7(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,DT,
     1     DELTA,TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     DELTA PRESSURE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2, DELTA
      REAL TIME, DT
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI, PERIOD
      REAL TMP, TMPT
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      PERIOD = 2 * PI / SIGMA
CC---STANDING WAVE
      DO I = 1, NXMOD
         IF ( TIME .LE. 0 ) THEN
C            TMP = AKAX / FR2 / SIGMA**2 * 1. / 2. / DELTA**2 * PI 
C     1           * SIN(PI/DELTA*TIME) * COS(NSWAVEX*XP(I)) 
            TMPT = ABS(TIME + DELTA)
            TMPT = TMPT / DT
            IF ( TMPT .GT. 2 ) THEN
               TMPT = 0
            ELSE IF ( TMPT .GT. 1 ) THEN
               TMPT = 2. / 3. * 1. / 4. * ( 2 - TMPT )**3
            ELSE
               TMPT = 2. / 3. 
     1              * ( 1 - 3. / 2. * TMPT**2 + 3. / 4. * TMPT**3 )
            END IF
            TMP = AKAX / FR2 / SIGMA**2 * TMPT * COS(NSWAVEX*XP(I)) / DT 
     1           / 2. 
         ELSE
            TMP = 0
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS7 END HERE

C==========================================================================
      SUBROUTINE PRESS8(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,N,
     1     DELTA,TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     DGD'S FIRST
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2, N, DELTA
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI, PERIOD
      REAL TMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      PERIOD = 2 * PI / SIGMA
C      WRITE(*,*) N, PERIOD, DELTA
CC---STANDING WAVE
      DO I = 1, NXMOD
         IF ( TIME .GE. 0 ) THEN
            TMP = P0 * SIN(NSWAVEX*XP(I)) 
     1           * ( 1 - EXP(-TIME**N/DELTA**N) )
C            TMP = P0 * SIN(NSWAVEX*XP(I))
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS8 END HERE

C==========================================================================
      SUBROUTINE PRESS9(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,
     1     BETA,TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     GRADUAL STANDING
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL      
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2, BETA
      REAL TIME
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI, PERIOD
      REAL TMP
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      PERIOD = 2 * PI / SIGMA
CC---STANDING WAVE
      DO I = 1, NXMOD
         TMP = P0 * BETA / SIGMA**2 / FR2 * EXP(-BETA*TIME)
     1        * ( - 2 * SIGMA * COS(SIGMA*TIME) 
     1        + BETA * SIN(SIGMA*TIME) ) * COS(NSWAVEX*XP(I))
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS9 END HERE

C==========================================================================
      SUBROUTINE PRESS10(NXMOD,NXMAX,NYMAX,P0,SIGMA,AKAX,NSWAVEX,FR2,
     1     BETA,DT,DELTA,TIME,XP,PA,NCPU)
CC--BY XIN GUO------------------------------------------------CC
C     GRADUAL TRAVELLING WAVE
CC------------------------------------------------------------CC
CC--IMPLICIT
      IMPLICIT NONE
CC--INCLUDE
      INCLUDE "mpif.h"
CC--PARALLEL
      INTEGER MYID, NUMPROCS, IERR
CC--INPUT
      INTEGER NXMOD, NXMAX, NYMAX, NCPU
      REAL P0, SIGMA, AKAX
      INTEGER NSWAVEX
      REAL FR2, BETA, DELTA
      REAL TIME, DT
      REAL XP(*)
CC--OUTPUT
      REAL PA(NXMAX,*)
CC--INTERNAL
      INTEGER I, J
      REAL PI, PERIOD
      REAL TMP, TMPT
CC--START WORKING
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
CC
      PI = ACOS(-1.)
      PERIOD = 2 * PI / SIGMA
CC---TRAVELING WAVE
      DO I = 1, NXMOD
         IF ( TIME .LE. 0 ) THEN
C            TMP = AKAX / FR2 / SIGMA**2 * 1. / 2. / DELTA 
C     1           * ( 1 - COS(PI/DELTA*TIME) ) * SIN(NSWAVEX*XP(I)) 
C     1           * BETA
            TMPT = ABS(TIME + DELTA)
            TMPT = TMPT / DT
            IF ( TMPT .GT. 2 ) THEN
               TMPT = 0
            ELSE IF ( TMPT .GT. 1 ) THEN
               TMPT = 2. / 3. * 1. / 4. * ( 2 - TMPT )**3
            ELSE
               TMPT = 2. / 3. 
     1              * ( 1 - 3. / 2. * TMPT**2 + 3. / 4. * TMPT**3 )
            END IF
            TMP = AKAX / FR2 / SIGMA**2 * TMPT * SIN(NSWAVEX*XP(I)) 
     1           * BETA                        
         ELSE
            TMP = AKAX / FR2 / SIGMA**2 * ( 2 * BETA * SIGMA 
     1           * EXP(-BETA*TIME) * COS(NSWAVEX*XP(I)+SIGMA*TIME) 
     1           - BETA**2 * EXP(-BETA*TIME) 
     1           * SIN(NSWAVEX*XP(I)+SIGMA*TIME) )
         END IF
         DO J = 1, NYMAX / NCPU
            PA(I,J) = TMP
         END DO
      END DO
C     
      RETURN
      END

C=====SUBROUTINE PRESS10 END HERE
