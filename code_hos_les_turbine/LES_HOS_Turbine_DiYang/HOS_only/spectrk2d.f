      SUBROUTINE SPECTK2D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,NMAX,DWK,
     1           WORK,TRIGSX,TRIGSY,IFAX,FT,F,SPF)
C
C Compute 2-D spectrum on x-y plane 
C
      REAL F(NXMAX,*),FT(NXMAX,*)
      REAL SPF(*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)
C
      NMAX=NXMOD/2-1
      DWK=PEX
C
C-Compute mean 
      FM=0.
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
 10   FM=FM+F(I,J)
C
      FM=FM/NXMOD/NYMOD
C     
C-Compute fluctuations
      DO 20 I=1,NXMOD
      DO 20 J=1,NYMOD
 20   FT(I,J)=F(I,J)-FM
C     
C-Compute spectrum
CC      CALL FFTXY_PARA(FT,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
CC     1     NYMAX,1,-1)
      CALL FFTXY(FT,WORK,TRIGSX,IFAX,NXMOD,NYMOD,NXMAX,-1)
C
      DO 30 N=1,NMAX
 30   SPF(N)=0.
C
      DO 40 L=1,NXMOD-1,2
       LP1=L+1
       WKX=(L-1)/2*PEX
       DO 40 M=1,NYMOD-1,2
        WKY=(M-1)/2*PEY
        WAV=SQRT(WKX**2+WKY**2)
        DO 45 N=1,NMAX
        IF(WAV.GE.((N-0.5)*DWK).AND.WAV.LT.((N+0.5)*DWK))THEN
           GOTO 50
        ENDIF
 45     CONTINUE
        GOTO 40
C
 50     CONTINUE
        MP1=M+1
        FTN=FT(L,M)**2+FT(L,MP1)**2+FT(LP1,M)**2+FT(LP1,MP1)**2
        CON=4.
        IF(L.EQ.1.OR.M.EQ.1) CON=2.
        IF(L.EQ.1.AND.M.EQ.1) CON=1.
        SPF(N)=SPF(N)+CON*FTN/DWK
 40   CONTINUE
C
      RETURN
      END


C----------------------------------------------------------------------C
      SUBROUTINE SPECTKX(NXMOD,NYMOD,NXMAX,NYMAX,PEX,NMAX,DWKX,
     1           WORK,TRIGSX,TRIGSY,IFAX,FT,F,SXF)
C
C Compute 1-D spectrum in x-direction
C Average in y-direction 
C
      REAL F(NXMAX,*),FT(NXMAX,*)
      REAL SXF(*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)
C
      NMAX=NXMOD/2-1
      DWKX=PEX
C
C-Compute mean 
      FM=0.
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
 10   FM=FM+F(I,J)
C
      FM=FM/NXMOD/NYMOD
C     
C-Compute fluctuations
      DO 20 I=1,NXMOD
      DO 20 J=1,NYMOD
 20   FT(I,J)=F(I,J)-FM
C     
C-Compute spectrum
CC      CALL FFTX_PARA(FT,WORK,TRIGSX,IFAX,NXMOD,NXMAX,NYMAX, 1,-1)
      CALL FFTX(FT,WORK,TRIGSX,IFAX,NXMOD,NYMOD,NXMAX,-1)
C
      DO 30 N=1,NMAX
 30   SXF(N)=0.
C
      DO 40 L=1,NXMOD-1,2
       LP1=L+1
       DO 40 M=1,NYMOD
        N=(L+1)/2  
        FTN=FT(L,M)**2+FT(LP1,M)**2
        CON=2.
        IF(L.EQ.1) CON=1.
        SXF(N)=SXF(N)+CON*FTN/DWKX/NYMOD
 40   CONTINUE
C
      RETURN
      END

C----------------------------------------------------------------------C
      SUBROUTINE SPECTKY(NXMOD,NYMOD,NXMAX,NYMAX,PEY,N2MAX,DWKY,
     1           WORK,TRIGSX,TRIGSY,IFAX,FT,F,SYF)
C
C Compute 1-D spectrum in y-direction
C Average in x-direction 
C
      REAL F(NXMAX,*),FT(NXMAX,*)
      REAL SYF(*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)
C
      N2MAX=NYMOD/2-1
      DWKY=PEY
C
C-Compute mean 
      FM=0.
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
 10   FM=FM+F(I,J)
C
      FM=FM/NXMOD/NYMOD
C     
C-Compute fluctuations
      DO 20 I=1,NXMOD
      DO 20 J=1,NYMOD
 20   FT(I,J)=F(I,J)-FM
C     
C-Compute spectrum
CC      CALL FFTY_PARA(FT,WORK,TRIGSY,IFAX,NYMOD,NXMAX,NYMAX,1,-1)
      CALL FFTY(FT,WORK,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,-1)
C
      DO 30 N=1,NMAX
 30   SYF(N)=0.
C
      DO 40 L=1,NXMOD-1
       DO 40 M=1,NYMOD-1,2
        MP1=M+1 
        WKY=(M-1)/2*PEY
        N=(M+1)/2  
        FTN=FT(L,M)**2+FT(L,MP1)**2
        CON=2.
        IF(M.EQ.1) CON=1.
        SYF(N)=SYF(N)+CON*FTN/DWKY/NXMOD
 40   CONTINUE
C
      RETURN
      END

CC----------------------------------------------------------
      SUBROUTINE SPECTK3D(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,DWK,
     1           WORK,TRIGSX,TRIGSY,IFAX,FT,F,SPF)
C
C Compute 2-D spectrum on x-y plane, KX AND KY
C
      REAL F(NXMAX,*),FT(NXMAX,*)
      REAL SPF(NXMAX,*)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)
C
C      NMAX=NXMOD/2-1
      DWK=PEX
C
C-Compute mean 
      FM=0.
      DO 10 I=1,NXMOD
      DO 10 J=1,NYMOD
 10   FM=FM+F(I,J)
C
      FM=FM/NXMOD/NYMOD
C     
C-Compute fluctuations
      DO 20 I=1,NXMOD
      DO 20 J=1,NYMOD
 20   FT(I,J)=F(I,J)-FM
C     
C-Compute spectrum

CC--TEST
C      WRITE(*,*) 'BEFORE FFT'
C      STOP
CC---@@@@@@@@@@@@@@@@@@@@@@2

CC      CALL FFTXY_PARA(FT,WORK,TRIGSX,TRIGSY,IFAX,NXMOD,NYMOD,NXMAX,
CC     1     NYMAX,1,-1)
      CALL FFTXY(FT,WORK,TRIGSX,IFAX,NXMOD,NYMOD,NXMAX,-1)
C
CC--TEST
C      WRITE(*,*) 'AFTER FFT'
C      STOP
CC---@@@@@@@@@@@@@@@@@@@@@@2

      DO 30 L=1,NXMAX/2-1
         DO 30 M=1,NYMAX/2-1
 30   SPF(L,M)=0.
C
      DO 40 L=1,NXMOD-1,2
         LP1=L+1
         WKX=(L-1)/2*PEX
         I = ( L + 1) / 2
         DO 40 M=1,NYMOD-1,2
            WKY=(M-1)/2*PEY
            J = ( M + 1) / 2
            MP1=M+1
            FTN=FT(L,M)**2+FT(L,MP1)**2+FT(LP1,M)**2+FT(LP1,MP1)**2
            CON=4.
            IF(L.EQ.1.OR.M.EQ.1) CON=2.
            IF(L.EQ.1.AND.M.EQ.1) CON=1.
            SPF(I,J)=SPF(I,J)+CON*FTN/DWK
 40   CONTINUE
C
      RETURN
      END

CC=================================================================





C===================================================================
      SUBROUTINE FFTXY(F,WORK,TRIGS,IFAX,NXMOD,NYMOD,NXMAX,ISGN)

C     SERIAL FFTXY

C
      REAL F(NXMAX,*),WORK(*),TRIGS(*)
      INTEGER IFAX(*)
C
      IF(ISGN.EQ.1) GO TO 50
C
      DO 10 J=NYMOD+1,NYMOD+2
      DO 10 I=1,NXMOD
 10   F(I,J)=0.
      DO 20 J=1,NYMOD+2
      DO 20 I=NXMOD+1,NXMOD+2
 20   F(I,J)=0.
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD+2,ISGN)
      GOTO 99
C
 50   CONTINUE
C
      CALL FFTFAX(NYMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,NXMAX,1,NYMOD,NXMOD+2,ISGN)
C
      CALL FFTFAX(NXMOD,IFAX,TRIGS)
      CALL FFT991(F,WORK,TRIGS,IFAX,1,NXMAX,NXMOD,NYMOD,ISGN)
C
 99   CONTINUE
      RETURN
      END

C=====SUBROUTINE FFTXY END HERE

C======================================================================
      SUBROUTINE SPECTK2Dt(NXMOD,NYMOD,NXMAX,NYMAX,PEX,PEY,F,WORK,IFAX,
     1     TRIGSX,TRIGSY,TIME,NCPU,IFILE)

C     COMPUTE THE 1D SPECTRUM OF 2D FUNCTION .VS. WAVENUMBER

      IMPLICIT NONE
      
      INTEGER I,J,K,JS,N,L,M,LP1,MP1
      INTEGER NXMOD,NYMOD,NXMAX,NYMAX,NCPU
      INTEGER NUMSEND,NUMRECV
      INTEGER ITAG,IFILE

      INTEGER NMAX,NXYMOD

      REAL TIME,PI,PEX,PEY,DWK,FM,FQ,WKX,WKY,WAV,FKN,CON
      
      INCLUDE "mpif.h"

      INTEGER MYID,NUMPROCS,IERR,STAT(MPI_STATUS_SIZE)

      REAL F(NXMAX,*),FSP(NXMAX),FK(NXMAX,NYMAX)
      INTEGER NWK(NXMAX)
      REAL WORK(*),TRIGSX(*),TRIGSY(*)
      INTEGER IFAX(*)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C-----------------------------------
C     TRANSFER DATA INTO MYID=0
C-----------------------------------

      NUMSEND=NXMAX*NYMAX/NCPU
      NUMRECV=NUMSEND
      
      IF(NCPU.EQ.1) THEN
         
         DO I=1,NXMOD
            DO J=1,NYMOD
               FK(I,J)=F(I,J)
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
                  FK(I,J)=F(I,J)
               ENDDO
            ENDDO
            
            DO K=1,NCPU-1
               JS=K*NYMAX/NCPU+1
               ITAG=K
               CALL MPI_RECV(FK(1,JS),NUMRECV,MPI_DOUBLE_PRECISION,K,
     1              ITAG,MPI_COMM_WORLD,STAT,IERR)
            ENDDO
      
         ENDIF

      ENDIF

C-----END HERE

C************************************************************
C     THE FOLLOWING COMPUTATIONS ARE DONE ONLY ON MYID=0
C************************************************************
      
      IF(MYID.EQ.0) THEN

C------------------------------------
C     COMPUTE PLANE AVERAGE OF F
C------------------------------------

         NMAX=NXMOD/2-1
         DWK=PEX
         PI=ACOS(-1.)
         NXYMOD=NXMOD*NYMOD
         
         FM=0.
         DO I=1,NXMOD
            DO J=1,NYMOD
               FM=FM+FK(I,J)
            ENDDO
         ENDDO
         FM=FM/NXYMOD

C-----END HERE

C-----------------------------------------------------
C     COMPUTE FLUCUATION OF F AND ITS MEAN SQUARE
C-----------------------------------------------------

         FQ=0.
         DO I=1,NXMOD
            DO J=1,NYMOD
               FK(I,J)=FK(I,J)-FM
               FQ=FQ+FK(I,J)**2
            ENDDO
         ENDDO
         FQ=FQ/NXYMOD
         
         CALL FFTXY(FK,WORK,TRIGSX,IFAX,NXMOD,NYMOD,NXMAX,-1)
         
         DO N=1,NMAX
            FSP(N)=0.
            NWK(N)=0
         ENDDO
         
         DO L=1,NXMOD-1,2
            LP1=L+1
            WKX=(L-1)/2*PEX
            DO M=1,NYMOD-1,2
               WKY=(M-1)/2*PEY
               WAV=SQRT(WKX**2+WKY**2)
               DO N=1,NMAX
                  IF(WAV.GE.((N-0.5)*DWK).AND.WAV.LT.((N+0.5)*DWK))THEN
                     NWK(N)=NWK(N)+1
                     GOTO 35
                  ENDIF
               ENDDO
 35            CONTINUE
               MP1=M+1
               FKN=FK(L,M)**2+FK(L,MP1)**2+FK(LP1,M)**2+FK(LP1,MP1)**2
               CON=4.
               IF(L.EQ.1.OR.M.EQ.1) CON=2.
               IF(L.EQ.1.AND.M.EQ.1) CON=1.
               FSP(N)=FSP(N)+0.5*CON*FKN
            ENDDO
         ENDDO
         
         DO N=1,NMAX
            FSP(N)=FSP(N)/DWK
         ENDDO
         
         CALL FFTXY(FK,WORK,TRIGSX,IFAX,NXMOD,NYMOD,NXMAX,1)
         
         WRITE(IFILE,900) TIME
 900     FORMAT(' ZONE T="',F8.3,'"')
         DO N=1,NMAX
            WAV=N*DWK  
            WRITE(IFILE,1001) WAV,FSP(N),FSP(N)/(0.5*FQ),NWK(N)
 1001       FORMAT(F11.5,2E12.4,I4)
         ENDDO
         
C-----END HERE

      ENDIF

C*****END HERE

      RETURN
      END

C=====SUBROUTINE SPECTK2D END HERE
