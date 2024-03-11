C     PERFORMANCE
C     0!
C     Lluís Torres 
C     TFM
C     FORTRAN 95

      PROGRAM PERFORMANCE

      USE MODEL

C-----(SYSTEM)------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     +1 -1 EDGES RATIO (1 => ALL +1), (0 => ALL -1)
      REAL*8 p
C     TEMPERATURE (TEMP := k_B·T)
      REAL*8 TEMP
C     p-LIST, TEMP-LIST, Γ-LIST
      INTEGER p_SIZE, TEMP_SIZE, H_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:), TEMP_LIST(:), H_LIST(:)
C-----(SIMULATION)---------------------------------------------
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SEED NUMBER
      INTEGER SEED
C     PERFORMANCE ESTIMATOR
      REAL*8 gMEAN, gERROR, g
C-----(DUMMY)-------------------------------------------------
      INTEGER ITEMP, Ip, iSEED
      CHARACTER(4) str
      CHARACTER(3) str1, str2, str3
      INTEGER zmax
      REAL*8 H !TRANSVERSE FIELD (EQUAL TO 0 IN CLASSIC MODEL)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER SC !NOT USED IN THIS PROGRAM
      INTEGER R !NOT USED IN THIS PROGRAM
      INTEGER TAU !NOT USED IN THIS PROGRAM
      INTEGER C !NOT USED IN THIS PROGRAM
      INTEGER zip_size !NOT USED IN THIS PROGRAM
      INTEGER MCINI !NOT USED IN THIS PROGRAM

C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      PRINT*, 'PERFORMANCE ESTIMATION'

C***********************************************************************
C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,MCINI,NSEEDS,SC,zip_size,TAU)
      M = z*N/2
      zmax = N-1
C***********************************************************************
      CALL SYSTEM('mkdir -p results/data/')
      OPEN(UNIT=1,FILE='results/data/gamma.dat')
300   FORMAT(A,4X,A,7X,A,7X,A,16X,A)
      WRITE(1,300) '#TEMP', 'Γ', 'p', 'γ', 'Δγ'
C***********************************************************************

C     FOR ALL TEMP VALUES
      DO ITEMP = 1,TEMP_SIZE
      TEMP = TEMP_LIST(ITEMP)
      WRITE(str,'(f4.2)') TEMP
      str1 = str(1:1)//str(3:4)

C     FOR ALL Γ VALUES      
      DO IH = 1,H_SIZE
      H = H_LIST(IH)
      WRITE(str,'(f4.2)') H
      str2 = str(1:1)//str(3:4)

C     FOR ALL p VALUES
      DO Ip = 1,p_SIZE
      p = p_LIST(IP)
      WRITE(str,'(f4.2)') p
      str3 = str(1:1)//str(3:4)

C***********************************************************************
C     GAMMA MEAN
      gMEAN = 0.d0
      OPEN(UNIT=2,FILE='results/accuracy/T'//str1//'_Γ'//str2//
     .'/g_'//str3//'.dat')
      DO ISEED = 1,NSEEDS
            READ(2,*) SEED, g
            gMEAN = gMEAN + g
      END DO !ISEED
      gMEAN = gMEAN / NSEEDS
      CLOSE(2)
C***********************************************************************

C***********************************************************************
C     GAMMA ERROR
      gERROR = 0.d0
      OPEN(UNIT=2,FILE='results/accuracy/T'//str1//'_Γ'//str2//
     .'/g_'//str3//'.dat')
      DO ISEED = 1,NSEEDS
            READ(2,*) SEED, g
            gERROR = gERROR + (g - gMEAN)**2
      END DO !ISEED
      gERROR = gERROR/SQRT(NSEEDS*(NSEEDS - 1.d0))
      CLOSE(2)
C***********************************************************************

C     SAVE RESULTS
200   FORMAT(F4.2,4X,F4.2,4X,F4.2,4X,F13.10,4X,F13.10)
      WRITE(1,200) TEMP, H, p, gMEAN, gERROR

      END DO !Ip
      END DO !IH
      END DO !ITEMP

      CLOSE(1)

      END PROGRAM PERFORMANCE