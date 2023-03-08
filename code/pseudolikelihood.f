C     RECONSTRUCTION
C     0!
C     Lluís Torres 
C     TFM
C     FORTRAN 77

      PROGRAM RECONSTRUCTION

      USE MODEL

C-----(SYSTEM)------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     TROTTER INDEX
      INTEGER R
C     +1 -1 EDGES RATIO (1 => ALL +1), (0 => ALL -1)
      REAL*8 p
C     TEMPERATURE (TEMP := k_B·T)
      REAL*8 TEMP
C     H (:= Γ = TRANSVERSE FIELD)
      REAL*8 H
C     p-LIST, TEMP-LIST, Γ-LIST
      INTEGER p_SIZE, TEMP_SIZE, H_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:), TEMP_LIST(:), H_LIST(:)
C     ORIGINAL SYSTEM ARRAYS
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ_0(:)
C-----(SIMULATION)---------------------------------------------
C     SAMPLE SIZE (# OF SPIN CONFIGURATIONS)
      INTEGER C
C     TOTAL MONTE-CARLO STEPS (MCS) FOR RECONSTRUCTION
      INTEGER TAU
C     FICTICIOUS TEMPERATURE
      REAL*8 TEMP_F
C     FICTICIOUS TEMPERATURE STEP
      REAL*8 TF_STEP
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SEED NUMBER, INITIAL SEED NUMBER
      INTEGER SEED, SEEDini
      PARAMETER(SEEDini = 100)
C     RANDOM NUMBER GENERATOR
      EXTERNAL r1279
C     ESTIMATE TIME VARIABLES
      REAL*4 TIME1, TIME2, time
C     SIMULATION VARIABLES
      INTEGER, ALLOCATABLE:: S_SET(:,:,:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      LOGICAL valid
      REAL*8 DPL
      REAL*8 PL
C-----(SPIN CONFIGURATION SAVING VARIABLES)-------------------
C     STORE SPIN CONFIGURATION AS N/zip_size INTEGERS
      INTEGER zip_size
      INTEGER, ALLOCATABLE:: bin(:)
      INTEGER, ALLOCATABLE:: decimal(:)
      INTEGER, ALLOCATABLE:: array(:)
C-----(DUMMY)-------------------------------------------------
      INTEGER ITEMP, IH, Ip, IC, IR
      INTEGER IMC
      CHARACTER(4) str
      CHARACTER(3) str1, str2, str3, str4
      !CHARACTER(2) strN
      INTEGER SC !NOT USED IN THIS PROGRAM
      
C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      PRINT*, 'SYSTEM RECONSTRUCTION'

C***********************************************************************
C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,NSEEDS,SC,zip_size,TAU)
      M = z*N/2
C     ALLOCATION
      ALLOCATE(decimal(1:N/zip_size))
      ALLOCATE(array(1:N))
      ALLOCATE(bin(1:N))
      ALLOCATE(S_SET(1:C,1:R,1:N))
C***********************************************************************
      CALL CPU_TIME(TIME1)
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

C     FOR ALL SEEDS
      DO SEED = SEEDini,SEEDini+NSEEDS-1
      WRITE(str4,'(i3)') SEED
      CALL setr1279(SEED)

C***********************************************************************
C     ORIGINAL RANDOM SYSTEM (GRAPH+COUPLINGS)
      CALL IRS(N,M,p,NBR_0,INBR_0,JJ_0)
C***********************************************************************
C     READ THE SAMPLE
      OPEN(UNIT=1,FILE='results/sample/T'//str1//'_Γ'//str2//
     .'/S_'//str3//'_'//str4//'.bin',FORM='UNFORMATTED')
      DO IC = 1,2*C
            DO IR=1,R
                  READ(1) decimal
                  CALL DEC2BIN(N,zip_size,bin,decimal)
                  CALL BIN2ARRAY(N,bin,array)
                  S_SET(IC,IR,:) = array
            END DO
      END DO
      CLOSE(1)
C***********************************************************************
C     INITIAL FICTICIOUS TEMPERATURE  
      TEMP_F = -LOG(0.5d0*(1+TANH(1.D0/TEMP)))/z
      TF_STEP = TEMP_F/TAU
C***********************************************************************
C     INITIAL RANDOM SYSTEM (GRAPH+COUPLINGS)
      CALL setr1279(999)
      CALL IRS(N,M,p,NBR,INBR,JJ)
C***********************************************************************
      !CALL CLASS_FUNCTION(zmax,TEMP,funct)
      !CALL CLASS_LAMBDA(N,C,S_SET,NBR,JJ,LAMBDA)
C***********************************************************************
C     INITIAL PSEUDOLIKELIHOOD
      PL = PSEUDO(N,R,C,S_SET,TEMP,NBR,JJ)
C***********************************************************************







      END DO !SEED
      END DO !Ip
      END DO !IH
      END DO !ITEMP

      END PROGRAM RECONSTRUCTION