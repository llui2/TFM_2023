C     OBSERVABLES
C     0!
C     Lluís Torres 
C     TFM
C     FORTRAN 77

      PROGRAM OBSERVABLES

      USE MODEL

C-----(SYSTEM)------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     NUMBER OF REPLICAS 
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
C     PRINCIPAL ARRAYS
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
C-----(SIMULATION)---------------------------------------------
C     SAMPLE SIZE (# OF SPIN CONFIGURATIONS)
      INTEGER C
C     TOTAL MONTE-CARLO STEPS (MCS)
      INTEGER MCTOT
C     MCS TILL WE CONSIDER EQUILIBRIUM
      INTEGER MCINI
C     SAVE SPIN CONFIGURATIONS EVERY SC (MCS)
      INTEGER SC
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SEED NUMBER, INITIAL SEED NUMBER
      INTEGER SEED, SEEDini
C     SIMULATION VARIABLES
      INTEGER, ALLOCATABLE:: S(:,:)
C-----(SPIN CONFIGURATION SAVING VARIABLES)-------------------
C     STORE SPIN CONFIGURATION AS N/zip_size INTEGERS
      INTEGER zip_size
      INTEGER, ALLOCATABLE:: bin(:)
      INTEGER, ALLOCATABLE:: decimal(:)
      INTEGER, ALLOCATABLE:: array(:)
C-----(OBSERVABLES)-------------------------------------------
C     MAGNETIZATION (Z COMPONENT)
      REAL*8 MZ
C     MAGNETIZATION (X COMPONENT)
      REAL*8 MX
C-----(DUMMY)-------------------------------------------------
      INTEGER ITEMP, IH, Ip, IC
      CHARACTER(4) str
      CHARACTER(3) str1, str2, str3, str4
      
C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      PRINT*, 'OBSERVABLES'

C***********************************************************************
C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,NSEEDS,SC,zip_size)
      M = z*N/2
      MCTOT = 3*C*SC/2
      MCINI = MCTOT/3
C     ALLOCATION
      ALLOCATE(decimal(1:N/zip_size))
      ALLOCATE(array(1:N))
      ALLOCATE(bin(1:N))
      ALLOCATE(S(1:R,1:N))
C***********************************************************************
C     INITIAL SEED NUMBER
      SEEDini = 100
      CALL CPU_TIME(TIME1)
C***********************************************************************
      CALL SYSTEM('mkdir -p results/observables/')
      OPEN(UNIT=1,FILE='results/observables/MZ.dat')
300   FORMAT(A,4X,A,7X,A,7X,A)
      WRITE(1,300) '#TEMP', 'H', 'p', '<|MZ|>'
      OPEN(UNIT=2,FILE='results/observables/MX.dat')
      WRITE(2,300) '#TEMP', 'H', 'p', '<MX>'
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
C     RESET OBSERVABLES
      MZ = 0.D0
      MX = 0.D0
C***********************************************************************

C     FOR ALL SEEDS
      DO SEED = SEEDini,SEEDini+NSEEDS-1
      WRITE(str4,'(i3)') SEED
      CALL setr1279(SEED)

C***********************************************************************
C     RECALL RANDOM SYSTEM (GRAPH+COUPLINGS)
      CALL IRS(N,M,p,NBR,INBR,JJ)
C***********************************************************************
C     READ SPIN CONFIGURATIONS FROM FILE
      OPEN(UNIT=3,FILE='results/sample/T'//str1//'_Γ'//str2//
     *'/S_'//str3//'_'//str4//'.bin',FORM='UNFORMATTED')
C***********************************************************************
      DO IC = 1,2*C
            DO i=1,R
                  READ(3) decimal
                  CALL DEC2BIN(N,zip_size,bin,decimal)
                  CALL BIN2ARRAY(N,bin,array)
                  S(i,:) = array
                  ! READ(3,*) bin
                  ! CALL BIN2ARRAY(N,bin,S(i,:))
            END DO
C           OBSERVABLES
            MZ = MZ + ABS(MAGNET_Z(N,R,S))
            MX = MX + MAGNET_X(N,R,S,TEMP,H)
      END DO !IC
C***********************************************************************
      CLOSE(3)
C***********************************************************************
C     DELLOCATE ARRAYS
      DO i = 1,N
            DEALLOCATE(NBR(i)%v)
            DEALLOCATE(INBR(i)%v)
            DEALLOCATE(JJ(i)%v)
      END DO
      DEALLOCATE(NBR)
      DEALLOCATE(INBR)
      DEALLOCATE(JJ)
C***********************************************************************

      END DO !SEED

C***********************************************************************
C     SAVE OBSERVABLES     
200   FORMAT(F4.2,4X,F4.2,4X,F4.2,4X,F13.10)
      WRITE(1,200) TEMP, H, p, MZ/(2*C*NSEEDS)
      WRITE(2,200) TEMP, H, p, MX/(2*C*NSEEDS)
C***********************************************************************

      END DO !Ip
      END DO !IH
      END DO !ITEMP

C***********************************************************************
      CLOSE(1)
      CLOSE(2)
C***********************************************************************

      END PROGRAM OBSERVABLES