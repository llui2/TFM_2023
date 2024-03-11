C     OBSERVABLES
C     0!
C     Lluís Torres 
C     TFM
C     FORTRAN 95

      PROGRAM OBSERVABLES

      USE MODEL

C-----(SYSTEM)------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     +1 -1 EDGES RATIO (1 => ALL +1), (0 => ALL -1)
      REAL*8 p
C     TEMPERATURE (TEMP := k_B·T)
      REAL*8 TEMP
C     p-LIST, TEMP-LIST
      INTEGER p_SIZE, TEMP_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:), TEMP_LIST(:)
C-----(SIMULATION)---------------------------------------------
C     SAMPLE SIZE (# OF SPIN CONFIGURATIONS)
      INTEGER C
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SEED NUMBER, INITIAL SEED NUMBER
      INTEGER SEED, SEEDini
      PARAMETER(SEEDini = 100)
C     SIMULATION VARIABLES
      INTEGER, ALLOCATABLE:: S(:)
C-----(SPIN CONFIGURATION SAVING VARIABLES)-------------------
C     STORE SPIN CONFIGURATION AS N/zip_size INTEGERS
      INTEGER zip_size
      INTEGER, ALLOCATABLE:: bin(:)
      INTEGER, ALLOCATABLE:: decimal(:)
      INTEGER, ALLOCATABLE:: array(:)
C-----(OBSERVABLES)-------------------------------------------
C     MAGNETIZATION (Z COMPONENT)
      REAL*8 MZ
C-----(DUMMY)-------------------------------------------------
      INTEGER ITEMP, Ip, IC
      CHARACTER(4) str
      CHARACTER(3) str1, str2, str3, str4
      REAL*8 H !TRANSVERSE FIELD (EQUAL TO 0 IN CLASSIC MODEL)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER TAU !NOT USED IN THIS PROGRAM
      INTEGER H_SIZE !NOT USED IN THIS PROGRAM
      REAL*8,ALLOCATABLE:: H_LIST(:) !NOT USED IN THIS PROGRAM
      INTEGER R !NOT USED IN THIS PROGRAM
      INTEGER MCINI !NOT USED IN THIS PROGRAM
      INTEGER SC !NOT USED IN THIS PROGRAM

C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      PRINT*, 'OBSERVABLES'

C***********************************************************************
C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,MCINI,NSEEDS,SC,zip_size,TAU)
      M = z*N/2
C     ALLOCATION
      ALLOCATE(decimal(1:N/zip_size))
      ALLOCATE(array(1:N))
      ALLOCATE(bin(1:N))
      ALLOCATE(S(1:N))
C***********************************************************************
      CALL SYSTEM('mkdir -p results/data/')
      OPEN(UNIT=1,FILE='results/data/mz.dat')
300   FORMAT(A,4X,A,7X,A,7X,A)
      WRITE(1,300) '#TEMP', 'Γ', 'p', '<|mz|>'
C***********************************************************************
C     THERE IS NO TRANSVERSE FIELD
      H = 0.00
      WRITE(str,'(f4.2)') H
      str2 = str(1:1)//str(3:4)
C***********************************************************************

C     FOR ALL TEMP VALUES
      DO ITEMP = 1,TEMP_SIZE
      TEMP = TEMP_LIST(ITEMP)
      WRITE(str,'(f4.2)') TEMP
      str1 = str(1:1)//str(3:4)

C     FOR ALL p VALUES
      DO Ip = 1,p_SIZE
      p = p_LIST(IP)
      WRITE(str,'(f4.2)') p
      str3 = str(1:1)//str(3:4)

C***********************************************************************
C     RESET OBSERVABLES
      MZ = 0.D0
C***********************************************************************

C     FOR ALL SEEDS
      DO SEED = SEEDini,SEEDini+NSEEDS-1
      WRITE(str4,'(i3)') SEED

C***********************************************************************
C     READ SPIN CONFIGURATIONS FROM FILE
      OPEN(UNIT=3,FILE='results/sample/T'//str1//'_Γ'//str2//
     .'/S_'//str3//'_'//str4//'.bin',FORM='UNFORMATTED')
C***********************************************************************
      DO IC = 1,C
            READ(3) decimal
            CALL DEC2BIN(N,zip_size,bin,decimal)
            CALL BIN2ARRAY(N,bin,array)
            S = array
C           OBSERVABLES
            MZ = MZ + ABS(MAGNET_Z(N,S))
      END DO !IC
C***********************************************************************
      CLOSE(3)
C***********************************************************************

      END DO !SEED

C***********************************************************************
C     SAVE OBSERVABLES     
200   FORMAT(F4.2,4X,F4.2,4X,F4.2,4X,F13.10)
      WRITE(1,200) TEMP, H, p, MZ/(C*NSEEDS)
C***********************************************************************

      END DO !Ip
      END DO !ITEMP

C***********************************************************************
      CLOSE(1)
C***********************************************************************

      END PROGRAM OBSERVABLES