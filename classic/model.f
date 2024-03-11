C     SPIN MODEL MODULE
C     0!
C     Lluís Torres 
C     TFG
C     FORTRAN 95

      MODULE MODEL 
C     VIANA-BRAY SPIN-GLASS MODEL WITH DISCRETE DISTRIBUTION OF THE COUPLING

C     MULTI ARRAY TYPE
      TYPE :: MULTI_ARRAY
      INTEGER,ALLOCATABLE :: v(:)
      END TYPE MULTI_ARRAY

      CONTAINS

C-----------------------------------------------------------------------
C     READ INPUT FILE
      SUBROUTINE READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,MCINI,NSEEDS,SC,zip_size,TAU)

      INTEGER N,z
      INTEGER R
      INTEGER TEMP_SIZE
      REAL*8,ALLOCATABLE:: TEMP_LIST(:)
      INTEGER H_SIZE
      REAL*8,ALLOCATABLE:: H_LIST(:)
      INTEGER p_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:)
      INTEGER C
      INTEGER MCINI
      INTEGER NSEEDS
      INTEGER SC
      INTEGER zip_size
      INTEGER TAU

      OPEN(UNIT=0,FILE="input.txt")
      
      READ(0,*)
      READ(0,*) N,z
      READ(0,*)
      READ(0,*) R
      READ(0,*)
      READ(0,*) TEMP_SIZE
      ALLOCATE(TEMP_LIST(1:TEMP_SIZE))
      READ(0,*) TEMP_LIST
      READ(0,*)
      READ(0,*) H_SIZE
      ALLOCATE(H_LIST(1:H_SIZE))
      READ(0,*) H_LIST
      READ(0,*)
      READ(0,*) p_SIZE
      ALLOCATE(p_LIST(1:p_SIZE))
      READ(0,*) p_LIST
      READ(0,*)
      READ(0,*) C
      READ(0,*)
      READ(0,*) MCINI
      READ(0,*)
      READ(0,*) NSEEDS
      READ(0,*)
      READ(0,*) SC
      READ(0,*)
      READ(0,*) zip_size
      READ(0,*)
      READ(0,*) TAU
      
      CLOSE(0)

      RETURN
      END SUBROUTINE READ_INPUT
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     RANDOM ERDÖS-RÉNYI GRAPH WITH RANDOM COUPLINGS GENERATOR
      SUBROUTINE IRS(N,M,p,NBR,INBR,JJ) ! Initial Random System (Graph and Couplings)
C     THIS SUBROUTINE GENERATES A RANDOM ERDÖS-RÉNYI GRAPH WITH p*M EDGES
C     WITH A WEIGHT OF 1 AND (1-p)*M EDGES WITH A WEIGHT OF 1 AND SAVES
C     IT IN THE NBR, INBR AND JJ ARRAYS.

      INTEGER N, M
      REAL*8 p
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      INTEGER i, j, k
      EXTERNAL r1279
      INTEGER edges_p, edges_n

      edges_p = INT(p*M)
      edges_n = INT((1-p)*M)

      IF (edges_p+edges_n.NE.M) THEN
            PRINT*, 'ERROR: p VALUE NOT VALID'
            STOP
      END IF
      IF (p.GT.1.OR.p.LT.0) THEN
            PRINT*, 'ERROR: p VALUE NOT VALID'
            STOP
      END IF

      ALLOCATE(NBR(N))
      ALLOCATE(INBR(N))
      ALLOCATE(JJ(N))

      DO i = 1,N
            ALLOCATE(NBR(i)%v(0))
            ALLOCATE(JJ(i)%v(0))
      END DO

C     GENERATE M/2 EDGES OF WEIGHT 1
      k = 0
      DO WHILE(k<edges_p)
            i = INT(r1279()*N) + 1 
            j = INT(r1279()*N) + 1 

            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,1)
                  CALL ADDTOLIST(JJ(j)%v,1)
                  k = k + 1 
            END IF
      END DO
      
C     GENERATE M/2 EDGES OF WEIGHT -1 
      k = 0
      DO WHILE(k<edges_n)
            i = INT(r1279()*N) + 1 
            j = INT(r1279()*N) + 1

            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN    
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,-1)
                  CALL ADDTOLIST(JJ(j)%v,-1)
                  k = k + 1 
            END IF
      END DO

C	INBR GENERATION
      DO i = 1,N
            ALLOCATE(INBR(i)%v(SIZE(NBR(i)%v)))
      END DO
      DO i = 1,N 
            DO j = 1,SIZE(NBR(i)%v)
                  DO k = 1,SIZE(NBR(NBR(i)%v(j))%v)
                        IF ((NBR(NBR(i)%v(j))%v(k).EQ.i)) THEN
                              INBR(i)%v(j) = k
                        END IF
                  END DO
            END DO 
      END DO

      RETURN
      END SUBROUTINE IRS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     RANDOM ERDÖS-RÉNYI GRAPH WITHOUT COUPLINGS GENERATOR
      SUBROUTINE IRG(N,M,NBR,INBR,JJ) ! Initial Random Graph
C     THIS SUBROUTINE GENERATES A RANDOM ERDÖS-RÉNYI GRAPH WITH M EDGES
C     AND SAVES IT IN THE NBR, INBR AND JJ ARRAYS.

      INTEGER N, M
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      
      INTEGER i, j, k
      EXTERNAL r1279

      ALLOCATE(NBR(N))
      ALLOCATE(INBR(N))
      ALLOCATE(JJ(N))

      DO i = 1,N
            ALLOCATE(NBR(i)%v(0))
            ALLOCATE(JJ(i)%v(0))
      ENDDO

C     GENERATE M EDGES OF WEIGHT 0
      k = 0
      DO WHILE(k<M)
            i = INT(r1279()*N) + 1 
            j = INT(r1279()*N) + 1 
            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,0)
                  CALL ADDTOLIST(JJ(j)%v,0)
                  k = k + 1
            END IF
      END DO

C	INBR GENERATION
      DO i = 1,N
            ALLOCATE(INBR(i)%v(SIZE(NBR(i)%v)))
      END DO
      DO i=1,N 
            DO j=1,SIZE(NBR(i)%v)
                  DO k=1,SIZE(NBR(NBR(i)%v(j))%v)
                        IF ((NBR(NBR(i)%v(j))%v(k).EQ.i)) THEN
                              INBR(i)%v(j) = k
                        END IF
                  END DO
            END DO
      END DO

      RETURN
      END SUBROUTINE IRG

C-----------------------------------------------------------------------

C------------------------------------------------------------------
C     RANDOM COUPLING ASSIGMENT
      SUBROUTINE RCA(N,M,p,NBR,INBR,JJ)
C     THIS SUBROUTINE ASSIGNS TO THE GRAPH p*M EDGES WITH
C     A WEIGHT OF 1 AND (1-p)*M EDGES WITH A WEIGHT OF -1
C     AND SAVES IT IN THE JJ ARRAY.

      INTEGER N, M
      REAL*8 p
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      INTEGER i, j, k
      EXTERNAL r1279
      INTEGER edges_p, edges_n

      edges_p = INT(p*M)
      edges_n = INT((1-p)*M)

      IF (edges_p+edges_n.NE.M) THEN
            WRITE(*,*)
            PRINT*, 'GRAPH PROBLEM: p value not possible'
            WRITE(*,*)
            STOP
      END IF

C     GENERATE M*p EDGES OF WEIGHT 1
      k = 0
      DO WHILE(k<edges_p)
1           i = INT(r1279()*N) + 1 
            j = INT(r1279()*SIZE(JJ(i)%v)) + 1 
            IF (SIZE(NBR(i)%v).EQ.0) GO TO 1
            IF (JJ(i)%v(j).EQ.0) THEN
                  JJ(i)%v(j) = 1
                  JJ(NBR(i)%v(j))%v(INBR(i)%v(j)) = 1
                  k = k + 1 
            END IF
      END DO

C     GENERATE M*(1-p) EDGES OF WEIGHT -1
      k = 0
      DO WHILE(k<edges_n)
2           i = INT(r1279()*N) + 1 
            j = INT(r1279()*SIZE(JJ(i)%v)) + 1 
            IF (SIZE(NBR(i)%v).EQ.0) GO TO 2
            IF (JJ(i)%v(j).EQ.0) THEN
                  JJ(i)%v(j) = -1
                  JJ(NBR(i)%v(j))%v(INBR(i)%v(j)) = -1
                  k = k + 1 
            END IF
      END DO

      END SUBROUTINE RCA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION ENERG(N,S,NBR,JJ)
C     THIS FUNCTION CALCULATES THE ENERGY OF THE SYSTEM GIVEN AN
C     S CONFIGURATION

      INTEGER N
      INTEGER S(1:N)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      
      INTEGER i
      REAL*8 ENE

      ENE = 0.0d0
      
      DO i = 1,N
            DO k= 1,SIZE(NBR(i)%v)
                  ENE = ENE - S(i)*S(NBR(i)%v(k))*JJ(i)%v(k)
            END DO
      END DO

      ENERG =  ENE/2
      
      RETURN
      END FUNCTION ENERG
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION MAGNET_Z(N,S)

      INTEGER N
      INTEGER S(1:N)

      INTEGER i
      REAL*8 SUM
      
      SUM = 0.D0
      DO i = 1,N
            SUM = SUM + S(i)
      END DO

      MAGNET_Z = SUM/N

      RETURN
      END FUNCTION MAGNET_Z
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     METROPOLIS ALGORITHM 
      SUBROUTINE METROPOLIS(N,S,valid,TEMP,DE,NBR,JJ)
C     THIS SUBROUTINE PROPOSES A CHANGE OF SPIN IN A RANDOM NODE, CALCULATES
C     THE ENERGY VARIATION OF THE SYSTEM (ΔH_eff) DUE TO IT, IF ΔH_eff < 0
C     THEN THE CHANGE IS ACCPETED, ELSE IF ΔH_eff > 0 THEN THE CHANGE IS 
C     ACCEPTED WITH A PROBABILITY OF EXP(-ΔH_eff/k_BT).

      INTEGER N
      INTEGER S(1:N)
      LOGICAL valid
      REAL*8 TEMP, DE
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      EXTERNAL r1279
      REAL*8 SUM

      valid = .FALSE.

C     RANDOM NODE SELECTION
      i = INT(r1279()*N) + 1 

C     CALCULATION OF ΔH
      SUM = 0
      DO k = 1,SIZE(NBR(i)%v)
            SUM = SUM + JJ(i)%v(k)*S(NBR(i)%v(k))
      END DO
      DE = 2*S(i)*SUM

C     CHECK IF THE CHANGE IS ACCEPTED
      IF (r1279().LT.MIN(1.d0,exp(-DE/TEMP))) THEN
            S(i) = -S(i)
            valid = .TRUE.
      END IF

      RETURN
      END SUBROUTINE METROPOLIS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     ARRAY TO BINARY
      SUBROUTINE ARRAY2BIN(N,binary,array)
C     THIS SUBRUTINE CONVERTS AN ARRAY OF -1 AND 1 TO A BINARY NUMBER

      INTEGER N
      INTEGER, ALLOCATABLE :: binary(:)
      INTEGER array(1:N)

      INTEGER i

      DO i = 1,N
            binary(i) = 0
      END DO
      DO i = 1,N
            IF (array(i) == 1) THEN
                  binary(i) = 1
            END IF
      END DO

      RETURN
      END SUBROUTINE ARRAY2BIN
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     BINARY TO DECIMALS
      SUBROUTINE BIN2DEC(N,zip_size,binary,decimal)
C     THIS SUBRUTINE CONVERTS A BINARY NUMBER TO N/zip_size DECIMAL NUMBERS

      INTEGER N, zip_size
      INTEGER, ALLOCATABLE :: binary(:)
      INTEGER decimal(1:N/zip_size)

      INTEGER i, j, scale

      DO j = 1,N/zip_size
            decimal(j) = 0
            scale = (j-1)*zip_size
            DO i = 1,zip_size
                  IF (binary(scale+i).EQ.1) THEN
                        decimal(j) = decimal(j) + 2**(zip_size-i)
                  END IF
            END DO
      END DO

      RETURN
      END SUBROUTINE BIN2DEC
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     DECIMALS TO BINARY
      SUBROUTINE DEC2BIN(N,zip_size,binary,decimal)
C     THIS SUBRUTINE CONVERTS N/zip_size DECIMAL NUMBERS TO A BINARY NUMBER

      INTEGER N, zip_size
      INTEGER, ALLOCATABLE :: binary(:)
      INTEGER decimal(1:N/zip_size)

      INTEGER i, j, scale
      INTEGER decimal_copy(1:N/zip_size)

      decimal_copy = decimal
      DO j=1,N/zip_size
      scale = (j-1)*zip_size
      DO i = zip_size,1,-1
            IF (MOD(decimal_copy(j),2)==1) THEN
                  binary(scale+i) = 1
            ELSE
                  binary(scale+i) = 0
            END IF
            decimal_copy(j) = decimal_copy(j)/2
      END DO
      END DO

      RETURN
      END SUBROUTINE DEC2BIN
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     BINARY TO ARRAY
      SUBROUTINE BIN2ARRAY(N,binary,array)
C     THIS SUBRUTINE CONVERTS A BINARY NUMBER TO AN ARRAY OF -1 AND 1

      INTEGER N
      INTEGER, ALLOCATABLE :: binary(:)
      INTEGER array(1:N)

      INTEGER i

      DO i = 1,N
            IF (binary(i).EQ.1) THEN
                  array(i) = 1
            ELSE
                  array(i) = -1
            END IF
      END DO

      RETURN
      END SUBROUTINE BIN2ARRAY
C-----------------------------------------------------------------------

C------------------------------------------------------------------
C     ADD ELEMENT TO LIST
      SUBROUTINE ADDTOLIST(list,element)

      INTEGER, DIMENSION(:), ALLOCATABLE:: list
      INTEGER element

      INTEGER i, isize
      INTEGER, DIMENSION(:), ALLOCATABLE :: clist

      IF (ALLOCATED(list)) THEN
            isize = SIZE(list)
            ALLOCATE(clist(isize+1))
            DO i = 1,isize          
                  clist(i) = list(i)
            END DO
            clist(isize+1) = element
            DEALLOCATE(list)
            CALL MOVE_ALLOC(clist, list)
      ELSE
            ALLOCATE(list(1))
            list(1) = element
      END IF

      RETURN
      END SUBROUTINE ADDTOLIST
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     REMOVE INDEX FROM LIST
      SUBROUTINE RMVOFLIST(list,index)

      INTEGER, DIMENSION(:), ALLOCATABLE:: list
      INTEGER index

      INTEGER i, isize
      INTEGER, DIMENSION(:), ALLOCATABLE:: clist

      IF (ALLOCATED(list)) THEN
            isize = SIZE(list)
            ALLOCATE(clist(isize-1))
            DO i = 1,index-1
                  clist(i) = list(i)
            END DO
            DO i = index,isize-1
                  clist(i) = list(i+1)
            END DO
            DEALLOCATE(list)
            CALL MOVE_ALLOC(clist, list)

      END IF

      RETURN
      END SUBROUTINE RMVOFLIST
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     TABULATE LAMBDA FOR A GIVEN D        
      SUBROUTINE CLASS_LAMBDA(N,C,D,NBR,JJ,LAMBDA)
C     THIS SUBROUTINE TABULATES THE VALUES OF THE SUM INSIDE THE
C     TANH OF THE PSEUDOLIKELIHOOD FOR A GIVEN D AND GRAPH
C     DEFINED BY NBR,INBR AND JJ.

      INTEGER N, C
      INTEGER D(1:C,1:N)
      TYPE(multi_array),ALLOCATABLE:: NBR(:)
      TYPE(multi_array),ALLOCATABLE:: JJ(:)
      INTEGER LAMBDA(1:C,1:N)

      INTEGER i, k, m
      INTEGER SUM

      SUM = 0
      DO m = 1,C
            DO i = 1,N
                  DO k = 1, SIZE(NBR(i)%v)
                        SUM = SUM + JJ(i)%v(k)*D(m,NBR(i)%v(k))
                  ENDDO
                  LAMBDA(m,i) = SUM
                  SUM = 0
            ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CLASS_LAMBDA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     TABULATE FUNCT FOR A GIVEN TEMP           
      SUBROUTINE CLASS_FUNCTION(zmax,TEMP,funct)
C     THIS SUBROUTINE TABULATES THE FUNCTION IN THE SUMATORIES OF
C     THE PSEUDOLIKELIHOOD.

      INTEGER zmax
      REAL*8 TEMP
      REAL*8 funct(-1:1,-zmax:zmax)

      INTEGER i, j

      DO i = -1,1,2
            DO j = -zmax,zmax
                  funct(i,j) = LOG(0.5d0*(1 + i*TANH(j/TEMP)))
            ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CLASS_FUNCTION
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION PSEUDO(N,C,D,TEMP,NBR,JJ)
C     THIS FUNCTION CALCULATES THE PSEUDOLIKELIHOOD FOR A GIVEN TEMP,
C     D AND GRAPH DEFINED BY NBR AND JJ

      INTEGER N, C
      INTEGER D(1:C,1:N)
      REAL*8 TEMP
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      INTEGER i, k, m
      REAL*8 L
      REAL*8 SUM

      PSEUDO = 0.0d0
      L = 0.0d0
      SUM = 0.0d0

      DO i = 1,N
      DO m = 1,C
            DO k = 1, SIZE(JJ(i)%v)
                  SUM = SUM + JJ(i)%v(k)*D(m,NBR(i)%v(k))
            ENDDO
            
            L = L + LOG(0.5d0*(1 + D(m,i)*TANH(SUM/TEMP)))
            SUM = 0.0d0
      ENDDO
      PSEUDO = PSEUDO + L/C
      L = 0.0d0
      ENDDO

      RETURN
      END FUNCTION PSEUDO
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     FIND LOCATION OF A VALUE IN AN ARRAY
      INTEGER FUNCTION FLOC(array,value)

      INTEGER,ALLOCATABLE:: array(:)
      INTEGER value
      
      INTEGER r
      LOGICAL out

      FLOC = 0
      DO r = 1, SIZE(array)
            IF (array(r).EQ.value) THEN
                  FLOC = r
                  out = .TRUE.
            END IF
      END DO

      RETURN
      END FUNCTION FLOC
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     PSEUDOLIKELIHOOD ALGORITHM 
      SUBROUTINE PSEUDOLIKELIHOOD(N,C,D,
     . valid,TEMP_F,DPL,NBR,INBR,JJ,zmax,funct,LAMBDA)
C     THIS SUBROUTINE PROPOSES A CHANGE OF PAIRWISE COUPLING,
C     CALCULATES THE PSEUDOLIKELIHOOD VARIATION OF THE SYSTEM (ΔPL) DUE TO IT,
C     IF ΔPL > 0 THEN THE CHANGE IS ACCPETED, ELSE IF ΔPL < 0 THEN THE
C     CHANGE IS ACCEPTED WITH A PROBABILITY OF EXP(-ΔPL/k_BT'), WHERE
C     T' IS THE FICTICIOUS TEMPERATURE

      INTEGER N, C
      INTEGER D(1:C,1:N)
      LOGICAL valid
      REAL*8 TEMP_F
      REAL*8  DPL
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      INTEGER zmax
      REAL*8 funct(-1:1,-zmax:zmax)
      INTEGER LAMBDA(1:C,1:N)

      LOGICAL change
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newNBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newINBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
      EXTERNAL r1279
      REAL*8 SUM
      INTEGER m,i,i1,i2,i3,i4
      INTEGER ip,ip2
      INTEGER newLAMBDA(1:C,1:N)
      
      valid = .FALSE.
      change = .FALSE. 

      ALLOCATE(newNBR(N))
      ALLOCATE(newINBR(N))
      ALLOCATE(newJJ(N))

      DO i = 1,N
            ALLOCATE(newNBR(i)%v(SIZE(NBR(i)%v)))
            ALLOCATE(newINBR(i)%v(SIZE(NBR(i)%v)))
            ALLOCATE(newJJ(i)%v(SIZE(NBR(i)%v)))
      END DO

      newNBR = NBR
      newINBR = INBR
      newJJ = JJ 

      newLAMBDA = LAMBDA

C     RANDOM PAIRWISE COUPLING CHANGE
      DO WHILE (change.EQV..FALSE.)
      CALL JJ_CHANGE(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .              i1,i2,i3,i4)
!       CALL JJ_CHANGE_2(N,NBR,INBR,JJ,newJJ,change,
!      .              i1,i2,i3,i4)
      END DO

C     CALCULATE newLAMBDA
      DO m = 1,C

      ip = FLOC(NBR(i1)%v,i3)
      ip2 = FLOC(newNBR(i1)%v,i3)
      IF (ip.EQ.0) THEN
      newLAMBDA(m,i1) = newLAMBDA(m,i1)+newJJ(i1)%v(ip2)*D(m,i3)
      ELSE IF (ip.NE.0) THEN
      newLAMBDA(m,i1) = newLAMBDA(m,i1)+2*newJJ(i1)%v(ip2)*D(m,i3)
      END IF

      ip = FLOC(NBR(i3)%v,i1)
      ip2 = FLOC(newNBR(i3)%v,i1)
      IF (ip.EQ.0) THEN
      newLAMBDA(m,i3) = newLAMBDA(m,i3)+newJJ(i3)%v(ip2)*D(m,i1)
      ELSE IF (ip.NE.0) THEN
      newLAMBDA(m,i3) = newLAMBDA(m,i3)+2*newJJ(i3)%v(ip2)*D(m,i1)
      END IF

      ip = FLOC(NBR(i2)%v,i4)
      ip2 = FLOC(newNBR(i2)%v,i4)
      IF (ip2.EQ.0) THEN
      newLAMBDA(m,i2) = newLAMBDA(m,i2)-JJ(i2)%v(ip)*D(m,i4)
      ELSE IF (ip2.NE.0) THEN
      newLAMBDA(m,i2) = newLAMBDA(m,i2)+2*newJJ(i2)%v(ip2)*D(m,i4)
      END IF 

      ip = FLOC(NBR(i4)%v,i2)
      ip2 = FLOC(newNBR(i4)%v,i2)
      IF (ip2.EQ.0) THEN
      newLAMBDA(m,i4) = newLAMBDA(m,i4)-JJ(i4)%v(ip)*D(m,i2)
      ELSE IF (ip2.NE.0) THEN
      newLAMBDA(m,i4) = newLAMBDA(m,i4)+2*newJJ(i4)%v(ip2)*D(m,i2)
      END IF 
      
      END DO      

C     CALCULATE DPL
      SUM = 0

      IF ((i1.NE.i4).AND.(i2.NE.i3).AND.(i3.NE.i4)) THEN
      DO m = 1,C
            SUM = SUM + 
     .      funct(D(m,i1),newLAMBDA(m,i1))
     .    - funct(D(m,i1),LAMBDA(m,i1)) +
     .      funct(D(m,i2),newLAMBDA(m,i2))
     .     -funct(D(m,i2),LAMBDA(m,i2)) + 
     .      funct(D(m,i3),newLAMBDA(m,i3)) 
     .     -funct(D(m,i3),LAMBDA(m,i3)) +
     .      funct(D(m,i4),newLAMBDA(m,i4))
     .     -funct(D(m,i4),LAMBDA(m,i4))
      END DO
      DPL = SUM/C
      END IF 

      IF (i1.EQ.i4) THEN
      do m = 1,C
            SUM = SUM + 
     .      funct(D(m,i1),newLAMBDA(m,i1))
     .    - funct(D(m,i1),LAMBDA(m,i1)) +
     .      funct(D(m,i2),newLAMBDA(m,i2))
     .     -funct(D(m,i2),LAMBDA(m,i2)) + 
     .      funct(D(m,i3),newLAMBDA(m,i3)) 
     .     -funct(D(m,i3),LAMBDA(m,i3))
      END DO
      DPL = SUM/C
      END IF

      IF ((i2.EQ.i3).OR.(i4.EQ.i3)) THEN
      DO m = 1,C
            SUM = SUM + 
     .      funct(D(m,i1),newLAMBDA(m,i1))
     .    - funct(D(m,i1),LAMBDA(m,i1)) +
     .      funct(D(m,i2),newLAMBDA(m,i2))
     .     -funct(D(m,i2),LAMBDA(m,i2)) + 
     .      funct(D(m,i4),newLAMBDA(m,i4))
     .     -funct(D(m,i4),LAMBDA(m,i4))
      END DO
      DPL = SUM/C
      END IF

C     CHECK IF THE CHANGE IS ACCEPTED
      IF (r1279().LT.MIN(1.d0,EXP(DPL/TEMP_F))) THEN
            valid = .TRUE.
            NBR = newNBR
            INBR = newINBR
            JJ = newJJ
            
            LAMBDA = newLAMBDA
      END IF

      RETURN
      END SUBROUTINE PSEUDOLIKELIHOOD
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     MAKE A JJ EXCHANGE OF TYPE 1 OR 2
      SUBROUTINE JJ_CHANGE(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .                     i1,i2,i3,i4)
C     THIS SUBROUTINE CHANGES PAIRWISE COUPLING OF i1-i3 FOR THE
C     PAIRWISE COUPLING OF i2-i4, AND GENERATES A THE NEW NBR, INBR
C     AND JJ ARRAYS.

      INTEGER N
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newNBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newINBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
      LOGICAL change
      INTEGER i1, i2, i3, i4

      INTEGER k1, k2, k3, k4
      EXTERNAL r1279
      INTEGER ip
      INTEGER i, j, k
      LOGICAL valid1, valid2

      change = .FALSE.
      valid1 = .FALSE.
      valid2 = .FALSE.

C     SELECTING A RANDOM NODE
      i1 = INT(r1279()*N) + 1 

C     SELECTING A RANDOM EDGE OF THAT NODE, THE EXTRA +1
C     IS BECAUSE WE CAN SELECT A NODE WITH VALUE 0
      k1 = int(r1279()*SIZE(NBR(i1)%v)) + 1 + 1 

C     IN CASE SIZE(NBR(i,k))=0 k1 SHOULD BE EQUAL TO 1
      IF (SIZE(NBR(i1)%v) == 0) THEN
            k1 = 1
      ENDIF

C     IF THE MAX EDGES PER NODE IS EXCEDED WE DENY THE CHANGE
      IF (k1 > N-1) THEN
            CHANGE = .FALSE.
            RETURN
      ENDIF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = 0 WE WANT A 0 <----> 1,-1 EXCHANGE
C-------------------------------------------------------------
      IF (k1 > SIZE(NBR(i1)%v)) THEN

C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid1.EQV..FALSE.)
      i2 = INT(r1279()*N) + 1
      IF ((i1.NE.i2).AND.(SIZE(NBR(i2)%v).GT.1)) THEN 
            valid1 = .TRUE.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     TO AVOID PROBLEMS i2 NEEDS TO HAVE MORE THAN 1 NEIGHBOR
C     OTHERWISE A NODE COULD END UP WITH 0 EDGES

C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      k2 = INT(r1279()*SIZE(NBR(i2)%v)) + 1

C     SELECTING A NODE i3 THAT IS NOT A NEIGBOR  OF i1 ALREADY
      DO WHILE (valid2.EQV..FALSE.)
      i3 = INT(r1279()*N) + 1
C     THIS ALLOW TO SEE IF i3 IS ALREADY A NEIGBOUR OF i1
      ip = FLOC(NBR(i1)%v,i3)
C     ip IS 0 IF i3 IS NOT IN NBR(i1)%v
      IF ((ip.EQ.0).AND.(i1.NE.i3).AND.(SIZE(NBR(i3)%v).LT.N-1)) THEN
            valid2 = .TRUE.
      END IF
      END DO 

      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

C     ADDING THE NEW EDGE
      CALL ADDTOLIST(newNBR(i1)%v,i3)
      CALL ADDTOLIST(newNBR(i3)%v,i1)
      CALL ADDTOLIST(newJJ(i1)%v,JJ(i2)%v(k2))
      CALL ADDTOLIST(newJJ(i3)%v,JJ(i4)%v(k4))

C     REMOVING THE OLD EDGE
      CALL RMVOFLIST(newNBR(i2)%v,k2)
      CALL RMVOFLIST(newNBR(i4)%v,k4)
      CALL RMVOFLIST(newJJ(i2)%v,k2)
      CALL RMVOFLIST(newJJ(i4)%v,k4)

C     UPDATING THE INBR     
      newINBR = newNBR
      DO i = 1,N 
            DO j = 1,SIZE(newNBR(i)%v)
                  DO k = 1,SIZE(newNBR(newNBR(i)%v(j))%v)
                        IF ((newNBR(newNBR(i)%v(j))%v(k).EQ.i)) THEN
                              newINBR(i)%v(j) = k
                        END IF
                  END DO
            END DO
      END DO

      change = .TRUE.
      RETURN
      END IF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = +1,-1 WE WANT a +-1 <----> -+1 EXCHANGE
C-------------------------------------------------------------
      IF (k1 <= SIZE(NBR(i1)%v)) THEN
      
C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid1.EQV..false.)
      i2 = INT(r1279()*N) + 1

      IF ((i1.NE.i2).AND.(SIZE(NBR(i2)%v).gt.1)) THEN
            valid1 = .TRUE.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      DO WHILE (valid2.EQV..false.)
      k2 = INT(r1279()*SIZE(NBR(i2)%v)) + 1
      IF (NBR(i2)%v(k2).NE.i1) THEN
            valid2 = .TRUE.
      END IF
      END DO
C     IF BOTH EDGES HAVE THE SAME VALUE THERE IS NO CHANGE
      IF (JJ(i1)%v(k1).EQ.JJ(i2)%v(k2)) THEN
            change = .FALSE.
            RETURN
      END IF
      
C     CHANGING THE EDGE VALUES
      i3 = NBR(i1)%v(k1)
      k3 = INBR(i1)%v(k1)
      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

      newJJ(i1)%v(k1) = JJ(i2)%v(k2)
      newJJ(i3)%v(k3) = JJ(i2)%v(k2)
      newJJ(i2)%v(k2) = JJ(i1)%v(k1)
      newJJ(i4)%v(k4) = JJ(i1)%v(k1)

      change = .TRUE.
      RETURN
      END IF 

      END SUBROUTINE JJ_CHANGE
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     MAKE A JJ EXCHANGE OF TYPE 2 (THIS IMPLIES THAT THE INITIAL GRAPH IS THE SAME AS THE ORIGINAL ONE BUT WITH A DIFFERENT JJ ARRAY)
      SUBROUTINE JJ_CHANGE_2(N,NBR,INBR,JJ,newJJ,change,
     .                     i1,i2,i3,i4)
C     THIS SUBROUTINE CHANGES PAIRWISE COUPLING OF i1-i3 FOR THE
C     PAIRWISE COUPLING OF i2-i4, AND GENERATES A THE NEW NBR, INBR
C     AND JJ ARRAYS.

      INTEGER N
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
      LOGICAL change
      INTEGER i1, i2, i3, i4

      INTEGER k1, k2, k3, k4
      EXTERNAL r1279
      LOGICAL valid1, valid2

      change = .FALSE.
      valid1 = .FALSE.
      valid2 = .FALSE.

C     SELECTING A RANDOM NODE
      i1 = INT(r1279()*N) + 1 

C     SELECTING A RANDOM EDGE OF THAT NODE, THE EXTRA +1
C     IS BECAUSE WE CAN SELECT A NODE WITH VALUE 0
      k1 = int(r1279()*SIZE(NBR(i1)%v)) + 1 + 1 

C     IN CASE SIZE(NBR(i,k))=0 k1 SHOULD BE EQUAL TO 1
      IF (SIZE(NBR(i1)%v) == 0) THEN
            k1 = 1
      ENDIF

C     IF THE MAX EDGES PER NODE IS EXCEDED WE DENY THE CHANGE
      IF (k1 > N-1) THEN
            CHANGE = .FALSE.
            RETURN
      ENDIF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = +1,-1 WE WANT a +-1 <----> -+1 EXCHANGE
C-------------------------------------------------------------
      IF (k1 <= SIZE(NBR(i1)%v)) THEN
      
C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid1.EQV..false.)
      i2 = INT(r1279()*N) + 1

      IF ((i1.NE.i2).AND.(SIZE(NBR(i2)%v).gt.1)) THEN
            valid1 = .TRUE.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      DO WHILE (valid2.EQV..false.)
      k2 = INT(r1279()*SIZE(NBR(i2)%v)) + 1
      IF (NBR(i2)%v(k2).NE.i1) THEN
            valid2 = .TRUE.
      END IF
      END DO
C     IF BOTH EDGES HAVE THE SAME VALUE THERE IS NO CHANGE
      IF (JJ(i1)%v(k1).EQ.JJ(i2)%v(k2)) THEN
            change = .FALSE.
            RETURN
      END IF
      
C     CHANGING THE EDGE VALUES
      i3 = NBR(i1)%v(k1)
      k3 = INBR(i1)%v(k1)
      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

      newJJ(i1)%v(k1) = JJ(i2)%v(k2)
      newJJ(i3)%v(k3) = JJ(i2)%v(k2)
      newJJ(i2)%v(k2) = JJ(i1)%v(k1)
      newJJ(i4)%v(k4) = JJ(i1)%v(k1)

      change = .TRUE.
      RETURN
      END IF 

      END SUBROUTINE JJ_CHANGE_2
C-----------------------------------------------------------------------

C------------------------------------------------------------------
      REAL*8 FUNCTION GAMMAA(N,M,NBR,JJ,NBR_0,JJ_0)
C     THIS FUNCTION CALCULATES THE GAMMA BETWEEN JJ AND JJ_0

      INTEGER N, M
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ_0(:)

      REAL*8 SUM
      INTEGER ip

      SUM = 0.d0
      DO i = 1,N
            DO k = 1,SIZE(NBR_0(i)%v)
                  IF (i.LT.NBR_0(i)%v(k)) THEN
                  ip = FLOC(NBR(i)%v,NBR_0(i)%v(k))
                  IF (ip.EQ.0) THEN
                        SUM = SUM + (JJ_0(i)%v(k))**2
                  END IF
                  IF (ip.ne.0) THEN
                        SUM = SUM + (JJ_0(i)%v(k)-JJ(i)%v(ip))**2
                  END IF
                  END IF
            END DO
      END DO

      GAMMAA = SUM/M
      RETURN
      END FUNCTION GAMMAA
C-----------------------------------------------------------------------

      END MODULE MODEL