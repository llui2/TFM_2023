C     SPIN MODEL MODULE
C     0!
C     Lluís Torres 
C     TFG
C     FORTRAN 95

      MODULE MODEL 
C     VIANA-BRAY SPIN-GLASS MODEL WITH DISCRETE DISTRIBUTION OF THE COUPLING IN A TRANSVERSE FIELD

C     MULTI ARRAY TYPE
      TYPE :: MULTI_ARRAY
      INTEGER,ALLOCATABLE :: v(:)
      END TYPE MULTI_ARRAY

      CONTAINS

C-----------------------------------------------------------------------
C     METROPOLIS.F
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     READ INPUT FILE
      SUBROUTINE READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,NSEEDS,SC,zip_size,TAU)

      INTEGER N,z
      INTEGER R
      INTEGER TEMP_SIZE
      REAL*8,ALLOCATABLE:: TEMP_LIST(:)
      INTEGER H_SIZE
      REAL*8,ALLOCATABLE:: H_LIST(:)
      INTEGER p_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:)
      INTEGER C
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
      SUBROUTINE IRS(N,M,p,NBR,INBR,JJ)
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
            i = MOD(INT(r1279()*N),N) + 1
            j = MOD(INT(r1279()*N),N) + 1

            IF ((i.gt.N).or.(j.gt.N)) THEN
                  PRINT*, 'ERROR: i OR j OUT OF RANGE'
                  STOP
            END IF

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
            i = MOD(INT(r1279()*N),N) + 1 
            j = MOD(INT(r1279()*N),N) + 1 

            IF ((i.gt.N).or.(j.gt.N)) THEN
                  PRINT*, 'ERROR: i OR j OUT OF RANGE'
                  STOP
            END IF

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
      REAL*8 FUNCTION ENERG(N,R,S,TEMP,H,NBR,JJ)
C     THIS FUNCTION CALCULATES THE ENERGY OF THE SYSTEM GIVEN A CONFIGURATION

      INTEGER N, R
      INTEGER S(1:R,1:N)
      REAL*8 TEMP, H
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      INTEGER i, j, k
      REAL*8 K2
      REAL*8 HD, V
      INTEGER ABOVE(1:R) !SPIN ABOVE i

      DO i=1,R-1
            ABOVE(i) = i+1
      END DO
      ABOVE(R) = 1

      K2 = -(TEMP/2)*LOG(TANH(H/(TEMP*R)))

      HD = 0.0d0 !DIAGONAL TERM
      V = 0.0d0 !TRANSVERSE FIELD TERM

      DO i = 1,R
            DO j = 1,N
                  DO k = 1,SIZE(NBR(j)%v)
                        HD = HD + JJ(j)%v(k)*S(i,j)*S(i,NBR(j)%v(k))
                  END DO
                  V = V + S(i,j)*S(ABOVE(i),j)
            END DO
      END DO

      ENERG =  - HD/(2*R) - K2*V
      
      RETURN
      END FUNCTION ENERG
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION MAGNET_Z(N,R,S)

      INTEGER N, R
      INTEGER S(1:R,1:N)

      INTEGER i, j
      REAL*8 SUM
      
      SUM = 0.D0
      DO i = 1,R
            DO j = 1,N
                  SUM = SUM + S(i,j)
            END DO
      END DO

      MAGNET_Z = SUM/(N*R)

      RETURN
      END FUNCTION MAGNET_Z
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION MAGNET_X(N,R,S,TEMP,H)

      INTEGER N, R
      INTEGER S(1:R,1:N)
      REAL*8 TEMP, H

      INTEGER i, j
      INTEGER ABOVE(1:R) !SPIN ABOVE i
      REAL*8 K, MAG

      DO i=1,R-1
            ABOVE(i) = i+1
      END DO
      ABOVE(R) = 1

      K = TAN(H/(TEMP*R))

      MAG = 0.D0
      DO i = 1,R
            DO j = 1,N
                  MAG = MAG + K**(S(i,j)*S(ABOVE(i),j))
            END DO
      END DO
      MAGNET_X = MAG/(N*R)

      RETURN
      END FUNCTION MAGNET_X
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     METROPOLIS ALGORITHM 
      SUBROUTINE METROPOLIS(N,R,S,valid,TEMP,H,DE,NBR,JJ)
C     THIS SUBROUTINE PROPOSES A CHANGE OF SPIN IN A RANDOM NODE, CALCULATES
C     THE ENERGY VARIATION OF THE SYSTEM (ΔH_eff) DUE TO IT, IF ΔH_eff < 0
C     THEN THE CHANGE IS ACCPETED, ELSE IF ΔH_eff > 0 THEN THE CHANGE IS 
C     ACCEPTED WITH A PROBABILITY OF EXP(-ΔH_eff/k_BT).

      INTEGER N, R
      INTEGER S(1:R,1:N)
      LOGICAL valid
      REAL*8 TEMP, H, DE
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      EXTERNAL r1279
      REAL*8 DHD, DV
      REAL*8 K2

      INTEGER PBC(0:R+1) !PBC IN THE TROTTER DIRECTION

      PBC(0)=R
      DO i=1,R
            PBC(i) = i
      END DO
      PBC(R+1) = 1

      valid = .FALSE.

C     RANDOM NODE SELECTION
      i = MOD(INT(r1279()*R),R) + 1
      j = MOD(INT(r1279()*N),N) + 1
      
C     CALCULATION OF ΔHD
      DHD = 0
      DO k=1,SIZE(NBR(j)%v)
            DHD = DHD + JJ(j)%v(k)*S(i,NBR(j)%v(k))
      END DO
      DHD = 2*DHD*S(i,j)/R

C     CALCULATION OF ΔV
      K2 = -(TEMP/2.)*LOG(TANH(H/(TEMP*R)))
      DV = K2*S(i,j)*(S(PBC(i+1),j)+S(PBC(i-1),j))

C     CALCULATION OF ΔH
      DE = DHD + DV

C     CHECK IF THE CHANGE IS ACCEPTED
      IF (r1279().LT.min(1.d0,exp(-DE/TEMP))) THEN
            S(i,j) = -S(i,j)
            valid = .true.
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
            isize = size(list)
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
      REAL*8 FUNCTION PSEUDO(N,R,C,D,TEMP,H,NBR,JJ)
C     THIS FUNCTION CALCULATES THE PSEUDOLIKELIHOOD FOR A GIVEN TEMP, H,
C     D AND GRAPH DEFINED BY NBR AND JJ

      INTEGER N, R, C
      INTEGER D(1:C,1:R,1:N)
      REAL*8 TEMP, H
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      INTEGER a, i, m, k
      REAL*8 L
      REAL*8 SUM
      REAL*8 K2
      INTEGER PBC(0:R+1) !PBC IN THE TROTTER DIRECTION

      PBC(0)=R
      DO i=1,R
            PBC(i) = i
      END DO
      PBC(R+1) = 1

      PSEUDO = 0.0d0
      L = 0.0d0
      SUM = 0.0d0

      K2 = -(TEMP/2.)*LOG(TANH(H/(TEMP*R)))

      DO a = 1,R
      DO i = 1,N
      DO m = 1,C
      DO k = 1,SIZE(JJ(i)%v)
            SUM = SUM + JJ(i)%V(k)*D(a,m,NBR(i)%V(k)) !COUPLING CONTRIBUTION
      END DO !k neighbours
      SUM = SUM + K2*D(a,m,i)*(D(PBC(a-1),m,i) + D(PBC(a+1),m,i)) !TRANSVERSE CONTRIBUTION
      L = L + LOG(0.5d0*(1 + D(a,m,i)*TANH(SUM/TEMP)))
      SUM = 0.0d0
      END DO !m configurations
      PSEUDO = PSEUDO + L/C
      L = 0.0d0
      END DO !i nodes
      END DO !a Trotter slices

      RETURN
      END FUNCTION PSEUDO
C-----------------------------------------------------------------------

      END MODULE MODEL