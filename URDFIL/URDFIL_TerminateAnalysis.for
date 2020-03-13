      SUBROUTINE URDFIL(LSTOP, LOVRWRT, KSTEP, KINC,
     1 DTIME,TIME)

      INCLUDE 'ABA_PARAM.INC'
      DIMENSION ARRAY(513), JRRAY(NPRECD, 513), TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
      PARAMETER(TOL=5.0D8)
      PARAMETER(DEFL=1.5D-1)
      LMISES=0
      LDEFL=0
      PRINT *,KINC,KSTEP,TIME,DTIME
      ! Assume that we do not mind if .fil file results are
      ! overwritten.

      LOVRWRT=1
      
      ! Find current increment

      CALL POSFIL(KSTEP, KINC, ARRAY, JRCD)

      ! Loop over all of the records

      DO K1=1,999999
        CALL DBFILE(0,ARRAY,JRCD)
        IF (JRCD .NE. 0) GO TO 110
          KEY=JRRAY(1,2)
          ! Record 1 contains element information for subsequent
          ! records
          IF (KEY .EQ. 1) THEN
            IELM = JRRAY(1, 3)
            IMATPT = JRRAY(1, 4)
          END IF
	      
          ! Record 12 contains values for SINV
	      
          IF (KEY .EQ. 12) THEN
            IF (ARRAY(3) .GT. TOL) THEN
              LMISES=1
              GOTO 210
            END IF
          END IF
	      
          ! Record 101 contains U
	      
          IF (KEY .EQ. 101) THEN
            IF (JRRAY(1, 3) .EQ. 63) THEN
              IF (ABS(ARRAY(5)) .GT. DEFL) THEN
                LDEFL=1
              GOTO 210
            END IF
          END IF
        END IF
      END DO
  110 CONTINUE

  210 IF (LMISES .EQ. 1) THEN
        WRITE(7, *)
        WRITE(7, 1023) TOL, IELEM, IMATPT
        WRITE(7, *)
        LSTOP=1
      END IF
      IF (LDEFL .EQ. 1) THEN
        WRITE(7, *)
        WRITE(7, 1024) DEFL
        WRITE(7, *)
        LSTOP=1
      END IF
      RETURN


 1023 FORMAT ('***NOTE: ANALYSIS TERMINATES MISES 
     &STRESS EXCEEDS', 2X, E9.3, 1X, 'IN ELEMENT', 1X, I6,
     &  1X, 'AT INT. PT.', 1X, I6)
 1024 FORMAT ('***NOTE: ANALYSIS TERMINATES AS DEFLECTION
     & OF NODE 63 EXCEEDS', 2X, E9.3)
      END