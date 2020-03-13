      MODULE XFEM_VARS
        IMPLICIT NONE
        PARAMETER(N=1000000)
        REAL*8 COORDS(N,3)
        REAL*8 XPHI(N),XPSI(N)
      END MODULE

C
C***********************************************************************
C
C  *User element, nodes=4, type=U12, properties=2, iproperties=5, coordinates=3, variables=100
C  1,2,3,4,5,6,7,11,12,13,14,15
C  
C  *Element, type=U12,
C
C  *Uel property, elset=ELEMTOPXU12
C  1.0e7, 0.333,1,5,7,5,2
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      !user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      
      DIMENSION F(NDOFEL)
      Parameter(Nint=1)

      AMATRX=0.0D0
      RHS(1:NDOFEL,1)=0.0D0
      SVARS(1:NSVARS)=0.0D0
      ENERGY=0.0D0

      ! For Element type U12
      IF(JTYPE.EQ.12) THEN
      DO I=1,Nint
        DO J=1,NDOFEL
          AMATRX(J,J)=1.0E6
        END DO
        F=MATMUL(AMATRX,U)
        RHS(1:NDOFEL,1)=RHS(1:NDOFEL,1)-F(1:NDOFEL)
      END DO
      PRINT *,KSTEP,KINC,TIME,DTIME
      PRINT *,JELEM,MCRD,NNODE,COORDS
      PRINT *,NPROPS,PROPS(1:NPROPS),NJPROP,JPROPS(1:NJPROP)
      PRINT *,NDOFEL,U(1:NDOFEL),V(1:NDOFEL),MLVARX,DU(1:MLVARX,1)
      PRINT *,NDOFEL,RHS(1:NDOFEL,1)!,NDOFEL,NDOFEL,AMATRX
      PRINT *,NSVARS,SVARS(1:NSVARS)
      END IF

      RETURN
      END
C
C***********************************************************************
C
C   用于计算裂纹形态和裂尖坐标
C
        SUBROUTINE KLEVELSET()

        RETURN
        END SUBROUTINE

C
C             v 1.0.0
C   Caculate Gauss Integrate Point
C   
C   the order of integration points is same with abaqus
C   search Node ordering and Numbering of integration for detail
C   
C   input
C       NNODE : number of nodes
C       NDIM: number of coordinate
C       nINTP: number of gauss integration point
C   return  :
C       GP  : coordinate of Guass Integration Point
C       GW  :   Weight of corespond point
      SUBROUTINE KGUASS(NNODE,NDIM,nINTP,GP,GW)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION GP(NDIM,nINTP),GW(nINTP)

        IF (NNODE.EQ.4 .AND. NDIM.EQ.2) THEN
          ! 
          IF (NINTP.EQ.2) THEN
            CN=0.5773502691896260D0
            GP(1,:)=(/-CN,CN/)
            GP(2,:)=0.D0
            GW=(/2.D0,2.D0/)
          ELSEIF (NINTP.EQ.1) THEN
            GP=0.D0
            GW(1)=4.D0
          ELSEIF (NINTP.EQ.4) THEN
            CN=0.5773502691896260D0
            GP(1,:)=(/-CN,CN,-CN,CN/)
            GP(2,:)=(/-CN,-CN,CN,CN/)
            GW=(/1.D0,1.D0,1.D0,1.D0/)
          ELSEIF (NINTP.EQ.3) THEN
            CN=0.774596669241483D0
            GP(1,:)=(/-CN,0.D0,CN/)
            GP(2,:)=0.D0
            GW=2D0*(/0.55555555556D0,0.88888888889D0,0.55555555556D0/)
          ENDIF
        ELSEIF(NNODE.EQ.8 .AND. NDIM.EQ.3) THEN
          IF (NINTP.EQ.2) THEN
            CN=0.5773502691896260D0
            GP(1,:)=(/-CN,CN/)
            GP(2:3,:)=0.D0
            GW=(/4.D0,4.D0/)
          ELSEIF (NINTP.EQ.8) THEN
            CN=0.5773502691896260D0
            GP(1,:)=(/-CN,CN,-CN,CN,-CN,CN,-CN,CN/)
            GP(2,:)=(/-CN,-CN,CN,CN,-CN,-CN,CN,CN/)
            GP(3,:)=(/-CN,-CN,-CN,-CN,CN,CN,CN,CN/)
            GW=(/1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0/)
          ELSEIF (NINTP.EQ.1) THEN
            GP=0.D0
            GW(1)=8.D0
          ELSEIF (NINTP.EQ.4) THEN
            CN=0.5773502691896260D0
            GP(1,:)=(/-CN,CN,-CN,CN/)
            GP(2,:)=(/-CN,-CN,CN,CN/)
            GP(3,:)=0.D0
            GW=(/2.D0,2.D0,2.D0,2.D0/)
          ELSEIF (NINTP.EQ.3) THEN
            CN=0.774596669241483D0
            GP(1,:)=(/-CN,0.D0,CN/)
            GP(2:3,:)=0.D0
            GW=4D0*(/0.55555555556D0,0.88888888889D0,0.55555555556D0/)
          ENDIF
        END IF
        RETURN
      END SUBROUTINE

C
C   Caculate Shpae Function value and it's derivation at point Coords
C       
C       dNdr[i,j]= \frac{d N_j}{d r_i}
C
C   input:
C       nNode  : number of nodes
C       nCoord : number of coordinate
C       rCoords : local coordinate of input point
C   return:
C       shapeN : shape Function Value 
C       dNdr   :  derivation of shape function at the local point
      SUBROUTINE KSHAPE(nNode,nCoord,rCoords,shapeN,dNdr)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION rCoords(nCoord),shapeN(nNode),dNdr(nCoord,nNode)

        DIMENSION rs(4),ss(4),ts(2)

        rs=(/-1.D0,1.D0,1.D0,-1.D0/)
        ss=(/-1.D0,-1.D0,1.D0,1.D0/)
        ts=(/-1.D0,1.D0/)

        IF (nNode.EQ.4 .AND. nCoord.EQ.2) THEN
          r=rCoords(1)
          s=rCoords(2)

          DO I=1,4
            ShapeN(I)=(1.D0+rs(I)*r)*(1.D0+ss(I)*s)/4.D0
            dNdr(1,I)= rs(I)*(1.D0+ss(I)*s)/4.D0
            dNdr(2,I)= ss(I)*(1.D0+rs(I)*r)/4.D0
          END DO

        ELSEIF (nNode.EQ.8 .AND. nCoord.EQ.3) THEN 
          r=rCoords(1)
          s=rCoords(2)
          t=rCoords(3)

          DO J=1,2
          DO I=1,4
            K=I+4*(J-1)
            ShapeN(K)=(1.D0+rs(I)*r)*(1.D0+ss(I)*s)*(1.D0+ts(J)*t)/8.D0
            dNdr(1,K)= rs(I)*(1.D0+ss(I)*s)*(1.D0+ts(j)*t)/8.D0
            dNdr(2,K)= ss(I)*(1.D0+rs(I)*r)*(1.D0+ts(j)*t)/8.D0
            dNdr(3,K)= ts(J)*(1.D0+rs(I)*r)*(1.D0+ss(I)*s)/8.D0
          END DO
          END DO
        END IF
        RETURN
      END SUBROUTINE

C
C   Caculate Jaccobi Matrix 
C       DJ[i,j]= \frac{d x_j}{d r_i}  i=1,nNode j=1,MCRD
C   input:
C       nNode : number of nodes
C       nCoord: number of coordinate
C       COORDS: COORDS[i,j] is the i-th coordinate of j-th node
C       rCoords: local coordinate of input point
C   return:
C       gCoords: global Coords of input point
C       DJ     :  Jacobi Matrix at input point
C       DET     : determinant of Jacobi Matrix
      SUBROUTINE KJACOBI(nNode,nCoord,COORDS,rCoords,DJ,DET)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION COORDS(nCoord,nNode),rCoords(nCoord)
        DIMENSION gCoords(nCoord),DJ(nCoord,nCoord)

        DIMENSION shapeN(1,nNode),dNdr(nCoord,nNode),xyz(nNode,nCoord)

        CALL KSHAPE(nNode,nCoord,rCoords,shapeN,dNdr)
        
        DJ=MATMUL(dNdr,TRANSPOSE(COORDS))

        CALL KDET(nCoord,DJ,DET)

        RETURN
      END SUBROUTINE

C   
C   CALCULATE DETERMINANT FOR SMALL MATRIX DJ (N<=3)
c
      SUBROUTINE KDET(N,DJ,DET)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION DJ(N,N)

        DET=0.D0
        IF(N.EQ.2) THEN
            DET=DJ(1,1)*DJ(2,2)-DJ(1,2)*DJ(2,1)
        ELSEIF(N.EQ.3) THEN
            DET=(DJ(1,1)*DJ(2,2)-DJ(1,2)*DJ(2,1))*DJ(3,3)
            DET=DET-(DJ(1,1)*DJ(2,3)-DJ(1,3)*DJ(2,1))*DJ(3,2)
            DET=DET+(DJ(1,2)*DJ(2,3)-DJ(1,3)*DJ(2,2))*DJ(3,1)
        ELSEIF(N.EQ.1) THEN
            DET=DJ(1,1)
        END IF
        RETURN 
      END SUBROUTINE


C
C     Normalize vector A and return its magnitude as AMAG
C
      SUBROUTINE KUNITV(N,A,AMAG)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION A(N)
        
        AMAG=0.D0
        DO I=1,N
          AMAG=AMAG+A(I)*A(I)
        END DO
        AMAG = DSQRT(AMAG)
        DO I=1,N
          A(I) = A(I)/AMAG
        END DO
        RETURN
      END SUBROUTINE

C
C   Cross Product for 3d vector A X B
C
      SUBROUTINE CROSS(A,B,AXB)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION A(3),B(3)
        
        DIMENSION AXB(3)

        AXB(1) = A(2)*B(3) - A(3)*B(2)
        AXB(2) = A(3)*B(1) - A(1)*B(3)
        AXB(3) = A(1)*B(2) - A(2)*B(1)
        RETURN
      END SUBROUTINE
C
C***********************************************************************
C
C     *DAMAGE INITIATION, CRITERION=USER, PROPERTIES=m, FAILURE MECHANISMS=n   
C
C
      SUBROUTINE UDMGINI(FINDEX,NFINDEX,FNORMAL,NDI,NSHR,NTENS,PROPS,
     1 NPROPS,STATEV,NSTATV,STRESS,STRAIN,STRAINEE,LXFEM,TIME,
     2 DTIME,TEMP,DTEMP,PREDEF,DPRED,NFIELD,COORDS,NOEL,NPT,LAYER,
     3 KSPT,KSTEP,KINC,KDIRCYC,KCYCLELCF,TIMECYC,SSE,SPD,SCD,SVD,
     4 SMD,JMAC,JMATYP,MATLAYO,LACCFLA,CELENT,DROT,ORI)
  
      INCLUDE 'ABA_PARAM.INC'
        
      DIMENSION FINDEX(NFINDEX),FNORMAL(NDI,NFINDEX),COORDS(*),
     1 STRESS(NTENS),STRAIN(NTENS),STRAINEE(NTENS),PROPS(NPROPS),
     2 STATEV(NSTATV),PREDEF(NFIELD),DPRED(NFIELD),TIME(2),JMAC(*),
     3 JMATYP(*),DROT(3,3),ORI(3,3)
      
      !user coding to define FINDEX,  and  FNORMAL
      DIMENSION PS(3), AN(3,3), WT(6)
      

      ! ROTATE THE STRESS TO GLOBAL SYSTEM IF THERE IS ORIENTATION
      CALL ROTSIG(STRESS,ORI,WT,1,NDI,NSHR)
      
      IF(MOD(NOEL,100).EQ.1) THEN
        PRINT *,'KSTEP,KINC,DTIME,TIME'
        PRINT *,KSTEP,KINC,DTIME,TIME
        PRINT *,'NFINDEX,NDI,NSHR,NTENS,NPROPS,PROPS'
        PRINT *,NFINDEX,NDI,NSHR,NTENS,NPROPS,PROPS
        PRINT *,'NOEL,NPT,LAYER,KSPT,CELENT,LXFEM'
        PRINT *,NOEL,NPT,LAYER,KSPT,CELENT,LXFEM
        PRINT *,'ORI,DROT,STRESS'
        PRINT *,ORI,DROT,STRESS
      END IF

      ! MAXIMUM PRINCIPAL STRESS CRITERION
      PS=0.0D0
      CALL SPRIND(WT,PS,AN,1,NDI,NSHR)
      SIG1 = PS(1)
      KMAX=1
      DO K1 = 2, NDI  
         IF(PS(K1).GT.SIG1) THEN
            SIG1 = PS(K1)
            KMAX = K1
         END IF
      END DO
      FINDEX(1) = SIG1/PROPS(1)
      DO K1=1, NDI
      	FNORMAL(K1,1) = AN(KMAX,K1)
      END DO

      ! QUADRATIC TRACTION-INTERACTION CRITERION
      FINDEX(2)=(STRESS(1)/PROPS( 2))**2.0+(STRESS(NDI+1)/
     $     PROPS(3))**2.0+(STRESS(NDI+2)/PROPS(4))**2.0

      FINDEX(2)=sqrt(FINDEX(2))

      DO K1=1, NDI
      	FNORMAL(K1,2)=ORI(K1,1)
      END DO
      RETURN
      END


C
C***********************************************************************
C
C   将会在分析步开始时、分析步结束时、增量步开始时、增量步结束时等被调用
C   用于
C     manage user-defined external databases 
C     calculate model-independent history information
C   参考    Abaqus/Standard subroutines> 1.1.31 UEXTERNALDB
C
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
C
      ! user coding to set up the FORTRAN environment, open files, close files, 
      ! calculate user-defined model-independent history information,
      ! write history information to external files,
      ! recover history information during restart analyses, etc.
      ! do not include calls to utility routine XIT
      CHARACTER*256 JOBNAME,CP,DIR
      INTEGER LJ,LD 

      IF(LOP.EQ.0.OR. LOP.EQ.4) THEN
        IF(LOP.EQ.0) THEN 
          ! Start of Analysis
          PRINT *,'Analysis Start!'
        ELSE
          ! Beginning of a restart analysis
          PRINT *,'Restart Analysis Start!'
        END IF

        CALL GETJOBNAME(JOBNAME, LJ )
        DIR='E:/UMAT-OUTPUT/'
        LD=LEN_TRIM(DIR)
        
        !OPEN(101,FILE='E:\dadelamination.TXT',STATUS='UNKNOWN')

      ELSE IF(LOP.EQ.3) THEN
        ! End of Analysis
        
        ! CLOSE(101)
        PRINT *,'Analysis End!'
      ELSE IF(LOP.EQ.1) THEN
        ! Start of the current analysis increment
        WRITE(*,'(A,I3,A,I3,\)') ' Start Step',KSTEP,' Increment ',KINC
        WRITE(*,'(A,F6.3,\)') ' Step Start Time:',TIME(1)
        WRITE(*,'(A,F6.3,\)') ' Step End Time:',TIME(1)+DTIME
        WRITE(*,'(A,F6.3,\)') ' Total Time:',TIME(2)
        WRITE(*,'(A,F6.3)') ' DTime:',DTIME
      
      ELSE IF(LOP.EQ.2) THEN
        ! End of the current analysis increment
        WRITE(*,'(A,I3,A,I3,\)') ' End Step',KSTEP,' Increment ',KINC
        WRITE(*,'(A,F6.3,\)') ' Step Time:',TIME(1)
        WRITE(*,'(A,F6.3,\)') ' Total Time:',TIME(2)
        WRITE(*,'(A,F6.3)') ' DTime:',DTIME
      
      ELSE IF(LOP.EQ.5) THEN
        ! Start of a step
        PRINT *,'Step',KSTEP,' Start!'

      ELSE IF(LOP.EQ.6) THEN
        ! End of a step
        PRINT *,'Step',KSTEP,' End!'

      END IF

      RETURN
      END

C
C***********************************************************************
C
C    在每个增量的结束时调用
C    检查是否有单元发生纤维破坏STATEV(14)==1，如果有纤维发生破坏，停止计算
C
C     要使URDFIL被调用，需要在inp文件step 段中添加 
C     *EL File
C     SDV
C     CAE界面里： Model->edit Keyword下添加
C     参见Analysis User Guide>4.1.2 Output to the data and results files
C     
      SUBROUTINE URDFIL(LSTOP, LOVRWRT, KSTEP, KINC,DTIME,TIME)

        INCLUDE 'ABA_PARAM.INC'
        DIMENSION ARRAY(513), JRRAY(NPRECD, 513), TIME(2)
        EQUIVALENCE (ARRAY(1),JRRAY(1,1)) ! 分别用于获取浮点数和整数
        
        DIMENSION STATEV(15)

        WRITE(*,'(A,I3,A,I3,\)') ' URDFIL Step:',KSTEP,' Increment:',KINC
        WRITE(*,'(A,F6.3,\)') ' Step Time:',TIME(1)
        WRITE(*,'(A,F6.3,\)') ' Total Time:',TIME(2)
        WRITE(*,'(A,F6.3)') ' DTime:',DTIME

        LSTOP=0   ! can't be termianated
        LOVRWRT=1 ! this inc can be over write

        CALL POSFIL(KSTEP, KINC, ARRAY, JRCD)
        DO K1=1,999999
          CALL DBFILE(0,ARRAY,JRCD)
          IF (JRCD .NE. 0) GO TO 110
          KEY=JRRAY(1,2)
          IF(KEY.EQ.1) THEN 
            IELM=JRRAY(1,3)
          END IF

          IF(KEY.EQ.5) THEN
            DO I=1,15
              STATEV(I)=ARRAY(I+2)
            END DO
            IF(STATEV(14).GE.1.0 ) THEN ! .OR. STATEV(15).GE.1.0
              PRINT '(I6,15F10.4)',IELM,STATEV
              ! 当不处于分析步的结尾时
              IF(TIME(1).NE.1.0) THEN
                LSTOP=1   ! terminate analysis 
                LOVRWRT=0 ! this inc can't be over write
              END IF
            END IF
          END IF
        END DO
  110   CONTINUE

        PRINT '(A,I3,A,I3)', ' LSTOP: ',LSTOP, ' LOVRWRT: ',LOVRWRT
      END SUBROUTINE