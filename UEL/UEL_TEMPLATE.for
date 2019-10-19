C
C
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
        
        PARAMETER(nINTP=2)
        DIMENSION GP(MCRD,nINTP),GW(nINTP),PT(MCRD)
        DIMENSION shapeN(nNODE),dNdr(MCRD,nNode),DJ(MCRD,MCRD)

        ! Initialize RHS and AMATRX
        RHS(1:NDOFEL,1)=0.D0
        AMATRX= 0.D0
        PRINT *,DTIME,NDOFEL,MLVARX
        ! set up gauss integration points GP and weights  GW
        CALL KGUASS(NNODE,MCRD,nINTP,GP,GW)

        DO Intp=1,nINTP
          ! operation on every gauss point
          W=GW(Intp)
          DO J=1,MCRD
            PT(J)=GP(J,Intp)
          END DO
          
        END DO
      RETURN
      END

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
