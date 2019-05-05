C  Node Ordering
C   4---5---6
C   |   |   |
C   1---2---3
C	2-5 is crack tip
C	
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,
     2     DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     3     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     4     JPROPS,NJPROP,PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0 )
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
C
      DX=COORDS(1,5)=COORDS(1,2)
	  DY=COORDS(2,5)=COORDS(2,2)
	  DL2=DX**2+DY**2
	  DL=SQRT(DL2)
      SKX=PROPS(1)/DL
	  SKY=PROPS(2)/DL
      DO K1=1,6
		RHS(K1,1)=ZERO
		DO K2=1,6
		  AMATRX(K1,K2)=ZERO
		END DO
	  END DO
	  
	  DO K1=1,5,2
		AMATRX(K1,K1)=SKX
		AMATRIX(K1+6,K1+6)=SKX
		AMATRIX(K1,K1+6)=-SKX
		AMATRIX(K1+6,K1)=-SKX
		
		AMATRX(K1+1,K1+1)=SKY
		AMATRIX(K1+7,K1+7)=SKY
		AMATRIX(K1+1,K1+7)=-SKY
		AMATRIX(K1+7,K1+1)=-SKY
		
	  END DO
	  
	  DO K1=1,3
		RHS(2*K1-1,1)=-SKX*(U(2*K1-1)-U*(2*K1+5))
		RHS(2*K1,1)=-SKY*(U(2*K1)-U*(2*K1+6))
		RHS(2*K1+5,1)=SKX*(U(2*K1-1)-U*(2*K1+5))
		RHS(2*K1+6,1)=SKY*(U(2*K1)-U*(2*K1+6))		
	  END DO
	  
	  DO K1=1,12
		PRINT(RHS(K1,1))
	  END DO
      RETURN
      END         
