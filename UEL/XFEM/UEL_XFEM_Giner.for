       SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     * PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     * KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     * LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
c
       INCLUDE 'ABA_PARAM.INC'
c
c      User subroutine for computation of element stiffness matrix and
c      the element force vector
c
c      ABAQUS defined variables:
c
       DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     * SVARS(NSVARS), ENERGY(8), COORDS(MCRD,NNODE), U(NDOFEL),
     * DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     * JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     * PREDEF(2,NPREDF,NNODE), LFLAGS(*), JPROPS(*)
c
c      Important variables (list not exhaustive)
c           JELEM  Current element number
c	      AMATRX Element stiffness matrix (element contribution to the stiffness 
c                  matrix of the overall system of equations)
c           RHS    Element residual force vector (element contribution to the right-hand-side 
c                  vector of the overall system of equations)
c           F      Element force vector (AMATRX times the updated solution for nodal dofs)
c           E      Young's modulus
c           Nu     Poisson's ratio
c           PSS    1 - Plane stress
c                  2 - Plane strain
c           orderQ Vector that stores the following quadrature orders:
c		  orderQ(1) = Quadrature order for quadrilaterals (in each direction)
c                   NOTE: only for non-subdivided enriched elements (quadrilaterals) 
c           orderQ(2) = Quadrature order for triangles (total points)
c                   NOTE: only for enriched elements subdivided into triangles 
c           orderQ(3) = Quadrature order for quadrilaterals (in each direction)
c             NOTE: only for enriched elements subdivided into 2 quadrilaterals (elemX = type 4)

c           dimens Dimension of the physical domain: 2=2D
c                  Actually, it should suffice with the ABAQUS variable MCRD, but ABAQUS automatically
c                  sets MCRD=3, even though for a 2D problem, since we are using the third and further 
c                  available dof for the enriched nodes.

c           NNODE  Number of nodes per element
c           NelmX  Number of enriched elements
c           NnodX  Number of nodes that belong to enriched elements
c           TypeX  Matrix that stores the key number to the type of node in an enriched element: 
c                   0 - Non-enriched (2 dof)
c                   1 - Heaviside enrichment (4 dof)
c                   2 - Crack tip enrichment (10 dof)
  
c           TypeXe Vector that stores the key in TypeX for the nodes of the current element
c           ix     Vector that stores the node numbers of the current element (connectivity)
c           Xe(8)	 X nodal coordinates of the current element
c                  (it is duplicated to ease the counting from the 4th to the 1st node)
c		  Ye(8)  Y nodal coordinates of the current element

c           NCracks Number of cracks
c           NCP    Number of crack path points (vertices)
c           maxNCP Maximum number of crack path points (vertices)
c           XYC    Matrix that stores the coordinates of crack path points
c           XYC0   Crack tip coordinates associated with the crack tip enriched element
c           XYCPrev Crack tip coordinates associated with the previous crack path point

             
c           gint   Total number of integration points (either with or without subdivision)
c		  flag   Subdivisión indicator (1 for subdivision)
c           mpg    Maximum expected number of integration points for an enriched element
c           sg     Matrix that stores the coordinates and weights of the integration points
c           xypg   Matrix that stores the coordinates of the integration points
c           Dist   Matrix that stores distances to crack from nodes of enriched elements
c                  (this information is previously preprocessed, for example in Matlab)  
c           ElemGG Matrix that stores information about the elements to be enriched, type of 
c                  crack intersection and points of interesection
c                  (this information is previously preprocessed, for example in Matlab)  
c
c           BatG   Matrix that stores the [B] matrix for each enriched element at Gauss points
c           DBatG  Matrix that stores the [D][B] matrix for each enriched element at Gauss points
c           JatG   Vector that stores the Jacobian = det([J]) for each enriched element at Gauss points

                
c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

 
c      Declaration of variables for XFEM user element
c
       CHARACTER*256 OUTDIR  ! to read the working directory
       INTEGER LENOUTDIR  ! working directory string length
       
       INTEGER i,j,k,PSS,orderQ(3),gint,flag,dimens
       INTEGER NCracks,maxNCP,NelmX,NnodX,TypeXe(NNODE),ix(NNODE)
       INTEGER,PARAMETER :: mpg=1650 ! up to more than 40x40 Gauss integration points per element
       INTEGER,ALLOCATABLE:: TypeX(:,:),NCP(:)

       REAL*8 E, Nu
       REAL*8 F(NDOFEL)
       REAL*8 sg(3,mpg),xypg(2,mpg),Xe(8),Ye(8),XYC0(2),XYCPrev(2)       
C
       REAL*8, ALLOCATABLE:: XYC(:,:,:),Dist(:,:),ElemGG(:,:)
       REAL*8, ALLOCATABLE:: BatG(:,:),DBatG(:,:),JatG(:)


c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


c23456789012345678901234567890123456789012345678901234567890123456789012 
 


c      Read real and integer properties set at the ABAQUS input file 
       E = PROPS(1)
       Nu = PROPS(2)
       PSS = JPROPS(1)
       orderQ(1) = JPROPS(2)
       orderQ(2) = JPROPS(3)
       orderQ(3) = JPROPS(4)
       dimens = JPROPS(5)

c      Read the working directory
       CALL GETOUTDIR(OUTDIR,LENOUTDIR)

c      *************************************************************************   
c      **** Read information previously preprocessed, for example in Matlab ****
c      *************************************************************************  
      
c      Read number of cracks, max number of crack path points, 
c           number of enriched elements and enriched nodes.   
       OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\files\GGInfoX')
       READ(68,*) NCracks,maxNCP,NelmX,NnodX
       CLOSE(68)
                     
c      Allocate dimensions
       ALLOCATE (TypeX(NnodX,2), NCP(NCracks))
       ALLOCATE (XYC(NCracks,maxNCP,2), Dist(NnodX,3), ElemGG(NelmX,10))

c      Read coordinates of path points for each crack
       OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\files\GGXYC')
       DO i=1,NCracks
         READ(68,*) NCP(i)
         DO j=1,NCP(i)
           READ(68,*) (XYC(i,j,k),k=1,2)  
         END DO
       END DO   	  
       CLOSE(68)

c      Read list of enriched nodes, type of enrichment and distances
       OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\files\GGnodeX')
       DO i=1,NnodX
         READ(68,*) (TypeX(i,j),j=1,2),(Dist(i,j),j=2,3)
         Dist(i,1)=TypeX(i,1)
       END DO   	  
       CLOSE(68)

c      Read list of enriched elements, type of enrichment and intersection points
       OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\files\GGelemX')
       DO i=1,NelmX
         READ(68,*) (ElemGG(i,j),j=1,10)
       END DO   	  
       CLOSE(68)


c      Call initializing routines for matrix and vectors
      CALL initializeM(RHS,NDOFEL,NRHS)
      CALL initializeM(AMATRX,NDOFEL,NDOFEL) 
      CALL initializeV(ENERGY,8)
      CALL initializeV(SVARS,NSVARS)


      
c      Verification of element type (type=12 for enriched element)  
c
       IF (JTYPE.EQ.12) THEN
c       **************************************
c       *    4 NODE ENRICHED ELEMENT WITH    *
c       *    UP TO 12 DOF/NODE FOR X-FEM     *			
c       **************************************
    
        IF (LFLAGS(1).EQ.71) THEN
c       Coupled thermal-stress, steady state analysis
         IF (LFLAGS(3).EQ.1) THEN
c        Normal implicit time incrementation procedure.
c        User subroutine UEL must define the residual vector in RHS 
c        and the stiffness matrix in AMATRX

c          Routine that defines the location of integration points according to
c          the appropriate subdivision. This enables to know the total number of 
c          integration points for the current element, stored in gint, and whether
c          the element is subdivided for integration (flag=1) or not.
           CALL int2d_X(JELEM,NelmX,ElemGG,MCRD,NNODE,COORDS,orderQ,
     *       	NCracks,maxNCP,NCP,XYC,gint,sg,Xe,Ye,flag,mpg,xypg,
     *        XYC0,XYCPrev)

c          Allocate dimensions once the total number of integration points gint is known
             ALLOCATE(BatG(3*gint,NDOFEL),DBatG(3*gint,NDOFEL),JatG(gint))
           CALL initializeM(BatG,3*gint,NDOFEL)
             CALL initializeM(DBatG,3*gint,NDOFEL) 
             CALL initializeV(JatG,gint)

c          Search of the enrichment type for the nodes of the current element.
c          The keys to the enrichment types are stored in the element vector TypeXe
           CALL TypeXelement(OUTDIR,LENOUTDIR,JELEM,NNODE,NelmX,
     *        		 ix,TypeXe)            

c		 Element stiffness matrix computation, stored in AMATRX	
           CALL K_U12(E,Nu,AMATRX,NDOFEL,NNODE,dimens,MCRD,
     *          COORDS,PSS,NnodX,ix,TypeXe,Dist,XYC0,XYCPrev,
     *          gint,sg,Xe,Ye,flag,BatG,DBatG,JatG)

c          Routine that multiplies AMATRX times U to obtain the force vector F
c          at the end of the current increment  
           CALL MULT_V(AMATRX,NDOFEL,NDOFEL,U,F,NDOFEL)

c          Compute the residual force vector
           DO I=1,NDOFEL
               RHS(I,1) = RHS(I,1) - F(I)
           END DO

c          Compute stresses at Gauss points for post-processing purposes
c          Store them as SVARS for output to the results file (.fil)
           CALL  SVARS_U12(JTYPE,JELEM,SVARS,NSVARS,U,NDOFEL,BatG,
     *                    DBatG,JatG,gint,mpg,xypg)
         END IF
        END IF       
       END IF


       RETURN
       END



C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                         C
C   Contiene subrutinas de calculo utilizadas por         C
C   todos los elementos                                   C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes


      SUBROUTINE CALC_D(PSS,D,E,nu)
      IMPLICIT NONE
      INTEGER PSS
      REAL*8 D(3,3), E, nu
C 
      IF (PSS.EQ.1) THEN
C         Tension plana
          D(1,1) = 1.0d0*E/(1.0d0-nu*nu)
          D(1,2) = nu*E/(1.0d0-nu*nu)
          D(1,3) = 0.0d0
          D(2,1) = nu*E/(1.0d0-nu*nu)
          D(2,2) = 1.0d0*E/(1.0d0-nu*nu)
          D(2,3) = 0.0d0
          D(3,1) = 0.0d0
          D(3,2) = 0.0d0
          D(3,3) = E/((1.0d0+nu)*2.0d0)
      ELSEIF (PSS.EQ.2) THEN
C         Deformacion plana
          D(1,1) = (1.0d0-nu)*E/((1.0d0-2.0d0*nu)*(1.0d0+nu))
          D(1,2) = nu*E/((1.0d0-2.0d0*nu)*(1.0d0+nu))
          D(1,3) = 0.0d0
          D(2,1) = nu*E/((1.0d0-2.0d0*nu)*(1.0d0+nu))
          D(2,2) = (1.0d0-nu)*E/((1.0d0-2.0d0*nu)*(1.0d0+nu))
          D(2,3) = 0.0d0
          D(3,1) = 0.0d0
          D(3,2) = 0.0d0
          D(3,3) = E/((1.0d0+nu)*2.0d0)
      END IF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Not used

      SUBROUTINE CALC_C(PSS,D,E,nu)
      IMPLICIT NONE
      INTEGER PSS
      REAL*8 D(3,3), E, nu

      IF (PSS.EQ.1) THEN
          D(1,1) = 1.0d0/E
          D(1,2) = -nu/E
          D(1,3) = 0.0d0
          D(2,1) = -nu/E
          D(2,2) = 1.0d0/E
          D(2,3) = 0.0d0
          D(3,1) = 0.0d0
          D(3,2) = 0.0d0
          D(3,3) = 2.0d0*(1.0d0+nu)/E
      ELSEIF (PSS.EQ.2) THEN
          D(1,1) = (1.0d0-nu*nu)/E
          D(1,2) = -nu*(1.0d0+nu)/E
          D(1,3) = 0.0d0
          D(2,1) = -nu*(1.0d0+nu)/E
          D(2,2) = (1.0d0-nu*nu)/E
          D(2,3) = 0.0d0
          D(3,1) = 0.0d0
          D(3,2) = 0.0d0
          D(3,3) = 2.0d0*(1.0d0+nu)/E
      END IF
      RETURN
      END    
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Not used 

      SUBROUTINE MULT_M(S,A,FA,CA,B,FB,CB)
      INTEGER FA, FB, CA, CB
      REAL*8 S(FA,CB), A(FA,CA), B(FB, CB)
      INTEGER I, J, K
      DO 10 I=1,FA
          DO 11 J=1,CB
              S(I,J) = 0.0D0
              DO 12 K=1,CA
                  S(I,J) = S(I,J) + A(I,K)*B(K,J)
 12           CONTINUE
 11       CONTINUE
 10   CONTINUE
      RETURN
      END                        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
      SUBROUTINE MULT_V(A,FA,CA,B,S,MS)
      INTEGER FA, CA
      REAL*8 S(MS), A(FA,CA), B(CA)
      INTEGER I, K
      DO 10 I=1,FA
          S(I) = 0.0D0
          DO 12 K=1,CA
              S(I) = S(I) + A(I,K)*B(K)
 12       CONTINUE
 10   CONTINUE
      RETURN
      END                        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
      SUBROUTINE INV_J(J,Ngdl)
      INTEGER Ngdl
      REAL*8 J(Ngdl,Ngdl)
C
      REAL*8 V, Det
      IF (Ngdl.EQ.2) THEN
          Det = J(1,1)*J(2,2) - J(2,1)*J(1,2)
          V = J(1,1)/Det
          J(1,1) = J(2,2)/Det
          J(2,2) = V
          J(1,2) = -J(1,2)/Det
          J(2,1) = -J(2,1)/Det
      END IF
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Not used 

      SUBROUTINE CALC_BT(B,BT,NF,NC)	
      IMPLICIT NONE
C     Subrutina de calculo de la matriz transpuesta
C
      INTEGER NF, NC
      REAL*8 B(NF,NC), BT(NC,NF)
      INTEGER I, J
C
      DO 10 I=1,NC
          DO 20 J=1,NF
              BT(I,J) = B(J,I)
 20       CONTINUE
 10   CONTINUE
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
      SUBROUTINE initializeM(M,A,B)	
      IMPLICIT NONE
C     initialize matrix
C
      INTEGER A,B,I,J
      REAL*8 M(A,B)
      
      DO  I=1,A
         DO  J=1,B
              M(I,J) = 0.0D0 
         END DO	     
      END DO
       
      
      RETURN
      END

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
      SUBROUTINE initializeV(V,A)	
      IMPLICIT NONE
C     initialize vector
C
      integer A,I
      Real*8 V(A)

      DO  I=1,A
         V(I) = 0.0D0      
      end do
       
      
      RETURN
      END


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

      subroutine int2d_X(n,NelmX,ElemGG,MCRD,nel,COORDS,orderQ,
     *	     NCracks,maxNCP,NCP,XYC,gint,sg,Xe,Ye,flag,mpg,xypg,
     *         XYC0,XYCPrev)

c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Calcular el número, coordenadas locales y pesos de los puntos de Gauss.
c		     Distingue entre subdivisión por triángulos y cuadriláteros 
c     --------------------------------------

c      Inputs:

c			ElemGG		Matriz de caracterización de elementos X
c			n			Número del elemento actual
c			COORDS(MCRD,nel)	Coordenadas Nodales del elemento:  hr(np(44))
c			l			Orden de cuadratura para cada direcc en cuadr (orderQ)
c                          y definido en UEL property
c			numpgt		Orden de cuadratura para triángulos (orderQ) y definido en UEL property
c			numpgc   	Orden de cuadratura para cada direcc en cuadr (orderQ)
c                          y definido en UEL property para elementos tipo 4

c			XYC(NCP,2)	Coordenadas de los puntos de Grieta	
c 		    nel         núm. de nodos del elemento (en principio 4 para cuadrilat)
                                    
c      Outputs:			
C             l     orden de cuadratura (nº ptos en cada dirección en cuadril)
c                   OJO, sólo para los elementos cuadriláteros XFEM no intersectados
c			gint		Número de Puntos de Gauss
c			sg(3,*)		Coordenadas y pesos
c			Xe(8)		Vector de coordenadas nodales x para el elemento
c			Ye(8)		Vector de coordenadas nodales y para el elemento
c						(Se utilizan duplicados para facilitar los contadores
c						 en el salto del cuarto nodo al primero)
c			flag        Indicador de que el elemento se ha subdividido
c						para integración (si flag=1)
c             mpg         Número máximo de ptos de Gauss esperado (de momento, mpg=50)
c             xypg        Coordenadas globales de los puntos de Gauss 

c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

c     include  'comblk.h'
c 	include  'eldata.h'
c 	include  'xfem2d.h'
c     include  'iofile.h'

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2

c	logical esx
      integer NCracks,maxNCP,NCP(NCracks)
      integer i,j,k,l,m,NumSub,gint,flag,NelmX,nel,n,MCRD,mpg
      integer numpgt,numpgc,orderQ(3)
      real*8  ElemGG(NelmX,10),ElemGGe(10),COORDS(MCRD,nel),
     &    	XYC(NCracks,maxNCP,2),XYC0(2),XYCPrev(2),xpg,ypg,
     &		Xe(8),Ye(8),SubXe(10,4),SubYe(10,4),sg(3,*),sgt(:,:),s,t,
     &		Ni(4),Area,dNl(4,2),Jac(2,2),sgc(:,:),xypg(2,mpg)
      allocatable :: sgt,sgc


c     Paso a los nombres de las variables locales para los órdenes de integreación
      
      l=orderQ(1)
      numpgt=orderQ(2)
      numpgc=orderQ(3)

c	Inicializacion de variables
      
      flag=0

      do i=1,2
        do j=1,mpg
          xypg(i,j)=0.0d0
        end do
      end do

c	Elemento estándar
c	if (.not.esx) then
c		call int2d(l,gint,sg)
c		return
c	endif

c	Lectura de características del elemento actual de ElemGG
c     Notar que ElemGGe es con tipo enteros           
c
      do i=1,NelmX
            if (int(ElemGG(i,1)).eq.n)then
                  ElemGGe=ElemGG(i,:)
              exit
            end if
      end do
      
c      write(dfich2,*) 'Nº elem=',n,'   ElemGGe(1)=',ElemGGe(1)
c      write(dfich2,*) '                ElemGGe(2)=',ElemGGe(2)


c	Coordenadas del Extremo de Grieta relacionado con el elemento

      if (ElemGGe(9).eq.0) then 
c       El elemento no es de extremo de grieta y no necesitará el extremo de grieta
        XYC0(1)=0
        XYC0(2)=0
        XYCPrev(1)=0
        XYCPrev(2)=0
      else 
        if (ElemGGe(10).eq.1) then
c       Si el extremo de grieta es el 1 (primeras coordenadas en la tabla) 
        XYC0=XYC(ElemGGe(9),1,1:2)
        XYCPrev=XYC(ElemGGe(9),2,1:2)
        else
c       El extremo de grieta es el último de la tabla
        XYC0=XYC(ElemGGe(9),NCP(ElemGGe(9)),1:2)
        XYCPrev=XYC(ElemGGe(9),NCP(ElemGGe(9))-1,1:2)
        endif	  
      endif
      
      


c	Vectores duplicados de coordenadas nodales
                  
      do j=0,1
            do i=1,4
            Xe(i+4*j)=COORDS(1,i)
            Ye(i+4*j)=COORDS(2,i)
            end do
      end do

c     write(dfich2,*)  'En xint2D_X:'
c	write(dfich2,*)
c	write(dfich2,*) 'Xe(i)=', (Xe(i),i=1,8)
c	write(dfich2,*)
c	write(dfich2,*) 'Ye(i)=', (Ye(i),i=1,8)
c	write(dfich2,*)
      

c	Se distinguen 3 posibilidades para el cálculo de los PGauss:
c	 División por triángulos, por cuadriláteros o integrar el elemento completo	

c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c	Elementos no intersectados (ElemX con parte entera tipo 1 o -1)
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    


      if (ElemGGe(2).le.1) then

            call int2d(l,gint,sg,nel) !Cuadratura de orderQ x orderQ = lxl, en pº,5x5

c          -- coordenadas x,y  globales de punto de gauss
            do i= 1,gint
c				  --- Funciones de forma
                          Ni(1)= (1-sg(1,i))*(1-sg(2,i))/4
                          Ni(2)= (1+sg(1,i))*(1-sg(2,i))/4
                          Ni(3)= (1+sg(1,i))*(1+sg(2,i))/4
                          Ni(4)= (1-sg(1,i))*(1+sg(2,i))/4

c                   -- coordenadas x,y  globales de punto de gauss
                          xpg= dot_product(Ni,Xe(1:4))
                      ypg= dot_product(Ni,Ye(1:4))

                      xypg(1,i)=xpg
                    xypg(2,i)=ypg
           end do



c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c	Elementos con triángulos.
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    


c	(El dominio del elemento se subdivide en triángulos alineados con la 
c	 geometría de la grieta sobre los cuales se calculan las coordenadas  
c	 y pesos de los puntos de Gauss)

      elseif (ElemGGe(2).ne.4)then ! Entonces la subdivisión es a base de triángulos

            flag=1

c		Número de puntos de integración en cada triángulo
c		numpg=-3 !Cuadratura de 3 puntos interiores
c		numpg=7  !Cuadratura de 7 puntos interiores
c		numpg=12 !Cuadratura de 12 puntos interiores

            allocate(sgt(4,abs(numpgt)))

            call tint2d(numpgt,gint,sgt) 

c		Dividir el elemento
            call subelem(ElemGGe,Xe,Ye,XYC0(1),XYC0(2),NumSub,SubXe,
     &                 SubYe,n,nel)

c         Coordenadas globales y locales  de punto de gauss
            do j=1,NumSub
            do i=1,gint
c			--- coordenadas x, y globales de punto de gauss
                xpg= dot_product(sgt(1:3,i),SubXe(j,1:3))
              ypg= dot_product(sgt(1:3,i),SubYe(j,1:3))		
                
              xypg(1,gint*(j-1)+i)=xpg
                xypg(2,gint*(j-1)+i)=ypg

c		    --- Area del subelemento			  
                Area= (SubXe(j,2)*SubYe(j,3)+SubXe(j,1)*SubYe(j,2)+
     &		       SubYe(j,1)*SubXe(j,3)-SubYe(j,1)*SubXe(j,2)-
     &               SubXe(j,1)*SubYe(j,3)-SubYe(j,2)*SubXe(j,3))/2.0d0
                  

c		    -- coordenadas s,t locales al elemento de pg 
              call invcuad4(Xe(1:4),Ye(1:4),xpg,ypg,s,t,Ni)
                sg(1,(i+(j-1)*gint))=s
                sg(2,(i+(j-1)*gint))=t
                sg(3,(i+(j-1)*gint))=sgt(4,i)*Area

              end do !i
          end do !j

c		Número de puntos de integración
          gint=NumSub*gint


c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c	Elementos subdivididos con cuadriláteros (tipo 4)
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

c	(El dominio del elemento se subdivide en dos cuadriláteros alineados con la 
c	 geometría de la grieta sobre los cuales se calculan las coordenadas  
c	 y pesos de los puntos de Gauss)

      elseif (ElemGGe(2).eq.4)then

            flag=1

c		Orden de la cuadratura (puntos de gauss en 1d)
c		numpg=l !La estándar de FEAP
c		numpg=5 !Cuadratura de 5x5
          

            allocate(sgc(3,numpgc*numpgc))

            call int2d(numpgc,gint,sgc,nel)

          if (gint.gt.numpgc*numpgc) then
            write(dfich2,*)
            write(dfich2,*) '*ERROR* El número de columnas de sgc ',
     &                   'en int2d_X.f debe aumentarse a ', gint
            stop
          endif

c		Dividir el elemento
            call subelem(ElemGGe,Xe,Ye,XYC0(1),XYC0(2),NumSub,SubXe
     &		,SubYe,n,nel)

c          -- coordenadas x,y  globales de punto de gauss
            do j= 1,NumSub
                  do i= 1,gint

c				  --- Funciones de forma
                          Ni(1)= (1-sgc(1,i))*(1-sgc(2,i))/4;
                          Ni(2)= (1+sgc(1,i))*(1-sgc(2,i))/4;
                          Ni(3)= (1+sgc(1,i))*(1+sgc(2,i))/4;
                          Ni(4)= (1-sgc(1,i))*(1+sgc(2,i))/4;

                          dNl(1,1)= -(1-sgc(2,i))/4;
                          dNl(1,2)= -(1-sgc(1,i))/4;
                          dNl(2,1)=  (1-sgc(2,i))/4;
                          dNl(2,2)= -(1+sgc(1,i))/4;
                          dNl(3,1)=  (1+sgc(2,i))/4;
                          dNl(3,2)=  (1+sgc(1,i))/4;
                          dNl(4,1)= -(1+sgc(2,i))/4;
                          dNl(4,2)=  (1-sgc(1,i))/4;

c			      --- Matriz Jacobiana
                          do k= 1,2
                                Jac(k,1)= 0.0;
                                Jac(k,2)= 0.0;
                                do m= 1,nel
                                      Jac(k,1)= Jac(k,1) + dNl(m,k)*SubXe(j,m);
                                      Jac(k,2)= Jac(k,2) + dNl(m,k)*SubYe(j,m);
                                end do
                          end do
                          
                          Area=Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)

c                   -- coordenadas x,y  globales de punto de gauss
                          xpg= dot_product(Ni,SubXe(j,1:4));
                      ypg= dot_product(Ni,SubYe(j,1:4));

                      xypg(1,gint*(j-1)+i)=xpg
                    xypg(2,gint*(j-1)+i)=ypg

                        
c				  -- coordenadas s,t locales al elemento de pg
                          call invcuad4(Xe(1:4),Ye(1:4),xpg,ypg,s,t,Ni)

                          sg(1,(i+(j-1)*gint))=s
                             sg(2,(i+(j-1)*gint))=t
                          sg(3,(i+(j-1)*gint))=sgc(3,i)*Area

                end do  !i
          end do !j

c		Número de puntos de integración
            gint=NumSub*gint

                  
      end if

      if (gint.gt.mpg) then
          write(dfich2,*) '*ERROR* El parámetro mpg=',mpg,
     &                   '  debe aumentarse a mpg=',gint
          stop
      end if


      return
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

      subroutine int2d(l,gint,sg,nel)


c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of points/direction
c         nel     - num. de nodos en el elemento      

c      Outputs:
c         gint    - Total number of points
c         sg(3,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2


c      include  'eldata.h'
c      include  'iofile.h'

      integer   i,j,k,l,gint, lr(9),lz(9),lw(9),nel
      real*8    g,h, third, sg(3,*),ss(5),ww(5)

      save

      data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data      lw/4*25,4*40,64/
      data      third / 0.3333333333333333d0 /

c     Set number of total points

      gint = l*l

c     5 pt. integration

      if(l.eq.0) then

        gint = 5
        g    = sqrt(0.6d0)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 5.d0/9.d0
        end do ! i

        sg(1,5) = 0.0d0
        sg(2,5) = 0.0d0
        sg(3,5) = 16.d0/9.d0

c     1x1 integration

      elseif(l.eq.1) then
        sg(1,1) = 0.d0
        sg(2,1) = 0.d0
        if(nel.eq.3) sg(2,1) = -third
        sg(3,1) = 4.d0

c     2x2 integration

      elseif(l.eq.2) then
        g = sqrt(third)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 1.d0
        end do ! i

c     3x3 integration

      elseif(l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.d0/81.d0
        do i = 1,9
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = h*lw(i)
        end do ! i

c     4x4 integration

      elseif(l.eq.4) then
        g     = sqrt(4.8d0)
        h     = third/g
        ss(1) = sqrt((3.d0+g)/7.d0)
        ss(4) = - ss(1)
        ss(2) = sqrt((3.d0-g)/7.d0)
        ss(3) = -ss(2)
        ww(1) = 0.5d0 - h
        ww(2) = 0.5d0 + h
        ww(3) = 0.5d0 + h
        ww(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! i

c     5x5 integration

      elseif(l.eq.5) then

        g     =  sqrt(1120.d0)
        ss(1) =  sqrt((70.d0 + g)/126.d0)
        ss(2) =  sqrt((70.d0 - g)/126.d0)
        ss(3) =  0.0d0
        ss(4) = -ss(2)
        ss(5) = -ss(1)

        ww(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
        ww(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
        ww(3) =  2.d0*(1.d0 - ww(1) - ww(2))
        ww(4) =  ww(2)
        ww(5) =  ww(1)

        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do ! k
        end do ! j

c     Error

      else
        write(dfich2,2000) l
        stop    
      endif

c     Format

2000  format(' *ERROR* INT2D: Illegal quadrature order =',i16)
      
      return
      end


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

      subroutine tint2d(l,gint,el)


c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set gauss points and weights for triangular elements

c      Inputs:
c         l       - Number of gauss points indicator

c      Outputs:
c         gint    - Total number of points
c         el(4,*) - Area coordinate points and weights for quadrature
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2


      integer   l, gint,i,j
      real*8    el(4,*), r0,r1,r2, ww, eta
      real*8    Lpt(73,3),W(73)



      save

      data ww, eta / 0.3333333333333333d0 , 0.1666666666666667d0 /

c     1-point gauss integration

      if(l.eq.1) then
        el(1,1) = ww
        el(2,1) = ww
        el(3,1) = ww
        el(4,1) = 1.d0
        gint    = 1

c     3-point integration: mid-edge points

      elseif(l.eq.3) then
        el(1,1) = 0.d0
        el(2,1) = 0.5d0
        el(3,1) = 0.5d0
        el(4,1) = ww

        el(1,2) = 0.5d0
        el(2,2) = 0.d0
        el(3,2) = 0.5d0
        el(4,2) = ww

        el(1,3) = 0.5d0
        el(2,3) = 0.5d0
        el(3,3) = 0.d0
        el(4,3) = ww

        gint    = 3

c     3-point integration: interior points

      elseif(l.eq.-3) then

        el(1,1) = 1.0d0 - ww
        el(2,1) = eta
        el(3,1) = eta
        el(4,1) = ww

        el(1,2) = eta
        el(2,2) = 1.0d0 - ww
        el(3,2) = eta
        el(4,2) = ww

        el(1,3) = eta
        el(2,3) = eta
        el(3,3) = 1.0d0 - ww
        el(4,3) = ww

        gint    = 3

c     4-point gauss integration: NOT RECOMMENDED DUE TO NEGATIVE WEIGHT

      elseif(l.eq.4) then
        el(1,1) =  ww
        el(2,1) =  ww
        el(3,1) =  ww
        el(4,1) = -27.d0/48.d0

        el(1,2) =  0.6d0
        el(2,2) =  0.2d0
        el(3,2) =  0.2d0
        el(4,2) =  25.d0/48.d0

        el(1,3) =  0.2d0
        el(2,3) =  0.6d0
        el(3,3) =  0.2d0
        el(4,3) =  el(4,2)

        el(1,4) =  0.2d0
        el(2,4) =  0.2d0
        el(3,4) =  0.6d0
        el(4,4) =  el(4,2)

        gint    =  4

c     6-point nodal integration

      elseif(l.eq.6) then

        el(1,1) =  1.0d0
        el(2,1) =  0.0d0
        el(3,1) =  0.0d0
        el(4,1) =  eta

        el(1,2) =  0.0d0
        el(2,2) =  1.0d0
        el(3,2) =  0.0d0
        el(4,2) =  eta

        el(1,3) =  0.0d0
        el(2,3) =  0.0d0
        el(3,3) =  1.0d0
        el(4,3) =  eta

        el(1,4) =  0.5d0
        el(2,4) =  0.5d0
        el(3,4) =  0.0d0
        el(4,4) =  eta

        el(1,5) =  0.0d0
        el(2,5) =  0.5d0
        el(3,5) =  0.5d0
        el(4,5) =  eta

        el(1,6) =  0.5d0
        el(2,6) =  0.0d0
        el(3,6) =  0.5d0
        el(4,6) =  eta

        gint    =  6

c     6-point order 4 formula

      elseif(l.eq.-6) then

        el(1,1) = 0.816847572980459d0
        el(2,1) = 0.091576213509771d0
        el(3,1) = 0.091576213509771d0
        el(4,1) = 0.109951743655322d0

        el(1,2) = 0.091576213509771d0
        el(2,2) = 0.816847572980459d0
        el(3,2) = 0.091576213509771d0
        el(4,2) = 0.109951743655322d0

        el(2,3) = 0.091576213509771d0
        el(1,3) = 0.091576213509771d0
        el(3,3) = 0.816847572980459d0
        el(4,3) = 0.109951743655322d0

        el(1,4) = 0.108103018168070d0
        el(2,4) = 0.445948490915965d0
        el(3,4) = 0.445948490915965d0
        el(4,4) = 0.223381589678011d0

        el(1,5) = 0.445948490915965d0
        el(2,5) = 0.108103018168070d0
        el(3,5) = 0.445948490915965d0
        el(4,5) = 0.223381589678011d0

        el(1,6) = 0.445948490915965d0
        el(2,6) = 0.445948490915965d0
        el(3,6) = 0.108103018168070d0
        el(4,6) = 0.223381589678011d0

        gint    = 6

c     7-point gauss integration

      elseif(l.eq.7) then
        r0      =  sqrt(15.0d0)
        r1      =  3.d0/7.d0
        r2      =  (r0 + r0)/21.d0

        el(1,1) =  ww
        el(2,1) =  el(1,1)
        el(3,1) =  el(1,1)
        el(4,1) =  0.225d0

        el(1,2) =  r1 + r2
        el(2,2) =  0.5d0 - 0.5d0*el(1,2)
        el(3,2) =  el(2,2)
        el(4,2) =  (155.d0 - r0)/1200.d0

        el(1,3) =  el(2,2)
        el(2,3) =  el(1,2)
        el(3,3) =  el(2,2)
        el(4,3) =  el(4,2)

        el(1,4) =  el(2,2)
        el(2,4) =  el(2,2)
        el(3,4) =  el(1,2)
        el(4,4) =  el(4,2)

        el(1,5) =  r1 - r2
        el(2,5) =  0.5d0 - 0.5d0*el(1,5)
        el(3,5) =  el(2,5)
        el(4,5) =  (155.d0 + r0)/1200.d0

        el(1,6) =  el(2,5)
        el(2,6) =  el(1,5)
        el(3,6) =  el(2,5)
        el(4,6) =  el(4,5)

        el(1,7) =  el(2,5)
        el(2,7) =  el(2,5)
        el(3,7) =  el(1,5)
        el(4,7) =  el(4,5)

        gint    =  7

c     12-point order 6 formula

      elseif(l.eq.12) then

        el(1, 1) = 0.873821971016996d0
        el(2, 1) = 0.063089014491502d0
        el(3, 1) = 0.063089014491502d0
        el(4, 1) = 0.050844906370207d0

        el(1, 2) = 0.063089014491502d0
        el(2, 2) = 0.873821971016996d0
        el(3, 2) = 0.063089014491502d0
        el(4, 2) = 0.050844906370207d0

        el(1, 3) = 0.063089014491502d0
        el(2, 3) = 0.063089014491502d0
        el(3, 3) = 0.873821971016996d0
        el(4, 3) = 0.050844906370207d0

        el(1, 4) = 0.501426509658179d0
        el(2, 4) = 0.249286745170910d0
        el(3, 4) = 0.249286745170910d0
        el(4, 4) = 0.116786275726379d0

        el(1, 5) = 0.249286745170910d0
        el(2, 5) = 0.501426509658179d0
        el(3, 5) = 0.249286745170910d0
        el(4, 5) = 0.116786275726379d0

        el(1, 6) = 0.249286745170910d0
        el(2, 6) = 0.249286745170910d0
        el(3, 6) = 0.501426509658179d0
        el(4, 6) = 0.116786275726379d0

        el(1, 7) = 0.636502499121399d0
        el(2, 7) = 0.310352451033785d0
        el(3, 7) = 0.053145049844816d0
        el(4, 7) = 0.082851075618374d0

        el(1, 8) = 0.636502499121399d0
        el(2, 8) = 0.053145049844816d0
        el(3, 8) = 0.310352451033785d0
        el(4, 8) = 0.082851075618374d0

        el(1, 9) = 0.310352451033785d0
        el(2, 9) = 0.636502499121399d0
        el(3, 9) = 0.053145049844816d0
        el(4, 9) = 0.082851075618374d0

        el(1,10) = 0.053145049844816d0
        el(2,10) = 0.636502499121399d0
        el(3,10) = 0.310352451033785d0
        el(4,10) = 0.082851075618374d0

        el(1,11) = 0.310352451033785d0
        el(2,11) = 0.053145049844816d0
        el(3,11) = 0.636502499121399d0
        el(4,11) = 0.082851075618374d0

        el(1,12) = 0.053145049844816d0
        el(2,12) = 0.310352451033785d0
        el(3,12) = 0.636502499121399d0
        el(4,12) = 0.082851075618374d0

        gint     = 12

c     13-point order 7 formula: NOT RECOMMENDED DUE TO NEGATIVE WEIGHT

      elseif(l.eq.13) then

        el(1, 1) =  ww
        el(2, 1) =  el(1,1)
        el(3, 1) =  el(1,1)
        el(4, 1) = -0.149570044467670d0

        el(1, 2) =  0.479308067841923d0
        el(2, 2) =  0.260345966079038d0
        el(3, 2) =  0.260345966079038d0
        el(4, 2) =  0.175615257433204d0

        el(1, 3) =  0.260345966079038d0
        el(2, 3) =  0.479308067841923d0
        el(3, 3) =  0.260345966079038d0
        el(4, 3) =  0.175615257433204d0

        el(1, 4) =  0.260345966079038d0
        el(2, 4) =  0.260345966079038d0
        el(3, 4) =  0.479308067841923d0
        el(4, 4) =  0.175615257433204d0

        el(1, 5) =  0.869739794195568d0
        el(2, 5) =  0.065130102902216d0
        el(3, 5) =  0.065130102902216d0
        el(4, 5) =  0.053347235608839d0

        el(1, 6) =  0.065130102902216d0
        el(2, 6) =  0.869739794195568d0
        el(3, 6) =  0.065130102902216d0
        el(4, 6) =  0.053347235608839d0

        el(1, 7) =  0.065130102902216d0
        el(2, 7) =  0.065130102902216d0
        el(3, 7) =  0.869739794195568d0
        el(4, 7) =  0.053347235608839d0

        el(1, 8) =  0.638444188569809d0
        el(2, 8) =  0.312865496004875d0
        el(3, 8) =  0.048690315425316d0
        el(4, 8) =  0.077113760890257d0

        el(1, 9) =  0.638444188569809d0
        el(2, 9) =  0.048690315425316d0
        el(3, 9) =  0.312865496004875d0
        el(4, 9) =  0.077113760890257d0

        el(1,10) =  0.312865496004875d0
        el(2,10) =  0.638444188569809d0
        el(3,10) =  0.048690315425316d0
        el(4,10) =  0.077113760890257d0

        el(1,11) =  0.048690315425316d0
        el(2,11) =  0.638444188569809d0
        el(3,11) =  0.312865496004875d0
        el(4,11) =  0.077113760890257d0

        el(1,12) =  0.312865496004875d0
        el(2,12) =  0.048690315425316d0
        el(3,12) =  0.638444188569809d0
        el(4,12) =  0.077113760890257d0

        el(1,13) =  0.048690315425316d0
        el(2,13) =  0.312865496004875d0
        el(3,13) =  0.638444188569809d0
        el(4,13) =  0.077113760890257d0

        gint     = 13


c     73-point order 19 formula

      elseif(l.eq.73) then

        Lpt(1,1)=3.33333333333333D-001
        Lpt(1,2)=3.33333333333333D-001
        Lpt(1,3)=3.33333333333333D-001
        W(1)=3.29063313889190D-002
        Lpt(2,1)=2.07800258539870D-002
        Lpt(2,2)=4.89609987073006D-001
        Lpt(2,3)=4.89609987073006D-001
        W(2)=1.03307318912720D-002
        Lpt(3,1)=4.89609987073006D-001
        Lpt(3,2)=4.89609987073006D-001
        Lpt(3,3)=2.07800258539870D-002
        W(3)=1.03307318912720D-002
        Lpt(4,1)=4.89609987073006D-001
        Lpt(4,2)=2.07800258539870D-002
        Lpt(4,3)=4.89609987073006D-001
        W(4)=1.03307318912720D-002
        Lpt(5,1)=9.09262146042150D-002
        Lpt(5,2)=4.54536892697893D-001
        Lpt(5,3)=4.54536892697893D-001
        W(5)=2.23872472630160D-002
        Lpt(6,1)=4.54536892697893D-001
        Lpt(6,2)=4.54536892697893D-001
        Lpt(6,3)=9.09262146042150D-002
        W(6)=2.23872472630160D-002
        Lpt(7,1)=4.54536892697893D-001
        Lpt(7,2)=9.09262146042150D-002
        Lpt(7,3)=4.54536892697893D-001
        W(7)=2.23872472630160D-002
        Lpt(8,1)=1.97166638701138D-001
        Lpt(8,2)=4.01416680649431D-001
        Lpt(8,3)=4.01416680649431D-001
        W(8)=3.02661258694680D-002
        Lpt(9,1)=4.01416680649431D-001
        Lpt(9,2)=4.01416680649431D-001
        Lpt(9,3)=1.97166638701138D-001
        W(9)=3.02661258694680D-002
        Lpt(10,1)=4.01416680649431D-001
        Lpt(10,2)=1.97166638701138D-001
        Lpt(10,3)=4.01416680649431D-001
        W(10)=3.02661258694680D-002
        Lpt(11,1)=4.88896691193805D-001
        Lpt(11,2)=2.55551654403098D-001
        Lpt(11,3)=2.55551654403098D-001
        W(11)=3.04909678021980D-002
        Lpt(12,1)=2.55551654403098D-001
        Lpt(12,2)=2.55551654403098D-001
        Lpt(12,3)=4.88896691193805D-001
        W(12)=3.04909678021980D-002
        Lpt(13,1)=2.55551654403098D-001
        Lpt(13,2)=4.88896691193805D-001
        Lpt(13,3)=2.55551654403098D-001
        W(13)=3.04909678021980D-002
        Lpt(14,1)=6.45844115695741D-001
        Lpt(14,2)=1.77077942152130D-001
        Lpt(14,3)=1.77077942152130D-001
        W(14)=2.41592127416410D-002
        Lpt(15,1)=1.77077942152130D-001
        Lpt(15,2)=1.77077942152130D-001
        Lpt(15,3)=6.45844115695741D-001
        W(15)=2.41592127416410D-002
        Lpt(16,1)=1.77077942152130D-001
        Lpt(16,2)=6.45844115695741D-001
        Lpt(16,3)=1.77077942152130D-001
        W(16)=2.41592127416410D-002
        Lpt(17,1)=7.79877893544096D-001
        Lpt(17,2)=1.10061053227952D-001
        Lpt(17,3)=1.10061053227952D-001
        W(17)=1.60508035868010D-002
        Lpt(18,1)=1.10061053227952D-001
        Lpt(18,2)=1.10061053227952D-001
        Lpt(18,3)=7.79877893544096D-001
        W(18)=1.60508035868010D-002
        Lpt(19,1)=1.10061053227952D-001
        Lpt(19,2)=7.79877893544096D-001
        Lpt(19,3)=1.10061053227952D-001
        W(19)=1.60508035868010D-002
        Lpt(20,1)=8.88942751496321D-001
        Lpt(20,2)=5.55286242518400D-002
        Lpt(20,3)=5.55286242518400D-002
        W(20)=8.08458026178400D-003
        Lpt(21,1)=5.55286242518400D-002
        Lpt(21,2)=5.55286242518400D-002
        Lpt(21,3)=8.88942751496321D-001
        W(21)=8.08458026178400D-003
        Lpt(22,1)=5.55286242518400D-002
        Lpt(22,2)=8.88942751496321D-001
        Lpt(22,3)=5.55286242518400D-002
        W(22)=8.08458026178400D-003
        Lpt(23,1)=9.74756272445543D-001
        Lpt(23,2)=1.26218637772290D-002
        Lpt(23,3)=1.26218637772290D-002
        W(23)=2.07936202748500D-003
        Lpt(24,1)=1.26218637772290D-002
        Lpt(24,2)=1.26218637772290D-002
        Lpt(24,3)=9.74756272445543D-001
        W(24)=2.07936202748500D-003
        Lpt(25,1)=1.26218637772290D-002
        Lpt(25,2)=9.74756272445543D-001
        Lpt(25,3)=1.26218637772290D-002
        W(25)=2.07936202748500D-003
        Lpt(26,1)=3.61141784841200D-003
        Lpt(26,2)=3.95754787356943D-001
        Lpt(26,3)=6.00633794794645D-001
        W(26)=3.88487690498100D-003
        Lpt(27,1)=3.61141784841200D-003
        Lpt(27,2)=6.00633794794645D-001
        Lpt(27,3)=3.95754787356943D-001
        W(27)=3.88487690498100D-003
        Lpt(28,1)=6.00633794794645D-001
        Lpt(28,2)=3.61141784841200D-003
        Lpt(28,3)=3.95754787356943D-001
        W(28)=3.88487690498100D-003
        Lpt(29,1)=3.95754787356943D-001
        Lpt(29,2)=3.61141784841200D-003
        Lpt(29,3)=6.00633794794645D-001
        W(29)=3.88487690498100D-003
        Lpt(30,1)=3.95754787356943D-001
        Lpt(30,2)=6.00633794794645D-001
        Lpt(30,3)=3.61141784841200D-003
        W(30)=3.88487690498100D-003
        Lpt(31,1)=6.00633794794645D-001
        Lpt(31,2)=3.95754787356943D-001
        Lpt(31,3)=3.61141784841200D-003
        W(31)=3.88487690498100D-003
        Lpt(32,1)=1.34466754530780D-001
        Lpt(32,2)=3.07929983880436D-001
        Lpt(32,3)=5.57603261588784D-001
        W(32)=2.55741606120220D-002
        Lpt(33,1)=1.34466754530780D-001
        Lpt(33,2)=5.57603261588784D-001
        Lpt(33,3)=3.07929983880436D-001
        W(33)=2.55741606120220D-002
        Lpt(34,1)=5.57603261588784D-001
        Lpt(34,2)=1.34466754530780D-001
        Lpt(34,3)=3.07929983880436D-001
        W(34)=2.55741606120220D-002
        Lpt(35,1)=3.07929983880436D-001
        Lpt(35,2)=1.34466754530780D-001
        Lpt(35,3)=5.57603261588784D-001
        W(35)=2.55741606120220D-002
        Lpt(36,1)=3.07929983880436D-001
        Lpt(36,2)=5.57603261588784D-001
        Lpt(36,3)=1.34466754530780D-001
        W(36)=2.55741606120220D-002
        Lpt(37,1)=5.57603261588784D-001
        Lpt(37,2)=3.07929983880436D-001
        Lpt(37,3)=1.34466754530780D-001
        W(37)=2.55741606120220D-002
        Lpt(38,1)=1.44460257761150D-002
        Lpt(38,2)=2.64566948406520D-001
        Lpt(38,3)=7.20987025817365D-001
        W(38)=8.88090357333800D-003
        Lpt(39,1)=1.44460257761150D-002
        Lpt(39,2)=7.20987025817365D-001
        Lpt(39,3)=2.64566948406520D-001
        W(39)=8.88090357333800D-003
        Lpt(40,1)=7.20987025817365D-001
        Lpt(40,2)=1.44460257761150D-002
        Lpt(40,3)=2.64566948406520D-001
        W(40)=8.88090357333800D-003
        Lpt(41,1)=2.64566948406520D-001
        Lpt(41,2)=1.44460257761150D-002
        Lpt(41,3)=7.20987025817365D-001
        W(41)=8.88090357333800D-003
        Lpt(42,1)=2.64566948406520D-001
        Lpt(42,2)=7.20987025817365D-001
        Lpt(42,3)=1.44460257761150D-002
        W(42)=8.88090357333800D-003
        Lpt(43,1)=7.20987025817365D-001
        Lpt(43,2)=2.64566948406520D-001
        Lpt(43,3)=1.44460257761150D-002
        W(43)=8.88090357333800D-003
        Lpt(44,1)=4.69335788381780D-002
        Lpt(44,2)=3.58539352205951D-001
        Lpt(44,3)=5.94527068955871D-001
        W(44)=1.61245467617310D-002
        Lpt(45,1)=4.69335788381780D-002
        Lpt(45,2)=5.94527068955871D-001
        Lpt(45,3)=3.58539352205951D-001
        W(45)=1.61245467617310D-002
        Lpt(46,1)=5.94527068955871D-001
        Lpt(46,2)=4.69335788381780D-002
        Lpt(46,3)=3.58539352205951D-001
        W(46)=1.61245467617310D-002
        Lpt(47,1)=3.58539352205951D-001
        Lpt(47,2)=4.69335788381780D-002
        Lpt(47,3)=5.94527068955871D-001
        W(47)=1.61245467617310D-002
        Lpt(48,1)=3.58539352205951D-001
        Lpt(48,2)=5.94527068955871D-001
        Lpt(48,3)=4.69335788381780D-002
        W(48)=1.61245467617310D-002
        Lpt(49,1)=5.94527068955871D-001
        Lpt(49,2)=3.58539352205951D-001
        Lpt(49,3)=4.69335788381780D-002
        W(49)=1.61245467617310D-002
        Lpt(50,1)=2.86112035056700D-003
        Lpt(50,2)=1.57807405968595D-001
        Lpt(50,3)=8.39331473680839D-001
        W(50)=2.49194181749100D-003
        Lpt(51,1)=2.86112035056700D-003
        Lpt(51,2)=8.39331473680839D-001
        Lpt(51,3)=1.57807405968595D-001
        W(51)=2.49194181749100D-003
        Lpt(52,1)=8.39331473680839D-001
        Lpt(52,2)=2.86112035056700D-003
        Lpt(52,3)=1.57807405968595D-001
        W(52)=2.49194181749100D-003
        Lpt(53,1)=1.57807405968595D-001
        Lpt(53,2)=2.86112035056700D-003
        Lpt(53,3)=8.39331473680839D-001
        W(53)=2.49194181749100D-003
        Lpt(54,1)=1.57807405968595D-001
        Lpt(54,2)=8.39331473680839D-001
        Lpt(54,3)=2.86112035056700D-003
        W(54)=2.49194181749100D-003
        Lpt(55,1)=8.39331473680839D-001
        Lpt(55,2)=1.57807405968595D-001
        Lpt(55,3)=2.86112035056700D-003
        W(55)=2.49194181749100D-003
        Lpt(56,1)=2.23861424097916D-001
        Lpt(56,2)=7.50505969759110D-002
        Lpt(56,3)=7.01087978926173D-001
        W(56)=1.82428401189510D-002
        Lpt(57,1)=2.23861424097916D-001
        Lpt(57,2)=7.01087978926173D-001
        Lpt(57,3)=7.50505969759110D-002
        W(57)=1.82428401189510D-002
        Lpt(58,1)=7.01087978926173D-001
        Lpt(58,2)=2.23861424097916D-001
        Lpt(58,3)=7.50505969759110D-002
        W(58)=1.82428401189510D-002
        Lpt(59,1)=7.50505969759110D-002
        Lpt(59,2)=2.23861424097916D-001
        Lpt(59,3)=7.01087978926173D-001
        W(59)=1.82428401189510D-002
        Lpt(60,1)=7.50505969759110D-002
        Lpt(60,2)=7.01087978926173D-001
        Lpt(60,3)=2.23861424097916D-001
        W(60)=1.82428401189510D-002
        Lpt(61,1)=7.01087978926173D-001
        Lpt(61,2)=7.50505969759110D-002
        Lpt(61,3)=2.23861424097916D-001
        W(61)=1.82428401189510D-002
        Lpt(62,1)=3.46470748167600D-002
        Lpt(62,2)=1.42421601113383D-001
        Lpt(62,3)=8.22931324069857D-001
        W(62)=1.02585637361990D-002
        Lpt(63,1)=3.46470748167600D-002
        Lpt(63,2)=8.22931324069857D-001
        Lpt(63,3)=1.42421601113383D-001
        W(63)=1.02585637361990D-002
        Lpt(64,1)=8.22931324069857D-001
        Lpt(64,2)=3.46470748167600D-002
        Lpt(64,3)=1.42421601113383D-001
        W(64)=1.02585637361990D-002
        Lpt(65,1)=1.42421601113383D-001
        Lpt(  65,2)=3.46470748167600D-002
        Lpt(65,3)=8.22931324069857D-001
        W(65)=1.02585637361990D-002
        Lpt(66,1)=1.42421601113383D-001
        Lpt(66,2)=8.22931324069857D-001
        Lpt(66,3)=3.46470748167600D-002
        W(66)=1.02585637361990D-002
        Lpt(67,1)=8.22931324069857D-001
        Lpt(67,2)=1.42421601113383D-001
        Lpt(67,3)=3.46470748167600D-002
        W(67)=1.02585637361990D-002
        Lpt(68,1)=1.01611192962780D-002
        Lpt(68,2)=6.54946280829380D-002
        Lpt(68,3)=9.24344252620784D-001
        W(68)=3.79992885530200D-003
        Lpt(69,1)=1.01611192962780D-002
        Lpt(69,2)=9.24344252620784D-001
        Lpt(69,3)=6.54946280829380D-002
        W(69)=3.79992885530200D-003
        Lpt(70,1)=9.24344252620784D-001
        Lpt(70,2)=1.01611192962780D-002
        Lpt(70,3)=6.54946280829380D-002
        W(70)=3.79992885530200D-003
        Lpt(71,1)=6.54946280829380D-002
        Lpt(71,2)=1.01611192962780D-002
        Lpt(71,3)=9.24344252620784D-001
        W(71)=3.79992885530200D-003
        Lpt(72,1)=6.54946280829380D-002
        Lpt(72,2)=9.24344252620784D-001
        Lpt(72,3)=1.01611192962780D-002
        W(72)=3.79992885530200D-003
        Lpt(73,1)=9.24344252620784D-001
        Lpt(73,2)=6.54946280829380D-002
        Lpt(73,3)=1.01611192962780D-002
        W(73)=3.79992885530200D-003

      do i=1,73
        do j=1,3
          el(j,i)=Lpt(i,j)
        enddo
        el(4,i)=W(i)
      enddo

      gint=73


c     Unspecified quadrature specified

      else
        write(dfich2,2000) l
        gint    = -1
        stop
      endif

c     Format

2000  format(' *ERROR* TINT2D: Wrong quadrature, l =',i3)

      return  
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes


      subroutine subelem(ElemGGe,Xe,Ye,XG,YG,NumSub,SubXe,SubYe,n,nel)


c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Descomposición del dominio del elemento en subelementos
c			 para la integracion
c     --------------------------------------

c      Inputs:
c			ElemGGe		Matriz de caracterización del elemento
c			Xe(8)		Coordenadas X de los nodos del elemento
c			Ye(8)		Coordenadas Y de los nodos del elemento
c			XG			Coordenada X del extremo de Grieta
c			YG			Coordenada Y del extremo de Grieta		

c      Outputs:
c			NumSub		Número de subelementos generados
c			SubXe(10,4)	Coordenadas X de los puntos de los subelementos			
c			SubYe(10,4)	Coordenadas X de los puntos de los subelementos
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none
      
c     include  'eldata.h'
c  	include  'iofile.h'

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2


      integer NumSub,n1,n2,l1,l2,lado,lado2,nl,i,tipo,lad,lad2,lads
      integer n,nel
      real*8  ElemGGe(8),XG,YG,xm,ym,x1,y1,x2,y2,xl,yl,xl2,yl2,xls,yls,
     &		Xe(8),Ye(8),SubXe(10,4),SubYe(10,4),tol
      logical comp

      
      tipo=ElemGGe(2)
      tol=1d-6


c      write(dfich2,*) 'Nº elem=',n,'   ElemGG(i,1)=',ElemGGe(1)
c      write(dfich2,*) 'ElemGG(i,2)=',ElemGGe(2),'   tipo=',tipo

c     write(dfich2,*)
c	write(dfich2,*)
c	write(dfich2,*) 'Xe(i)=', (Xe(i),i=1,8)
c	write(dfich2,*)
c	write(dfich2,*) 'Ye(i)=', (Ye(i),i=1,8)
c	write(dfich2,*)


      
c	Interseccion en 2 nodos opuestos
c
c		2									
c		3
c		4
c		5
c		11
c		12
c		20
c		22
c		30
c		31
c		32
c		40
c		41
c		50
c		51
c		80 


c	!!!Convierto a enteros con la asignacion en variable entera:preguntar 



      select case(tipo)

      case(2)
        NumSub= 2;
        SubXe= 0
        SubYe= 0
        n1= min(ElemGGe(3),ElemGGe(4));
        SubXe(1,1)= Xe(n1);      SubYe(1,1)= Ye(n1);
        SubXe(1,2)= Xe(n1+1);    SubYe(1,2)= Ye(n1+1);
        SubXe(1,3)= Xe(n1+2);    SubYe(1,3)= Ye(n1+2);
        SubXe(2,1)= Xe(n1);      SubYe(2,1)= Ye(n1);
        SubXe(2,2)= Xe(n1+2);    SubYe(2,2)= Ye(n1+2);
        SubXe(2,3)= Xe(n1+3);    SubYe(2,3)= Ye(n1+3);

      

c	 === Interseccion en nodo y lado
      case(3)
        NumSub= 5;
        SubXe=0
        SubYe=0
        n1= ElemGGe(3);
        l1= ElemGGe(4);
        x1= ElemGGe(5);
        y1= ElemGGe(6);

        if ((abs(n1-l1).eq.1).or.((n1.eq.4).and.(l1.eq.1))) then

c	    % - triangulo + cuadrilatero
    
              SubXe(1,1)= Xe(n1);      SubYe(1,1)= Ye(n1);
          SubXe(1,2)= Xe(n1+1);    SubYe(1,2)= Ye(n1+1);
          SubXe(1,3)= x1;          SubYe(1,3)= y1;
          
            xm= (Xe(n1)+x1+Xe(n1+2)+Xe(n1+3))/4;
          ym= (Ye(n1)+y1+Ye(n1+2)+Ye(n1+3))/4;
    
          SubXe(2,1)= Xe(n1);      SubYe(2,1)= Ye(n1);
          SubXe(2,2)= x1;          SubYe(2,2)= y1;
          SubXe(2,3)= xm;          SubYe(2,3)= ym;
    
          SubXe(3,1)= x1;			 SubYe(3,1)= y1;
          SubXe(3,2)= Xe(n1+2);    SubYe(3,2)= Ye(n1+2);
          SubXe(3,3)= xm;			 SubYe(3,3)= ym;
    
          SubXe(4,1)= Xe(n1+2);    SubYe(4,1)= Ye(n1+2);
          SubXe(4,2)= Xe(n1+3);    SubYe(4,2)= Ye(n1+3);
          SubXe(4,3)= xm;			 SubYe(4,3)= ym;
    
          SubXe(5,1)= Xe(n1+3);    SubYe(5,1)= Ye(n1+3);
          SubXe(5,2)= Xe(n1+4);    SubYe(5,2)= Ye(n1+4);
          SubXe(5,3)= xm;			 SubYe(5,3)= ym;
    
        else
c	    % - cuadrilatero + triangulo
    
          SubXe(1,1)= Xe(n1);      SubYe(1,1)= Ye(n1);
          SubXe(1,2)= x1;          SubYe(1,2)= y1;
          SubXe(1,3)= Xe(n1+3);    SubYe(1,3)= Ye(n1+3);

          xm= (Xe(n1)+Xe(n1+1)+Xe(n1+2)+x1)/4;
          ym= (Ye(n1)+Ye(n1+1)+Ye(n1+2)+y1)/4;
    
          SubXe(2,1)= Xe(n1);      SubYe(2,1)= Ye(n1);
          SubXe(2,2)= Xe(n1+1);    SubYe(2,2)= Ye(n1+1);
          SubXe(2,3)= xm;          SubYe(2,3)= ym;
    
          SubXe(3,1)= Xe(n1+1);    SubYe(3,1)= Ye(n1+1);
          SubXe(3,2)= Xe(n1+2);    SubYe(3,2)= Ye(n1+2);
          SubXe(3,3)= xm;          SubYe(3,3)= ym;
    
          SubXe(4,1)= Xe(n1+2);    SubYe(4,1)= Ye(n1+2);
          SubXe(4,2)= x1;          SubYe(4,2)= y1;
          SubXe(4,3)= xm;          SubYe(4,3)= ym;
    
          SubXe(5,1)= x1;          SubYe(5,1)= y1;
          SubXe(5,2)= Xe(n1);      SubYe(5,2)= Ye(n1);
          SubXe(5,3)= xm;          SubYe(5,3)= ym;
        end if
  


c		% === Interseccion en 2 lados opuestos
      case(4)

          NumSub=2;
            SubXe(1,:)=
     &	(/ElemGGe(5),Xe(ElemGGe(3)+1),Xe(ElemGGe(3)+2),ElemGGe(7)/);
          SubXe(2,:)=
     &	(/ElemGGe(7),Xe(ElemGGe(3)+3),Xe(ElemGGe(3)+4),ElemGGe(5)/);
          SubYe(1,:)=
     &	(/ElemGGe(6),Ye(ElemGGe(3)+1),Ye(ElemGGe(3)+2),ElemGGe(8)/);
          SubYe(2,:)=
     &	(/ElemGGe(8),Ye(ElemGGe(3)+3),Ye(ElemGGe(3)+4),ElemGGe(6)/);



c	 === Interseccion en 2 lados contiguos
      case(5)
        l1= ElemGGe(3);
        l2= ElemGGe(4);
        if ((l1.gt.l2).or.((l1.eq.1).and.(l2.eq.4)))then
         l2= ElemGGe(3);
         l1= ElemGGe(4);
         x2= ElemGGe(5);
         y2= ElemGGe(6);
         x1= ElemGGe(7);
         y1= ElemGGe(8);
        else
         l1= ElemGGe(3)
         l2= ElemGGe(4);
         x1= ElemGGe(5);
         y1= ElemGGe(6);
         x2= ElemGGe(7);
         y2= ElemGGe(8);
        end if
        NumSub= 6;  
        SubXe=0
        SubYe=0
  
        SubXe(1,1)= x1;          SubYe(1,1)= y1;
        SubXe(1,2)= Xe(l1+1);    SubYe(1,2)= Ye(l1+1);
        SubXe(1,3)= x2;          SubYe(1,3)= y2;
  
        xm= (x1+x2+Xe(l1+2)+Xe(l1+3)+Xe(l1+4))/5;
        ym= (y1+y2+Ye(l1+2)+Ye(l1+3)+Ye(l1+4))/5;
  
        SubXe(2,1)= x1;          SubYe(2,1)= y1;
        SubXe(2,2)= x2;          SubYe(2,2)= y2;
        SubXe(2,3)= xm;          SubYe(2,3)= ym;
  
        SubXe(3,1)= x2;          SubYe(3,1)= y2;
        SubXe(3,2)= Xe(l1+2);    SubYe(3,2)= Ye(l1+2);
        SubXe(3,3)= xm;          SubYe(3,3)= ym;
  
        SubXe(4,1)= Xe(l1+2);    SubYe(4,1)= Ye(l1+2);
        SubXe(4,2)= Xe(l1+3);    SubYe(4,2)= Ye(l1+3);
        SubXe(4,3)= xm;          SubYe(4,3)= ym;
  
        SubXe(5,1)= Xe(l1+3);    SubYe(5,1)= Ye(l1+3);
        SubXe(5,2)= Xe(l1+4);    SubYe(5,2)= Ye(l1+4);
        SubXe(5,3)= xm;          SubYe(5,3)= ym;
  
        SubXe(6,1)= Xe(l1+4);    SubYe(6,1)= Ye(l1+4);
        SubXe(6,2)= x1;          SubYe(6,2)= y1;
        SubXe(6,3)= xm;          SubYe(6,3)= ym;



c	% == Tipo 11
      case(11)
          NumSub=3;
          nl=ElemGGe(3);
          SubXe= 0
          SubYe= 0
          do i= 1,3
            SubXe(i,1)=XG;           SubYe(i,1)= YG;
            SubXe(i,2)=Xe(nl+i);     SubYe(i,2)= Ye(nl+i);
            SubXe(i,3)=Xe(nl+i+1);   SubYe(i,3)= Ye(nl+i+1);
          end do 
          

c	% == Tipo 12
      case(12)
          NumSub=2;
            n1=ElemGGe(3);

            if ((abs(XG-Xe(n1)).le.tol).and.(abs(YG-Ye(n1)).le.tol)) then
              n1=ElemGGe(3);
          elseif ((abs(XG-Xe(n1+1)).le.tol).and.
     &     		(abs(YG-Ye(n1+1)).le.tol)) then
              n1=n1+1
            else
            write(dfich2,*) 'Paso por aquí'
            write(dfich2,*) '*ERROR* El elemento ',n,' de tipo ',tipo,
     &                   ' no ha podido subdividirse en subelem.f'
            stop

          end if   

          SubXe= 0;
          SubYe= 0;
          do i= 1,2
            SubXe(i,1)=XG;           SubYe(i,1)= YG;
            SubXe(i,2)=Xe(n1+i);     SubYe(i,2)= Ye(n1+i);
            SubXe(i,3)=Xe(n1+i+1);   SubYe(i,3)= Ye(n1+i+1);
          end do
       

c	% == Tipo 20
      case(20)
          NumSub=4;
          n1=ElemGGe(3);
          SubXe= 0
          SubYe= 0
          do i= 1,4 
            SubXe(i,1)=XG;           SubYe(i,1)= YG;
            SubXe(i,2)=Xe(n1+i-1);   SubYe(i,2)= Ye(n1+i-1);
            SubXe(i,3)=Xe(n1+i);     SubYe(i,3)= Ye(n1+i);
          end  do
      

c	% == Tipo 22
      case(22)
          NumSub=2;
          n1=ElemGGe(3);
          SubXe= 0
          SubYe= 0
          do i= 1,2
            SubXe(i,1)=Xe(n1+2*(i-1));       SubYe(i,1)= Ye(n1+2*(i-1));
            SubXe(i,2)=Xe(n1+1+2*(i-1));     SubYe(i,2)= Ye(n1+1+2*(i-1));
            SubXe(i,3)=Xe(n1+2+2*(i-1));     SubYe(i,3)= Ye(n1+2+2*(i-1));
          end  do
      

c	% == Tipo 30
      case(30)
          NumSub=5;
          lado=ElemGGe(4);
          xl=ElemGGe(5);
          yl=ElemGGe(6);
          SubXe= 0
          SubYe= 0
    
          SubXe(1,1)=XG;              SubYe(1,1)= YG;
          SubXe(1,2)=xl;              SubYe(1,2)= yl;
          SubXe(1,3)=Xe(lado+1);      SubYe(1,3)= Ye(lado+1);
          do i= 2,4
            SubXe(i,1)=XG;              SubYe(i,1)= YG;
            SubXe(i,2)=Xe(lado+i-1);    SubYe(i,2)= Ye(lado+i-1);
            SubXe(i,3)=Xe(lado+i);      SubYe(i,3)= Ye(lado+i);
          end do
          SubXe(5,1)=XG;              SubYe(5,1)= YG;
          SubXe(5,2)=Xe(lado+4);      SubYe(5,2)= Ye(lado+4);
          SubXe(5,3)=xl;              SubYe(5,3)= yl;
      

c	% == Tipo 31
      case(31)
          NumSub=3;
          lado=ElemGGe(4);
c	    %xl=ElemGGe(5);   es el mismo XG   
c	    %yl=ElemGGe(6);
          SubXe= 0
          SubYe= 0
    
          do i= 1,3
            SubXe(i,1)=XG;              SubYe(i,1)= YG;
            SubXe(i,2)=Xe(lado+i);      SubYe(i,2)= Ye(lado+i);
            SubXe(i,3)=Xe(lado+i+1);    SubYe(i,3)= Ye(lado+i+1);
          end do
      

c	% == Tipo 32
      case(32)
          NumSub=3;
          lado=ElemGGe(4);
          xl=ElemGGe(5);
          yl=ElemGGe(6);
          SubXe= 0;
          SubYe= 0;
    
          SubXe(1,1)=XG;              SubYe(1,1)= YG;
          SubXe(1,2)=Xe(lado);        SubYe(1,2)= Ye(lado);
          SubXe(1,3)=xl;              SubYe(1,3)= yl;

          SubXe(2,1)=XG;              SubYe(2,1)= YG;
          SubXe(2,2)=xl;              SubYe(2,2)= yl;
          SubXe(2,3)=Xe(lado+1);      SubYe(2,3)= Ye(lado+1);
    
          if ((lado-n).eq.1) then
              SubXe(3,1)=XG;              SubYe(3,1)= YG;
              SubXe(3,2)=Xe(lado+1);      SubYe(3,2)= Ye(lado+1);
              SubXe(3,3)=Xe(lado+2);      SubYe(3,3)= Ye(lado+2);
          else 
              SubXe(3,1)=XG;              SubYe(3,1)= YG;
              SubXe(3,2)=Xe(lado+3);      SubYe(3,2)= Ye(lado+3);
              SubXe(3,3)=Xe(lado+4);      SubYe(3,3)= Ye(lado+4);
          end  if
      

c	% == Tipo 40
      case(40)
          NumSub=6;
          lad=ElemGGe(3);
          lad2=ElemGGe(4);
          xl=ElemGGe(5);
          yl=ElemGGe(6);    
          xl2=ElemGGe(7);
          yl2=ElemGGe(8);
          SubXe= 0;
          SubYe= 0;
    
          SubXe(1,1)=XG;              SubYe(1,1)= YG;
          SubXe(1,2)=Xe(lad);         SubYe(1,2)= Ye(lad);
          SubXe(1,3)=xl;              SubYe(1,3)= yl;

          SubXe(2,1)=XG;              SubYe(2,1)= YG;
          SubXe(2,2)=xl;              SubYe(2,2)= yl;
          SubXe(2,3)=Xe(lad+1);       SubYe(2,3)= Ye(lad+1);
    
          SubXe(3,1)=XG;              SubYe(3,1)= YG;
          SubXe(3,2)=Xe(lad+1);       SubYe(3,2)= Ye(lad+1);
          SubXe(3,3)=Xe(lad+2);       SubYe(3,3)= Ye(lad+2);

          SubXe(4,1)=XG;              SubYe(4,1)= YG;
          SubXe(4,2)=Xe(lad2);        SubYe(4,2)= Ye(lad2);
          SubXe(4,3)=xl2;             SubYe(4,3)= yl2;
    
          SubXe(5,1)=XG;              SubYe(5,1)= YG;
          SubXe(5,2)=xl2;             SubYe(5,2)= yl2;
          SubXe(5,3)=Xe(lad2+1);      SubYe(5,3)= Ye(lad2+1);

          SubXe(6,1)=XG;              SubYe(6,1)= YG;
          SubXe(6,2)=Xe(lad2+1);      SubYe(6,2)= Ye(lad2+1);
          SubXe(6,3)=Xe(lad2+2);      SubYe(6,3)= Ye(lad2+2);
      

c	% == Tipo 41
      case(41)
          NumSub=4;
          lad=ElemGGe(3);
          lad2=ElemGGe(4);
          xl=ElemGGe(5);
          yl=ElemGGe(6);    
          xl2=ElemGGe(7);
          yl2=ElemGGe(8);


c		comp=
    
c	    if ((XG.eq.xl).and.(YG.eq.yl)) then  !!!!!!!!!!reales
            if ((abs(XG-xl).lt.tol).and.(abs(YG-yl).lt.tol)) then
              lads=lad;
              xls=xl2;
              yls=yl2;
            elseif ((abs(XG-xl2).lt.tol).and.(abs(YG-yl2).lt.tol)) then
c	    elseif ((XG.eq.xl2).and.(YG .eq.yl2)) then !!!!!!!!!!reales
              lads=lad2;
              xls=xl;
              yls=yl;
          else
            write(dfich2,*) '*ERROR* El elemento ',n,' de tipo ',tipo,
     &                   ' no ha podido subdividirse en subelem.f'
            stop
          end if   
    
          SubXe= 0;
          SubYe= 0;
    
          SubXe(1,1)=XG;              SubYe(1,1)= YG;
          SubXe(1,2)=Xe(lads+1);      SubYe(1,2)= Ye(lads+1);
          SubXe(1,3)=Xe(lads+2);      SubYe(1,3)= Ye(lads+2);

          SubXe(2,1)=XG;              SubYe(2,1)= YG;
          SubXe(2,2)=Xe(lads+2);      SubYe(2,2)= Ye(lads+2);
          SubXe(2,3)=xls;             SubYe(2,3)= yls;
    
          SubXe(3,1)=XG;              SubYe(3,1)= YG;
          SubXe(3,2)=xls;             SubYe(3,2)= yls;
          SubXe(3,3)=Xe(lads+3);      SubYe(3,3)= Ye(lads+3);

          SubXe(4,1)=XG;              SubYe(4,1)= YG;
          SubXe(4,2)=Xe(lads+3);      SubYe(4,2)= Ye(lads+3);
          SubXe(4,3)=Xe(lads+4);      SubYe(4,3)= Ye(lads+4);
    
      

c	% === Tipo 50
      case(50)
          NumSub=6;
          lad=ElemGGe(3);
            lad2=ElemGGe(4)
          xl=ElemGGe(5);
          yl=ElemGGe(6);    
          xl2=ElemGGe(7);
          yl2=ElemGGe(8);
          if ((lad.eq.1) .and. (lad2.eq.4)) then
              lad=ElemGGe(4);
              xl=ElemGGe(7);
              yl=ElemGGe(8);    
              xl2=ElemGGe(5);
              yl2=ElemGGe(6);      
          end if

          SubXe= 0;
          SubYe= 0;
    
          SubXe(1,1)=XG;              SubYe(1,1)= YG;
          SubXe(1,2)=Xe(lad);         SubYe(1,2)= Ye(lad);
          SubXe(1,3)=xl;              SubYe(1,3)= yl;

          SubXe(2,1)=XG;              SubYe(2,1)= YG;
          SubXe(2,2)=xl;              SubYe(2,2)= yl;
          SubXe(2,3)=Xe(lad+1);       SubYe(2,3)= Ye(lad+1);
    
          SubXe(3,1)=XG;              SubYe(3,1)= YG;
          SubXe(3,2)=Xe(lad+1);       SubYe(3,2)= Ye(lad+1);
          SubXe(3,3)=xl2;             SubYe(3,3)= yl2;

          SubXe(4,1)=XG;              SubYe(4,1)= YG;
          SubXe(4,2)=xl2;             SubYe(4,2)= yl2;
          SubXe(4,3)=Xe(lad+2);       SubYe(4,3)= Ye(lad+2);
    
          SubXe(5,1)=XG;              SubYe(5,1)= YG;
          SubXe(5,2)=Xe(lad+2);       SubYe(5,2)= Ye(lad+2);
          SubXe(5,3)=Xe(lad+3);       SubYe(5,3)= Ye(lad+3);

          SubXe(6,1)=XG;              SubYe(6,1)= YG;
          SubXe(6,2)=Xe(lad+3);       SubYe(6,2)= Ye(lad+3);
          SubXe(6,3)=Xe(lad+4);       SubYe(6,3)= Ye(lad+4);
      

c	% === Tipo 51
      case(51)
          NumSub=4;
          lad=ElemGGe(3);
          lad2=ElemGGe(4);
          xl=ElemGGe(5);
          yl=ElemGGe(6);    
          xl2=ElemGGe(7);
          yl2=ElemGGe(8);
      
          SubXe= 0;
          SubYe= 0;

c		comp=((abs(XG-xl).lt.tol).and.(abs(YG-yl).lt.tol))
            
          if ((XG.eq.xl).and.(YG.eq.yl))  then !Lado con fin de grieta

c		if(comp) then

              lads=lad;
                                

              SubXe(1,1)=XG;              SubYe(1,1)= YG;
              SubXe(1,2)=Xe(lads+1);      SubYe(1,2)= Ye(lads+1);
              SubXe(1,3)=xl2;             SubYe(1,3)= yl2;

              SubXe(2,1)=XG;              SubYe(2,1)= YG;
              SubXe(2,2)=xl2;             SubYe(2,2)= yl2;
              SubXe(2,3)=Xe(lads+2);      SubYe(2,3)= Ye(lads+2);
    
              SubXe(3,1)=XG;              SubYe(3,1)= YG;
              SubXe(3,2)=Xe(lads+2);      SubYe(3,2)= Ye(lads+2);
              SubXe(3,3)=Xe(lads+3);      SubYe(3,3)= Ye(lads+3);

              SubXe(4,1)=XG;              SubYe(4,1)= YG;
              SubXe(4,2)=Xe(lads+3);      SubYe(4,2)= Ye(lads+3);
              SubXe(4,3)=Xe(lads+4);      SubYe(4,3)= Ye(lads+4);
    
           elseif ((XG.eq.xl2).and.(YG.eq.yl2))  then !!!!!reales
c	     elseif ((abs(XG-xl2).lt.tol).and.(abs(YG-yl2).lt.tol))then 
 
              lads=lad2;
              SubXe(1,1)=XG;              SubYe(1,1)= YG;
              SubXe(1,2)=Xe(lads+1);      SubYe(1,2)= Ye(lads+1);
              SubXe(1,3)=Xe(lads+2);      SubYe(1,3)= Ye(lads+2);

              SubXe(2,1)=XG;              SubYe(2,1)= YG;
              SubXe(2,2)=Xe(lads+2);      SubYe(2,2)= Ye(lads+2);
              SubXe(2,3)=Xe(lads+3);      SubYe(2,3)= Ye(lads+3);
    
              SubXe(3,1)=XG;              SubYe(3,1)= YG;
              SubXe(3,2)=Xe(lads+3);      SubYe(3,2)= Ye(lads+3);
              SubXe(3,3)=xl;              SubYe(3,3)= yl;

              SubXe(4,1)=XG;              SubYe(4,1)= YG;
              SubXe(4,2)=xl;              SubYe(4,2)= yl;
              SubXe(4,3)=Xe(lads+4);      SubYe(4,3)= Ye(lads+4);
           else
            write(dfich2,*) '*ERROR* El elemento ',n,' de tipo ',tipo,
     &                   ' no ha podido subdividirse en subelem.f'
            stop
          end if
    
      

c    	% == Tipo 80
        case(80)
          n1=0
          do i=1,nel
            if ((XG.eq.Xe(i)).and.(YG.eq.Ye(i))) then
                n1=i
                endif
             end do
        if (n1.eq.0) then
       
            write(dfich2,*) '*ERROR* El elemento ',n,' de tipo ',tipo,
     &                   ' no ha podido subdividirse en subelem.f'
            stop
        endif
        NumSub= 2;  
        SubXe= 0;
        SubYe= 0;
  
        do i= 1,2
          SubXe(i,1)= XG;          SubYe(i,1)= YG;
          SubXe(i,2)= Xe(n1+i);    SubYe(i,2)= Ye(n1+i);
          SubXe(i,3)= Xe(n1+i+1);  SubYe(i,3)= Ye(n1+i+1);
        end do


        end select
 



       end


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

      subroutine invcuad4(Xe,Ye,x,y,s,t,N)


c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Calcular las coordenadas locales en el elemento de los puntos de gauss
c				a partir de sus coordenadas globales
c     --------------------------------------

c      Inputs:
c			Xe		Coordenadas X de los nodos del elemento
c			Ye		Coordenadas Y de los nodos del elemento
c			x		Coordenada x global del punto de Gauss
c			y	    Coordenada y global del punto de Gauss

c      Outputs:
c			s		Coordenadas x local del punto de Gauss
c			t		Coordenadas y local del punto de Gauss			
c			N		Funciones de Forma
      
c-----[--.----+----.----+----.-----------------------------------------]


      implicit  none


      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2



      integer iter
      real*8	Xe(4),Ye(4),x,y,s,t,N(4),eps,error,XY(2),dNl(4,2),M(2,2),
     &		delta(2)			

      s= 0;
      t= 0;
      eps= 1.0d-10;
      Iter= 0;
      N= (/0.25, 0.25, 0.25, 0.25/);
      error = 1.0d0;

c	Iterativamente calcula las coordenadas locales del punto
      do iter=1,100
        if (error.lt.eps) exit

        XY(1)= dot_product(N,Xe);
        XY(2)= dot_product(N,Ye);
        
        XY(1)= x-XY(1);
        XY(2)= y-XY(2);
        
        error= sqrt( XY(1)**2 + XY(2)**2);
        
        dNl(1,1)= -(1-t)/4;
        dNl(1,2)= -(1-s)/4;
        dNl(2,1)=  (1-t)/4;
        dNl(2,2)= -(1+s)/4;
        dNl(3,1)=  (1+t)/4;
        dNl(3,2)=  (1+s)/4;
        dNl(4,1)= -(1+t)/4;
        dNl(4,2)=  (1-s)/4;
        M(1,1)= dot_product(dNl(:,1),Xe);
        M(1,2)= dot_product(dNl(:,2),Xe);
        M(2,1)= dot_product(dNl(:,1),Ye);
        M(2,2)= dot_product(dNl(:,2),Ye);
        

        call invertEug(M,2,2)
        delta= matmul(M,XY);               
        s= s + delta(1);
        t= t + delta(2);
      
        N(1)= (1-s)*(1-t)/4;
        N(2)= (1+s)*(1-t)/4;
        N(3)= (1+s)*(1+t)/4;
        N(4)= (1-s)*(1+t)/4;
      
      end do

      if (error.gt.eps) then
        write(dfich2,*) ' *WARNING* No convergence in INVCUAD4'
        stop
      endif

      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes


      subroutine invertEug(a,nmax,ndm)



c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Invert small square matrix

c      Inputs:
c         a(ndm,*) - Matrix to be inverted
c         nmax     - Size of upper submatrix to invert
c         ndm      - Dimension of array

c      Outputs:
c         a(ndm,*) - Submatrix replaces original terms, others not
c                    changed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2


      integer   i,j,n,ndm,nmax
      real*8    d, a(ndm,*)
      logical   flag

      save


      flag=.false.

      do n = 1,nmax
        if(a(n,n).ne.0.0d0) then
          d = 1.d0/a(n,n)
          do j = 1,nmax
            a(n,j) = -a(n,j)*d
          end do ! j

          do i = 1,nmax
            if(n.ne.i) then
              do j = 1,nmax
                if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
              end do ! j
            endif
            a(i,n) = a(i,n)*d
          end do ! i
          a(n,n) = d
        else
c          write(dfich2,*) ' *WARNING* Zero pivot in INVERTEUG'
c          write(dfich2,*) ' Utilizo la otra rutina para invertir'
          flag=.true. 
            exit         
        endif
      end do ! n

      if (flag) then !Utilizo la otra rutina para invertir
        call INV_J(a,ndm) 
      endif

      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
C
      SUBROUTINE TypeXelement(OUTDIR,LENOUTDIR,JELEM,NNODE,NelmX,
     *			 ix,TypeXe)            

      IMPLICIT NONE

C       Rutina que busca los num. de nodos del elemto en cuestión (conectividad)
c       Tb. busca para cada uno de esos nodos del elemento en cuestión
c       qué TypeX tiene y los guarda en TypeXe

c       OUTPUTS: ix  :  conectividad (num. de nodos del elemento)
c                TypeXe  : TypeX para cada nodo del elemento


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2

      character*256 OUTDIR
      integer LENOUTDIR,i,j,k
      integer JELEM,NNODE,NelmX,numelem,valores(NNODE)
      integer TypeXe(NNODE)
 
      integer ne,nn(NNODE),ix(NNODE)
      character*1 basu !para quitar las comas de la conectividad
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\files\TopXTypeX')
      OPEN(69,FILE=OUTDIR(1:LENOUTDIR)//'\files\TopX')
      
      do i=1,NelmX
         READ(68,*) numelem,(valores(j),j=1,NNODE)
         READ(69,1000) ne,basu,nn(1),basu,nn(2),basu,nn(3),basu,nn(4)
         if (numelem.eq.JELEM) then
            do k=1,NNODE
              TypeXe(k)=valores(k)
              ix(k)=nn(k)
              end do
            exit
         end if  
      end do   	  
      

      CLOSE(68)
           CLOSE(69)


c	write(dfich2,*) (TypeXe(k),k=1,NNODE)
c	write(dfich2,*) (ix(k),k=1,NNODE)

1000  format(I6,A1,I6,A1,I6,A1,I6,A1,I6)

      RETURN

      END


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes




      SUBROUTINE K_U12(E,Nu,AMATRX,NDOFEL,NNODE,dimens,MCRD,
     *	         COORDS,PSS,NnodX,ix,TypeXe,Dist,XYC0,XYCPrev,         	
     *	         gint,sg,Xe,Ye,flag,BatG,DBatG,JatG)

      IMPLICIT NONE

      INTEGER NDOFEL,NNODE,dimens,MCRD,PSS,NnodX,gint,flag,pos
      INTEGER l,i,j,kk,TypeXe(NNODE),ix(NNODE)

      REAL*8 E,Nu,Dist(NnodX,3),sg(3,*)
      REAL*8 AMATRX(NDOFEL,NDOFEL),XYC0(2),XYCPrev(2)
      REAL*8 Xe(2*NNODE),Ye(2*NNODE),COORDS(MCRD,NNODE),xl(dimens,NNODE)
      REAL*8 xsj(gint),shp(3,4)
      REAL*8 dNF(NNODE,2,4),Fnode(NNODE,4),H,Hnode(NNODE)
      REAL*8 B(3,NDOFEL), DB(3,NDOFEL), BT(NDOFEL,3), D(3,3)
      REAL*8 BatG(3*gint,NDOFEL),DBatG(3*gint,NDOFEL),JatG(gint)

      LOGICAL NodeType1,NodeType2

c     NOTES:   
c      Routine shapef2D is called to compute standard shape functions,
c      derivatives and jacobian at integration points. This routine outputs: 
c  
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c
c      Local coordinates of integration points are passed in sg(1,*), sg(2,*)
c      Integration weights are passed in sg(3,*)
c      



C23456789012345678901234567890123456789012345678901234567890123456789012

c	initialize AMATRX and logical variables
      CALL initializeM(AMATRX,NDOFEL,NDOFEL)          
      NodeType1=.false.
      NodeType2=.false.

c     Reduce info passed thru COORDS (3D) to xl (2D)
      DO i=1,dimens
       DO j=1,NNODE
         xl(i,j)=COORDS(i,j)
       END DO
      END DO

c     Define constitutive stress-strain elastic matrix
      CALL CALC_D(PSS,D,E,Nu)

c	Specify the type of nodal enrichment
      DO i=1,NNODE
          IF (TypeXe(i).eq.1) THEN
            NodeType1=.true.
          ELSEIF (TypeXe(i).eq.2) THEN
            NodeType2=.true.
          END IF
      END DO

c     Numerical integration loop over gint integration points

      DO l = 1,gint

c		Compute shape functions, derivatives and jacobian at integration point
          CALL shapef2D(sg(1,l),xl,shp,xsj(l),dimens,NNODE,ix,.false.)
            IF (flag.eq.1) THEN !Element is subdivided for integration
            xsj(l) = sg(3,l) !The integration weight includes the jacobian
            ELSE !Element is not subdivided. Standard integration
              xsj(l) = xsj(l)*sg(3,l) 
            ENDIF

c		Value of the Heaviside function at integration point
c         (This call is also used to store the values of H 
c          at nodes of the element for modified enrichment)
            IF (NodeType1) THEN
              CALL heaviside(NnodX,Dist,NNODE,ix,shp,H,Hnode)
            ENDIF

c	    Derivatives of shape functions Ni times enrichment functions Fj at integration point
c         (This call is also used to compute the derivatives of shape functions Ni times 
c          enrichment functions Fj at nodes of the element for modified enrichment)
          IF (NodeType2) THEN
              CALL fCrackTip(XYC0,XYCPrev,shp,Xe,Ye,dNF,Fnode)	
            ENDIF


c         STIFFNESS MATRIX COMPUTATION: 
c         Assembly of element matrix B (denoted as B) at integration point

            CALL initializeM(B,3,NDOFEL)
              Pos=1			
c		  Loop over nodes
              DO i= 1,NNODE

c			Contribution to B of derivatives of standard shape functions
                  B(1,Pos)  = shp(1,i)
                  B(2,Pos+1)= shp(2,i)
                  B(3,Pos)  = shp(2,i)
                  B(3,Pos+1)= shp(1,i)		
                                    
c			Contribution to B of derivatives of shape functions times Heaviside function 
                  IF (TypeXe(i).eq.1) THEN
                    B(1,2+Pos)  = shp(1,i)*(H-Hnode(i))
                    B(2,3+Pos)  = shp(2,i)*(H-Hnode(i))
                    B(3,2+Pos)  = shp(2,i)*(H-Hnode(i))
                    B(3,3+Pos)  = shp(1,i)*(H-Hnode(i))		

c			Contribution to B of derivatives of shape functions times crack tip functions 
                  ELSEIF(TypeXe(i).eq.2) THEN
                    DO kk= 1,4
                      B(1,2*kk+2+Pos)= dNF(i,1,kk)-shp(1,i)*Fnode(i,kk)
                    B(2,2*kk+3+Pos)= dNF(i,2,kk)-shp(2,i)*Fnode(i,kk)
                      B(3,2*kk+2+Pos)= dNF(i,2,kk)-shp(2,i)*Fnode(i,kk)
                      B(3,2*kk+3+Pos)= dNF(i,1,kk)-shp(1,i)*Fnode(i,kk)
                    END DO					
                  END IF 
                  Pos=Pos+12 !Each node has 12 dof

              END DO ! i = End loop over element nodes

              DB=matmul(D,B) ! Matrix D*B
              BT=transpose(B)  ! B transpose 
c		  Integration of Bt*D*B
              AMATRX= AMATRX + matmul(BT,DB)*xsj(l)

c       Store information at each integration point for further post-processing
            DO  i=1,3
              DO j=1,NDOFEL
               BatG(3*(l-1)+i,j)=B(i,j)
               DBatG(3*(l-1)+i,j)=DB(i,j)
              END DO
              END DO
            JatG(l)=xsj(l)
      
        END DO ! l = End loop for each integration point
 
      RETURN
      END



C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes


      subroutine shapef2D(ss,xl,shp,xsj,ndm,nel,ix,flg)


c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Computes shape function and derivatives for
c               quadrilateral elements

c      Inputs:
c         ss(2)     - Natural coordinates for point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element
c         ix(*)     - Nodes attached to element. 
c         flg       - Flag, compute global x/y derivatives if false,
c                           else derivatives are w/r natural coords.

c      Outputs:
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer dfich,dfich2   
      common /debugfich/ dfich,dfich2



      logical   flg
      integer   ndm,nel, i,j,k, ix(*)
      real*8    xsj, temp
      real*8    shp(3,*),xl(ndm,*), s(4),t(4),xs(3,2),sx(2,2),ss(2)

      save

c     Set values of half natural coords at nodes

      data s/-0.5d0,0.5d0,0.5d0,-0.5d0/,t/-0.5d0,-0.5d0,0.5d0,0.5d0/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        OJO  Comento las rutinas que no vamos a usar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      write(dfich2,*) ' ********************* '
c      write(dfich2,*) ((xl(i,j),j=1,nel),i=1,ndm)
c      write(dfich2,*) ' ********************* '



c     Form 4-node quadrilateral shape functions

      if(nel.eq.4 .and. .not.flg) then
         call shapef(ss(1),ss(2),xl,shp,xsj,ndm,flg)
      elseif(nel.eq.16) then
c        call shp2dc(ss,xl,shp,temp,xsj,1, flg)
         write(dfich2,*) '**ERROR** Comentada la rutina en shapef2D.f'
         stop
      elseif(nel.eq.12) then
c        call shp2ds(ss,xl,shp,temp,xsj,1, flg)
         write(dfich2,*) '**ERROR** Comentada la rutina en shapef2D.f'
         stop
      else
        do i = 1,4
          shp(3,i) = (0.5d0+s(i)*ss(1))*(0.5d0+t(i)*ss(2))
          shp(1,i) = s(i)*(0.5d0+t(i)*ss(2))
          shp(2,i) = t(i)*(0.5d0+s(i)*ss(1))
        end do ! i

c       Form triangle by adding third and fourth together

        if(nel.eq.3) then
          do i = 1,3
            shp(i,3) = shp(i,3)+shp(i,4)
          end do ! i
        end if

c       Add quadratic terms if necessary

        if(nel.gt.4) then
c          call shap2(ss(1),ss(2),shp,ix,nel)
           write(dfich2,*) '**ERROR** Comentada la rutina en shapef2D.f'
           stop
        endif 
         
c       Construct jacobian and its inverse

        do i = 1,max(3,ndm)
          do j = 1,2
            xs(i,j) = 0.0d0
            do k = 1,nel
              xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
            end do ! k
          end do ! j
        end do ! i
        if(ndm.eq.2) then
          xsj = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
        elseif(ndm.eq.3) then
          xsj = sqrt((xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))**2
     &             + (xs(3,1)*xs(1,2)-xs(3,2)*xs(1,1))**2
     &             + (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))**2)
        endif
        if(.not.flg) then
          if(xsj.eq.0.0d0) then
            temp = 1.0d0
          else
            temp = 1.d0/xsj
          endif
          sx(1,1) = xs(2,2)*temp
          sx(2,2) = xs(1,1)*temp
          sx(1,2) =-xs(1,2)*temp
          sx(2,1) =-xs(2,1)*temp

c         Form global derivatives

          do i = 1,nel
            temp     = shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
            shp(2,i) = shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
            shp(1,i) = temp
          end do ! i

c         Return center node in hierarchical form for 8-nodes

          if(nel.eq.8) then
            temp     = shp(1,9)*sx(1,1)+shp(2,9)*sx(2,1)
            shp(2,9) = shp(1,9)*sx(1,2)+shp(2,9)*sx(2,2)
            shp(1,9) = temp
          endif
        endif
      endif

      end

C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

 
      subroutine shapef(s,t,xl,shp,xsj,ndm,flg)


c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Shape function routine for 4-node isoparametric
c               quadrilaterals

c      Inputs:
c         s,t       - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         flg       - Flag, Compute global derivatives if true,
c                           else compute derivatives w/r natural coords.

c      Outputs:
c         shp(3,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx  or dN_i/dxi_1
c                     shp(2,i) = dN_i/dy  or dN_i/dxi_2
c                     shp(3,i) = N_i
c         xsj       - Jacobian determinant at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flg
      integer   ndm
      real*8    xo,xs,xt, yo,ys,yt, xsm,xsp,xtm,xtp, ysm,ysp,ytm,ytp
      real*8    s,t, xsj,xsj1, sh,th,sp,tp,sm,tm, xl(ndm,4),shp(3,4)

      save

c     Set up interpolations

      sh = 0.5d0*s
      th = 0.5d0*t
      sp = 0.5d0 + sh
      tp = 0.5d0 + th
      sm = 0.5d0 - sh
      tm = 0.5d0 - th
      shp(3,1) =   sm*tm
      shp(3,2) =   sp*tm
      shp(3,3) =   sp*tp
      shp(3,4) =   sm*tp

c     Set up natural coordinate functions (times 4)

      xo =  xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4)
      xs = -xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4) + xo*t
      xt = -xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4) + xo*s
      yo =  xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4)
      ys = -xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4) + yo*t
      yt = -xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4) + yo*s

c     Compute jacobian (times 16)

      xsj1 = xs*yt - xt*ys

c     Divide jacobian by 16 (multiply by .0625)

      xsj = 0.0625d0*xsj1
      if(.not.flg) then
        if(xsj1.eq.0.0d0) then
          xsj1 = 1.0d0
        else
          xsj1 = 1.0d0/xsj1
        endif

c       Divide functions by jacobian

        xs  = (xs+xs)*xsj1
        xt  = (xt+xt)*xsj1
        ys  = (ys+ys)*xsj1
        yt  = (yt+yt)*xsj1

c       Multiply by interpolations

        ytm =  yt*tm
        ysm =  ys*sm
        ytp =  yt*tp
        ysp =  ys*sp
        xtm =  xt*tm
        xsm =  xs*sm
        xtp =  xt*tp
        xsp =  xs*sp

c       Compute shape functions

        shp(1,1) = - ytm+ysm
        shp(1,2) =   ytm+ysp
        shp(1,3) =   ytp-ysp
        shp(1,4) = - ytp-ysm
        shp(2,1) =   xtm-xsm
        shp(2,2) = - xtm-xsp
        shp(2,3) = - xtp+xsp
        shp(2,4) =   xtp+xsm
      endif

      end


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes
 
      subroutine fCrackTip(XYC0,XYCPrev,shp,Xe,Ye,dNF,Fnode)



c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Calcular las funciones de extremo de grieta en un punto.  
c			 Se hace una transformación de coordenadas a polares con origen  
c			 en extremo de grieta. Se calculan las funciones de extremo de 
c			 grieta y luego sus derivadas respecto a coordenadas globales.
c     --------------------------------------

c      Inputs:
c			X		Coordenada X punto de Gauss
c			Y		Coordenada Y punto de Gauss
c			XYC		Coordenadas de los puntos de Grieta
c			shp		Funciones de forma y sus derivadas en el punto
c       

c      Outputs:
c			dNF(nel,2,4)	Derivadas de N*funciones de extremo de grieta en PG
c			Fnode(nel,2,4) Funciones de extremo de grieta en nodos del elem

c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none



c      integer dfich,dfich2   
c      common /debugfich/ dfich,dfich2



      integer	i,j,k
      real*8	X,Y,XYC0(2),XYCPrev(2),shp(3,4),longG
      real*8  tx,ty,rt,rn,r,r2
      real*8  rr,theta,dthx,dthy,drrx,drry,F(4),dF(4,2),dNF(4,2,4)	
      real*8  Xe(8),Ye(8),betta,Fnode(4,4)
                 

c	XYC0=XYC(NCP,1:2)

c	Coordenadas globales del punto de Gauss
      X= dot_product(shp(3,1:4),Xe(1:4))
      Y= dot_product(shp(3,1:4),Ye(1:4))


c	Calculo de vectores unitarios tangente y normal en extremo de grieta
      longG= sqrt((XYCPrev(1)-XYC0(1))**2+(XYCPrev(2)-XYC0(2))**2);
      tx= (XYC0(1)-XYCPrev(1))/longG;
      ty= (XYC0(2)-XYCPrev(2))/longG;

c	% --- Coordenadas polares en extremo de grieta 
      rt=  (X-XYC0(1))*tx + (Y-XYC0(2))*ty;
      rn= -(X-XYC0(1))*ty + (Y-XYC0(2))*tx;
      r2= (X-XYC0(1))**2 + (Y-XYC0(2))**2;
      r= sqrt(r2);
      rr= sqrt(r);
      theta= atan2(rn,rt);
      betta= atan2(Y-XYC0(2),X-XYC0(1))

c      % --- Derivadas de coordenadas de extremo de grieta
c      dthx= -(ty*rt+tx*rn)/r2;
c      dthy=  (tx*rt-ty*rn)/r2;
      dthx=-(Y-XYC0(2))/r2
      dthy=(X-XYC0(1))/r2
      drrx=(X-XYC0(1))/2/r/rr;
      drry=(Y-XYC0(2))/2/r/rr;

c      % --- Funciones singulares extremo de grieta
c			Fi(x)
      F(1)= rr*sin(theta/2);
      F(2)= rr*cos(theta/2);
      F(3)= rr*sin(theta/2)*sin(theta);
      F(4)= rr*cos(theta/2)*sin(theta);

c      % --- Derivadas de funciones singulares extremo de grieta
c			dFi/dx, dFi/dy
c      dF(1,1)= sin(theta/2)*drrx + rr*cos(theta/2)*dthx/2;
c      dF(1,2)= sin(theta/2)*drry + rr*cos(theta/2)*dthy/2;
c      dF(2,1)= cos(theta/2)*drrx - rr*sin(theta/2)*dthx/2;
c      dF(2,2)= cos(theta/2)*drry - rr*sin(theta/2)*dthy/2;
      dF(1,1)= sin(theta/2-betta)/2/rr
      dF(1,2)= cos(theta/2-betta)/2/rr
      dF(2,1)= dF(1,2)
      dF(2,2)= -dF(1,1)
      dF(3,1)= sin(theta/2)*sin(theta)*drrx + rr*( cos(theta/2)*
     &	     sin(theta)/2 + sin(theta/2)*cos(theta) )*dthx;
      dF(3,2)= sin(theta/2)*sin(theta)*drry + rr*( cos(theta/2)*
     &	     sin(theta)/2 + sin(theta/2)*cos(theta) )*dthy;
      dF(4,1)= cos(theta/2)*sin(theta)*drrx + rr*(-sin(theta/2)*
     &	     sin(theta)/2 + cos(theta/2)*cos(theta) )*dthx;
      dF(4,2)= cos(theta/2)*sin(theta)*drry + rr*(-sin(theta/2)*
     &	     sin(theta)/2 + cos(theta/2)*cos(theta) )*dthy;

c      %--- Derivadas de interpolacion
      do i=1,4
        do j=1,2
          do k=1,4
            dNF(i,j,k)= shp(j,i)*F(k) + shp(3,i)*dF(k,j);
          end do
        end do
      end do


c	Para conocer las deriv funciones de extremo de grieta en los nodos del elemento (Eugenio)
c       y poder tener los gdl físicos en los dos primeros gdl. 

c       ******* Bucle para las posiciones nodales *********

      do i=1,4

c	  Coordenadas globales del nodo i
        X= Xe(i)
        Y= Ye(i)

c	   --- Coordenadas polares en extremo de grieta 
        rt=  (X-XYC0(1))*tx + (Y-XYC0(2))*ty;
        rn= -(X-XYC0(1))*ty + (Y-XYC0(2))*tx;
        r2= (X-XYC0(1))**2 + (Y-XYC0(2))**2;
        r= sqrt(r2);
        rr= sqrt(r);
        theta= atan2(rn,rt);
        betta= atan2(Y-XYC0(2),X-XYC0(1))


C      write(dfich2,*) 'Angulo theta=',theta
C      write(dfich2,*) 'Angulo beta=',betta

c       --- Derivadas de coordenadas de extremo de grieta
        dthx=-(Y-XYC0(2))/r2
        dthy=(X-XYC0(1))/r2
        drrx=(X-XYC0(1))/2/r/rr;
        drry=(Y-XYC0(2))/2/r/rr;

c       --- Funciones singulares extremo de grieta
c			Fi(x)
        Fnode(i,1)= rr*sin(theta/2);
        Fnode(i,2)= rr*cos(theta/2);
        Fnode(i,3)= rr*sin(theta/2)*sin(theta);
        Fnode(i,4)= rr*cos(theta/2)*sin(theta);

c        --- Derivadas de funciones singulares extremo de grieta
c			dFi/dx, dFi/dy
        
c        No hacen falta las derivadas de F porque el valor de F en los nodos 
c        es una constante que no depende de x,y en la expresión de la aprox. de desplaz. 

c      %--- Derivadas de interpolacion  /// ESTO YA NO LO HAGO AQUÍ. LO HAGO EN LA RUTINA K_U12
c        do j=1,2
c          do k=1,4
c            dFInode(i,j,k)= shp(j,i)*F(k);
c          end do
c        end do     

      end do



 

      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes


      subroutine heaviside(NnodX,Dist,nel,ix,shp,H,Hnode)


c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: Calcular la función de Heaviside en un punto x. 
c			 Evalua la distancia del punto a la grieta interpolando
c			 la distancia de los nodos a la grieta con las funciones de
c			 forma en el punto x.
c     --------------------------------------

c      Inputs:
c			Dist(n,3) = Dist(nel,3)	Matriz de distancia de nodos a grieta
c			ix(*) = ix(nel)		Conectividad del elemento
c			shp(3,nel)	Funciones de forma y sus derivadas en el punto

c      Outputs:
c			H 			Función heaviside en Punto x,y
c             Hnode   	Función heaviside en nodos del elemento
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

c      integer dfich,dfich2   
c      common /debugfich/ dfich,dfich2
      
      integer i,nel,ix(nel),j,NnodX
      real*8  VDist(4),Dist(NnodX,3),DisPG,shp(3,nel),H,Hnode(4)
      
c	Distancias de los nodos del elemento a la grieta extraídas de la matriz [Dist]
      do i=1,4
        do j=1,NnodX
          if (Dist(j,1).eq.ix(i)) then
            VDist(i)=Dist(j,2)
            exit
          endif
        enddo
      end do

c	Proyectar la distancia en los nodos al punto x
      DisPG=dot_product(VDist,shp(3,1:4))

c	Función de Heaviside en el punto. 
c		(Unitaria y de signo igual al de la distancia evaluada en x)
      H=sign(1.0d0,DisPG)
c      write(dfich2,*) 'Heaviside H=',H

c	Para conocer la función de Heaviside en los nodos del elemento (Eugenio)
c       y poder tener los gdl físicos en los dos primeros gdl. 
c		(Unitaria y de signo igual al de la distancia evaluada en x)
      do i=1,4
        Hnode(i)=sign(1.0d0,VDist(i))
      end do

c      write(dfich2,*) 'Heaviside en nodos Hnode=',(Hnode(i),i=1,4)

      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Yes

      SUBROUTINE SVARS_U12(JTYPE,JELEM,SVARS,NSVARS,U,Dof,BatG,DBatG,
     *                     JatG,gint,mpg,xypg)

C      Calculates and/or stores the following magnitudes at the element integration points, 
C      storing them in SVARS: strains, stresses, strain energy density, dv/dx, du/dy, jacobian, 
C      dNi/dx, dNi/dy, global coordinates of integration points.



      IMPLICIT NONE
      
      INTEGER i,j,k,NSVARS, Dof, gint, JTYPE,JELEM,mpg
      REAL*8  SVARS(NSVARS), U(Dof),BatG(3*gint,Dof),DBatG(3*gint,Dof)
      REAL*8  JatG(gint),B(3,Dof),DB(3,Dof),Bdvdx(3,Dof),Bdudy(3,Dof)
      REAL*8  EPS(3),SIG(3),W,dvdx(3),dudy(3),JAC,xypg(2,mpg)

C23456789012345678901234567890123456789012345678901234567890123456789012   

C     First value stored in SVARS is the total number of integration points
C     of the enriched element 
      SVARS(1)=gint

        DO i=1,gint
           JAC=JatG(i)
           DO k=1,3
             DO j=1,Dof    
               B(k,j)=BatG(3*(i-1)+k,j)
               Bdvdx(k,j)=B(k,j) ! For computation of dv/dx
               Bdudy(k,j)=B(k,j) ! For computation of du/dy
               DB(k,j)=DBatG(3*(i-1)+k,j) 
             END DO
           END DO
           CALL MULT_V(B,3,Dof,U,EPS,3) ! Compute strains EPS
           CALL MULT_V(DB,3,Dof,U,SIG,3) ! Compute stresses SIG
           W=0.5d0*(EPS(1)*SIG(1)+EPS(2)*SIG(2)+EPS(3)*SIG(3))

C          Computation of dv/dx & du/dy
C          Set to zero positions in the 3rd row of B associated with dN/dy 
           DO j=1,Dof,2
               Bdvdx(3,j)=0.0d0
           END DO
           CALL MULT_V(Bdvdx,3,Dof,U,dvdx,3) !compute dv/dx, stored in dvdx(3)                                 

C          Set to zero positions in the 3rd row of B associated with dN/dx 
           DO j=2,Dof,2
               Bdudy(3,j)=0.0d0
           END DO
           CALL MULT_V(Bdudy,3,Dof,U,dudy,3) !compute du/dy, stored in dudy(3)                                                                  

C          Store in SVARS the following information at integration points
           SVARS(1+20*(i-1)+1)=EPS(1)
           SVARS(1+20*(i-1)+2)=EPS(2)
           SVARS(1+20*(i-1)+3)=EPS(3)
           SVARS(1+20*(i-1)+4)=SIG(1)
           SVARS(1+20*(i-1)+5)=SIG(2)
           SVARS(1+20*(i-1)+6)=SIG(3)
           SVARS(1+20*(i-1)+7)=W
           SVARS(1+20*(i-1)+8)=dvdx(3)
           SVARS(1+20*(i-1)+9)=dudy(3)
           SVARS(1+20*(i-1)+10)=JAC   ! Jacobian includes integration weight

C          Store in SVARS the shape functions derivatives dNi/dx, dNi/dy for external computation
C          of dq/dx, dq/dy (used in domain interaction integrals). 
C          (we take them from the positions associated with the standard dofs)
           SVARS(1+20*(i-1)+11)=B(1,1)
           SVARS(1+20*(i-1)+12)=B(1,13)
           SVARS(1+20*(i-1)+13)=B(1,25)
           SVARS(1+20*(i-1)+14)=B(1,37)
           SVARS(1+20*(i-1)+15)=B(2,2)
           SVARS(1+20*(i-1)+16)=B(2,14)
           SVARS(1+20*(i-1)+17)=B(2,26)
           SVARS(1+20*(i-1)+18)=B(2,38)

C          Store in SVARS the global coordinates of integration points
           SVARS(1+20*(i-1)+19)=xypg(1,i)
           SVARS(1+20*(i-1)+20)=xypg(2,i)
        
        END DO !i loop over all integration points of the element

        RETURN
      END



C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C