      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      
        EMOD=PROPS(1)
        ENU=PROPS(2)
        EPS0=PROPS(3)
        EPS1=PROPS(4)

        EBULK3=EMOD/(1-2*ENU)
        EG2=EMOD/(1+ENU)
        EG=EG2/2
        ELAM=(EBULK3-EG2)/3
        DO I=1,NDI
          DO J=1,NDI
            DDSDDE(I,J)=ELAM
          END DO
          DDSDDE(I,I)=EG2
        END DO
        
        DO I=NDI+1,NTENS
          DDSDDE(I,I)=EG
        END DO      
      RETURN
      END
