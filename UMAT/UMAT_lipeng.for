C                   V2.1
C     V2.1 说明:
C       添加基体层-基体材料的UMAT支持，注意复合材料用UMAT-Composite 基体材料用UMAT-Matrix
C			注意
C     1. 在UMAT的计算中中间变量都是用传统的应力应变顺序
C     FIBER,RESIN中模量的顺序为E11,E22,E33,v23,v13,v12,G23,G13,G12    
C			SSM SSF顺序为传统顺序 11 22 33 23 13 12  
C     2. 在子过程中，注意默认类型要与输入参数匹配  一般所有的BUG都是这个问题
C       不要声明了与变量名不符合的类型后 
C
C     
C                   PROPS范例
C     #基体的九个弹性常数加上三个拉压剪模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
C     4100., 4100., 4100., 0.46, 0.46,   0.46,  1400.,  1400., 1400.,121.,210.,76.,
C     #纤维的九个弹性常数加上两个拉压模量 E11,E22,E33,V12,V13,V23,G12,G13,G23
C     276000.,19000., 19000., 0.2, 0.2, 0.36, 27000.,27000., 6980., 4850., 3000.,
C     # VF ALPHA BETA
C     0.575,  0.3,    0.3,
C     #MSEG
C     6., 
C     #ETM 基体塑性阶段折线强度和模量
C     42.7,   53.7,   76.5,  101.5,  111.3,   125.,  
C     4100.,  2500.,  2000.,  1400.,   800.,   410.,   
C     #STATEV(13): SSF(6) SSM(6) 衰减次数
C     1., 2.,3.,4.,5.,6.,7.,8.,9., 10.,11.,12.,
C     #KT22  KTT22 KC22 K12 K23 
C     2.322,  5.083,  1.656,  1.447,1.87,1.72,1.69 ,
C     LMISE RF 破坏判据的选择
C     54.4,   0.01,  1
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)

        IF (CMNAME(1:14) .EQ. 'UMAT-COMPOSITE') THEN
          CALL UMAT_MAT1(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

        ELSE IF(CMNAME(1:11) .EQ. 'UMAT-MATRIX') THEN
          CALL UMAT_MAT2(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
        END IF
      END SUBROUTINE
      
      SUBROUTINE UMAT_MAT1(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)

      DIMENSION SCFS(7) !KT22,KTT22,KT22BI,KC22BI,KC22,K23,K12
      
      REAL*8 LMISE,FAILCHS
      
      DIMENSION RESIN(12),FIBER(11),EGM(9),SUUM(3),EGF(9),SUUF(2),
     1SF(6,6),SM(6,6),SL(6,6),A(6,6),B(6,6),ETM(2,20),SSF(6),SSM(6),
     2DSTRESS(6)

      CHARACTER*256 JOBNAME,CP,DIR
      INTEGER LJ,LD 
      
      CALL GETJOBNAME(JOBNAME, LJ )
      DIR='E:/UMAT-OUTPUT/'
      LD=LEN_TRIM(DIR)
      CP(1:LD)=DIR
      CP((LD+1):(LJ+LD))=JOBNAME(1:LJ)
      !OPEN(101,FILE='E:\dadelamination.TXT',STATUS='UNKNOWN')
      OPEN(102,FILE=CP(1:LJ+LD)//'-matrix.txt',STATUS='UNKNOWN')
      OPEN(103,FILE=CP(1:LJ+LD)//'-fiber.txt',STATUS='UNKNOWN')
      
      CALL KADJUST(PROPS,NPROPS,STATEV,NSTATV,RESIN,FIBER,VF,ALFA,
     1 BETA,MSEG,ETM,SSF,SSM,RF,SCFS,LMISE,FAILCHS,NOEL,KINC)
            
      CALL KFIBER(FIBER,EGF,SUUF,SF,NOEL)
      
      CALL KMATRIX(RESIN,EGM,SM,SSM,ETM,MSEG,STATEV(13),RF,
     1 SUUM,NOEL,KINC)

      ! 将判断强度破坏放在这里 这里判断上一increment的纤维基体应力是否可能出现破坏
      ! 不放在后面的原因是如果这一步开始时没破坏但是结束时在某个积分点处发生破坏 
      ! 那么当前步的计算结果可能被废弃
      CALL KSTRTH(SSF,SSM,SUUM,SUUF,STATEV,NSTATV,FAILCHS,NOEL,KINC)
      !IF(STATEV(14).EQ.1)THEN
      !  CALL XIT
      !END IF

      CALL KBRIGE(VF,EGF,SF,EGM,SM,ALFA,BETA,A,B,SL,NOEL)
      
      CALL ABQ2NORM_1D(STRESS)
      ! 注意DSTRAN的应变分量顺序
      CALL KSTRESS(VF,SF,SM,A,B,SL,DDSDDE,STRESS,DSTRAN,STATEV,
     1 NSTATV,RF,NOEL)

      CALL KINTER(DDSDDE,DSTRAN,DSTRESS,A,B,VF,SSF,SSM,SCFS,LMISE,
     1 STATEV,NSTATV,NOEL)
      
      CALL NORM2ABQ_1D(STRESS)
      CALL NORM2ABQ_2D(DDSDDE)

      IF(NOEL.EQ.1) THEN
        WRITE(*,'(A6,I7,\)') 'KINC',KINC,'KSTEP',KSTEP,'NOEL',NOEL
        PRINT *,''
        !PRINT *,NSTATV
        PRINT *,"STATEV"
        WRITE(*,"(12F20.5)") (STATEV(I),I=1,12)
        !PRINT *,"SCFS"
        !WRITE(*,*) SCFS
        PRINT *,"SF"
        WRITE(*,'(6F20.15)') ((SF(I,J),J=1,6),I=1,6)
        PRINT *,"SM"
        WRITE(*,'(6F20.15)') ((SM(I,J),J=1,6),I=1,6)
        PRINT *,"A"
        WRITE(*,'(6F20.15)') ((A(I,J),J=1,6),I=1,6)
        !PRINT *,"B"
        !WRITE(*,'(6F20.15)') ((B(I,J),J=1,6),I=1,6)      
        !PRINT *,"SL"
        !WRITE(*,'(6F20.15)') ((SL(I,J),J=1,6),I=1,6)
        PRINT *,"DSTRAN"
        WRITE(*,'(6F20.15)') (DSTRAN(J),J=1,6)      
        PRINT *,"STRESS"
        WRITE(*,'(6F20.5)') (STRESS(J),J=1,6)   
        PRINT *,"SSF"
        WRITE(*,'(6F20.5)') (SSF(J),J=1,6)
        PRINT *,"SSM"
        WRITE(*,'(6F20.5)') (SSM(J),J=1,6)
        PRINT *,"DDSDDE"
        WRITE(*,'(6F20.5)') ((DDSDDE(I,J),J=1,6),I=1,6)
      END IF

      RETURN 
      END
C
C***********************************************************************
C
C     从PROPS中读取材料参数 和 从STATEV中读取纤维和基体的应力状态
C     注意PROPS中设计的顺序是12 13 23 所以要经过转换
C 输入：
C     PROPS: 在UMAT定义的时候输入的数组
C     NPROPS: PROPS的大小
C     STATEV: 求解过程中存入的数据（纤维和基体平均应力）
C     NSTATV: STATEV数组的大小
C     NOEL: 单元号
C     KINC: 增量步编号
C 输出:
C     RESIN: 基体的12个材料常数数组
C     FIBER: 纤维的11个材料常数数组
C     VF:   纤维体积含量
C     ALFA: 桥联参数ALPHA
C     BETA: 桥梁参数 BETA
C     MSEG: 基体性能段数
C     ETM: 基体的分段性能参数 数组
C     SSF: 纤维的平均应力
C     SSM: 基体的平均应力
C     RF: 基体的衰减参数
C     KT22,KC22,K12,K13: 应力集中系数
C     LMISE:  界面开裂强度
C
      SUBROUTINE KADJUST(PROPS,NPROPS,STATEV,NSTATV,RESIN,FIBER,VF,ALFA,
     1BETA,MSEG,ETM,SSF,SSM,RF,SCFS,LMISE,FAILCHS,NOEL,KINC)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION PROPS(NPROPS),STATEV(NSTATV),RESIN(12),FIBER(11),
     1  ETM(2,20),SSF(6),SSM(6),SCFS(7)
        REAL*8 LMISE,FAILCHS

        REAL*8 KT22,KTT22,KT22BI,KC22BI,KC22,K23,K12
        
        N=0
        DO I=1,12
          RESIN(I)=PROPS(I)
        ENDDO
        TMP=RESIN(4)
        RESIN(4)=RESIN(6)
        RESIN(6)=TMP
        TMP=RESIN(7)
        RESIN(7)=RESIN(9)
        RESIN(9)=TMP
        N=N+12

        DO I=1,11
          FIBER(I)=PROPS(I+N)
        ENDDO
        TMP=FIBER(4)
        FIBER(4)=FIBER(6)
        FIBER(6)=TMP
        TMP=FIBER(7)
        FIBER(7)=FIBER(9)
        FIBER(9)=TMP
        N=N+11

        VF=PROPS(N+1)
        ALFA=PROPS(N+2)
        BETA=PROPS(N+3)
        N=N+3

        MSEG=PROPS(N+1)
        N=N+1

        DO I=1,MSEG
        	ETM(1,I)=PROPS(N+I)
        END DO
        N=N+MSEG

        DO I=1,MSEG
        	ETM(2,I)=PROPS(N+I)
        END DO
        N=N+MSEG

        IF(KINC.EQ.1)THEN
          DO I=1,12
            STATEV(I)=PROPS(I+N)
          END DO
        ENDIF
        N=N+12

        DO I=1,6
          SSF(I)=STATEV(I)
          SSM(I)=STATEV(I+6)
        END DO

        !!KT22,KTT22,KC22,K23,K12
        !DO I=1,5
        !  SCFS(I)=PROPS(N+I)
        !END DO
        !N=N+5
        !KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI
        KT22=PROPS(N+1)
        KTT22=PROPS(N+2)
        KC22=PROPS(N+3)
        K12=PROPS(N+4)
        K23=PROPS(N+5)
        KT22BI=PROPS(N+6)! KT22
        KC22BI=PROPS(N+7)! KC22
        N=N+7
        !KT22,KTT22,KT22BI,KC22BI,KC22,K23,K12
        SCFS(1)=KT22
        SCFS(2)=KTT22
        SCFS(3)=KT22BI
        SCFS(4)=KC22BI
        SCFS(5)=KC22
        SCFS(6)=K23
        SCFS(7)=K12

        LMISE=PROPS(N+1)
        RF=PROPS(N+2)
        FAILCHS=PROPS(N+3)

        RETURN
      END
C
C***********************************************************************
C
C     计算纤维的柔度矩阵 注意 返回为传统顺序
C 输入参数
C     FIBER: 纤维的11个参数,
C 输出：
C     EGF: 纤维的9个弹性常数
C     SUUF: 纤维的拉压强度 
C     SF： 纤维的柔度矩阵
      SUBROUTINE KFIBER(FIBER,EGF,SUUF,SF,NOEL)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION FIBER(11)
        DIMENSION EGF(9),SUUF(2),SF(6,6)

        DO I=1,9
          EGF(I)=FIBER(I)
        END DO
        DO I=1,2
          SUUF(I)=FIBER(I+9)
        END DO
        CALL KELASC(EGF(1),EGF(2),EGF(3),EGF(4),EGF(5),EGF(6),EGF(7),
     1    EGF(8),EGF(9),SF)

        RETURN
      END
C
C***********************************************************************
C
C     计算基体的柔度矩阵 注意 返回为传统顺序
C
      SUBROUTINE KMATRIX(RESIN,EGM,SM,SSM,ETM,MSEG,STATEV13,RF,
     1SUUM,NOEL,KINC)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION EGM(9),SUUM(3),RESIN(12),SM(6,6),SSM(6),ETM(2,20),
     1SP(6,6)

      EGM(1)=RESIN(1)*RF**STATEV13
      EGM(2)=RESIN(2)*RF**STATEV13
      EGM(3)=RESIN(3)*RF**STATEV13
      EGM(4)=RESIN(4)
      EGM(5)=RESIN(5)
      EGM(6)=RESIN(6)
      EGM(7)=RESIN(7)*RF**STATEV13
      EGM(8)=RESIN(8)*RF**STATEV13
      EGM(9)=RESIN(9)*RF**STATEV13
      SUUM(1)=RESIN(10)
      SUUM(2)=RESIN(11)
      SUUM(3)=RESIN(12)

      CALL KELASC(EGM(1),EGM(2),EGM(3),EGM(4),EGM(5),EGM(6),EGM(7),
     1EGM(8),EGM(9),SM)

      CALL KPLASC(SSM,ETM,MSEG,SP,EGM,ID,STATEV13,RF,NOEL)
      IF(ID.EQ.0)RETURN
      IF(STATEV13.EQ.1.)THEN
      DO 20 I=1,6
      DO 20 J=1,6
  20  SM(I,J)=SM(I,J)
      ENDIF
      IF(STATEV13.EQ.0.)THEN
      DO 10 I=1,6
      DO 10 J=1,6
  10  SM(I,J)=SM(I,J)+SP(I,J)
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
C     计算弹性柔度矩阵 注意 返回为传统顺序
C
      SUBROUTINE KELASC(E11,E22,E33,V23,V13,V12,G23,G13,G12,SE)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SE(6,6)
        
        INTEGER I,J
        DO I=1,6
          DO J=1,6
            SE(I,J)=0.
          END DO
        END DO
        SE(1,1)=1/E11
        SE(2,2)=1/E22
        SE(3,3)=1/E33
        SE(1,2)=-V12/E11
        SE(1,3)=-V13/E11
        SE(2,3)=-V23/E22
        SE(4,4)=1/G23
        SE(5,5)=1/G13
        SE(6,6)=1/G12
        DO I=1,3
          DO J=I,3
            SE(J,I)=SE(I,J)
          END DO
        END DO
        RETURN
      END
C
C***********************************************************************
C
C     计算柔度矩阵的塑性分量 注意 返回为传统顺序
C	输入:
C			SSM: 基体的平均应力
C			ETM: 基体的分段性能
C     MSEG: ETM分段数 
C			RF: 衰减系数
C 输出:
C			SP: 塑性分量
C			EGM: 基体塑性阶段的等效性能
C			ID:	是否达到塑性, 是为1 否为0
      SUBROUTINE KPLASC(SSM,ETM,MSEG,SP,EGM,ID,STATEV13,RF,NOEL)

        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SA(6),SSM(6),ETM(2,20),SP(6,6),EGM(9)

        SA(1)=SSM(1)-(SSM(1)+SSM(2)+SSM(3))/3
        SA(2)=SSM(2)-(SSM(1)+SSM(2)+SSM(3))/3
        SA(3)=SSM(3)-(SSM(1)+SSM(2)+SSM(3))/3
        SA(4)=SSM(4)
        SA(5)=SSM(5)
        SA(6)=SSM(6)

        ID=0
        S1=0.0
        DO I=1,3
          S1=S1+SA(I)**2+2*SA(I+3)**2
        END DO
        S1=SQRT(S1/6)*3
      
        IF(S1.LE.ETM(1,1)) RETURN
        ID=1
        DO 10 I=1,MSEG-1
            IF(S1.GT.ETM(1,I).AND.S1.LE.ETM(1,I+1))THEN
                ET=ETM(2,I+1)
                GOTO 20
            ENDIF
10    	CONTINUE
        ET=ETM(2,MSEG)
  
20    	E=ETM(2,1)
        C=9*(E-ET)/(4*E*ET*S1*S1)
        ! 防止溢出造成不对称
        DO I=1,3
          DO J=1,3
            SP(I,J)=C*(SA(I)*SA(J))
            SP(I,J+3)=2*C*(SA(I)*SA(J+3))
            SP(I+3,J)=2*C*(SA(J)*SA(I+3))
            SP(I+3,J+3)=4*C*(SA(I+3)*SA(J+3))
          END DO
        END DO

        ! 基体属性
        ISFAIL=STATEV13
        IF(ISFAIL.EQ.0) THEN
          EGM(1)=ET
          EGM(2)=ET
          EGM(3)=ET
          EGM(4)=0.5
          EGM(5)=0.5
          EGM(6)=0.5
          EGM(7)=ET/3
          EGM(8)=ET/3
          EGM(9)=ET/3
        ELSE
          EGM(1)=ET*RF
          EGM(2)=ET*RF
          EGM(3)=ET*RF
          EGM(4)=0.5
          EGM(5)=0.5
          EGM(6)=0.5
          EGM(7)=ET/(2*(1+EGM(4)))*RF
          EGM(8)=EGM(7)
          EGM(9)=EGM(7)
        END IF

        RETURN
      END SUBROUTINE
C
C***********************************************************************
C
C 利用 (SM-SF)*(VM*A+VF*I)^{-1} 是对称矩阵且A是对角线已知的上三角矩阵来计算A
C DIAG_A(6)是A对角线上的六个元素
C 输入参数：
C     VF: 纤维体积分数
C     SF: 纤维柔度矩阵
C     SM: 基体柔度矩阵
C		  EGM: 基体的等效材料参数
C		  EGF: 纤维的材料常数
C		  ALFA: 桥联参数ALPHA
C		  BETA: 桥联参数BETA
C 返回参数:
C     A:    桥联矩阵A
C     B:    桥联矩阵B
C     SL:  向复合材料的柔度矩阵
C
      SUBROUTINE KBRIGE(VF,EGF,SF,EGM,SM,ALFA,BETA,A,B,SL,NOEL)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION EGF(9),SF(6,6),EGM(9),SM(6,6)
        DIMENSION A(6,6),B(6,6),SL(6,6)

        !临时变量 
        DIMENSION DIAG_A(6),S(6,6),F(15,15),X(15),Y(15),TMP(6,6)
        INTEGER I,J,K,N
        !REAL*8 VM,EM,E11F,E22F,GM,G12F
        
        E11F=EGF(1)
        E22F=EGF(2)
        G12F=EGF(9) 
        EM=EGM(1)
        GM=EGM(7)

        DIAG_A(1)=EM/E11F
        DIAG_A(2)=BETA+(1-BETA)*EM/E22F
        DIAG_A(3)=DIAG_A(2)
        DIAG_A(4)=DIAG_A(2)
        DIAG_A(5)=ALFA+(1-ALFA)*GM/G12F
        DIAG_A(6)=DIAG_A(5)          
        
        !PRINT *,'VF',VF
        !WRITE(*,'(A,6F15.8)') 'DIAG_A:',DIAG_A
        ! 初始化临时变量
        VM=1.0-VF
        MAXS=0.0
        DO I=1,6
          DO J=1,6
            A(J,I)=0
            B(J,I)=0
            S(J,I)=SM(J,I)-SF(J,I)
          END DO
        END DO

        DO I=1,15
          DO J=1,15
            F(I,J)=0.0
          ENDDO
          Y(I)=0.0
        ENDDO
        
        DO I =1,6
          B(I,I)=1.0/(VM*DIAG_A(I)+VF)
        ENDDO

        DO I=2,6
          N=(I-1)*(I-2)/2
          DO J=1,I-1
            DO K=1,J-1
              F(N+J,(J-1)*(J-2)/2+K)=F(N+J,(J-1)*(J-2)/2+K)+S(I,K)
            ENDDO

            DO K=1,I-1
              F(N+J,N+K)=F(N+J,N+K)-S(J,K)
            ENDDO
            Y(N+J)=S(I,J)*(B(I,I)-B(J,J))
          ENDDO
        ENDDO
        
        !WRITE(*,*) 'F:'
        !WRITE(*,'(15F15.8)') F
        !WRITE(*,*) 'Y:'
        !WRITE(*,'(15F15.8)') Y
        CALL KGAUSS(F,Y,15,X,I)
        !WRITE(*,*) 'X:'
        !WRITE(*,'(15F15.8)') X

        N=1
        DO I=2,6
          DO J=1,I-1
            B(J,I)=X(N)
            N=N+1
          ENDDO
        ENDDO

        SL=MATMUL(S,B)
        DO I=1,6
          DO J=1,6
            SL(J,I)=SM(J,I)-VF*SL(J,I)
          END DO
        END DO

        ! 计算A矩阵
        DO I=1,6
          DO J=1,6
            TMP(J,I)=0.0
          ENDDO
        ENDDO
        
        CALL KINVER2(B,TMP,6,I)
        
        DO I=1,6
          A(I,I)=(TMP(I,I)-VF)/(1.0-VF)
        END DO
        DO I=2,6
          DO J=1,I-1
            A(J,I)=TMP(J,I)/(1.0-VF)
          ENDDO
        ENDDO

        !WRITE(*,*) 'SM:'
        !WRITE(*,'(6F15.8,2x)') ((SM(I,J),J=1,6),I=1,6)  
        !WRITE(*,*) 'SF:'
        !WRITE(*,'(6F15.8,2x)') ((SF(I,J),J=1,6),I=1,6)  
        !WRITE(*,*) 'A:'
        !WRITE(*,'(6F15.8,2x)') ((B(I,J),J=1,6),I=1,6)  
        !WRITE(*,*) 'A:'
        !WRITE(*,'(6F15.8,2x)') ((A(I,J),J=1,6),I=1,6)  
        !WRITE(*,*) 'SL:'
        !WRITE(*,'(6F20.15,2x)') ((SL(I,J),J=1,6),I=1,6)  
        RETURN
      END SUBROUTINE

C
C*** SUBROUTINE TO FIND THE INVERSION S OF THE POSITIVE DEFINITE MATRIX C
C***  N<7, OTHERWISE, A & B SHOULD BE MODIFIED
C
      SUBROUTINE KINVER2(D,A,N,L)
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION D(N,N),ME(6),B(6),C(6),A(N,N)
      REAL*8 T
      T=1.0
      DO 5 I=1,N
          DO 5 J=1,N
              A(I,J)=D(I,J)
  5   CONTINUE
      L=1
      DO 10 J=1,N
  10  ME(J)=J
      DO 20 I=1,N
          Y=0.
          DO 30 J=I,N
              IF(ABS(A(I,J)).LE.ABS(Y))GOTO 30
              K=J
              Y=A(I,J)
  30      CONTINUE
          IF((ABS(Y)+1.).EQ.1.)THEN
              L=0
              STOP
          ENDIF
          Y=1./Y
          DO 40 J=1,N
              C(J)=A(J,K)
              A(J,K)=A(J,I)
              A(J,I)=-C(J)*Y
              B(J)=A(I,J)*Y
  40      A(I,J)=A(I,J)*Y
          A(I,I)=Y
          J=ME(I)
          ME(I)=ME(K)
          ME(K)=J
          DO 11 K=1,N
              IF(K.EQ.I)GOTO 11
              DO 12 J=1,N
                  IF(J.EQ.I)GOTO 12
                  A(K,J)=A(K,J)-B(J)*C(K)
  12          CONTINUE
  11      CONTINUE
  20  CONTINUE
      DO 33 I=1,N
          DO 44 K=1,N
              IF (ME(K).EQ.I) GOTO 55
  44      CONTINUE
  55      IF(K.EQ.I)GOTO 33
          DO 66 J=1,N
              W=A(I,J)
              A(I,J)=A(K,J)
   66     A(K,J)=W
          IW=ME(I)
          ME(I)=ME(K)
          ME(K)=IW
  33  CONTINUE
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE KGAUSS1(A,B,N,X,L)
C
      INCLUDE 'ABA_PARAM.INC'
C
C*** SUBROUTINE TO SOLVE [A]{X}={B}, SOLUTION IS BACKED IN {X}
C
      DIMENSION A(N,N),X(N),B(N),JS(N)
C
      L=1
      DO 50 K=1,N-1
      D=0.
      DO 210 I=K,N
      DO 210 J=K,N
      IF(ABS(A(I,J)).GT.D) THEN
      D=ABS(A(I,J))
      JS(K)=J
      IS=I
      ENDIF
 210  CONTINUE
      IF((D+1.0).EQ.1.) THEN
      L=0
      ELSE
      IF(JS(K).NE.K) THEN
      DO 220 I=1,N
      T=A(I,K)
      A(I,K)=A(I,JS(K))
      A(I,JS(K))=T
 220  CONTINUE
      ENDIF
      IF(IS.NE.K) THEN
      DO 230 J=K,N
      T=A(K,J)
      A(K,J)=A(IS,J)
      A(IS,J)=T
 230  CONTINUE
      T=B(K)
      B(K)=B(IS)
      B(IS)=T
      ENDIF
      ENDIF
      IF(L.EQ.0) THEN
      WRITE(*,100)
      RETURN
      ENDIF
      DO 10 J=K+1,N
      A(K,J)=A(K,J)/A(K,K)
 10   CONTINUE
      B(K)=B(K)/A(K,K)
      DO 30 I=K+1,N
      DO 20 J=K+1,N
      A(I,J)=A(I,J)-A(I,K)*A(K,J)
 20   CONTINUE
      B(I)=B(I)-A(I,K)*B(K)
 30   CONTINUE
 50   CONTINUE
      IF(ABS(A(N,N))+1.EQ.1.) THEN
      L=0
      WRITE(*,100)
      RETURN
      ENDIF
      X(N)=B(N)/A(N,N)
      DO 70 I=N-1,1,-1
      T=0.
      DO 60 J=I+1,N
      T=T+A(I,J)*X(J)
 60   CONTINUE
      X(I)=B(I)-T
 70   CONTINUE
 100  FORMAT(/10X,'LINEAR FAIL')
      JS(N)=N
      DO 150 K=N,1,-1
      IF(JS(K).NE.K) THEN
      T=X(K)
      X(K)=X(JS(K))
      X(JS(K))=T
      ENDIF
 150  CONTINUE
      RETURN
      END

C 计算 [A]{X}={B}
C **注意输入A,B都会被修改 所以调用KGAUSS的参数不能直接传入，最好用临时数组传入**
C 输入:
C   A: 
C   B: 
C   N: 矩阵维度
C 输出:
C   X:  
C   L: 本次求解的状态信息 0 表示正确求解 1表示出现溢出错误 2表示解不唯一
      SUBROUTINE KGAUSS(A, B, N, X, L)
        INCLUDE 'ABA_PARAM.INC'
        
        INTEGER N, L
        DIMENSION A(N, N), X(N), B(N)

        INTEGER I, J,K, I_MAX
        REAL*8 V_MAX,T
        
        DO K=1,N-1
          I_MAX=K
          V_MAX=ABS(A(K,K))
          ! 寻找合适主元
          DO J=K+1,N
            IF(ABS(A(J,K)).GT.V_MAX) THEN
              V_MAX=A(J,K)
            END IF
          END DO
          
          ! 交换行
          IF(I_MAX.NE.K) THEN
            DO J=K,N
              T=A(K,J)
              A(K,J)=A(I_MAX,J)
              A(I_MAX,J)=T
            END DO
            T=B(K)
            B(K)=B(I_MAX)
            B(I_MAX)=T
          END IF
          
          ! 高斯消去 上三角化
          IF (V_MAX+1.0.NE.1.0) THEN
            DO I=K+1,N
              T=A(I,K)/A(K,K)
              DO J=K,N
                  A(I,J)=A(I,J)-A(K,J)*T
              END DO
              B(I)=B(I)-B(K)*T
            END DO
          END IF
        END DO
        ! 反向替换求解上三角形式方程
        L=0
        DO K=N,1,-1
          IF(A(K,K)+1.0.EQ.1.0) THEN
            IF(B(K)+1.0.EQ.1.0) THEN
              ! make NAN to 0
              X(K)=0.0
              L=2
            ELSE
              ! make INFTY error
              WRITE(*,*) 'KGAUSS ERROR : INFINITY',A(K,K),B(K)
              L=1
              STOP
            END IF
          ELSE
            X(K)=B(K)/A(K,K)
          END IF
          
          DO I=1,K-1
            B(I)=B(I)-A(I,K)*X(K)
          END DO
        END DO
        RETURN
      END SUBROUTINE
C
C***********************************************************************
C
C     计算单元的刚度矩阵并更新单元的应力
      SUBROUTINE KSTRESS(VF,SF,SM,A,B,SL,DDSDDE,STRESS,DSTRAN,STATEV,
     1NSTATV,RF,NOEL)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SF(6,6),SM(6,6),A(6,6),B(6,6),SL(6,6),DSTRAN(6)
         
        DIMENSION DDSDDE(6,6),STRESS(6),STATEV(NSTATV)
      
        DIMENSION SK(6,6),TMPD(6)
        
        CALL KINVER2(SL,SK,6,L)
        DO I=1,6
          DO J=1,6
            DDSDDE(I,J)=SK(I,J)*(RF**STATEV(14))
          END DO
        END DO

        DO I=1,6
          TMPD(I)=DSTRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD)
        DO I=1,6
          DO J=1,6
            STRESS(I)=STRESS(I)*(RF**STATEV(14))+DDSDDE(I,J)*TMPD(J)
          END DO
        END DO
        RETURN
      END SUBROUTINE
C
C***********************************************************************
C
! 计算基体和纤维中内应力SSF和SSM
      SUBROUTINE KINTER1(DDSDDE,DSTRAN,DSTRESS,A,B,VF,SSF,SSM,SCFS,LMISE,
     1 STATEV,NSTATV,NOEL)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION DSTRAN(6),A(6,6),B(6,6),SCFS(7),STATEV(NSTATV)
        
        DIMENSION DDSDDE(6,6),DSTRESS(6),SSF(6),SSM(6)

        DIMENSION C(6,6),DSSF(6),DSSM(6),TMPD(6),SS(3),VECS(3,3)
        REAL*8 LMISE,I1,I2,I3,P,Q,T,THETA,MISE
        REAL*8 K12,K23,K22,KT22,KTT22,KC22,KT22BI,KC22BI
        
        ! 应力集中系数
        KT22=SCFS(1)
        KTT22=SCFS(2)
        KT22BI=SCFS(3)
        KC22BI=SCFS(4)
        KC22=SCFS(5)
        K23=SCFS(6)
        K12=SCFS(7)

        !单元应力增量
        DO I=1,6
          TMPD(I)=DSTRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD)
        DO I=1,6
            DSTRESS(I)=0
        END DO
        DO I=1,6
          DO J=1,6
            DSTRESS(I)=DSTRESS(I)+DDSDDE(I,J)*TMPD(J)
          END DO
        END DO
        
        !纤维内应力SSF
        DO I=1,6
          SUM=0.
          DO  J=1,6
            SUM=SUM+B(I,J)*DSTRESS(J)
          END DO
          DSSF(I)=SUM
          SSF(I)=SSF(I)+DSSF(I)
        END DO
 
        !基体内应力SSM
        DO  I=1,6
          SUM=0.
          DO J=1,6
            SUM=SUM+A(I,J)*DSSF(J)
          END DO
          DSSM(I)=SUM
        END DO

        !C 计算基体的主应力S1：
        CALL KPRINC1D(SSM,SIGMA1,S2,S3)
        !C 计算基体的Mises应力MISE：
        MISE=SQRT(0.5*((SSM(1)-SSM(2))**2+(SSM(1)-SSM(3))**2
     1  +(SSM(2)-SSM(3))**2+6*(SSM(4)**2+SSM(5)**2+SSM(6)**2)))
        IF(DSSM(2).GT.0.AND.MISE.GT.LMISE.AND.SIGMA1.GT.0)THEN
            K22=KTT22
        ELSE IF(DSSM(2).LT.0)THEN
            K22=KC22
        ELSE
            K22=KT22
        ENDIF
        IF(DSSM(3).GT.0.AND.MISE.GT.LMISE.AND.SIGMA1.GT.0)THEN
            K33=KTT22
        ELSE IF(DSSM(2).LT.0)THEN
            K33=KC22
        ELSE
            K33=KT22
        ENDIF
        SSM(1)=SSM(1)+DSSM(1)
        SSM(2)=SSM(2)+DSSM(2)*K22
        SSM(3)=SSM(3)+DSSM(3)*K33
        SSM(4)=SSM(4)+DSSM(4)*K23
        SSM(5)=SSM(5)+DSSM(5)*K12
        SSM(6)=SSM(6)+DSSM(6)*K12
        !
        DO I=1,6
          STATEV(I)=SSF(I)
          STATEV(I+6)=SSM(I)
        END DO
        RETURN
      END
C     计算纤维和基体的内应力
      SUBROUTINE KINTER(DDSDDE,DSTRAN,DSTRESS,A,B,VF,SSF,SSM,SCFS,LMISE,
     1 STATEV,NSTATV,NOEL)
        INCLUDE 'ABA_PARAM.INC'
        DIMENSION DSTRAN(6),A(6,6),B(6,6),SCFS(7),STATEV(NSTATV)
        
        DIMENSION DDSDDE(6,6),DSTRESS(6),SSF(6),SSM(6)

        DIMENSION C(6,6),DSSF(6),DSSM(6),TMPD(6),SS(3),VECS(3,3)
        REAL*8 LMISE,I1,I2,I3,P,Q,T,THETA,SIGMA1,MISE
        REAL*8 K12,K23,K22,KT22,KTT22,KC22,KT22BI,KC22BI
        
        ! 应力集中系数
        KT22=SCFS(1)
        KTT22=SCFS(2)
        KT22BI=SCFS(3)
        KC22BI=SCFS(4)
        KC22=SCFS(5)
        K23=SCFS(6)
        K12=SCFS(7)

        !单元应力增量
        DO I=1,6
          TMPD(I)=DSTRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD)
        DO I=1,6
            DSTRESS(I)=0
        END DO
        DO I=1,6
          DO J=1,6
            DSTRESS(I)=DSTRESS(I)+DDSDDE(I,J)*TMPD(J)
          END DO
        END DO
        
        !纤维内应力SSF
        DO I=1,6
          SUM=0.
          DO  J=1,6
            SUM=SUM+B(I,J)*DSTRESS(J)
          END DO
          DSSF(I)=SUM
          SSF(I)=SSF(I)+DSSF(I)
        END DO
 
        !基体内应力SSM
        DO  I=1,6
          SUM=0.
          DO J=1,6
            SUM=SUM+A(I,J)*DSSF(J)
          END DO
          DSSM(I)=SUM
        END DO

        !C 计算基体的主应力S1：
        CALL KPRINC1D(SSM,S1,S2,S3)
        !C 计算基体的Mises应力MISE：
        MISE=SQRT(0.5*((SSM(1)-SSM(2))**2+(SSM(1)-SSM(3))**2
     1  +(SSM(2)-SSM(3))**2+6*(SSM(4)**2+SSM(5)**2+SSM(6)**2)))

        ! 更新基体的内应力SSM(6):
        ! 将基体的应力状态分成3种情况:
        ! Y=1表示基体单轴拉伸，
        ! Y=2表示基体双轴拉伸且2方向的应力大于3方向的应力，
        ! Y=3表示基体双轴拉伸且3方向的应力大于2方向的应力)
        IF(DSSM(3).LT.0.01)THEN
            Y=1
            !C 判断2方向的应力集中系数K22:
            IF(DSSM(2).GT.0.AND.MISE.GT.LMISE)THEN
                K22=KTT22
            ELSE IF(DSSM(2).LT.0)THEN
                K22=KC22
            ELSE
                K22=KT22
            ENDIF
            !C 判断3方向的应力集中系数K33:
            IF(DSSM(3).GT.0.AND.MISE.GT.LMISE)THEN
                K33=KTT22
            ELSE IF(DSSM(3).LT.0)THEN
                K33=KC22
            ELSE
                K33=KT22
            ENDIF

            !C 更新基体的内应力SSM(6):
            SSM(1)=SSM(1)+DSSM(1)
            SSM(2)=SSM(2)+DSSM(2)*K22
            SSM(3)=SSM(3)+DSSM(3)*K33
            SSM(4)=SSM(4)+DSSM(4)*K23
            SSM(5)=SSM(5)+DSSM(5)*K12
            SSM(6)=SSM(6)+DSSM(6)*K12
        ELSE IF(DSSM(2).GE.DSSM(3))THEN
            Y=2
            !C 判断单轴的应力集中系数K22：
            IF(S1.GT.0.AND.MISE.GT.LMISE)THEN
                K22=KTT22
            ELSE
                K22=KT22
            ENDIF
            !C 判断双轴应力集中系数KBI：
            IF(DSSM(3).GT.0)THEN
                KBI=KT22BI
            ELSE
                KBI=KC22BI
            ENDIF
            !C 更新基体的内应力SSM(6):
            SSM(1)=SSM(1)+DSSM(1)
            SSM(2)=SSM(2)+DSSM(3)*KBI+K22*(DSSM(2)-DSSM(3))
            SSM(3)=SSM(3)+DSSM(3)*KBI
            SSM(4)=SSM(4)+DSSM(4)*K23
            SSM(5)=SSM(5)+DSSM(5)*K12
            SSM(6)=SSM(6)+DSSM(6)*K12
        ELSE
            Y=3
            !C 判断单轴的应力集中系数K22：
            IF(S1.GT.0.AND.MISE.GT.LMISE)THEN
                K22=KTT22
            ELSE
                K22=KT22
            ENDIF
            !C 判断双轴应力集中系数KBI：
            IF(DSSM(2).GT.0)THEN
                KBI=KT22BI
            ELSE
                KBI=KC22BI
            ENDIF
            !C 更新基体的内应力SSM(6):
            SSM(1)=SSM(1)+DSSM(1)
            SSM(2)=SSM(2)+DSSM(2)*KBI
            SSM(3)=SSM(3)+DSSM(2)*KBI+K22*(DSSM(3)-DSSM(2))
            SSM(4)=SSM(4)+DSSM(4)*K23
            SSM(5)=SSM(5)+DSSM(5)*K12
            SSM(6)=SSM(6)+DSSM(6)*K12
        ENDIF

        DO  I=1,6
          STATEV(I)=SSF(I)
          STATEV(I+6)=SSM(I)
        END DO
        RETURN
      END
C
C***********************************************************************
C      
      SUBROUTINE KSTRTH(SSF,SSM,SUUM,SUUF,STATEV,NSTATV,FAILCHS,NOEL,
     1KINC)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SSF(6),SSM(6),SUUM(3),SUUF(2),STATEV(NSTATV)
        REAL*8 FAILCHS
        INTEGER ID,FAILTYPE

        ! 纤维破坏判据
        CALL KPRINC1D(SSF,S1,S2,S3)
        IF(S3.LT.0)THEN
            SE=S1
        ELSE 
            SE=(S1**3+S2**3+S3**3)**(1./3.)
        ENDIF
        IF(SE.GT.SUUF(1))THEN
          STATEV(14)=1
          WRITE(103,*) 'FIBER IS BROKEN',NOEL,KINC
        ENDIF

        ! 基体破坏判据
        IF(FAILCHS.EQ.1)THEN
          ID=0
          CALL KTSAIWU(SSM,SUUM,ID,FAILTYPE,KINC,NOEL)
          IF(ID.EQ.1)THEN
              IF(STATEV(13).NE.1) THEN
                WRITE(102,*) 'MATRIX IS BROKEN',NOEL,KINC
              END IF
              STATEV(13)=1
          ENDIF
        ELSE IF(FAILCHS.EQ.2)THEN
          ID=0
          CALL KMOHR(SSM,SUUM,ID,FAILTYPE,ANZ)
          IF(ID.EQ.1)THEN
            IF(STATEV(1).NE.1) THEN
              WRITE(102,*) 'MATRIX IS BROKEN',ANZ,NOEL,KINC
            END IF
              STATEV(13)=1
          ENDIF
        ENDIF

        RETURN
      END
C
C***********************************************************************
C
C     TSAI-WU判据
      SUBROUTINE KTSAIWU(SSM,SUUM,ISFAIL,FAILTYPE,KINC,NOEL)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SSM(6),SUUM(3)
 
        INTEGER ISFAIL,FAILTYPE,KINC,NOEL

        REAL*8 F1,F11,F66,MISE,DELAMI
  
        F1=(SUUM(2)-SUUM(1))/(SUUM(1)*SUUM(2))
        F11=1/(SUUM(1)*SUUM(2))
        F66=1/(SUUM(3)**2)
        !C
        !C      MISE=F1*(SSM(1)+SSM(2)+SSM(3))+F11*(SSM(1)**2+SSM(2)**2+SSM(3)**2-
        !C     1SSM(1)*SSM(2)-SSM(2)*SSM(3)-SSM(3)*SSM(1))+F66*(SSM(4)**2+
        !C     2SSM(5)**2+SSM(6)**2)
        !C
        MISE=F1*(SSM(1)+SSM(2))+F66*SSM(6)**2
        MISE=MISE+F11*(SSM(1)**2+SSM(2)**2-SSM(1)*SSM(2))

        DELAMI=F11*SSM(3)**2+F66*(SSM(5)**2+SSM(4)**2)+F1*SSM(3)

        ISFAIL=0
        FAILTYPE=0
        IF(MISE.GT.1)THEN
            ISFAIL=1
            WRITE(102,*)'MATRIX IS BROKEN',NOEL,KINC
            FAILTYPE=FAILTYPE+1
        ENDIF
        IF(DELAMI.GT.1)THEN
            ISFAIL=1
            FAILTYPE=FAILTYPE+2
            WRITE(102,*)'DELAMINATION',NOEL,KINC
        ENDIF
        RETURN
      END
C
C***********************************************************************
C
C 利用Mohr判据判断破坏
C 如果破坏 返回ID=1  否则 ID=0
C 输入参数:
C     SIGMA: 六个应力分量数组 sigma11 sigma22 sigma33 sigma23 sigma13 sigma12
C     USIGMA: 三个强度分量  ut uc us
C 输出参数:
C     ISFAIL: 整型  是否破坏，破坏为1 未破坏为0
C     FAILTYPE: 整型 破坏类型， 0表示拉伸破坏 1 表示压缩破坏 2表示剪切破坏
C     VEC:    浮点型  破环面法向
      SUBROUTINE KMOHR(SIGMA,USIGMA,ISFAIL,FAILTYPE,VEC)        
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SIGMA(6),USIGMA(3),VEC(3)
        INTEGER ISFAIL,FAILTYPE

        REAL*8 S(3,3),PRINCS(3),VECS(3,3)
        REAL*8 DELTA,ST,SC,SS,AT,AC,BT,BC,A,B,T,SN,R,THETA
        INTEGER I,J

        ST=USIGMA(1)
        SC=USIGMA(2)
        SS=USIGMA(3)

        CALL KPRIND1D(SIGMA,PRINCS,VECS)
        !WRITE(*,'(3F15.8)') (PRINCS(I),I=1,3)
        !WRITE(*,'(3F15.8)')((VECS(I,J),J=1,3),I=1,3)
        SN=(PRINCS(1)+PRINCS(3))/2
        R=ABS(PRINCS(1)-PRINCS(3))/2

        DELTA=0.0
        ISFAIL=0
        IF(ABS(SN).LT.1E-6) THEN
          !剪切
          DELTA=SS/R
          FAILTYPE=2
        ELSE 
          AT=2*ST/(ST**2-4*SS**2);
          BT=-AT*SS**2-1/(4*AT);
          
          AC=2*SC/(4*SS**2-SC**2);
          BC=-AC*SS**2-1/(4*AC);

          A=AC
          B=BC
          FAILTYPE=1
          IF(SN.GT.0.0) THEN
            A=AT
            B=BT
            FAILTYPE=0
          ENDIF
          DELTA=(SN/A+SQRT((SN/A)**2+4*(R**2)*(SS**2)))/(2*R**2)
          THETA=ATAN2(SQRT((DELTA*R)**2-1./(4*A**2)),-0.5/A)/2.0
          !THETA=ACOS(1.0/(-2*A*R*DELTA))/2
          
          ! 当应力圆与抛物线相切于外部时，应当考虑应力圆于拉压应力圆的相切
          IF (SN.GT.0.0) THEN
            IF((DELTA*SN).GT.(ST/2.0)) THEN
              DELTA=ST/(SN+R)
              THETA=0.0
            ENDIF
          ELSE
            IF((DELTA*ABS(SN)).GT.(SC/2.0)) THEN
              DELTA=SC/(ABS(SN)+R)
              THETA=3.141592653/2.0
            ENDIF
          ENDIF

          DO I=1,3
            VEC(I)=COS(THETA)*VECS(I,1)+SIN(THETA)*VECS(I,3)
          ENDDO
        ENDIF

        ISFAIL=0
        IF((DELTA.GE.0).AND.(DELTA.LE.1.0)) THEN
          ISFAIL=1
        ENDIF
      END SUBROUTINE    
C
C***********************************************************************
C
C     从传统顺序的应力向量(6维)求解主应力
C     注意: 应力分量时传统顺序 不是ABAQUS中的顺序
C 输入
C   SS: 传统的应力分量形式 11 22 33 23 13 12 注意不是ABAQUS中的顺序
C 输出:
C   S1,S2,S3: 从大到小的主应力
      SUBROUTINE KPRINC1D(SS,S1,S2,S3)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SS(6)

        DIMENSION CONV_SS(6),PRINCS(3)
        DO I=1,6
          CONV_SS(I)=SS(I)
        END DO
        CALL NORM2ABQ_1D(CONV_SS)
        
        CALL SPRINC(CONV_SS,PRINCS,1,3,3)
        S1=PRINCS(1)
        S2=PRINCS(2)
        S3=PRINCS(3)

        RETURN 
      END SUBROUTINE
C       调用SPRINC 求解3x3应力矩阵的特征值和特征向量
C 输入:
C   S: 3x3应力矩阵
C 输出:
C   S1,S2,S3: 从大到小的主应力
      SUBROUTINE KPRINC2D(S,S1,S2,S3)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION S(3,3)        
        DIMENSION STRESS(6),PRINCS(3)

        DO I=1,3
          STRESS(I)=S(I,I)
        END DO
        STRESS(4)=S(2,3)
        STRESS(5)=S(1,3)
        STRESS(6)=S(1,2)

        CALL KPRINC1D(STRESS,S1,S2,S3)
      END SUBROUTINE
C
C***********************************************************************
C
C     调用SPRINC 求解3x3对称矩阵的特征值和特征向量
C 输入:
C   S: 6维向量
C 输出:
C   PRINCS: 三个特征值(从大到小排序)
C   VECS: 列为LBDS对应的单位特征向量
      SUBROUTINE KPRIND1D(S,PRINCS,VECS)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION S(6)
        DIMENSION PRINCS(3),VECS(3,3)

        DIMENSION ABQ_SS(6),PV(3)
        
        DO I=1,6
          ABQ_SS(I)=S(I)
        END DO
        CALL NORM2ABQ_1D(ABQ_SS)
        CALL SPRIND(ABQ_SS,PRINCS,VECS,1,3,3)
        !冒泡排序
        DO I=2,3
          PIVOT=PRINCS(I)
          PV(1)=VECS(1,I)
          PV(2)=VECS(2,I)
          PV(3)=VECS(3,I)
          J=I-1
          DO WHILE(J.GE.1)
            IF(PRINCS(J).GE.PIVOT) THEN
              EXIT
            END IF
            PRINCS(J+1)=PRINCS(J)
            VECS(1,J+1)=VECS(1,J)
            VECS(2,J+1)=VECS(2,J)
            VECS(3,J+1)=VECS(3,J)
            J=J-1
          ENDDO
          PRINCS(J+1)=PIVOT
          VECS(1,J+1)=PV(1)
          VECS(2,J+1)=PV(2)
          VECS(3,J+1)=PV(3)
        ENDDO
      END SUBROUTINE
C
C***********************************************************************
C
C     调用SPRIND 求解3x3对称矩阵的特征值和特征向量
C 输入:
C   S: 3x3应力矩阵
C 输出:
C   PRINCS: 三个特征值(从大到小排序)
C   VECS: 列为LBDS对应的单位特征向量
      SUBROUTINE KPRIND2D(S,PRINCS,VECS)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION S(3,3)
        DIMENSION PRINCS(3),VECS(3,3)

        DIMENSION STRESS(6),STRAIN(6)

        DO I=1,3
          STRESS(I)=S(I,I)
        END DO
        STRESS(4)=S(2,3)
        STRESS(5)=S(1,3)
        STRESS(6)=S(1,2)

        CALL KPRIND1D(STRESS,PRINCS,VECS)
      END SUBROUTINE
C
C***********************************************************************
C
C     应力应变向量从一般顺序11 22 33 23 13 12转换为Abaqus的顺序 11 22 33 12 13 23
C
      SUBROUTINE NORM2ABQ_1D(VEC)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION VEC(6)

        TMP=VEC(4)
        VEC(4)=VEC(6)
        VEC(6)=TMP
      END SUBROUTINE
C
C***********************************************************************
C
C     应力应变向量 从Abaqus的顺序 11 22 33 12 13 23转换为一般顺序11 22 33 23 13 12
C
      SUBROUTINE ABQ2NORM_1D(VEC)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION VEC(6)

        TMP=VEC(4)
        VEC(4)=VEC(6)
        VEC(6)=TMP
      END SUBROUTINE      
C
C***********************************************************************
C
C     柔度矩阵和刚度矩阵 从一般顺序转换为Abaqus的顺序
C
      SUBROUTINE NORM2ABQ_2D(DM)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION DM(6,6)

        DO I=1,6
          TMP=DM(4,I)
          DM(4,I)=DM(6,I)
          DM(6,I)=TMP
        END DO

        DO I=1,6
          TMP=DM(I,4)
          DM(I,4)=DM(I,6)
          DM(I,6)=TMP
        END DO
      END SUBROUTINE
C
C***********************************************************************
C
C     柔度矩阵和刚度矩阵 从Abaqus的顺序转换为一般顺序
C
      SUBROUTINE ABQ2NORM_2D(DM)
        INCLUDE 'ABA_PARAM.INC'
        
        DIMENSION DM(6,6)

        DO I=1,6
          TMP=DM(4,I)
          DM(4,I)=DM(6,I)
          DM(6,I)=TMP
        END DO

        DO I=1,6
          TMP=DM(I,4)
          DM(I,4)=DM(I,6)
          DM(I,6)=TMP
        END DO
      END SUBROUTINE
C
C***********************************************************************
C
      SUBROUTINE UMAT_MAT2(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
      
      DIMENSION RESIN(12),EGM(9),SUUM(3),ETM(2,20),SSM(6)
      DIMENSION SM(6,6),SK(6,6),DSTRESS(6),TMPD(6)

      ! Initialize form PROPS
      I_P=1
      E=PROPS(I_P)
      V=PROPS(I_P+1)
      DO I=1,3
        RESIN(I)=E
        RESIN(3+I)=V
        RESIN(6+I)=E/(2+2*V)
        RESIN(9+I)=PROPS(I+3)
      END DO
      I_P=I_P+5
      
      NMSEG=PROPS(I_P)
      DO I=1,NMSEG
        ETM(1,I)=PROPS(I_P+I)
        ETM(2,I)=PROPS(I_P+NMSEG+I)
      END DO
      I_P=I_P+2*NMSEG+1
      RF=PROPS(I_P)

      ! Initialize from STATEV
      SSM(1:6)=STATEV(1:6)
      CALL KMATRIX(RESIN,EGM,SM,SSM,ETM,NMSEG,STATEV(13),RF,
     1  SUUM,NOEL,KINC)
      
      CALL KINVER2(SM,SK,6,L)
      CALL ABQ2NORM_1D(STRESS)
      DO I=1,6
        TMPD(I)=DSTRAN(I)
      END DO
      CALL ABQ2NORM_1D(TMPD)
      DO I=1,6
        SUM=0D0
        DO J=1,6
          DDSDDE(I,J)=SK(I,J)
          SUM=SUM+SK(I,J)*TMPD(J)
        END DO
        !STRESS(I)=STRESS(I)+SUM !增量式
        STRESS(I)=SUM !全量式
      END DO

      STATEV(1:6)=STRESS(1:6)
      IF(STATEV(13).EQ.0)THEN
        F1=(SUUM(2)-SUUM(1))/(SUUM(1)*SUUM(2))
        F11=1/(SUUM(1)*SUUM(2))
        F44=1/(SUUM(3)**2)

        SMISE=F1*(STRESS(1)+STRESS(2)+STRESS(3))+F11*
     1  (STRESS(1)**2+STRESS(2)**2+STRESS(3)**2-STRESS(1)*STRESS(2)
     2  -STRESS(2)*STRESS(3)-STRESS(3)*STRESS(1))
     3  +F44*(STRESS(4)**2+STRESS(5)**2+STRESS(6)**2)
          
        IF(SMISE.GT.1.)THEN
          STATEV(13)=1.0D0
          WRITE(*,*) 'Resin Layer Broken: ','NOEL=',NOEL,'KINC=',KINC
        ENDIF
      ENDIF

      CALL NORM2ABQ_1D(STRESS)
      CALL NORM2ABQ_2D(DDSDDE)
      RETURN 
      END SUBROUTINE