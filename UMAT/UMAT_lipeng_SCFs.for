C                   V2.0
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
C     1., 2.,3.,4.,5.,6.,7.,8.,9., 10.,11.,12.,   #残余应力
C     #KT22  KTT22 KC22 K12 K23      更改后的不需输入SCFs，由子程序自动计算，
C     2.322,  5.083,  1.656,  1.447,1.87,1.72,1.69 ,
C      !更改后不需输入开裂临界MISES应力（LMISE），改成输入单向板横向拉伸强度YT_UD
C     LMISE RF 破坏判据的选择(1-Tsai-Wu,2-Mohr)    
C     54.4,   0.01,  1
C----------------------------------------------------------------
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      include 'aba_param.inc'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
      
C*******************
       CHARACTER*256 JOBNAME,CP,DIR
      INTEGER LJ,LD 
      
      CALL GETJOBNAME(JOBNAME, LJ )
      DIR='E:\1830900\ABAQUS_simulation\m=4,n=4\'
      LD=LEN_TRIM(DIR)
      CP(1:LD)=DIR
      CP((LD+1):(LJ+LD))=JOBNAME(1:LJ)
!C      !OPEN(101,FILE='E:\dadelamination.TXT',STATUS='UNKNOWN')
!      OPEN(102,FILE=CP(1:LJ+LD)//'-matrix.txt',STATUS='UNKNOWN')
!      OPEN(103,FILE=CP(1:LJ+LD)//'-fiber.txt',STATUS='UNKNOWN')
      OPEN(101,FILE=CP(1:LJ+LD)//'-SHUCHU.TXT',STATUS='UNKNOWN')        
C***************
      IF(CMNAME.EQ.'ULAM1'.OR.CMNAME.EQ.'ULAM2')THEN      
          CALL UMAT1(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      ELSE IF(CMNAME.EQ.'URES')THEN
          CALL UMAT2(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      ENDIF
C
      RETURN
      END
C--------------------------------------------------------------
      
C
      SUBROUTINE UMAT1(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
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

     
c      PRINT *,KINC,NOEL    !输出在log文件里
      
          CALL KADJUST(PROPS,NPROPS,STATEV,NSTATV,RESIN,FIBER,VF,ALFA,
     1 BETA,MSEG,ETM,SSF,SSM,RF,SCFS,LMISE,FAILCHS,NOEL,KINC)

      CALL KFIBER(FIBER,EGF,SUUF,SF,NOEL)
      
      CALL KMATRIX(RESIN,EGM,SM,SSM,ETM,MSEG,STATEV,NSTATV,RF,
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
     1 NSTATV,RF,NOEL,STRAN)

      CALL KINTER(DDSDDE,DSTRAN,DSTRESS,A,B,VF,SSF,SSM,SCFS,LMISE,
     1 STATEV,NSTATV,NOEL)
      
      CALL NORM2ABQ_1D(STRESS)
      CALL NORM2ABQ_2D(DDSDDE)

      IF(NOEL.EQ.1) THEN
!        WRITE(*,'(A6,I7,\)') 'KINC',KINC,'KSTEP',KSTEP,'NOEL',NOEL
!        PRINT *,''
!        !PRINT *,NSTATV
!        PRINT *,"STATEV"
!        WRITE(*,"(12F20.5)") (STATEV(I),I=1,12)
!        !PRINT *,"SCFS"
!        !WRITE(*,*) SCFS
!        PRINT *,"SF"
!        WRITE(*,'(6F20.15)') ((SF(I,J),J=1,6),I=1,6)
!        PRINT *,"SM"
!        WRITE(*,'(6F20.15)') ((SM(I,J),J=1,6),I=1,6)
!        PRINT *,"A"
!        WRITE(*,'(6F20.15)') ((A(I,J),J=1,6),I=1,6)
!        !PRINT *,"B"
!        !WRITE(*,'(6F20.15)') ((B(I,J),J=1,6),I=1,6)      
!        !PRINT *,"SL"
!        !WRITE(*,'(6F20.15)') ((SL(I,J),J=1,6),I=1,6)
!        PRINT *,"DSTRAN"
!        WRITE(*,'(6F20.15)') (DSTRAN(J),J=1,6)      
!        PRINT *,"STRESS"
!        WRITE(*,'(6F20.5)') (STRESS(J),J=1,6)   
!        PRINT *,"SSF"
!        WRITE(*,'(6F20.5)') (SSF(J),J=1,6)
!        PRINT *,"SSM"
!        WRITE(*,'(6F20.5)') (SSM(J),J=1,6)
!        PRINT *,"DDSDDE"
!        WRITE(*,'(6F20.5)') ((DDSDDE(I,J),J=1,6),I=1,6)
!!C--------------------------------------------
!        WRITE(101,*)'KINC=',KINC
!        WRITE(101,*)'SCFS0001=',SCFS   !已验证SCFS
!        write(101,*)'LMISE=',LMISE

C---------------------------------------------
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
        REAL*8 LMISE,FAILCHS,YT_UD

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
!C******************************残余应力
!        IF(KINC.EQ.1)THEN
!          DO I=1,12
!            STATEV(I)=PROPS(I+N)
!          END DO
!        ENDIF
!        N=N+12
!C*******************************
        DO I=1,6
          SSF(I)=STATEV(I)
          SSM(I)=STATEV(I+6)
        END DO
C
!C------------------------------------------------  原来手动输入的SCFs改为自动计算
!        !!KT22,KTT22,KC22,K23,K12
!        !DO I=1,5
!        !  SCFS(I)=PROPS(N+I)
!        !END DO
!        !N=N+5
!        !KT22,KTT22,KC22,K12,K23,KT22BI,KC22BI
!        KT22=PROPS(N+1)
!        KTT22=PROPS(N+2)
!        KC22=PROPS(N+3)
!        K12=PROPS(N+4)
!        K23=PROPS(N+5)
!        KT22BI=PROPS(N+6)! KT22
!        KC22BI=PROPS(N+7)! KC22
!        N=N+7
!        !KT22,KTT22,KT22BI,KC22BI,KC22,K23,K12
!        SCFS(1)=KT22
!        SCFS(2)=KTT22
!        SCFS(3)=KT22BI
!        SCFS(4)=KC22BI
!        SCFS(5)=KC22
!        SCFS(6)=K23
!        SCFS(7)=K12
C---------------------------------------------------下面是调用子程序计算SCFs
        CALL KGETSCFS(RESIN,FIBER,SCFS,VF)
C        
        LMISE=PROPS(N+1)                      ! 原来输入开裂的临界MISES应力
!!C+++++++++++++++++++++++++++++++++++++++++++++++现在改为输入单向板横向拉伸强度，程序自动计算开裂临界MISES应力
!        YT_UD=PROPS(N+1)
!        CALL KGETLMISE(RESIN,FIBER,SCFS,VF,ALFA,BETA,YT_UD,LMISE)
!!C++++++++++++++++++++++++++++++++++++++++++++++++
        RF=PROPS(N+2)
        FAILCHS=PROPS(N+3)          !FAILCHS:基体的破坏判据选取，1-Tsai-Wu，2-Mohr
        N=N+3
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
      SUBROUTINE KMATRIX(RESIN,EGM,SM,SSM,ETM,MSEG,STATEV,NSTATV,RF,
     1SUUM,NOEL,KINC)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION EGM(9),SUUM(3),RESIN(12),SM(6,6),SSM(6),ETM(2,20),
     1SP(6,6),STATEV(NSTATV)

      EGM(1)=RESIN(1)*RF**STATEV(13)
      EGM(2)=RESIN(2)*RF**STATEV(13)
      EGM(3)=RESIN(3)*RF**STATEV(13)
      EGM(4)=RESIN(4)
      EGM(5)=RESIN(5)
      EGM(6)=RESIN(6)
      EGM(7)=RESIN(7)*RF**STATEV(13)
      EGM(8)=RESIN(8)*RF**STATEV(13)
      EGM(9)=RESIN(9)*RF**STATEV(13)
      SUUM(1)=RESIN(10)
      SUUM(2)=RESIN(11)
      SUUM(3)=RESIN(12)

      CALL KELASC(EGM(1),EGM(2),EGM(3),EGM(4),EGM(5),EGM(6),EGM(7),
     1EGM(8),EGM(9),SM)

      CALL KPLASC(SSM,ETM,MSEG,SP,EGM,ID,STATEV,NSTATV,RF,NOEL)
      IF(ID.EQ.0)RETURN
      IF(STATEV(13).EQ.1.)THEN
      DO 20 I=1,6
      DO 20 J=1,6
  20  SM(I,J)=SM(I,J)
      ENDIF
      IF(STATEV(13).EQ.0.)THEN
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
      SUBROUTINE KPLASC(SSM,ETM,MSEG,SP,EGM,ID,STATEV,NSTATV,RF,NOEL)

        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SA(6),SSM(6),ETM(2,20),SP(6,6),EGM(9),STATEV(NSTATV)

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
        ISFAIL=STATEV(13)
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
          !WRITE(101,*)'KINC=',KINC
          !WRITE(101,*)((A(I,J),J=1,6),I=1,6) 
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
     1NSTATV,RF,NOEL,STRAN)
        INCLUDE 'ABA_PARAM.INC'

        DIMENSION SF(6,6),SM(6,6),A(6,6),B(6,6),SL(6,6),DSTRAN(6)
         
        DIMENSION DDSDDE(6,6),STRESS(6),STATEV(NSTATV),STRAN(6)
      
        DIMENSION SK(6,6),TMPD(6),TMPD1(6)
        
        CALL KINVER2(SL,SK,6,L)
        DO I=1,6
          DO J=1,6
            DDSDDE(I,J)=SK(I,J)*(RF**STATEV(14))
          END DO
        END DO
C     
        
c        PRINT '(6F12.4)',((DDSDDE(I,J),J=1,6),I=1,6)       !输出在log文件里
         DO I=1,6
          TMPD(I)=DSTRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD)
!C
!        DO I=1,6
!          DO J=1,6
!            STRESS(I)=STRESS(I)*(RF**STATEV(14))+DDSDDE(I,J)*TMPD(J)  !########应力折减了，会出现不收敛##########
!          END DO
!        END DO
C      下面以增量更新：
        DO I =1,6
          DO J=1,6
              STRESS(I)=STRESS(I)+DDSDDE(I,J)*TMPD(J)
          ENDDO
      ENDDO
      
C      下面以全量形式更新：+++++++++++++++++++++++++++++++++++++++++++++++++
!C       
!        
!        DO I=1,6
!          TMPD1(I)=STRAN(I)
!        END DO
!        CALL ABQ2NORM_1D(TMPD1)
!      STRESS=0.
!      DO K1=1,6
!      DO K2=1,6
!      STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*(TMPD1(K1)+TMPD(K1))   !应力全量更新
!      END DO
!      END DO
!C     这种全量更新没有遇到不收敛的问题++++++++++++++++++++++++++++++++++++
        RETURN
      END SUBROUTINE
C
C***********************************************************************
C
! 计算基体和纤维中内应力SSF和SSM
C 不过貌似这个KINTRR1（未考虑双轴应力集中系数）没有被调用，调用的KINTER，其中有双轴应力集中系数KT22BI等
      SUBROUTINE KINTER1(DDSDDE,DSTRAN,DSTRESS,A,B,VF,SSF,SSM,SCFS,LMISE
     1 ,STATEV,NSTATV,NOEL)
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
C----------------------------------------------------------以下新加的，STATEV取到20
          STATEV(17)=1
          STATEV(18)=2
          WRITE(101,*) 'FIBER IS BROKEN',NOEL,KINC
C----------------------------------------------------------
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
C----------------------------------------------------------
          STATEV(16)=1
          STATEV(18)=1
          WRITE(101,*)  'MATRIX IS BROKEN',NOEL,KINC
C----------------------------------------------------------
          ENDIF
        ELSE IF(FAILCHS.EQ.2)THEN
          ID=0
          CALL KMOHR(SSM,SUUM,ID,FAILTYPE,ANZ)
          IF(ID.EQ.1)THEN
            IF(STATEV(1).NE.1) THEN
              WRITE(102,*) 'MATRIX IS BROKEN',ANZ,NOEL,KINC
            END IF
              STATEV(13)=1
C----------------------------------------------------------
          STATEV(16)=1
          STATEV(18)=1
          WRITE(101,*)  'MATRIX IS BROKEN',ANZ,NOEL,KINC
C----------------------------------------------------------
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
C       冒泡排序       
        DO I=2,3
          PIVOT=PRINCS(I)
          J=I-1
          DO WHILE(J.GE.1)
            IF(PRINCS(J).GE.PIVOT) THEN
              EXIT
            END IF
            PRINCS(J+1)=PRINCS(J)
            J=J-1
          ENDDO
          PRINCS(J+1)=PIVOT
        ENDDO
        
        S1=PRINCS(1)
        S2=PRINCS(2)
        S3=PRINCS(3)
!C-------------------------------------------
!        WRITE(101,*)'ABQ_SS=',CONV_SS
!        WRITE(101,*)'PRINCS=',PRINCS
!C-------------------------------------------
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
      SUBROUTINE UMAT_RESIN(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
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

      I_P=1
      E=PROPS(1)
      V=PROPS(2)
      DO I=1,3
        RESIN(I)=E
        RESIN(3+I)=V
        RESIN(6+I)=E/(2+2*V)
        RESIN(9+I)=PROPS(3+I)
      END DO
      I_P=I_P+5
      
      NMSEG=PROPS(I_P)
      I_P=I_P+1

      DO I=1,NMSEG
        ETM(1,I)=PROPS(I_P+I)
        ETM(2,I)=PROPS(I_P+NMSEG+I)
      END DO
      I_P=I_P+2*NMSEG
      RF=PROPS(I_P)
      YMISE=PROPS(I_P+1)

      CALL KMATRIX(RESIN,EGM,SM,SSM,ETM,MSEG,STATEV,NSTATV,RF,
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
        STRESS(I)=STRESS(I)+SUM
      END DO

      CALL NORM2ABQ_1D(STRESS)
      CALL NORM2ABQ_2D(DDSDDE)
      RETURN 
      END SUBROUTINE
C
C***********************************************************************
C
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
C
C INITIALIZE THE OLD MAXIMUM. FOR A JOB THAT IS BEING RESTARTED
C THIS VALUE SHOULD BE SET TO THE MAXIMUM MISES STRESS IN THE
C ORIGINAL ANALYSIS.
C
      DATA OLDMAX/-1.D0/
C
      CURRMAX = 0.D0
C
C FIND CURRENT INCREMENT.
C
      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)
C
C SEARCH FOR THE HIGHEST VALUE OF MISES STRESS
C AND STORE THIS IN CURRMAX
C
      DO K1=1,999999
         CALL DBFILE(0,ARRAY,JRCD)
         IF (JRCD.NE.0) GO TO 110
         KEY=JRRAY(1,2)
         IF (KEY.EQ.12) THEN
            IF (ARRAY(3).GT.CURRMAX) CURRMAX=ARRAY(3)
         END IF
      END DO
 110  CONTINUE
C
C COMPLETED READING OF CURRENT INCREMENT. NOW CHECK TO
C SEE IF VALUE OF MISES STRESS HAS INCREASED SINCE 
C LAST INCREMENT
C
      IF (CURRMAX.LE.OLDMAX) LSTOP=1
      OLDMAX=CURRMAX
      LOVRWRT=1
C
      RETURN
      END

C******************************************************************计算应力集中系数
      SUBROUTINE  KGETSCFS(RESIN,FIBER,SCFS,VF)
      INCLUDE 'ABA_PARAM.INC'
C      implicit none
      complex::z1,z2,X,N,N1,N2,N3,M,Rei,ans1,n11,n12,n13
     &,wz,ans5,U,Z3,VZ,N4,ans3,r2z,r3z,r4z
      real*8 pusi,r,k,k1,k2,b,Vf,V,E11f,E22f,G12f,v12f,v23f,Em,vm,Gm,
     &G23f,fi,v21f,C,D,ans2,F,H,xz,ans4,r2,xita2,r3,xita3,r4,xita4,
     &wzr,xitaw,wzx,wzy,x0,x1,G0,j1,j2,j3,gama,i,eq0,eq1,yp,ans,A,B0
     &,K23,K12,WVf,K22bit,K33bit,K22bic,K33bic,K22t,K22c,Sucm,Sutm,Susm
      DIMENSION RESIN(12),FIBER(11),SCFS(7)
c r表示λ,pusi表示开裂角 yp表示伊普瑟隆
!C
!      E11f=276000
!      E22f=19000
!      v12f=0.2
!      G12f=27000
!      v23f=0.36
!      Em=4.1*1000
!      vm=0.46
!      Sutm=121
!      Sucm=210
!      Susm=76
!      Vf=0.58974     
C
C
      E11F=FIBER(1)
      E22F=FIBER(2)
      V12F=FIBER(6)
      G12F=FIBER(9)
      V23F=FIBER(4)
      EM=RESIN(1)
      VM=RESIN(4)
      SUTM=RESIN(10)
      SUCM=RESIN(11)
      SUSM=RESIN(12)
C      VF=VF
!C
!      WRITE(101,*)'E11F=',e11F,'E22F=',E22F,'V12F=',V12F,'V23F=',V23F,
!     1'G12F=',G12F,'EM=',EM,'NUM=',VM,SUTM,SUCM,SUSM
      
      V=1-Vf
      v21f=v12f*E22f/E11f
      Gm=Em/(2*(1+vm))
      G23f=0.5*E22f/(1+v23f)
      k1=3-4*vm
      k2=(3-V23f-4*v12f*v21f)/(1+V23f)
      yp=(G23f+k2*Gm)/(Gm+k1*G23f)
      k=(Gm*(1+k2))/((1+yp)*(Gm+k1*G23f))
      b=1./sqrt(Vf)
      r=-log(yp)/(2*3.1415926)
      A=(2*E22f*Em*v12f**2+E11f*(Em*(v23f-1)-E22f*(2*vm**2+vm-1)))/(E11f
     &*(E22f+Em*(1-v23f)+E22f*vm)-2*E22f*Em*v12f**2)
      B0=(Em*(1+v23f)-E22f*(1+vm))/(E22f*(vm+4*vm**2-3)-Em*(1+v23f))
      K22t=(1+0.5*sqrt(Vf)*A+0.5*sqrt(Vf)*(3-Vf-sqrt(Vf))*B0)*(((Vf+0.3*
     &V)*E22f+0.7*V*Em)/(0.3*E22f+0.7*Em))
      K22c=(1-sqrt(Vf)*A*(Sucm-Sutm)/(4*Sucm)+(B0/(2*(1-sqrt(Vf))))*(
     &-Vf**2*(1-2*((Sucm-Sutm)/(2*Sucm))**2)+((Sucm+Sutm)*Vf/Sucm)*
     &(1+(Sucm-Sutm)/Sucm)-sqrt(Vf)*((Sucm-Sutm)/Sucm+1-2*((Sucm-Sutm)/(
     &2*Sucm))**2)))*(((Vf+0.3*V)*E22f+0.7*V*Em)/(0.3*E22f+0.7*Em))
c 由于纤维半径小b，故此处b0代表参数大B
c 横向压缩系数
      K23=2*Susm*sqrt((K22t*K22c)/(Sutm*Sucm))
c 横向剪切系数
      WVf=3.14*sqrt(Vf)*(1/(4*Vf)-0.03125-0.0039*Vf-0.00122*Vf**2)
      K12=(1-Vf*(G12f-Gm)/(G12f+Gm)*(WVf-(1./3.)))*((Vf+0.3*V)*G12f+0.7
     &*V*Gm)/(0.3*G12f+0.7*Gm)
     
c 轴向剪切系数
      K22Bit=((Vf+0.3*V)*E22f+0.7*V*Em)*(1+A*sqrt(Vf))/(0.3*E22f+
     &0.7*Em)
      K33Bit=K22Bit
c 横向双轴等值拉伸系数
      K22Bic=((Vf+0.3*V)*E22f+0.7*V*Em)/(0.3*E22f+0.7*Em)*(1-A*sqrt(Vf
     &)*(Sucm-Sutm)/(2*Sucm))
      K33Bic=K22Bic
!c 横向双轴等值压缩系数
!      write(*,*)K12
!       write(*,*)K23
!        write(*,*)K22t
!         write(*,*)K22c
!          write(*,*)K22Bit
!           write(*,*)K22Bic
!c 下面是二分法求解开裂角pusi
      
      x0=0.1
      x1=1.57
100   pusi=x0
      G0=(1-(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.14-pusi))+(1-k)*(1+4*r*
     &*2)*(sin(pusi))**2)/(2-k-k*(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.14
     &-pusi)))
      j1=k*G0-1-2*(1-k)*yp*exp(2*r*pusi)*cos(pusi)
      j2=2*(1-k)*yp*exp(2*r*pusi)*sin(pusi)
      j3=2*(1-k)*yp*exp(2*r*pusi)*(j1*cos(pusi)-j2*sin(pusi))/j2
      if (yp.LT.1) then
      gama=2*r*(j1**2+j2**2)/(j1**2+j2**2-2*j2*j3)
      end if
      if (yp .GE.1) then
      gama=-2*r*(j1**2+j2**2)/(j1**2+j2**2-2*j2*j3)
      end if
      fi=pusi-gama
      Rei=(exp(cmplx(0,fi))-exp(cmplx(0,pusi)))**cmplx(0.5,r)*
     &(exp(cmplx(0,fi))-exp(cmplx(0,-pusi)))**cmplx(0.5,-r)*
     &exp(cmplx(0,-fi))
      eq0=REAL((G0-1./k-2*(1-k)/(k*exp(cmplx(0,fi)))*exp(2*r*(pusi-3.14
     &15926)))*Rei)
      pusi=(x1+x0)/2
      G0=(1-(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.14-pusi))+(1-k)*(1+4*r*
     &*2)*(sin(pusi))**2)/(2-k-k*(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.14
     &-pusi)))
      j1=k*G0-1-2*(1-k)*yp*exp(2*r*pusi)*cos(pusi)
      j2=2*(1-k)*yp*exp(2*r*pusi)*sin(pusi)
      j3=2*(1-k)*yp*exp(2*r*pusi)*(j1*cos(pusi)-j2*sin(pusi))/j2
      if (yp .LT. 1) then
      gama=2*r*(j1**2+j2**2)/(j1**2+j2**2-2*j2*j3)
      end if
      if (yp .GE.1) then
      gama=-2*r*(j1**2+j2**2)/(j1**2+j2**2-2*j2*j3)
      end if
      fi=pusi-gama
      Rei=(exp(cmplx(0,fi))-exp(cmplx(0,pusi)))**cmplx(0.5,r)*
     &(exp(cmplx(0,fi))-exp(cmplx(0,-pusi)))**cmplx(0.5,-r)*
     &exp(cmplx(0,-fi))
      eq1=REAL((G0-1./k-2*(1-k)/(k*exp(cmplx(0,fi)))*exp(2*r*(pusi-3.14
     &15926)))*Rei)
      ans=eq0*eq1
      if (ans.GT. 0) then
         x0=(x0+x1)/2.
          end if
      if (ans .LT. 0) then
          x1=(x0+x1)/2.  
      end if
      if ((x1-x0) .GE.0.000001) THEN
      GO TO 100
      END IF
      z1=b*exp(cmplx(0,pusi))
      z2=exp(cmplx(0,-pusi))/b
      Z3=exp(cmplx(0,pusi))/b
c 总公式代入的是z2，但是excel表格用的是z3才能得到正确结果
      X=(z1-exp(cmplx(0,pusi)))**cmplx(-0.5,r)*(z1-exp(cmplx(0,-pusi)))*
     &*cmplx(-0.5,-r)
      D=(1-k)*exp(2*r*(pusi-3.1415926))
      C=-D*(cos(pusi)-2*r*sin(pusi))
      F=(1-(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.14-pusi))+(1-k)*(1+4*r**
     &2)*(sin(pusi))**2)/(4/k-2-2*(cos(pusi)+2*r*sin(pusi))*exp(2*r*(3.1
     &4-pusi)))
      N=F*z1+k/z1-(z1-exp(cmplx(0,pusi)))**cmplx(0.5,r)*(z1-exp(cmplx(0,
     &-pusi)))**cmplx(0.5,-r)*(F-0.5-D/z1)
c      N1=F*z3+k/z3+(z3-exp(cmplx(0,pusi)))**cmplx(0.5,r)*(z3-exp(cmplx(0
c     &,-pusi)))**cmplx(0.5,-r)*(F-0.5-D/z3)/yp
      r2z=z3-exp(cmplx(0,pusi))
      r2=sqrt(real(r2z)**2+imag(r2z)**2)
      xita2=atan(imag(r2z)/real(r2z))
      r3z=z3-exp(cmplx(0,-pusi))
      r3=sqrt(real(r3z)**2+imag(r3z)**2)
      xita3=atan(imag(r3z)/real(r3z))
      r4z=F-0.5-D/z3
      r4=sqrt(real(r4z)**2+imag(r4z)**2)
      xita4=atan(imag(r4z)/real(r4z))
      wzr=sqrt(r2*r3)*r4*exp(r*(xita3-xita2))/yp
      xitaw=0.5*(xita2+xita3)+xita4+r*(log(r2)-log(r3))
      wzx=wzr*cos(xitaw)
      wzy=wzr*sin(xitaw)
      N1=cmplx(wzx,wzy)+F*z3+k/z3
C n计算值正确，n1计算值错误，因此采取excel表格中实部与虚部分开计算的方法进行。
c 
      N2=F*exp(cmplx(0,-pusi))+k*exp(cmplx(0,pusi))
      N3=F*exp(cmplx(0,pusi))+k*exp(cmplx(0,-pusi))
      H=(cos(pusi)+2*r*sin(pusi))*(0.5-F)
      M=F-k/z1**2-((F-0.5)*z1+H+C/z1+D/z1**2)*X
      U=F-k/b**2*exp(cmplx(0,-2*pusi))-((F-0.5)*Z1+H+C*Z2+D*Z2**2)*X
      VZ=-cmplx(real(N3),-imag(N3))+cmplx(real(N1),-imag(N1))-exp(cmplx(
     &0,-2*pusi))*(N3-N)
      XZ=2*real((N3-N)*EXP(CMPLX(0,-pusi)))
      ans2=((Vf+0.3*V)*E22f+0.7*V*Em)/((0.3*E22f+0.7*Em)*2*(b-1))
c      ans1=(exp(cmplx(0,-2*pusi))*M*(1./b-b)-exp(cmplx(0,-pusi))*(N2-N1)
c     &+exp(cmplx(0,-pusi))*(2+exp(cmplx(0,-2*pusi)))*(N-N3))
c      ans3=real(ans1)*ans2
c      ans3=(z3-cmplx(cos(pusi),sin(pusi)))**cmplx(0.5,r)*(z2-
c     &cmplx(cos(pusi),-sin(pusi)))**cmplx(0.5,-r)*(F-0.5-D/z2)/yp+F*z3+k
c     &/z3
      ans5=z2*b
      ans1=(exp(cmplx(0,-pusi))*(1./b-b)*U+VZ)
      ans4=real(ans5*ans1-xz)*ans2
c 此处采取和excel中所列公式，并未采用word公式 ,ANS4=KTT22
C      write(*,*) ans4
C      pause
C      
      SCFS(1)=K22T
      SCFS(2)=ANS4
      SCFS(3)=K22BIT
      SCFS(4)=K22BIC
      SCFS(5)=K22C
      SCFS(6)=K23
      SCFS(7)=K12
C
C      WRITE(101,*) 'VF=',VF,'SCFS=',SCFS
      RETURN
      end
C
C**************************************************************计算开裂临界MISES应力
      SUBROUTINE  KGETLMISE(RESIN,FIBER,SCFS,VF,ALPHA,BETA,YT_UD,LMISE)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RESIN(12),FIBER(11),SCFS(7)
      DIMENSION BRIDGE_A(3,3),BRIDGE_B(3,3)   !二维桥联矩阵A、B，用于计算开裂临界MISES应力
      DIMENSION Vf_I(3,3),Vm_A(3,3),TEMP(3,3)
      DIMENSION SIGMA(3),SIGMA_M(3),SIGMA_F(3)
      REAL*8 ALPHA,BETA,LMISE,SIGMA22,YT_UD          !SIGMA22:开裂时的横向外载，YT_UD:单向板横向拉伸强度
      REAL*8 KT22,KTT22,KT22BI,KC22BI,KC22,K23,K12
      REAL*8 E11F,E22F,G12F,EM,NUM,NU12F,Sutm,Vf           !非常重要，否则传值错误
      
!C
!      E11f=276000
!      E22f=19000
!      NU12f=0.2
!      G12f=27000
!      NU23f=0.36
!      Em=4100
!      NUm=0.46
!      Sutm=121
!      Sucm=210
!      Susm=76
!      Vf=0.58974
!      YT_UD=73          !瞎写的，记得改
!      ALPHA=0.3
!      BETA=0.3
C
C
      E11F=FIBER(1)
      E22F=FIBER(2)
      NU12F=FIBER(6)
      G12F=FIBER(9)
C      NU23F=FIBER(4)
      EM=RESIN(1)
      NUM=RESIN(4)
      SUTM=RESIN(10)
C      SUCM=RESIN(11)
C      SUSM=RESIN(12)
      Gm=Em/(2.*(1.+NUm))
      KT22=SCFS(1)
      KTT22=SCFS(2)
!C
!      KT22BI=SCFS(3)
!      KC22BI=SCFS(4)
!      KC22=SCFS(5)
!      K23=SCFS(6)
!      K12=SCFS(7)
C
C      下面计算开裂时的临界MISES等效应力（LMISE）
C      下面求Aij
      BRIDGE_A=0.
      BRIDGE_A(1,1)=Em/E11f
      BRIDGE_A(2,2)=BETA+(1-BETA)*Em/E22f
      BRIDGE_A(3,3)=ALPHA+(1-ALPHA)*Gm/G12f
C      write(*,*)BRIDGE_A(3,3),'GM=',GM   !*******************************已验证
      BRIDGE_A(1,2)=(E11f*NUm-Em*NU12f)*(BRIDGE_A(2,2)-BRIDGE_A(1,1))
     &/(E11f-Em)
C      下面求Bij
      BRIDGE_B=0.
      Vf_I=0.
      Vm_A=0.
      TEMP=0.
      DO I=1,3
          Vf_I(I,I)=Vf
      ENDDO
      DO I=1,3
          DO J=1,3
              Vm_A(I,J)=(1-Vf)*BRIDGE_A(I,J)
          ENDDO
      ENDDO
      DO I=1,3
          DO J=1,3
              TEMP(I,J)=Vf_I(I,J)+Vm_A(I,J)
          ENDDO
      ENDDO
      
      CALL KINVER2(TEMP,BRIDGE_B,3,L)
C
      SIGMA=0.
      SIGMA_M=0.
      SIGMA_F=0.
      SIGMA22=Ktt22*YT_UD/(Ktt22-Kt22)-((Vf+0.3*(1-Vf))*E22f+0.7*
     &(1-Vf)*Em)*Sutm/((0.3*E22f+0.7*Em)*(Ktt22-Kt22))
      SIGMA(2)=SIGMA22
      DO I=1,3
          DO J=1,3
              SIGMA_F(I)=SIGMA_F(I)+BRIDGE_B(I,J)*SIGMA(J)
          ENDDO
      ENDDO
      DO I=1,3
          DO J=1,3
              SIGMA_M(I)=SIGMA_M(I)+BRIDGE_A(I,J)*SIGMA_F(J)
          ENDDO
      ENDDO
      LMISE=SQRT(SIGMA_M(1)**2+(KT22*SIGMA_M(2))**2-
     &KT22*SIGMA_M(1)*SIGMA_M(2))
C
C      write(*,*)'LMISE=',LMISE,'SIGMA_M(2)=',SIGMA_M(2),'KT22=',KT22
C     2,'sigma22=',sigma22,'SIGMAM11=',SIGMA_M(1)
C      write(*,'(3F20.15)') ((bridge_a(i,j),j=1,3),i=1,3)
C
      RETURN
      END
      
C**********************************************************************************基体层
      SUBROUTINE UMAT2(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     3 PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      include 'aba_param.inc'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
C ----------------------------------------------------------------
C ----------------------------------------------------------------
      REAL*8 KT22,KTT22,KC22,K12,K23,LMISE,K22,ZY,DELT,SK
      REAL*8 E11,E22,E33,V12,V13,V23,G12,G13,G23,V21,V31,V32,DELTA
      REAL*8 I1,I2,I3,P,Q,T,THETA,SIGMA1,SIGMA2,SIGMA3,MISE,
     1S1,S2,S3,F1,F11,F44    
      INTEGER N
      DIMENSION RESIN(12),FIBER(11),EGM(9),SUUM(3),EGF(9),SUUF(2),
     1SF(6,6),SM(6,6),SP(6,6),A(6,6),ETM(2,20),SSF(6),SSM(6),DSTRESS(6),
     2TMPD(6),TMPD1(6)

      N=0
      DO I=1,12
        RESIN(I)=PROPS(I)
      ENDDO
      N=N+12
C
      MSEG=PROPS(N+1)
      N=N+1
      DO 10 I=1,MSEG
10    ETM(1,I)=PROPS(N+I)
      N=N+MSEG
      DO 20 I=1,MSEG
20    ETM(2,I)=PROPS(N+I)
      N=N+MSEG
      DELT=PROPS(N+1)    !DELT，界面层衰减系数
      N=N+1
      zy=props(N+1)
      N=N+1
c 
C**********************************最大应变不超过一特定值ZY
      IF(STATEV(19).EQ.0)THEN
      DO I=1,6
          IF(STRAN(I).GT.ZY)THEN
              STATEV(19)=1
            write(101,*)'应变>',ZY
            write(101,*)KINC,NOEL,'STRANI=',I,stran(I)  !增量步，单元，应变
            DO K=1,6
                WRITE(101,*)K,STRAN(K)
            ENDDO
          !STOP
          ENDIF
      ENDDO
      ENDIF
C**********************************
      E11=RESIN(1)
      E22=RESIN(2)
      E33=RESIN(3)
      V12=RESIN(4)
      V13=RESIN(5)
      V23=RESIN(6)
      G12=RESIN(7)
      G13=RESIN(8)
      G23=RESIN(9)
      SUUM(1)=RESIN(10)
      SUUM(2)=RESIN(11)
      SUUM(3)=RESIN(12)
      V21=E22*V12/E11
      V31=E33*V13/E11
      V32=E33*V23/E22
C
      DDSDDE=0.
      DELTA=1-V12*V21-V23*V32-V31*V13-2*V21*V32*V13
      DDSDDE(1,1)=E11*(1-V23*V32)/DELTA
      DDSDDE(2,2)=E22*(1-V13*V31)/DELTA
      DDSDDE(3,3)=E33*(1-V12*V21)/DELTA
      DDSDDE(1,2)=E11*(V21+V31*V23)/DELTA
      DDSDDE(1,3)=E11*(V31+V21*V32)/DELTA
      DDSDDE(2,3)=E22*(V32+V12*V31)/DELTA
      DDSDDE(4,4)=G23
      DDSDDE(5,5)=G13
      DDSDDE(6,6)=G12               !改成传统顺序了
      DO I=1,3
        DO J=I,3
          DDSDDE(J,I)=DDSDDE(I,J)
        END DO
      END DO
      
      DO I=1,6
        DO J=1,6
          DDSDDE(I,J)=DDSDDE(I,J)*(DELT**STATEV(19))   !界面层用STATEV(19)
        END DO
      END DO
C
      CALL ABQ2NORM_1D(STRESS)     !应力顺序换成传统顺序
      DO I=1,6
          TMPD(I)=DSTRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD) 
        DO I=1,6
          TMPD1(I)=STRAN(I)
        END DO
        CALL ABQ2NORM_1D(TMPD1)
        
        STRESS=0.
      DO K1=1,6
      DO K2=1,6
      STRESS(K1)=STRESS(K1)+DDSDDE(K1,K2)*(TMPD1(K2)+TMPD(K2))   !全量更新
      END DO
      END DO
C------------------------------
      IF(STATEV(19).EQ.0)THEN
C------------------------------
      F1=(SUUM(2)-SUUM(1))/(SUUM(1)*SUUM(2))
      F11=1/(SUUM(1)*SUUM(2))
      F44=1/(SUUM(3)**2)

      MISE=F1*(STRESS(1)+STRESS(2)+STRESS(3))+F11*(STRESS(1)**2+
     1STRESS(2)**2+STRESS(3)**2-STRESS(1)*STRESS(2)-STRESS(2)*STRESS(3)
     2-STRESS(3)*STRESS(1))+F44*(STRESS(4)**2+STRESS(5)**2+STRESS(6)**2)
          
C*********************************TSAI-WU
      
        IF(MISE.GT.1.)THEN

          STATEV(19)=1
          
         WRITE(101,*)'界面层破坏','NOEL=',NOEL,'KINC=',KINC
       ENDIF
      ENDIF
c
      CALL NORM2ABQ_1D(STRESS)   !应力换成ABAQUS顺序
      CALL NORM2ABQ_2D(DDSDDE)   !刚度矩阵换成ABAQUS顺序
      RETURN
      END
     
