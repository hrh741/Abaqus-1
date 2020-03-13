      ! Motamedi D. Nonlinear XFEM modeling of delamination in fiber reinforced composites considering uncertain fracture properties and effect of fiber bridging[D]. University of British Columbia, 2013.
	  Subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
      props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, time, dtime,
      kstep, kinc, jelem, params, ndload, jdltyp, adlmag, predef, npredf,
      lflags, mlvarx, ddlmag, mdload, pnewdt, jprops, njprop, period)
      
      include 'ABA_PARAM.INC'
      parameter(zero=0.D0, half=0.5D0, one=1.D0, seven=7.0D0, eight=8.0D0)
      
      Dimension rhs(ndofel, *), amatrx(ndofel, ndofel), props(*),
      svars(*), energy(8), coords(mcrd, nnode), u(ndofel),
      du(ndofel, *), v(ndofel), a(ndofel), time(2), params(*),
      jdltyp(mdload, *), adlmag(mdload, *), ddlmag(mdload, *),
      predef(2, npredf, nnode), lflags(*), jprops(*)
      
      
      Real *8 gausspoint(140, mcrd), xcr(mcrd+1, mcrd), dndx(nnode, mcrd+1), u_midpnt(nnode, mcrd), b(mcrd*mcrd, ndofel), dndxi(nnode, mcrd), c_coords(nnode, mcrd), gs(2*mcrd, mcrd*mcrd), u_ms(mcrd*mcrd, mcrd*mcrd), c(2*mcrd, 2*mcrd), dep(mcrd, mcrd), qt(mcrd, mcrd), qtt(mcrd, mcrd), qr(mcrd, mcrd), qrt(mcrd, mcrd), dgep(mcrd, mcrd), fb(mcrd, mcrd), fbt(mcrd, mcrd), strglb(mcrd, mcrd), strglob(mcrd, mcrd), dudx(mcrd, mcrd), nx(1, nnode), gweight(140), ff(mcrd*mcrd), strsg(mcrd*mcrd), coor(mcrd), dltu(mcrd), dv(mcrd)
      ! user - elemnt arrays
      ! general element values
      ! gauss integration variables(3 integ point)
      ! arrays for 3D  element
      Real *8 u_phi(nnode)
      Real *8 uneg(mcrd), upos(mcrd), ung(mcrd), ups(mcrd), u1(mcrd), ubas(ndofel), ec(2*mcrd), strs(2*mcrd), btstrs(ndofel), trct(mcrd), lhsc(ndofel), rhsc(ndofel)
      Real *8 bt(ndofel, mcrd*mcrd), gst(mcrd*mcrd, 2*mcrd)
      Real *8 jacb(mcrd, mcrd), invjacb(mcrd, mcrd), invfb(mcrd, mcrd), invfbt(mcrd, mcrd), iunit(mcrd, mcrd), dgbc(mcrd, ndofel), btgt(ndofel, 2*mcrd), gb(2*mcrd, ndofel), msb(mcrd*mcrd, ndofel), btmsb(ndofel, ndofel), bctdgbc(ndofel, ndofel), kmat(ndofel, ndofel), kgem(ndofel, ndofel), kcon(ndofel, ndofel)
      Real *8 bc(mcrd, ndofel), bct(ndofel, mcrd), cstrn(mcrd, mcrd), eec(mcrd, mcrd), bcres(ndofel, mcrd)
      Real *8 xi, yi, zi, weight, enrcoh, enrjc, e11, e22, e33, g12, g23, g31, nu12, nu21, nu23, nu32, nu31, nu13, detfb, detj, nu, uacrt, jsdv, hpoint, enrelm, nnint, kpen, tmax, jcrt, dmg1, dmg2, dmg3
      
      Integer count1, count2, intp, iintp, i, j, k, l, m, n
      
      !*************************************************************************************
      !***initialisation: important !! FORTRAN DOES NOT PUT ZEROS IN THERE
      ! automatically***
      !*************************************************************************************
      ! open level - set values from pre - processing
      open(16, file='C:/ABAQUS-Matlab/IOstat /philvlst.dat ',STATUS=' old ')
      
      nintp=24 ! Number of integration points
      uacrt=props(1) ! Crack Lenght
      dv(1)=props(2) ! Material Initiation Opening Failure
      dv(2)=props(3) ! Material Final Opening Failure
      e11=props(4) ! Material Constant
      e22=props(5) ! Material Constant
      e33=props(6) ! Material Constant
      g12=props(7) ! Material Constant
      g23=props(8) ! Material Constant
      g31=props(9) ! Material Constant
      nu12=props(10) ! Material Constant
      nu23=props(11) ! Material Constant
      nu31=props(12) ! Material Constant
	  
      kpen=props(13) ! Material Constant
      jsdv=props(14) ! G Standard Deviation
      
      ! calculating the traction-separation constant
      tmax=kpen*dv(1)
      jcrt=tmax*dv(2)/2
      width=1 ! Element Tickness
      !**********************************
      !***zero the required matrices***
      !**********************************
      enrelm=0.0
      
      Do i=1, nnode
        nx(1,i)=0.0
        u_phi(i)=0.0
        Do j=1, mcrd+1
          dndx(i,j)=0.0
        End Do
      End Do
      
      Do j=1, mcrd
        coor(j)=0.0
      End Do
      
      Do i=1, 2*mcrd
        Do j=1, 2*mcrd
          c(i,j)=0.0
        End Do
      End Do
      
      Do i=1, mcrd*mcrd
        Do j=1, ndofel
          b(i,j)=0.0
          bt(j,i)=0.0
        End Do
      End Do
      
      Do i=1, 2*mcrd
        Do j=1, mcrd*mcrd
          gs(i,j)=0.0
          gst(j,i)=0.0
        End Do
      End Do
      
      Do i=1, nnode
        Do j=1, mcrd
          c_coords(i,j)=0.0
        End Do
      End Do
      
      Do i=1, mcrd
        Do j=1, mcrd
          iunit(i,j)=0.0
          jacb(i,j)=0.0
          invjacb(i,j)=0.0
          dep(i,j)=0.0
          dgep(i,j)=0.0
          
          qt(i,j)=0.0
          qtt(i,j)=0.0
          qr(i,j)=0.0
          qrt(i,j)=0.0
        End Do
      End Do
      
      Do i=1, nnode
        Do j=1, mcrd
          u_midpnt(i,j)=0.0
        End Do
      End Do
      
      Do i=1, nsvars
        svars(i)=0.0
      End Do
      
      Do i=1, ndofel
        Do j=1, ndofel
          amatrx(i,j)=0.0
          kgem(i,j)=0.0
          kmat(i,j)=0.0
          kcon(i,j)=0.0
        End Do
      End Do
      
      Do i=1, ndofel
        rhs(i,1)=0.0
        rhsc(i)=0.0
        lhsc(i)=0.0
        ubas(i)=0.0
      End Do
      !********************
      !***dummy arrays***
      !********************
      Do i=1, mcrd
        iunit(i,i)=1.0
      End Do
      ints=1 ! Integration point scheme (1: gauss)
      stype=1 ! Element type (1: B8, 2: T4)
      !****************
      !***hook law***
      !****************
      nu21=nu12*e22/e11
      nu32=nu23*e33/e22
      nu13=nu31*e11/e33
      delt=(1-nu12*nu21-nu23*nu32-nu13*nu31-2*nu12*nu23*nu31)
      
      c(1,1)=e11*(1-nu32*nu23)/delt
      c(1,2)=e11*(nu21+nu31*nu23)/delt
      c(1,3)=e11*(nu31+nu21*nu32)/delt
      c(2,1)=e22*(nu12+nu13*nu23)/delt
      c(2,2)=e22*(1-nu13*nu31)/delt
      c(2,3)=e22*(nu32+nu31*nu12)/delt
      c(3,1)=e33*(nu13+nu12*nu23)/delt
      c(3,2)=e33*(nu23+nu13*nu21)/delt
      c(3,3)=e33*(1-nu12*nu21)/delt
      c(4,4)=g12/2
      c(5,5)=g23/2
      
      c(6,6)=g31/2
      !*******************************************
      !***finding deform shape of coordinates***
      !*******************************************
      Do i=1, nnode
        Do j=1, mcrd
          c_coords(i,j)=coords(j,i)+u(2*mcrd*(i-1)+j)
        End Do
      End Do
      !************************************************
      !***finding reference coordinate deformation***
      !************************************************
      Do i=1, ndofel
        u(i)=u(i)-du(i,1)
      End Do
      !*************************
      !***calling level set***
      !*************************
      Call lvlsetrdr(u_phi,jelem,nnode) ! READ THE Level-Set FROM MATLAB
      !*******************************************************
      !***calling the local crack 'S PLANE IN EACH ELEMENT ***
      !*******************************************************
      Call midplnfind(nnint,u_midpnt,qt,u_phi,coords,nnode,mcrd)
      
      Do i=1, mcrd
        Do j=1, mcrd
          qtt(i,j)=qt(j,i)
        End Do
      End Do
      !*********************************
      !***calling gauss coordinates***
      !*********************************
      Call subtgauss(gausspoint,gweight,nnode,mcrd)
      !****************************************
      !***looking for opening displacement***
      !****************************************
      count1=0
      count2=0
      
      Do i=1, mcrd
        upos(i)=0.0
        uneg(i)=0.0
        ups(i)=0.0
        ung(i)=0.0
      End Do
      !************************************
      !***loop over integration points***
      !************************************
      Do intp=1, nintp
        hpoint=0.0
        Do i=1, mcrd
          u1(i)=0.0
          dltu(i)=0.0
        End Do
        
        coor(1)=gausspoint(intp,1)
        coor(2)=gausspoint(intp,2)
        coor(3)=gausspoint(intp,3)
        weight=gweight(intp)
        
        !***********************************
        !***calling the shape functions***
        !***********************************
        Call lagrangebasis(coor,nx,dndx,nnode,mcrd)
        !*********************************
        !***calculating gpt-level set***
        !*********************************
        Do i=1, nnode
          hpoint=hpoint+u_phi(i)*dndx(i,4)/abs(u_phi(i))
        End Do
        If(hpoint>0.0) Then
          hpoint=1.0
        Else If(hpoint<0.0) Then
          hpoint=-1.0
        Else
          hpoint=0.01
          weight=0.0
        End If
        !********************************************
        !***checking for gpt within contactbound***
        !********************************************
        Do i=1, nnode
          Do j=1, mcrd
            u1(j)=u1(j)+dndx(i,4)*u(6*(i-1)+j) +(hpoint-abs(u_phi(i))/u_phi(i))* (dndx(i,4)*u(6*(i)+j-3))
          End Do
        End Do
        
        If(hpoint<0.0) Then
          count1=count1+1
          uneg(1)=uneg(1)+u1(1)
          uneg(2)=uneg(2)+u1(2)
          uneg(3)=uneg(3)+u1(3)
        Else
          count2=count2+1
          upos(1)=upos(1)+u1(1)
          upos(2)=upos(2)+u1(2)
          upos(3)=upos(3)+u1(3)
        End If
        
      End Do
      Do i=1, mcrd
        Do j=1, mcrd
          ups(i)=ups(i)+qt(i,j)*upos(j)/count2
          ung(i)=ung(i)+qt(i,j)*uneg(j)/count1
        End Do
      End Do
      
      dltu(1)=ups(1)-ung(1)
      dltu(2)=ups(2)-ung(2)
      dltu(3)=ups(3)-ung(3)
      
      !*********************************************
      !***forming stiffness and residual matrix***
      !*********************************************
      
      !*******************************
      !***calling elast-plast relationship***
      
      !*******************************
      Do i=1, mcrd
        Do j=1, mcrd
          dgep(i,j)=0.0
        End Do
      End Do
      
      Call elsplc(dep,dltu,dv,uacrt,c,mcrd,enrelm,jtype,kpen,
      dmg1,dmg2,dmg3)
      
      ! assigning damage indeces and crack opening displacement to user-variables
      
      svars(1)=dmg1
      svars(2)=dmg2
      svars(3)=dmg3
      svars(4)=dltu(1)
      svars(5)=dltu(2)
      svars(6)=dltu(3)
      
      ! transforimg the elastic-plasti! relationship into local crack plane
      
      Do i=1, mcrd
        Do j=1, mcrd
          Do k=1, mcrd
            Do l=1, mcrd
              dgep(i,j)=dgep(i,j)+qtt(i,k)*dep(k,l)*qt(l,j)
            End Do
          End Do
        End Do
      End Do
      
      !************************************
      !***loop over integration points***
      !************************************
      Do i=1, ndofel
        Do j=1, ndofel
          btmsb(i,j)=0.0
        End Do
      End Do
      
      Do iintp=1, nintp
        coor(1)=gausspoint(iintp,1)
        coor(2)=gausspoint(iintp,2)
        coor(3)=gausspoint(iintp,3)
        weight=gweight(iintp)
        !***********************************
        !***calling the shape functions***
        !***********************************
        Call lagrangebasis(coor,nx,dndx,nnode,mcrd)
        xi=coor(1)
        yi=coor(2)
        zi=coor(3)
        Do i=1, mcrd
          Do j=1, mcrd
            jacb(i,j)=0.0
          End Do
        End Do
        !*******************************
        !***forming jacobian matrix***
        
        !*******************************
        Do i=1, mcrd
          Do j=1, mcrd
            Do k=1, nnode
              jacb(i,j)=jacb(i,j)+coords(i,k)*dndx(k,j)
            End Do
          End Do
        End Do
        
        detj=jacb(1,1)*(jacb(2,2)*jacb(3,3)-jacb(2,3)*jacb(3,2))- jacb(1,2)*(jacb(2,1)*jacb(3,3)-jacb(3,1)*jacb(2,3))+ jacb(1,3)*(jacb(2,1)*jacb(3,2)-jacb(2,2)*jacb(3,1))
        If(detj<0.0) Then
          detj=(-1)*detj
        End If
        invjacb(1,1)=(jacb(2,2)*jacb(3,3)-jacb(2,3)*jacb(3,2))/detj
        invjacb(1,2)=(jacb(1,3)*jacb(3,2)-jacb(1,2)*jacb(3,3))/detj
        invjacb(1,3)=(jacb(1,2)*jacb(2,3)-jacb(1,3)*jacb(2,2))/detj
        invjacb(2,1)=(jacb(2,3)*jacb(3,1)-jacb(2,1)*jacb(3,3))/detj
        invjacb(2,2)=(jacb(1,1)*jacb(3,3)-jacb(1,3)*jacb(3,1))/detj
        invjacb(2,3)=(jacb(1,3)*jacb(2,1)-jacb(1,1)*jacb(2,3))/detj
        invjacb(3,1)=(jacb(2,1)*jacb(3,2)-jacb(2,3)*jacb(3,1))/detj
        invjacb(3,2)=(jacb(1,2)*jacb(3,1)-jacb(1,1)*jacb(3,2))/detj
        invjacb(3,3)=(jacb(2,2)*jacb(1,1)-jacb(2,1)*jacb(1,2))/detj
        !********************************
        !***forming derivaties dn/dx***
        !********************************
        Do i=1, nnode
          Do j=1, mcrd
            dndxi(i,j)=0.0
          End Do
        End Do
        
        Do i=1, nnode
          Do j=1, mcrd
            Do k=1, mcrd
              dndxi(i,j)=dndxi(i,j)+dndx(i,k)*invjacb(k,j)
            End Do
          End Do
        End Do
        !*******************************
        !***calculating b-level set***
        !*******************************
        hpoint=0.0
        Do i=1, mcrd
          Do k=1, mcrd
            dudx(i,k)=0.0
          End Do
        End Do
        
        Do i=1, mcrd
          Do k=1, mcrd
            Do j=1, nnode
              dudx(i,k)=dudx(i,k)+dndxi(j,k)*u(6*(j-1)+i)
            End Do
          End Do
        End Do
        !*********************************
        !***calculating gpt-level set***
        
        !*********************************
        Do i=1, nnode
          hpoint=hpoint+u_phi(i)*dndx(i,4)/abs(u_phi(i))
        End Do
        If(hpoint>0.0) Then
          hpoint=1.0
        Else If(hpoint<0.0) Then
          hpoint=-1.0
        Else
          hpoint=0.01
          weight=0.0
        End If
        Do i=1, mcrd
          Do k=1, mcrd
            Do j=1, nnode
              dudx(i,k)=dudx(i,k)+(hpoint-abs(u_phi(i))/u_phi(i))* dndxi(j,k)*u(6*(j-1)+i+mcrd)
            End Do
          End Do
        End Do
        !************************
        !***forming b matrix***
        !************************
        Do i=1, mcrd*mcrd
          Do j=1, ndofel
            b(i,j)=0.0
            bt(j,i)=0.0
          End Do
        End Do
        
        Call gradb(b,dndxi,hpoint,u_phi,nnode,mcrd,ndofel)
        
        bt=transpose(b)
        
        Do i=1, mcrd*mcrd
          ff(i)=0.0
        End Do
        Do i=1, 2*mcrd
          ec(i)=0.0
        End Do
        
        Do i=1, mcrd
          Do j=1, mcrd
            fb(i,j)=0.0
            fbt(i,j)=0.0
          End Do
        End Do
        
        ff=matmul(b,u)
        
        dudx(1,1)=ff(1)
        dudx(2,1)=ff(2)
        dudx(3,1)=ff(3)
        dudx(1,2)=ff(4)
        dudx(2,2)=ff(5)
        dudx(3,2)=ff(6)
        dudx(1,3)=ff(7)
        dudx(2,3)=ff(8)
        dudx(3,3)=ff(9)
        
        
        !**************************************************
        !***forming green-lagrangian large strains***
        !**************************************************
        
        Do i=1, mcrd*mcrd
          If(i==1 .Or. i==5 .Or. i==9) Then
            ff(i)=ff(i)+1.0
          End If
        End Do
        
        fb(1,1)=ff(1)
        fb(2,1)=ff(2)
        fb(3,1)=ff(3)
        fb(1,2)=ff(4)
        fb(2,2)=ff(5)
        fb(3,2)=ff(6)
        fb(1,3)=ff(7)
        fb(2,3)=ff(8)
        fb(3,3)=ff(9)
        !*******************************
        !***forming invf gradf***
        !*******************************
        Do i=1, mcrd
          Do j=1, mcrd
            invfb(i,j)=0.0
          End Do
        End Do
        detfb=fb(1,1)*(fb(2,2)*fb(3,3)-fb(2,3)*fb(3,2))- fb(1,2)*(fb(2,1)*fb(3,3)-fb(3,1)*fb(2,3))+ fb(1,3)*(fb(2,1)*fb(3,2)-fb(2,2)*fb(3,1))
        If(detfb<0.0) Then
          detfb=(-1)*detfb
        End If
        invfb(1,1)=(fb(2,2)*fb(3,3)-fb(2,3)*fb(3,2))/detfb
        invfb(1,2)=(fb(1,3)*fb(3,2)-fb(1,2)*fb(3,3))/detfb
        invfb(1,3)=(fb(1,2)*fb(2,3)-fb(1,3)*fb(2,2))/detfb
        invfb(2,1)=(fb(2,3)*fb(3,1)-fb(2,1)*fb(3,3))/detfb
        invfb(2,2)=(fb(1,1)*fb(3,3)-fb(1,3)*fb(3,1))/detfb
        invfb(2,3)=(fb(1,3)*fb(2,1)-fb(1,1)*fb(2,3))/detfb
        invfb(3,1)=(fb(2,1)*fb(3,2)-fb(2,3)*fb(3,1))/detfb
        invfb(3,2)=(fb(1,2)*fb(3,1)-fb(1,1)*fb(3,2))/detfb
        invfb(3,3)=(fb(2,2)*fb(1,1)-fb(2,1)*fb(1,2))/detfb
        
        fbt=transpose(fb)
        invfbt=transpose(invfb)
        
        Do i=1, 2*mcrd
          Do j=1, mcrd*mcrd
            gs(i,j)=0.0
            gst(j,i)=0.0
          End Do
        End Do
        
        Call tangsms(gs,ff,mcrd)
        
        gst=transpose(gs)
        gb=matmul(gs,b)
        ec=matmul(gb,u)
        
        !**************************
        !***stress calculation***
        !**************************
        Do i=1, 2*mcrd
          strs(i)=0.0
        End Do
        Do i=1, mcrd*mcrd
          strsg(i)=0.0
        End Do
        Do i=1, mcrd
          Do j=1, mcrd
            strglb(i,j)=0.0
            strglob(i,j)=0.0
          End Do
        End Do
        
        strs=matmul(c,ec)
        
        strglob(1,1)=strs(1)
        strglob(2,2)=strs(2)
        strglob(3,3)=strs(3)
        strglob(2,1)=strs(4)
        strglob(1,2)=strs(4)
        strglob(1,3)=strs(5)
        strglob(3,1)=strs(5)
        strglob(2,3)=strs(6)
        strglob(3,2)=strs(6)
        
        strglb=matmul(strglob,iunit)
        strsg(1)=strglb(1,1)
        strsg(2)=strglb(1,2)
        strsg(3)=strglb(1,3)
        strsg(4)=strglb(2,1)
        strsg(5)=strglb(2,2)
        strsg(6)=strglb(2,3)
        strsg(7)=strglb(3,1)
        strsg(8)=strglb(3,2)
        strsg(9)=strglb(3,3)
        !************************************
        !****forming gs ms matrices****
        !************************************
        Do i=1, mcrd**2
          Do j=1, mcrd**2
            u_ms(i,j)=0.0
          End Do
        End Do
        
        Do i=1, mcrd
          u_ms(i,i)=strsg(1)
          u_ms(i+mcrd,i+mcrd)=strsg(5)
          u_ms(i+2*mcrd,i+2*mcrd)=strsg(9)
          u_ms(i+mcrd,i)=strsg(4)
          u_ms(i,i+mcrd)=strsg(2)
          u_ms(i+2*mcrd,i)=strsg(7)
          u_ms(i,i+2*mcrd)=strsg(3)
          u_ms(i+2*mcrd,i+mcrd)=strsg(8)
          u_ms(i+mcrd,i+2*mcrd)=strsg(6)
        End Do
        !********************************************
        
        !***forming cohesive contact b! matrix***
        !********************************************
        Do i=1, mcrd
          Do j=1, ndofel
            bc(i,j)=0.0
            bct(j,i)=0.0
            bcres(j,i)=0.0
          End Do
        End Do
        Do i=1, nnode
          bcres(6*i-5,1)=0.0
          bcres(6*i-4,2)=0.0
          bcres(6*i-3,3)=0.0
          bcres(6*i-2,1)=dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
          bcres(6*i-1,2)=dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
          bcres(6*i,3)=dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
          bc(1,6*i-2)=-2*dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
          bc(2,6*i-1)=-2*dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
          bc(3,6*i)=-2*dndx(i,4)*(hpoint-abs(u_phi(i))/u_phi(i))
        End Do
        
        bct=transpose(bc)
        !********************************
        !***forming stiffness matrix***
        !********************************
        Do i=1, 2*mcrd
          Do j=1, ndofel
            gb(i,j)=0.0
            btgt(j,i)=0.0
          End Do
        End Do
        
        gb=matmul(gs,b)
        btgt=transpose(gb)
        
        Do i=1, ndofel
          Do j=1, ndofel
            Do k=1, 2*mcrd
              Do l=1, 2*mcrd
                kmat(i,j)=kmat(i,j)+ btgt(i,k)*c(k,l)*gb(l,j)*weight*detj*width
              End Do
            End Do
          End Do
        End Do
        
        msb=matmul(u_ms,b)
        btmsb=matmul(bt,msb)
        
        Do i=1, ndofel
          Do j=1, ndofel
            kgem(i,j)=kgem(i,j)+btmsb(i,j)*weight*detj*width
          End Do
        End Do
        
        dgbc=matmul(dgep,bc)
        bctdgbc=matmul(bct,dgbc)
        
        Do i=1, ndofel
          
          Do j=1, ndofel
            kcon(i,j)=kcon(i,j)+bctdgbc(i,j)*weight*detj*width
          End Do
        End Do
        
        Do i=1, mcrd
          trct(i)=0.0
        End Do
        
        Do i=1, mcrd
          Do j=1, ndofel
            trct(i)=trct(i)+dgbc(i,j)*(u(j)+du(j,1))
          End Do
        End Do
        !****************************
        !***residual calculation***
        !****************************
        Do i=1, mcrd*mcrd
          ff(i)=0.0
        End Do
        Do i=1, 2*mcrd
          ec(i)=0.0
        End Do
        
        Do i=1, mcrd*mcrd
          Do j=1, ndofel
            ff(i)=ff(i)+b(i,j)*(u(j)+du(j,1))
          End Do
        End Do
        !*************************************************
        !***forming green-lagrangian large strain***
        !*************************************************
        Do i=1, mcrd*mcrd
          If(i==1 .Or. i==5 .Or. i==9) Then
            ff(i)=ff(i)+1.0
          End If
        End Do
        
        fb(1,1)=ff(1)
        fb(2,1)=ff(2)
        fb(3,1)=ff(3)
        fb(1,2)=ff(4)
        fb(2,2)=ff(5)
        fb(3,2)=ff(6)
        fb(1,3)=ff(7)
        fb(2,3)=ff(8)
        fb(3,3)=ff(9)
        
        Do i=1, 2*mcrd
          Do j=1, mcrd*mcrd
            gs(i,j)=0.0
          End Do
        End Do
        
        Call tangsms(gs,ff,mcrd)
        
        gb=matmul(gs,b)
        btgt=transpose(gb)
        
        
        Do i=1, mcrd
          Do j=1, mcrd
            cstrn(i,j)=0.0
            eec(i,j)=0.0
          End Do
        End Do
        
        Do i=1, mcrd
          Do j=1, mcrd
            Do k=1, mcrd
              cstrn(i,j)=fb(k,i)*fb(k,j)+cstrn(i,j)
            End Do
          End Do
        End Do
        Do i=1, mcrd
          Do j=1, mcrd
            If(i==j) Then
              eec(i,j)=0.5*cstrn(i,j)-0.5
            Else
              eec(i,j)=0.5*cstrn(i,j)
            End If
          End Do
        End Do
        
        ec(1)=eec(1,1)
        ec(2)=eec(2,2)
        ec(3)=eec(3,3)
        ec(4)=eec(1,2)
        ec(5)=eec(2,3)
        ec(6)=eec(3,1)
        
        strs=matmul(c,ec)
        
        ! xfem elements residual forces due to large deformation
        
        Do i=1, ndofel
          Do j=1, 2*mcrd
            rhs(i,1)=rhs(i,1)- btgt(i,j)*strs(j)*weight*detj*width
          End Do
        End Do
        
        ! xfem residual forces due to cohesive region or contact interface
        
        Do i=1, ndofel
          Do j=1, mcrd
            rhs(i,1)=rhs(i,1)+ bcres(i,j)*trct(j)*weight*detj*width
            rhsc(i)=rhsc(i)+ bcres(i,j)*trct(j)*weight*detj*width
          End Do
        End Do
      End Do ! END IF LOOP OVER GAUSS POINTS
      !**************************************
      !***forming right-hand-side matrix***
      !**************************************
      amatrx=kmat+kgem+kcon
      Return
      End Subroutine uel
