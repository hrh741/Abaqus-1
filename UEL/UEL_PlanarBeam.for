c
c simple 2-d linear beam element with generalized section properties
c
      subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
      1 props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, time, dtime,
      2 kstep, kinc, jelem, params, ndload, jdltyp, adlmag, predef, npredf,
      3 lflags, mlvarx, ddlmag, mdload, pnewdt, jprops, njprop, period)
c
      include ’aba_param.inc’
c
      dimension rhs(mlvarx, *), amatrx(ndofel, ndofel), svars(*), props(*),
      1 energy(7), coords(mcrd,nnode), u(ndofel), du(mlvarx, *), v(ndofel),
      2 a(ndofel), time(2), params(*), jdltyp(mdload, *), adlmag(mdload, *),
      3 ddlmag(mdload, *),predef(2, npredf, nnode), lflags(4), jprops(*)
c
      dimension b(2, 7), gauss(2)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0,
      1 six=6.d0, eight=8.d0, twelve=12.d0)
      data gauss/.211324865d0, .788675135d0/
c
c calculate length and direction cosines
c
      dx=coords(1, 2)-coords(1, 1)
      dy=coords(2, 2)-coords(2, 1)
      dl2=dx**2+dy**2
      dl=sqrt(dl2)
      hdl=dl/two
      acos=dx/dl
      asin=dy/dl
c
c initialize rhs and lhs
c
      do k1=1, 7
      rhs(k1, 1)= zero
      do k2=1, 7
      amatrx(k1, k2)= zero
      end do
      end do
c
      nsvint=nsvars/2
c
c loop over integration points
c
      do kintk=1, 2
      g=gauss(kintk)
c
c make b-matrix
c
      b(1, 1)=(-three+four*g)*acos/dl
      b(1, 2)=(-three+four*g)*asin/dl
      b(1, 3)=zero
      b(1, 4)=(-one+four*g)*acos/dl
      b(1, 5)=(-one+four*g)*asin/dl
      b(1, 6)=zero
      b(1, 7)=(four-eight*g)/dl
      b(2, 1)=(-six+twelve*g)*-asin/dl2
      b(2, 2)=(-six+twelve*g)* acos/dl2
      b(2, 3)=(-four+six*g)/dl
      b(2, 4)= (six-twelve*g)*-asin/dl2
      b(2, 5)= (six-twelve*g)* acos/dl2
      b(2, 6)= (-two+six*g)/dl
      b(2, 7)=zero
c calculate (incremental) strains and curvatures
c
      eps=zero
      deps=zero
      cap=zero
      dcap=zero
      do k=1, 7
      eps=eps+b(1, k)*u(k)
      deps=deps+b(1, k)*du(k, 1)
      cap=cap+b(2, k)*u(k)
      dcap=dcap+b(2, k)*du(k, 1)
      end do
c
c call constitutive routine ugenb
c
      isvint=1+(kintk-1)*nsvint
      bn=zero
      bm=zero
      daxial=zero
      dbend=zero
      dcoupl=zero
      call ugenb(bn, bm, daxial, dbend, dcoupl, eps, deps, cap, dcap,
      1 svars(isvint), nsvint, props, nprops)
c
c assemble rhs and lhs
c
      do k1=1, 7
      rhs(k1, 1)=rhs(k1, 1)-hdl*(bn*b(1, k1)+bm*b(2, k1))
      bd1=hdl*(daxial*b(1, k1)+dcoupl*b(2, k1))
      bd2=hdl*(dcoupl*b(1, k1)+dbend *b(2, k1))
      do k2=1, 7
      amatrx(k1, k2)=amatrx(k1, k2)+bd1*b(1, k2)+bd2*b(2, k2)
      end do
      end do
      end do
c
      return
      end
	  
	  subroutine ugenb(bn,bm,daxial,dbend,dcoupl,eps,deps,cap,dcap,
      1 svint,nsvint,props,nprops)
c
      include ’aba_param.inc’
c
      parameter(zero=0.d0,twelve=12.d0)
c
      dimension svint(*),props(*)
c
c variables to be defined by the user
c
c bn - axial force
c bm - bending moment
c daxial - current tangent axial stiffness
c dbend - current tangent bending stiffness
c dcoupl - tangent coupling term
c
c variables that may be updated
c
c svint - state variables for this integration point
c
c variables passed in for information
c
c eps - axial strain
c deps - incremental axial strain
c cap - curvature change
c dcap - incremental curvature change
c props - element properties
c nprops - # element properties
c nsvint - # state variables
c
c current assumption
c
c props(1) - section height
c props(2) - section width
c props(3) - Young’s modulus
c
      h=props(1)
      w=props(2)
      E=props(3)
c
c formulate linear stiffness
c
      daxial=E*h*w
      dbend=E*w*h**3/twelve
      dcoupl=zero
c
c calculate axial force and moment
c
      bn=svint(1)+daxial*deps
      bm=svint(2)+dbend*dcap
c
c store internal variables
c
      svint(1)=bn
      svint(2)=bm
      svint(3)=eps
      svint(4)=cap
c
      return
      end