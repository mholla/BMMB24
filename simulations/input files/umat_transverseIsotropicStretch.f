c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0    ! lambda_p, the amount of pre-stretch of current step
      statev(2)=1.0d0    ! J_e, elastic jacobian (elastic tangent)   
      statev(3)=0.0d0    ! Ce11
      statev(4)=0.0d0    ! Ce22 
      statev(5)=0.0d0    ! Ce33
      statev(6)=0.0d0    ! Ce12
      statev(7)=0.0d0    ! Ce13
      statev(8)=0.0d0    ! Ce23

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     #rpl,ddsddt,drplde,drpldt,
     #stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     #ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     #celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     # ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     # stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     # props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_area_morph(stress,statev,ddsdde,sse,
     #                     time,dtime,coords,props,dfgrd1,
     #                     ntens,ndi,nshr,nstatv,nprops,
     #                     noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_area_morph(stress,statev,ddsdde,sse,
     #                           time,dtime,coords,props,dfgrd1,
     #                           ntens,ndi,nshr,nstatv,nprops,
     #                           noel,npt,kstep,kinc)

c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It implements transverse isotropic areal pre-stretch, in the 
c * " plane, orthogonal to the undeformed unit normal (xn0). 
c * " The stretching is applied at the integration point of each
c * " element.
c * " It is not suitable for use with local material directions.
c * " This code had been modified for simulation of stretching
c * " within a eighth elliptical model of murine cranial dura mater. The 
c * " unit normals are with respect to an ellipsoid. More specifically,
c * " this UMAT subroutine for stretching is to simulate the
c * " "pre-stretched" state of dura mater in-vivo. The UMAT is applied
c * " to two models, one which is the pre-cut model and a second which 
c * " which is the post-cut model. The dura mater has a cut simulated by
c * " removal of tie constraints along the surface of two touching parts.
c * " Finally, the cut retraction dimensions are compared to experimental
c * " cut retraction dimensions, so that an estimate of experimental
c * " pre-stretch may be made.

c * " Written by Maria Holland in Aug 2012
c * " Modified by Jack Consolini on 2023-03-22

c...  ------------------------------------------------------------------
      implicit none

c...  STEP 1: INITIALIZATION OF ABAQUS VARIABLES & MATERIAL PARAMETERS

c...  variables to be defined
      real*8  stress(ntens),statev(nstatv),ddsdde(ntens,ntens),sse

c...  variables passed in for information
      real*8  time(2),dtime,coords(3),props(nprops),dfgrd1(3,3)
      integer ntens,ndi,nshr,nstatv,nprops,noel,npt,kstep,kinc

c...  material properties
      real*8  lam, mu, a, b, c, lpF, tsF

c...  local variables      
      integer i,j
      real*8  norm, xn0(3), xn(3)
      real*8  lambdaP, detf, detfe, lnJe, alpha
      real*8  fe(3,3), Ce(6), stretch(3), be(6), xi(6)
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/

c...  initialize material parameters 
      lam   = props(1)                ! lame constant
      mu    = props(2)                ! shear modulus
      a     = props(3)                ! ellipsoid axis in x
      b     = props(4)                ! ellipsoid axis in y
      c     = props(5)                ! ellipsoid axis in z
      lpF   = props(6)                ! Final stretch value (lambda_p final)
      tsF   = props(7)                ! Total step time for stretching

c...  STEP 2: CALCULATE UNDEFORMED ELASTIC NORMAL (normalized)
      xn0(1) = 2*(coords(1)/a**2.0)
      xn0(2) = 2*(coords(2)/b**2.0)
      xn0(3) = 2*(coords(3)/c**2.0)
      norm = sqrt(xn0(1)*xn0(1) + xn0(2)*xn0(2) + xn0(3)*xn0(3))
      xn0 = xn0/norm
            
c...  STEP 3: CALCULATE DEFORMED ELASTIC NORMAL (normalized)
      xn(1) = (dfgrd1(1,1)*xn0(1) + dfgrd1(1,2)*xn0(2) + dfgrd1(1,3)*xn0(3))
      xn(2) = (dfgrd1(2,1)*xn0(1) + dfgrd1(2,2)*xn0(2) + dfgrd1(2,3)*xn0(3))
      xn(3) = (dfgrd1(3,1)*xn0(1) + dfgrd1(3,2)*xn0(2) + dfgrd1(3,3)*xn0(3))

c...  STEP 4: UPDATE GROWTH FACTOR
      alpha = (lpF-1.0)/tsF
      lambdaP = alpha*(time(2)+dtime) + 1.0d0
      statev(1) = lambdaP 
      
c...  STEP 5: CALCULATE ELASTIC TENSOR (Fe = F\cdotF^p = F\cdot[lambda^p_i*[I-(n \otimes n_0)] 
c...  +1/(lambda^p_i)^2*(n \otimes n_0)]
      fe(1,1) = lambdaP*dfgrd1(1,1) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(1)*xn0(1)
      fe(1,2) = lambdaP*dfgrd1(1,2) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(1)*xn0(2)
      fe(1,3) = lambdaP*dfgrd1(1,3) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(1)*xn0(3)
      fe(2,1) = lambdaP*dfgrd1(2,1) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(2)*xn0(1)
      fe(2,2) = lambdaP*dfgrd1(2,2) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(2)*xn0(2)
      fe(2,3) = lambdaP*dfgrd1(2,3) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(2)*xn0(3)
      fe(3,1) = lambdaP*dfgrd1(3,1) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(3)*xn0(1)
      fe(3,2) = lambdaP*dfgrd1(3,2) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(3)*xn0(2)
      fe(3,3) = lambdaP*dfgrd1(3,3) + ((1/(lambdaP*lambdaP))-lambdaP)*xn(3)*xn0(3)

c...  STEP 6: CALCULATE DETERMINATE OF ELASTIC TENSOR (J_e = det(Fe))
      detf=+dfgrd1(1,1) * (dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #     -dfgrd1(1,2) * (dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #     +dfgrd1(1,3) * (dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))
     
      detfe = +fe(1,1) * (fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2) * (fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3) * (fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
     
      statev(2) = detfe
      lnJe = dlog(detfe)

c...  STEP 7: CALCULATE ELASTIC RIGHT-CAUCHY GREEN DEFORMATION TENSOR (Ce = fe^t * fe)
      Ce(1) = fe(1,1)*fe(1,1) + fe(2,1)*fe(2,1) + fe(3,1)*fe(3,1) ! Ce11
      Ce(2) = fe(1,2)*fe(1,2) + fe(2,2)*fe(2,2) + fe(3,2)*fe(3,2) ! Ce22
      Ce(3) = fe(1,3)*fe(1,3) + fe(2,3)*fe(2,3) + fe(3,3)*fe(3,3) ! Ce33
      Ce(4) = fe(1,1)*fe(1,2) + fe(2,1)*fe(2,2) + fe(3,1)*fe(3,2) ! Ce12, Ce21
      Ce(5) = fe(1,1)*fe(1,3) + fe(2,1)*fe(2,3) + fe(3,1)*fe(3,3) ! Ce13, Ce31
      Ce(6) = fe(1,2)*fe(1,3) + fe(2,2)*fe(2,3) + fe(3,2)*fe(3,3) ! Ce23, Ce32

      statev(3) = Ce(1)
      statev(4) = Ce(2)
      statev(5) = Ce(3)
      statev(6) = Ce(4)
      statev(7) = Ce(5)
      statev(8) = Ce(6)
     
c...  STEP 8: CALCULATE ELASTIC LEFT-CAUCHY GREEN DEFORMATION TENSOR (be = fe * fe^t)
      be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3) ! be11
      be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3) ! be22
      be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3) ! be33
      be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3) ! be12, be21
      be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3) ! be13, be31
      be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3) ! be23, be32

c...  STEP 9: CALCULATE CAUCHY STRESS
      do i = 1,ntens
        stress(i) = ((lam*lnJe-mu)*xi(i) + mu*be(i))/detfe
      enddo     

c...  STEP 10: CALCULATE ELASTIC AND GEOMETRIC TANGENT
      ddsdde(1,1) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(1)
      ddsdde(2,2) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(2)
      ddsdde(3,3) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(3)
      ddsdde(1,2) = (lam)/detfe
      ddsdde(1,3) = (lam)/detfe
      ddsdde(2,3) = (lam)/detfe
      ddsdde(1,4) = stress(4)
      ddsdde(2,4) = stress(4)
      ddsdde(3,4) = 0.d0
      ddsdde(4,4) = -(lam*lnJe - mu)/detfe + (stress(1) + stress(2))/2.d0
      
      if (ntens.eq.6) then
        ddsdde(1,5) = stress(5)
        ddsdde(2,5) = 0.d0
        ddsdde(3,5) = stress(5)
        ddsdde(1,6) = 0.d0
        ddsdde(2,6) = stress(6)
        ddsdde(3,6) = stress(6)
        ddsdde(5,5) = -(lam*lnJe - mu)/detfe + (stress(1) + stress(3))/2.d0
        ddsdde(6,6) = -(lam*lnJe - mu)/detfe + (stress(2) + stress(3))/2.d0
        ddsdde(4,5) = stress(6)/2.d0
        ddsdde(4,6) = stress(5)/2.d0
        ddsdde(5,6) = stress(4)/2.d0
      endif
    
c...  STEP 11: USE SYMMETRY TO FILL IN REMAINDER OF TANGENT TENSORS
      do i = 2,ntens
        do j = 1,i-1
          ddsdde(i,j) = ddsdde(j,i)
        enddo
      enddo    

c...  STEP 12: CALCULATE STRAIN ENERGY
      sse = (lam*lnJe**2.d0 + mu*(be(1)+be(2)+be(3) - 3.d0 - 2.d0*lnJe))/2.d0

      return
      end