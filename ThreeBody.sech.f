c     234567890
      program ThreeBody

      integer LegPoints,xNumPoints,yNumPoints,NLeft,NRight
      integer Type,J,Parity,NumStates,NumBound,NumStateInc,NumPrintStates,Order,Left,Right,Bottom,Top
      integer RSteps,RadialCouplingFlag,scount,nsol
      double precision m1,m2,m3
      double precision DcoefAB,DcoefAA,RscaleAB,RscaleAA
      double precision RDerivDelt,RFirst,RLast,XFirst,XLast,StepX,Rchange,r0xfactor,r0yfactor
      double precision xMin,xMax,yMin,yMax
      double precision, allocatable :: R(:)
      double precision, allocatable :: xPoints(:),yPoints(:)
      double precision Shift,Shift2

      integer k,l
      integer LeadDim,MatrixDim,HalfBandWidth
      integer xDim,yDim
      integer, allocatable :: xBounds(:),yBounds(:)
      double precision TotalMemory
      double precision M,mu,d31,phi31,d23,phi23,d12,phi12
      double precision RLeft,RRight
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),ux(:,:,:),v(:,:,:),vy(:,:,:)
      double precision, allocatable :: up(:,:),vp(:,:)
      double precision, allocatable :: P(:,:),Q(:,:)
      double precision, allocatable :: S(:,:),H(:,:)
      double precision, allocatable :: LUFac(:,:)
      double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:)
      double precision, allocatable :: PlotPsi(:,:),Energies(:,:)

      character*64 LegendreFile

      double precision ur(1:50000),acoef,bcoef,diff
      double precision sec,time,Rinitial,secp,timep
      integer ntrunc,ncount
      common /Rvalue/ Rvalue

      ncount = 0
      ntrunc = 0

C     TIME INICIALIZATION
      call cpu_time(time)
      SEC = time
      secp = time

c     read in number of energies and states to print
      read(5,*)
      read(5,*) NumStates,NumBound,NumStateInc,NumPrintStates,RadialCouplingFlag

      read(5,*)
      read(5,*)
      read(5,*) Type,J,Parity
      read(5,*)
      read(5,*)
      read(5,*) NLeft,NRight

      if ((Parity.ne.1).or.(J.gt.1)) then
         write(*,*)'Wrong Symmetry selected'
         write(*,*)'Program stoped.'
         stop
      endif


c     read in Gauss-Legendre info
      read(5,*)
      read(5,*)
      read(5,1002) LegendreFile

      read(5,*)
      read(5,*)
      read(5,*) LegPoints

c     read in boundary conditions
      read(5,*)
      read(5,*)
      read(5,*) Shift,Shift2,Order

c     read in masses
      read(5,*)
      read(5,*)
      read(5,*) m1,m2,m3

c     read in potential coefficients
      read(5,*)
      read(5,*)
      read(5,*)DcoefAB,DcoefAA,RscaleAB,RscaleAA

c     read in grid information
      read(5,*)
      read(5,*)
      read(5,*) xNumPoints,yNumPoints,r0xfactor,r0yfactor

      xMin = 0.0d0              ! units of Pi
      xMax = 1.0d0              ! units of Pi
      Left = NLeft              ! phi (x) boundary condition at 0
      Right = NRight            ! phi (x) boundary condition at pi

      yMin = 0.0d0              ! units of Pi
      yMax = 0.5d0              ! units of Pi
      Bottom = 2                ! theta (y) boundary condition at 0 -> no BC imposed (function is finite)
      Top = 2                   ! theta (y) boundary condition at pi/2 -> no BC imposed (function is finite)

      read(5,*)
      read(5,*)
      read(5,*) RSteps,RDerivDelt,RFirst,RLast,Rchange
      allocate(R(RSteps))
      Rinitial = RFirst
      if(RSteps .eq. 1) then
         R(1) = RFirst
         goto 55
      endif
      XFirst = RFirst**(1.0d0/3.0d0)
      XLast = RLast**(1.0d0/3.0d0)
      StepX = (XLast-XFirst)/dfloat(RSteps-1)
      do iR = 1,RSteps
         R(iR) = (XFirst+(iR-1)*StepX)**3
      enddo
 55   continue


      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

      xDim = xNumPoints+Order-3
      if (Left .eq. 2) xDim = xDim + 1
      if (Right .eq. 2) xDim = xDim + 1
      yDim = yNumPoints+Order-3
      if (Top .eq. 2) yDim = yDim + 1
      if (Bottom .eq. 2) yDim = yDim + 1

      MatrixDim = xDim*yDim
      HalfBandWidth = yDim*Order+Order
      LeadDim = 3*HalfBandWidth+1

      TotalMemory = 2.0d0*(HalfBandWidth+1)*MatrixDim ! S, H
      TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*xDim ! x splines
      TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*yDim ! y splines
      TotalMemory = TotalMemory + 2.0d0*NumStates*NumStates ! P and Q matrices
      TotalMemory = TotalMemory + LeadDim*MatrixDim ! LUFac
      TotalMemory = TotalMemory + 4.0d0*NumStates*MatrixDim ! channel functions
      TotalMemory = TotalMemory + 4*xDim*xDim ! (CalcHamiltonian)
      TotalMemory = TotalMemory + 4*yDim*yDim ! (CalcHamiltonian)
      TotalMemory = TotalMemory + LegPoints**2*yNumPoints*xNumPoints ! (CalcHamiltonian)
      TotalMemory = 8.0d0*TotalMemory/(1024.0d0*1024.0d0)

!write(6,*)
!write(6,*) 'MatrixDim ',MatrixDim
!write(6,*) 'HalfBandWidth ',HalfBandWidth
!write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
!write(6,*)

      M = m1+m2+m3
      mu = dsqrt(m1*m2*m3/M)

      phi12 =  2.0d0*datan(m3/mu)
      phi23 =  0.0d0
      phi31 = -2.0d0*datan(m2/mu)

      d12   = dsqrt(m3/mu*(1.0d0-m3/M))
      d23   = dsqrt(m1/mu*(1.0d0-m1/M))
      d31   = dsqrt(m2/mu*(1.0d0-m2/M))

!write(6,*) 'mu =',mu

      allocate(xPoints(xNumPoints),yPoints(yNumPoints))
      allocate(xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order))
      allocate(u(LegPoints,Order+2,xDim),ux(LegPoints,Order+2,xDim))
      allocate(v(LegPoints,Order+2,yDim),vy(LegPoints,Order+2,yDim))
      allocate(up(xDim,xNumPoints),vp(yDim,yNumPoints))
      allocate(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,MatrixDim))
      allocate(P(NumStates,NumStates),Q(NumStates,NumStates))

      allocate(LUFac(LeadDim,MatrixDim))
      allocate(lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates),Energies(NumStates,2))

      open(100,File='fort.100')!,STATUS='NEW') !,BLOCKSIZE=512)
      open(99, File='fort.99')!, STATUS='NEW') !,BLOCKSIZE=512)
      open(98, File='fort.98')!, STATUS='NEW') !,BLOCKSIZE=512)

      do iR = 1,RSteps

         Rvalue = R(iR)

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Choosing Rchange (jpdincao)  

c     if (R(iR) .lt. Rchange) then
c     NumFirst = NumStates
c     else
c     NumFirst = NumBound
c     endif
c     !write(6,*) R(iR)
c     !write(6,*) 'Shift ',Shift

         NumFirst = NumStates
c     if (iR.ge.2) then
         if (R(iR).ge.Rchange) then
c     if (ncount.eq.0) then   
c     Rchange = R(iR) 
c     ncount = ncount+1
c     endif
            NumFirst = NumBound
         endif
c     endif   

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     call GridMaker(R(iR),xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints,d31,RscaleAB,RscaleAA,r0xfactor,r0yfactor)

c         write(6,*) 'Calling GridMakerNew'
         call GridMakerNEW(R(iR),xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints,d12,d31,d23,
     .        phi12,phi31,phi23,RscaleAB,RscaleAA,r0xfactor,r0yfactor)
c         write(6,*) 'Calculating Basis Functions'
         call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
         call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,1,ux)
         call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,yDim,yBounds,yNumPoints,0,v)
         call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,yDim,yBounds,yNumPoints,1,vy)

c         write(6,*) 'Calculating the overlap S ....'
         call CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,
     >        xNumPoints,yNumPoints,u,v,xBounds,yBounds,MatrixDim,HalfBandWidth,S)
c         write(6,*) 'Done.'

         if (RadialCouplingFlag .ne. 0) then

            RLeft = 0.99d0*R(iR) !R(iR)-RDerivDelt
            call CalcHamiltonian(RLeft,J,mu,d12,phi12,d23,phi23,
     >           d31,phi31,DcoefAB,DcoefAA,RscaleAB,RscaleAA,Order,xPoints,yPoints,LegPoints,xLeg,
     >           wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,ux,vy,xBounds,yBounds,
     >           MatrixDim,HalfBandWidth,H)
            call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,lPsi,Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
            if (iR .gt. 1) call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,mPsi,lPsi)
            do k = 1,NumStates
               Energies(k,1) = Energies(k,1) + 1.875d0/(mu*RLeft*RLeft)
c     Energies(k,1) = Energies(k,1) - 0.125d0/(mu*RLeft*RLeft)
            enddo
!write(6,*)
!write(6,*) RLeft
            do k = 1,NumStates
!write(6,*) k,Energies(k,1),Energies(k,2)
            enddo

            RRight = 1.01d0*R(iR) !R(iR)+RDerivDelt
            call CalcHamiltonian(RRight,J,mu,d12,phi12,d23,phi23,
     >           d31,phi31,DcoefAB,DcoefAA,RscaleAB,RscaleAA,Order,xPoints,yPoints,LegPoints,xLeg,
     >           wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,ux,vy,xBounds,yBounds,
     >           MatrixDim,HalfBandWidth,H)
            call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,rPsi,Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,lPsi,rPsi)
            do k = 1,NumStates
               Energies(k,1) = Energies(k,1) + 1.875d0/(mu*RRight*RRight)
c     Energies(k,1) = Energies(k,1) - 0.125d0/(mu*RRight*RRight)
            enddo
!write(6,*)
!write(6,*) RRight
            do k = 1,NumStates
!write(6,*) k,Energies(k,1),Energies(k,2)
            enddo
         endif

c     mu = 0.5d0

!write(6,*) 'H ....'
         call CalcHamiltonian(R(iR),J,mu,d12,phi12,d23,phi23,
     >        d31,phi31,DcoefAB,DcoefAA,RscaleAB,RscaleAA,Order,xPoints,yPoints,LegPoints,xLeg,
     >        wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,ux,vy,xBounds,yBounds,
     >        MatrixDim,HalfBandWidth,H)
!write(6,*) 'Done.'

         call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,mPsi,Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
         if (RadialCouplingFlag .ne. 0) call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,rPsi,mPsi)
         do k = 1,NumStates
            Energies(k,1) = Energies(k,1) + 1.875d0/(mu*R(iR)*R(iR))
         enddo
!write(6,*)
!write(6,*) R(iR)
         do k = 1,NumStates
!write(6,100) k,Energies(k,1),Energies(k,2)
         enddo

         write(100,11)R(iR),(Energies(k,1),k=1,NumStates)
         

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Adjusting Shift
         ur(iR) = Energies(1,1)
         if (iR.ge.2) then
            if ((R(iR+1).lt.Rchange).and.((ur(iR)-ur(iR-1)).gt.0.d0)) then
               Shift = Energies(1,1)
               if (Energies(1,1).gt.0.d0) Shift = 0.1d0*Shift
               if (Energies(1,1).lt.0.d0) Shift = 10.d0*Shift
            endif
            if (R(iR+1).gt.Rchange) then
c     First group of energies
               Shift = Energies(1,1)
               if (Energies(1,1).gt.0.d0) then 
                  if (Energies(1,1).lt.1.d-10) then
                     Shift = -1.d-9
                  else
                     Shift = -10.d0*Shift  
c     if (Shift.lt.1.d-10) Shift = -1.d-10
                  endif
               endif   
               if (Energies(1,1).lt.0.d0) Shift = 10.d0*Shift
c     if ((Energies(1,1).lt.0.d0)) then
c     if ((R(iR).ge.300.d0).and.(ncount.le.10)) then
c     Shift = -1.d-9
c     ncount = ncount+1
c     else
c     Shift = 10.d0*Shift
c     endif
c     endif
               
c     Second group of energies
               Shift2 = Energies(NumBound+1,1)
               if (Energies(NumBound+1,1).gt.0.d0) Shift2 = 0.01d0*Shift2
               if (Energies(NumBound+1,1).lt.0.d0) Shift2 = 100.0d0*Shift2
               if (Shift2.le.Energies(NumBound,1)) Shift2 = Energies(NumBound,1)
            endif
         endif   
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     energies = energies-15.d0/4.d0/2.d0/mu/R(iR)**2.d0 !jpdincao
c     energies = energies*2.d0*mu

         call cpu_time(timep)
         write(*,11)R(iR),(Energies(k,1),k=1,NumStates),Shift,Shift2,timep-secp
         call cpu_time(timep)
         secp = timep 

         if (RadialCouplingFlag .ne. 0) then
            RDerivDelt = (RRight-RLeft)/2.d0
            call CalcRadialCoupling(NumStates,HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,P,Q)
            write(98,*) R(iR)
            write(99,*) R(iR)
            do k = 1,NumStates
               write(98,20) (P(k,l), l = 1,NumStates)
               write(99,20) (Q(k,l), l = 1,NumStates)
            enddo
         endif

         if (NumPrintStates .ne. 0) then
            do k = 1,xNumPoints
               write(9997,*) xPoints(k)
            enddo
            do k = 1,yNumPoints
               write(9998,*) yPoints(k)
            enddo
            do k = 1,MatrixDim
c     write(9999+iR,20) (mPsi(k,Index(l)), l = 1,NumPrintStates)
               write(9999+iR,20) (mPsi(k,(l)), l = 1,NumPrintStates)
            enddo
            close(unit=9999+iR)
         endif

      enddo

      deallocate(S,H)
      deallocate(lPsi,mPsi,rPsi,Energies)
      deallocate(LUFac)
      deallocate(xPoints,yPoints)
      deallocate(xLeg,wLeg)
      deallocate(xBounds,yBounds)
      deallocate(u,ux)
      deallocate(v,vy)
      deallocate(P,Q)

 10   format(f9.3,1P,100e25.15)
 11   format(f12.3,1P,100e20.12)
 111  format(f12.3,1P,100I5)
 20   format(1P,100e20.12)
 30   format(1P,100e25.15)
 100  format(i4,1P,10e25.15)
 1002 format(a64)

      call cpu_time(time)
      WRITE(*,*)'>>> EXECUTION TIME =',time-SEC

      close(100)
      close(99)
      close(98)

      stop
      end

      subroutine CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,
     >     xNumPoints,yNumPoints,u,v,xBounds,yBounds,MatrixDim,HalfBandWidth,S)

      integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints
      integer xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order)
      integer MatrixDim,HalfBandWidth
      double precision xPoints(*),yPoints(*),xLeg(LegPoints),wLeg(LegPoints)
      double precision u(LegPoints,Order+2,xDim)
      double precision v(LegPoints,Order+2,yDim)
      double precision S(HalfBandWidth+1,MatrixDim)

      integer ix,iy,ixp,iyp,kx,ky,lx,ly
      integer i1,i1p,i2,i2p,k0,kp0
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)
      double precision ax,bx,xTempS
      double precision y,ay,by,yTempS,yScaledZero
      double precision, allocatable :: xIntScale(:),xS(:,:)
      double precision, allocatable :: yIntScale(:),yS(:,:),sin2y(:,:)

      allocate(xIntScale(xNumPoints),xS(xDim,xDim))
      allocate(yIntScale(yNumPoints),yS(yDim,yDim))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(kyMin(yDim,yDim),kyMax(yDim,yDim),sin2y(LegPoints,yNumPoints))

      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
      enddo
      do ky = 1,yNumPoints-1
         ay = yPoints(ky)
         by = yPoints(ky+1)
         yIntScale(ky) = 0.5d0*(by-ay)
         yScaledZero = 0.5d0*(by+ay)
         do ly = 1,LegPoints
            y = yIntScale(ky)*xLeg(ly)+yScaledZero
            sin2y(ly,ky) = dsin(2.0d0*y)
         enddo
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      do iy = 1,yDim
         do iyp = 1,yDim
            kyMin(iyp,iy) = max(yBounds(iy),yBounds(iyp))
            kyMax(iyp,iy) = min(yBounds(iy+Order+1),yBounds(iyp+Order+1))-1
         enddo
      enddo

      xS = 0.0d0
      do ix = 1,xDim
         k0 = xBounds(ix)-1
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            kp0 = xBounds(ixp)-1
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               xTempS = 0.0d0
               do lx = 1,LegPoints
                  xTempS = xTempS + wLeg(lx)*u(lx,kx-k0,ix)*u(lx,kx-kp0,ixp)
               enddo
               xS(ixp,ix) = xS(ixp,ix) + xIntScale(kx)*xTempS
            enddo
         enddo
      enddo

      yS = 0.0d0
      do iy = 1,yDim
         k0 = yBounds(iy)-1
         do iyp = max(1,iy-Order),min(yDim,iy+Order)
            kp0 = yBounds(iyp)-1
            do ky = kyMin(iyp,iy),kyMax(iyp,iy)
               yTempS = 0.0d0
               do ly = 1,LegPoints
                  yTempS = yTempS + wLeg(ly)*v(ly,ky-k0,iy)*v(ly,ky-kp0,iyp)*sin2y(ly,ky)
               enddo
               yS(iyp,iy) = yS(iyp,iy) + yIntScale(ky)*yTempS
            enddo
         enddo
      enddo

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            do iy = 1,yDim
               i2 = i1+iy-1
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  i2p = i1p+iyp-1
                  Row = i2+1
                  Col = i2p+1
                  if (Col .ge. Row) then
                     NewRow = HalfBandWidth+1+Row-Col
                     S(NewRow,Col) = xS(ixp,ix)*yS(iyp,iy)
                  endif
               enddo
            enddo
         enddo
      enddo

      deallocate(xIntScale,xS)
      deallocate(yIntScale,yS,sin2y)
      deallocate(kxMin,kxMax)
      deallocate(kyMin,kyMax)

      return
      end

      subroutine CalcHamiltonian(R,J,mu,d12,phi12,d23,phi23,
     >     d31,phi31,DcoefAB,DcoefAA,RscaleAB,RscaleAA,Order,xPoints,yPoints,LegPoints,xLeg,wLeg,
     >     xDim,yDim,xNumPoints,yNumPoints,u,v,ux,vy,xBounds,yBounds,
     >     MatrixDim,HalfBandWidth,H)

      integer J,Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints
      integer xBounds(yNumPoints+2*Order),yBounds(yNumPoints+2*Order),MatrixDim,HalfBandWidth
      double precision R,mu,d12,phi12,d23,phi23,d31,phi31,DcoefAB,DcoefAA,RscaleAB,RscaleAA
      double precision xPoints(xNumPoints),yPoints(yNumPoints),xLeg(LegPoints),wLeg(LegPoints)
      double precision u(LegPoints,Order+2,xDim),ux(LegPoints,Order+2,xDim)
      double precision v(LegPoints,Order+2,yDim),vy(LegPoints,Order+2,yDim)
      double precision H(HalfBandWidth+1,MatrixDim)

      integer ix,iy,ixp,iyp,kx,ky,lx,ly
      integer i1,i1p,i2,i2p,kx0,kxp0,ky0,kyp0,k,kp
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)
      double precision Pi
      double precision a,b,c,r12,r23,r31
      double precision m,Rall,V12,V23,V31
      double precision VInt,VTempInt
      double precision x,ax,bx,xScaledZero,xTempS,xTempD1
      double precision y,ay,by,yScaledZero,yTempD0,yTempD1,yTempD2
      double precision Vpot
      double precision, allocatable :: TempPot(:,:,:,:)
      double precision, allocatable :: xIntScale(:),xTempV(:),xS(:,:),xD1(:,:)
      double precision, allocatable :: sinx12(:,:),sinx23(:,:),sinx31(:,:)
      double precision, allocatable :: yIntScale(:),yTempV(:,:),yD0(:,:),yD1(:,:),yD2(:,:)
      double precision, allocatable :: siny(:,:),sin2y(:,:),siny2(:,:),cosy(:,:),cosy2(:,:),cosy3(:,:)

      allocate(xIntScale(xNumPoints),xTempV(LegPoints),xS(xDim,xDim),xD1(xDim,xDim))
      allocate(sinx12(LegPoints,xNumPoints),sinx23(LegPoints,xNumPoints),sinx31(LegPoints,xNumPoints))
      allocate(yIntScale(yNumPoints),yTempV(LegPoints,yNumPoints))
      allocate(yD0(yDim,yDim),yD1(yDim,yDim),yD2(yDim,yDim))
      allocate(siny(LegPoints,yNumPoints),sin2y(LegPoints,yNumPoints),siny2(LegPoints,yNumPoints))
      allocate(cosy(LegPoints,yNumPoints),cosy2(LegPoints,yNumPoints),cosy3(LegPoints,yNumPoints))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(kyMin(yDim,yDim),kyMax(yDim,yDim))
      allocate(TempPot(LegPoints,LegPoints,xNumPoints,yNumPoints))

      Pi = 3.14159265358979323846d0
      H = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            x = xIntScale(kx)*xLeg(lx)+xScaledZero
c     sinx12(lx,kx) = dsin(x-Pi/6.0d0+phi12)
c     sinx23(lx,kx) = dsin(x-Pi/6.0d0+phi23)    ! phi23 will usually be zero
c     sinx31(lx,kx) = dsin(x-Pi/6.0d0+phi31)
            sinx12(lx,kx) = dsin(x+Pi/2.0d0+phi12)
            sinx23(lx,kx) = dsin(x+Pi/2.0d0+phi23) ! phi23 will usually be zero
            sinx31(lx,kx) = dsin(x+Pi/2.0d0+phi31)
         enddo
      enddo

      do ky = 1,yNumPoints-1
         ay = yPoints(ky)
         by = yPoints(ky+1)
         yIntScale(ky) = 0.5d0*(by-ay)
         yScaledZero = 0.5d0*(by+ay)
         do ly = 1,LegPoints
            y = yIntScale(ky)*xLeg(ly)+yScaledZero
            siny(ly,ky) = dsin(y)
            sin2y(ly,ky) = dsin(2.0d0*y)
            siny2(ly,ky) = siny(ly,ky)*siny(ly,ky)
            cosy(ly,ky) = dcos(y)
            cosy2(ly,ky) = cosy(ly,ky)*cosy(ly,ky)
            cosy3(ly,ky) = cosy2(ly,ky)*cosy(ly,ky)
         enddo
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      do iy = 1,yDim
         do iyp = 1,yDim
            kyMin(iyp,iy) = max(yBounds(iy),yBounds(iyp))
            kyMax(iyp,iy) = min(yBounds(iy+Order+1),yBounds(iyp+Order+1))-1
         enddo
      enddo

      xS = 0.0d0
      xD1 = 0.0d0
      do ix = 1,xDim
         kx0 = xBounds(ix)-1
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            kxp0 = xBounds(ixp)-1
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               k = kx-kx0
               kp = kx-kxp0
               xTempS = 0.0d0
               xTempD1 = 0.0d0
               do lx = 1,LegPoints
                  a = wLeg(lx)
                  b = a*u(lx,k,ix)*u(lx,kp,ixp)
                  c = a*ux(lx,kp,ixp)
                  xTempS = xTempS + b
                  xTempD1 = xTempD1 + c*ux(lx,k,ix)
               enddo
               xS(ixp,ix) = xS(ixp,ix) + xIntScale(kx)*xTempS
               xD1(ixp,ix) = xD1(ixp,ix) + xIntScale(kx)*xTempD1
            enddo
         enddo
      enddo

      yD0 = 0.0d0
      yD1 = 0.0d0
      yD2 = 0.0d0
      do iy = 1,yDim
         ky0 = yBounds(iy)-1
         do iyp = max(1,iy-Order),min(yDim,iy+Order)
            kyp0 = yBounds(iyp)-1
            do ky = kyMin(iyp,iy),kyMax(iyp,iy)
               k = ky-ky0
               kp = ky-kyp0
               yTempD0 = 0.0d0
               yTempD1 = 0.0d0
               yTempD2 = 0.0d0
               do ly = 1,LegPoints
                  a = wLeg(ly)
                  b = a*v(ly,k,iy)*v(ly,kp,iyp)
                  yTempD0 = yTempD0 + a*vy(ly,k,iy)*vy(ly,kp,iyp)*sin2y(ly,ky)
                  yTempD1 = yTempD1 + b*cosy(ly,ky)/siny(ly,ky)
                  yTempD2 = yTempD2 + b*siny(ly,ky)/cosy(ly,ky)
               enddo
               yD0(iyp,iy) = yD0(iyp,iy) +       yIntScale(ky)*yTempD0
               yD1(iyp,iy) = yD1(iyp,iy) + 2.0d0*yIntScale(ky)*yTempD1
               yD2(iyp,iy) = yD2(iyp,iy) + 2.0d0*yIntScale(ky)*yTempD2
            enddo
         enddo
      enddo

      m = 1.0d0/(mu*R*R)

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            do iy = 1,yDim
               i2 = i1+iy-1
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  i2p = i1p+iyp-1
                  Row = i2+1
                  Col = i2p+1
                  if (Col .ge. Row) then
                     NewRow = HalfBandWidth+1+Row-Col
                     H(NewRow,Col) = 2.0d0*m*(xS(ixp,ix)*yD0(iyp,iy)+xD1(ixp,ix)*yD1(iyp,iy))
     >                    +m*xS(ixp,ix)*dfloat(J*(J+1))*yD2(iyp,iy)
                  endif
               enddo
            enddo
         enddo
      enddo

c     if potential integral is not separable, use the following code section
c     to do 2D integrals
c     assume mu=sqrt(mu1*mu2)

      Rall = R/dsqrt(2.0d0)

      do ky = 1,yNumPoints-1
         do kx = 1,xNumPoints-1
            do ly = 1,LegPoints
               do lx = 1,LegPoints
                  r12 = d12*Rall*dsqrt(1.0d0+siny(ly,ky)*sinx12(lx,kx))
                  V12 = Vpot(r12,DcoefAB,RscaleAB)
                  r23 = d23*Rall*dsqrt(1.0d0+siny(ly,ky)*sinx23(lx,kx))
                  V23 = Vpot(r23,DcoefAA,RscaleAA)
                  r31 = d31*Rall*dsqrt(1.0d0+siny(ly,ky)*sinx31(lx,kx))
                  V31 = Vpot(r31,DcoefAB,RscaleAB)
                  TempPot(lx,ly,kx,ky) = V12+V23+V31
               enddo
            enddo
         enddo
      enddo

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         kx0 = xBounds(ix)-1
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            kxp0 = xBounds(ixp)-1
            do iy = 1,yDim
               i2 = i1+iy-1
               ky0 = yBounds(iy)-1
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  kyp0 = yBounds(iyp)-1
                  i2p = i1p+iyp-1
                  Row = i2+1
                  Col = i2p+1
                  if (Col .ge. Row) then
                     NewRow = HalfBandWidth+1+Row-Col
                     VInt = 0.0d0
                     do ky = kyMin(iyp,iy),kyMax(iyp,iy)
                        do ly = 1,LegPoints
                           yTempV(ly,ky) = wLeg(ly)*sin2y(ly,ky)*v(ly,ky-ky0,iy)*v(ly,ky-kyp0,iyp)
                        enddo
                     enddo
                     do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                        do lx = 1,LegPoints
                           xTempV(lx) = wLeg(lx)*u(lx,kx-kx0,ix)*u(lx,kx-kxp0,ixp)
                        enddo
                        do ky = kyMin(iyp,iy),kyMax(iyp,iy)
                           VTempInt = 0.0d0
                           do ly = 1,LegPoints
                              do lx = 1,LegPoints
                                 VTempInt = VTempInt + xTempV(lx)*yTempV(ly,ky)*TempPot(lx,ly,kx,ky)
                              enddo
                           enddo
                           VInt = VInt + xIntScale(kx)*yIntScale(ky)*VTempInt
                        enddo
                     enddo

                     H(NewRow,Col) = H(NewRow,Col)+VInt

                  endif
               enddo
            enddo
         enddo
      enddo

      deallocate(TempPot)
      deallocate(xIntScale,xTempV,xS,xD1)
      deallocate(sinx12,sinx23,sinx31)
      deallocate(yIntScale,yTempV,yD0,yD1)
      deallocate(siny,sin2y,siny2,cosy,cosy2,cosy3)
      deallocate(kxMin,kxMax)
      deallocate(kyMin,kyMax)

      return
      end

      double precision function Vpot(r,Dcoef,Rscale)
      double precision r,Dcoef,Rscale,x,cutoff

      x = r/Rscale
      Vpot = Dcoef/dcosh(r/Rscale)**2

      return
      end

      subroutine CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,MatrixDim,xBounds,xNumPoints,Deriv,u)

      integer Left,Right,Order,LegPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
      double precision xPoints(*),xLeg(*)
      double precision u(LegPoints,Order+2,MatrixDim)

      integer i,k,l,Count
      integer, allocatable :: t(:)
      double precision BSpline
      double precision ax,bx,xIntScale,xScaledZero
      double precision, allocatable :: x(:,:),tx(:),b(:)

      allocate(t(xNumPoints+2*Order),tx(xNumPoints+2*Order),b(xNumPoints+Order))
      allocate(x(LegPoints,xNumPoints))

      do i = 1,Order
         t(i) = 1
         tx(i) = xPoints(1)
      enddo
      do i = 1,xNumPoints
         t(i+Order) = i
         tx(i+Order) = xPoints(i)
      enddo
      do i = 1,Order
         t(i+Order+xNumPoints) = xNumPoints
         tx(i+Order+xNumPoints) = xPoints(xNumPoints)
      enddo

      select case (Left)
      case (0:1)
         select case (Right)
         case (0:1)
            do i = 2,xNumPoints+2*Order-1
               xBounds(i-1) = t(i)
            enddo
         case (2)
            do i = 2,xNumPoints+2*Order
               xBounds(i-1) = t(i)
            enddo
         end select
      case (2)
         select case (Right)
         case (0:1)
            do i = 1,xNumPoints+2*Order-1
               xBounds(i) = t(i)
            enddo
         case (2)
            do i = 1,xNumPoints+2*Order
               xBounds(i) = t(i)
            enddo
         end select
      end select

      deallocate(t)

      do k = 1,xNumPoints-1
         ax = xPoints(k)
         bx = xPoints(k+1)
         xIntScale = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do l = 1,LegPoints
            x(l,k) = xIntScale*xLeg(l)+xScaledZero
         enddo
      enddo

      u = 0.0d0

      Count = 1
      select case (Left)
      case (0)
         do k = xBounds(Count),xBounds(Count+Order+1)-1
            do l = 1,LegPoints
               u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,2,tx,b,x(l,k))
            enddo
         enddo
         Count = Count + 1
      case (1)
         do k = xBounds(Count),xBounds(Count+Order+1)-1
            do l = 1,LegPoints
               u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,1,tx,b,x(l,k))+
     >              BSpline(Order,Deriv,xNumPoints,2,tx,b,x(l,k))
            enddo
         enddo
         Count = Count + 1
      case(2)
         do i = 1,2
            do k = xBounds(Count),xBounds(Count+Order+1)-1
               do l = 1,LegPoints
                  u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,i,tx,b,x(l,k))
               enddo
            enddo
            Count = Count + 1
         enddo
      end select

      do i = 3,xNumPoints+Order-3
         do k = xBounds(Count),xBounds(Count+Order+1)-1
            do l = 1,LegPoints
               u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,i,tx,b,x(l,k))
            enddo
         enddo
         Count = Count + 1
      enddo

      select case (Right)
      case (0)
         do k = xBounds(Count),xBounds(Count+Order+1)-1
            do l = 1,LegPoints
               u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,xNumPoints+Order-2,tx,b,x(l,k))
            enddo
         enddo
      case (1)
         do k = xBounds(Count),xBounds(Count+Order+1)-1
            do l = 1,LegPoints
               u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,xNumPoints+Order-2,tx,b,x(l,k))+
     >              BSpline(Order,Deriv,xNumPoints,xNumPoints+Order-1,tx,b,x(l,k))
            enddo
         enddo
      case(2)
         do i = xNumPoints+Order-2,xNumPoints+Order-1
            do k = xBounds(Count),xBounds(Count+Order+1)-1
               do l = 1,LegPoints
                  u(l,k-xBounds(Count)+1,Count) = BSpline(Order,Deriv,xNumPoints,i,tx,b,x(l,k))
               enddo
            enddo
            Count = Count + 1
         enddo
      end select

      deallocate(x,tx,b)

      return
      end

      double precision function BSpline(Order,Deriv,xNumPoints,n,t,b,x)

      integer Order,Deriv,xNumPoints,n
      double precision t(*),b(xNumPoints+Order),x

      double precision bvalue

      b = 0.0d0
      b(n) = 1.0d0

      BSpline = bvalue(t,b,n,Order+1,x,Deriv)

      return
      end

      double precision function bvalue ( t, bcoef, n, k, x, jderiv )
c     from  * a practical guide to splines *  by c. de boor    
c     alls  interv
c     
c     alculates value at  x  of  jderiv-th derivative of spline from b-repr.
c     the spline is taken to be continuous from the right, EXCEPT at the
c     rightmost knot, where it is taken to be continuous from the left.
c     
c******i n p u t ******
c     t, bcoef, n, k......forms the b-representation of the spline  f  to
c     be evaluated. specifically,
c     t.....knot sequence, of length  n+k, assumed nondecreasing.
c     bcoef.....b-coefficient sequence, of length  n .
c     n.....length of  bcoef  and dimension of spline(k,t),
c     a s s u m e d  positive .
c     k.....order of the spline .
c     
c     w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c     arbitrarily by the dimension statement for  aj, dl, dr  below,
c     but is  n o w h e r e  c h e c k e d  for.
c     
c     x.....the point at which to evaluate .
c     jderiv.....integer giving the order of the derivative to be evaluated
c     a s s u m e d  to be zero or positive.
c     
c******o u t p u t  ******
c     bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c     
c******m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c     cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c     this interval are then obtained from  bcoef (or taken to be zero if
c     not explicitly available) and are then differenced  jderiv  times to
c     obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c     Precisely, with  j = jderiv, we have from x.(12) of the text that
c     
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c     
c     where
c     / bcoef(.),                     ,  j .eq. 0
c     /
c     bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c     / ----------------------------- ,  j .gt. 0
c     /    (t(.+k-j) - t(.))/(k-j)
c     
c     Then, we use repeatedly the fact that
c     
c     sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c     with
c     (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c     a(.,x)  =    ---------------------------------------
c     (x - t(.))      + (t(.+m-1) - x)
c     
c     to write  (d**j)f(x)  eventually as a linear combination of b-splines
c     of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c     desired number  (d**j)f(x). (see x.(17)-(19) of text).
c     
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     *     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
C     double precision bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      double precision bcoef(*),t(*),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
c     dimension t(n+k)
c     former fortran standard made it impossible to specify the length of  t
c     precisely without the introduction of otherwise superfluous addition-
c     al arguments.
      bvalue = 0.0d0
      if (jderiv .ge. k)                go to 99
c     
c     *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c     t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c     outside the support of  the spline  f , hence  bvalue = 0.
c     (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c     at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c     *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
      go to 99
c     
c     *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do j=1,i
         dl(j) = x - t(i+1-j)
      enddo
      do j=i,km1
         aj(k-j) = 0.0d0
         dl(j) = dl(i)
      enddo
      go to 10
    8 do j=1,km1
         dl(j) = x - t(i+1-j)
      enddo
c     
 10   jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do j=1,jcmax
         dr(j) = t(i+j) - x
      enddo
      do j=jcmax,km1
         aj(j+1) = 0.0d0
         dr(j) = dr(jcmax)
      enddo
      go to 20
 18   do j=1,km1
         dr(j) = t(i+j) - x
      enddo
c     
 20   do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
      enddo
c     
c     *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
         enddo
      enddo
c     
c     *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
 30   if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
            ilo = ilo - 1
         enddo
      enddo
 39   bvalue = aj(1)
c     
 99   return
      end
      subroutine interv ( xt, lxt, x, left, mflag )
c     from  * a practical guide to splines *  by C. de Boor    
c     omputes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c     
c******i n p u t  ******
c     xt.....a double precision sequence, of length  lxt , assumed to be nondecreasing
c     lxt.....number of terms in the sequence  xt .
c     x.....the point whose location with respect to the sequence  xt  is
c     to be determined.
c     
c******o u t p u t  ******
c     left, mflag.....both integers, whose value is
c     
c     1     -1      if               x .lt.  xt(1)
c     i      0      if   xt(i)  .le. x .lt. xt(i+1)
c     i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c     i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c     
c     In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c     indicates that  x  lies outside the CLOSED interval
c     xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c     intervals is due to the decision to make all pp functions cont-
c     inuous from the right, but, by returning  mflag = 0  even if
C     x = xt(lxt), there is the option of having the computed pp function
c     continuous from the left at  xt(lxt) .
c     
c******m e t h o d  ******
c     The program is designed to be efficient in the common situation that
c     it is called repeatedly, with  x  taken from an increasing or decrea-
c     sing sequence. This will happen, e.g., when a pp function is to be
c     graphed. The first guess for  left  is therefore taken to be the val-
c     ue returned at the previous call and stored in the  l o c a l  varia-
c     ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c     essary since the present call may have nothing to do with the previ-
c     ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c     ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c     while also moving  ilo  and  ihi  in the direction of  x , until
c     xt(ilo) .le. x .lt. xt(ihi) ,
c     after which we use bisection to get, in addition, ilo+1 = ihi .
c     left = ilo  is then returned.
c     
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
      if (x .ge. xt(lxt))            go to 110
      if (lxt .le. 1)                go to 90
      ilo = lxt - 1
      ihi = lxt
c     
 20   if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c     
c     **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
 31   ihi = ilo
      ilo = ihi - istep
      if (ilo .le. 1)                go to 35
      if (x .ge. xt(ilo))            go to 50
      istep = istep*2
      go to 31
 35   ilo = 1
      if (x .lt. xt(1))                 go to 90
      go to 50
c     **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
 40   istep = 1
 41   ilo = ihi
      ihi = ilo + istep
      if (ihi .ge. lxt)              go to 45
      if (x .lt. xt(ihi))            go to 50
      istep = istep*2
      go to 41
 45   if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c     
c     **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
 50   middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
      ilo = middle
      go to 50
 53   ihi = middle
      go to 50
c**** set output and return.
 90   mflag = -1
      left = 1
      return
 100  mflag = 0
      left = ilo
      return
 110  mflag = 1
      if (x .eq. xt(lxt)) mflag = 0
      left = lxt
 111  if (left .eq. 1)                  return
      left = left - 1
      if (xt(left) .lt. xt(lxt))        return
      go to 111
      end

c     This subroutine returns the converged approximations to eigenvalues
c     of A*z = lambda*B*z and (optionally):
c     
c     (1) The corresponding approximate eigenvectors;
c     
c     (2) An orthonormal (Lanczos) basis for the associated approximate
c     invariant subspace;
c     
c     (3) Both.
c     
c     Matrices A and B are stored in LAPACK-style band form.
c     
c     There is negligible additional cost to obtain eigenvectors.  An orthonormal
c     (Lanczos) basis is always computed.  There is an additional storage cost 
c     of n*nev if both are requested (in this case a separate array Z must be 
c     supplied).
c     
c     The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c     are called Ritz values and Ritz vectors respectively.  They are referred 
c     to as such in the comments that follow.  The computed orthonormal basis 
c     for the invariant subspace corresponding to these Ritz values is referred 
c     to as a Lanczos basis.
c     
c     \Usage
c     call dsband
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA, 
c     RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V, 
c     LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
c     
c     \Arguments
c     
c     RVEC    Logical (INPUT)
c     Specifies whether Ritz vectors corresponding to the Ritz value 
c     approximations to the eigenproblem A*z = lambda*B*z are computed.
c     
c     RVEC = .FALSE.     Compute Ritz values only.
c     
c     RVEC = .TRUE.      Compute the associated Ritz vectors. 
c     
c     HOWMNY  Character*1  (INPUT) 
c     Specifies how many Ritz vectors are wanted and the form of Z
c     the matrix of Ritz vectors. See remark 1 below.
c     = 'A': compute all Ritz vectors;
c     = 'S': compute some of the Ritz vectors, specified
c     by the logical array SELECT.
c     
c     SELECT  Logical array of dimension NCV.  (INPUT)
c     If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c     computed. To select the Ritz vector corresponding to a
c     Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c     If HOWMNY = 'A' , SELECT is not referenced.
c     
c     D       Double precision array of dimension NEV.  (OUTPUT)
c     On exit, D contains the Ritz value approximations to the
c     eigenvalues of A*z = lambda*B*z. The values are returned
c     in ascending order. If IPARAM(7) = 3,4,5 then D represents
c     the Ritz values of OP computed by dsaupd transformed to
c     those of the original eigensystem A*z = lambda*B*z. If 
c     IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c     as the those of A*z = lambda*B*z.
c     
c     Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c     On exit, Z contains the B-orthonormal Ritz vectors of the
c     eigensystem A*z = lambda*B*z corresponding to the Ritz
c     value approximations.
c     
c     If  RVEC = .FALSE. then Z is not referenced.
c     NOTE: The array Z may be set equal to first NEV columns of the 
c     Lanczos basis array V computed by DSAUPD.
c     
c     LDZ     Integer.  (INPUT) 
c     The leading dimension of the array Z.  If Ritz vectors are
c     desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c     
c     SIGMA   Double precision  (INPUT)
c     If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c     IPARAM(7) = 1 or 2.
c     
c     N       Integer.  (INPUT) 
c     Dimension of the eigenproblem.  
c     
c     AB      Double precision array of dimension LDA by N. (INPUT)
c     The matrix A in band storage, in rows KL+1 to
c     2*KL+KU+1; rows 1 to KL of the array need not be set.
c     The j-th column of A is stored in the j-th column of the
c     array AB as follows:
c     AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c     
c     MB      Double precision array of dimension LDA by N. (INPUT)
c     The matrix M in band storage, in rows KL+1 to
c     2*KL+KU+1; rows 1 to KL of the array need not be set. 
c     The j-th column of M is stored in the j-th column of the
c     array AB as follows:
c     MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c     Not referenced if IPARAM(7) = 1
c     
c     LDA     Integer. (INPUT)
c     Leading dimension of AB, MB, RFAC.
c     
c     RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c     RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c     is invoked.  It is used to store the LU factors of
c     (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c     It is not referenced when IPARAM(7) = 1.
c     
c     KL      Integer. (INPUT)
c     Max(number of subdiagonals of A, number of subdiagonals of M)
c     
c     KU      Integer. (OUTPUT)
c     Max(number of superdiagonals of A, number of superdiagonals of M)
c     
c     WHICH   Character*2.  (INPUT)
c     When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c     the following.
c     
c     'LM' -> want the NEV eigenvalues of largest magnitude.
c     'SM' -> want the NEV eigenvalues of smallest magnitude.
c     'LA' -> want the NEV eigenvalues of largest REAL part.
c     'SA' -> want the NEV eigenvalues of smallest REAL part.
c     'BE' -> Compute NEV eigenvalues, half from each end of the 
c     spectrum.  When NEV is odd, compute one more from 
c     the high end than from the low end. 
c     
c     When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
c     
c     BMAT    Character*1.  (INPUT)
c     BMAT specifies the type of the matrix B that defines the
c     semi-inner product for the operator OP.
c     BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c     BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c     NEV     Integer. (INPUT)
c     Number of eigenvalues of OP to be computed.
c     
c     TOL     Double precision scalar.  (INPUT)
c     Stopping criterion: the relative accuracy of the Ritz value 
c     is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c     If TOL .LE. 0. is passed a default is set:
c     DEFAULT = DLAMCH('EPS')  (machine precision as computed
c     by the LAPACK auxiliary subroutine DLAMCH).
c     
c     RESID   Double precision array of length N.  (INPUT/OUTPUT)
c     On INPUT:
c     If INFO .EQ. 0, a random initial residual vector is used.
c     If INFO .NE. 0, RESID contains the initial residual vector,
c     possibly from a previous run.
c     On OUTPUT:
c     RESID contains the final residual vector.
c     
c     NCV     Integer.  (INPUT)
c     Number of columns of the matrix V (less than or equal to N).
c     Represents the dimension of the Lanczos basis constructed
c     by dsaupd for OP.
c     
c     V       Double precision array N by NCV.  (OUTPUT)
c     Upon INPUT: the NCV columns of V contain the Lanczos basis 
c     vectors as constructed by dsaupd for OP.
c     Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c     represent the Ritz vectors that span the desired 
c     invariant subspace.
c     NOTE: The array Z may be set equal to first NEV columns of the 
c     Lanczos basis vector array V computed by dsaupd. In this case
c     if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c     the desired Ritz vectors. 
c     
c     LDV     Integer.  (INPUT)
c     Leading dimension of V exactly as declared in the calling
c     program.
c     
c     IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c     IPARAM(1) = ISHIFT: 
c     The shifts selected at each iteration are used to restart
c     the Arnoldi iteration in an implicit fashion.
c     It is set to 1 in this subroutine.  The user do not need
c     to set this parameter.
c     ------------------------------------------------------------
c     ISHIFT = 1: exact shifts with respect to the reduced 
c     tridiagonal matrix T.  This is equivalent to 
c     restarting the iteration with a starting vector 
c     that is a linear combination of Ritz vectors 
c     associated with the "wanted" Ritz values.
c     -------------------------------------------------------------
c     
c     IPARAM(2) = No longer referenced. 
c     
c     IPARAM(3) = MXITER
c     On INPUT:  max number of Arnoldi update iterations allowed.
c     On OUTPUT: actual number of Arnoldi update iterations taken.
c     
c     IPARAM(4) = NB: blocksize to be used in the recurrence.
c     The code currently works only for NB = 1.
c     
c     IPARAM(5) = NCONV: number of "converged" eigenvalues.
c     This represents the number of Ritz values that satisfy
c     the convergence criterion.
c     
c     IPARAM(6) = IUPD
c     No longer referenced. Implicit restarting is ALWAYS used. 
c     
c     IPARAM(7) = MODE
c     On INPUT determines what type of eigenproblem is being solved.
c     Must be 1,2,3,4,5; See under \Description of dsband for the 
c     five modes available.
c     
c     IPARAM(8) = NP
c     Not referenced.
c     
c     IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c     OUTPUT: NUMOP  = total number of OP*x operations,
c     NUMOPB = total number of B*x operations if BMAT='G',
c     NUMREO = total number of steps of re-orthogonalization.
c     
c     WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
c     
c     WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
c     
c     LWORKL   Integer.  (INPUT)
c     LWORKL must be at least NCV**2 + 8*NCV.
c     
c     IWORK    Integer array of dimension at least N. (WORKSPACE)
c     Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
c     factorization of M or (A-SIGMA*M).
c     
c     INFO     Integer.  (INPUT/OUTPUT)
c     Error flag on output.
c     =  0: Normal exit.
c     =  1: Maximum number of iterations taken.
c     All possible eigenvalues of OP has been found. IPARAM(5)  
c     returns the number of wanted converged Ritz values.
c     =  3: No shifts could be applied during a cycle of the 
c     Implicitly restarted Arnoldi iteration. One possibility 
c     is to increase the size of NCV relative to NEV. 
c     See remark 4 in DSAUPD.
c     
c     = -1: N must be positive.
c     = -2: NEV must be positive.
c     = -3: NCV-NEV >= 2 and less than or equal to N.
c     = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c     = -6: BMAT must be one of 'I' or 'G'.
c     = -7: Length of private work WORKL array is not sufficient.
c     = -8: Error return from trid. eigenvalue calculation;
c     Informational error from LAPACK routine dsteqr.
c     = -9: Starting vector is zero.
c     = -10: IPARAM(7) must be 1,2,3,4,5.
c     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c     = -12: NEV and WHICH = 'BE' are incompatible.
c     = -13: HOWMNY must be one of 'A' or 'P'
c     = -14: DSAUPD did not find any eigenvalues to sufficient
c     accuracy.
c     = -9999: Could not build an Arnoldi factorization.
c     IPARAM(5) returns the size of the current
c     Arnoldi factorization.
c     
c     \Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c     Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine. 
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c     
c     \Remarks
c     1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c     
c     2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c     
c     \Author    
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c     
      subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,
     >     tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)

      character        which*2, bmat, howmny
      integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info  
      Double precision tol, sigma
      logical          rvec

      integer          iparam(*), iwork(*)
      logical          select(*)
      Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*)

      integer          ipntr(14)

      integer          ido, i, j, Row, Col, type, ierr

      Double precision one, zero
      parameter        (one = 1.0d0, zero = 0.0d0)

      Double precision ddot, dnrm2, dlapy2
      external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy

c     iparam(3) : Max number of Arnoldi iterations
      iparam(3) = 300
      iparam(7) = 3
      rvec = .TRUE.
      howmny = 'A'
      which = 'LA'
      bmat = 'G'
      type = 4 
      ido = 0
      iparam(1) = 1

      rfac = 0.0d0
      do i = 1,n
         do j = i,min(i+k,n)
            Row = k+1+i-j
            Col = j
            rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
         enddo
         do j = max(1,i-k),i-1
            Row = 2*k+1
            Col = j
            rfac(Row+i-j,j) = rfac(Row+j-i,i)
         enddo
      enddo

      call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
      if ( ierr .ne. 0 )  then
         print*, ' '
         print*, '_SBAND: Error with _gbtrf:',ierr
         print*, ' '
         go to 9000
      end if

 90   continue 

      call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)

      if (ido .eq. -1) then
         call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
         call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
         if (ierr .ne. 0) then
            print*, ' ' 
            print*, '_SBAND: Error with _gbtrs.'
            print*, ' ' 
            go to 9000
         end if
      else if (ido .eq. 1) then
         call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
         call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
         if (ierr .ne. 0) then 
            print*, ' '
            print*, '_SBAND: Error with _gbtrs.' 
            print*, ' '
            go to 9000
         end if
      else if (ido .eq. 2) then
         call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
      else 
         if ( info .lt. 0) then
            print *, ' '
            print *, ' Error with _saupd info = ',info
            print *, ' Check the documentation of _saupd '
            print *, ' '
            go to 9000
         else 
            if ( info .eq. 1) then
c     print *, ' '
c     print *, ' Maximum number of iterations reached.'
c     print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               write(6,*) "No shifts could be applied during"
               write(6,*) "implicit Arnoldi update, try increasing NCV."
               print *, ' '
            end if
            if (iparam(5) .gt. 0) then
               call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
               if ( info .ne. 0) then
                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
               endif
            endif
         endif
         go to 9000
      endif

      go to 90 

 9000 continue

      end
      subroutine GetGaussFactors(File,Points,x,w)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     This subroutine retrieves the points and weights for
c     Gaussian Quadrature from a given file
c     
c     Variables:
c     File		name of file containing points and 
c     weights
c     Points		number of points to be used for 
c     quadrature
c     x		array of nodes
c     w		array of weights
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      integer Points
      double precision x(Points),w(Points)
      character*64 File
      
      integer i,tempPoints
      
      open(unit=7,file=File(1:index(File,' ')-1))

      do i = 1,18
         read(7,*)
      enddo
      read(7,*) tempPoints
      do while (tempPoints .ne. Points)
         do i = 1,tempPoints
            read(7,*)
         enddo
         read(7,*) tempPoints
      enddo
      
      do i = 1,Points
         read(7,*) x(i),w(i)
      enddo
      
      close(unit=7)

      return
      end

      subroutine GridMaker(R,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints,d31,Rscale,r0xfactor,r0yfactor)

      integer xNumPoints,yNumPoints,ncount
      double precision R,xMin,xMax,yMin,yMax,d31,Rscale,r0xfactor,r0yfactor
      double precision xPoints(xNumPoints),yPoints(yNumPoints)

      integer xNumPoints1,xNumPoints2,xNumPoints3,xNumPoints4,yNumPoints1,yNumPoints2
      integer i1,i2,j1,j2,k,ix,iy,i3,i4
      double precision Pi
      double precision xC1,xC2,xC3,xDelt1,xDelt2,xDelt3,xDelt4
      double precision yC,yDelt1,yDelt2
      double precision ax,axc,ay,ayc,r0,r0x,r0y 

      Pi = 3.1415926535897932385d0

      r0x = r0xfactor*Rscale

c     Three sectors grid (for two identical particles with v12=0)
      ax = (dsqrt(2.d0)*r0x/(d31*R))**2-1.d0
      axc = -dsqrt(3.d0)/2.d0
      if(ax.ge.axc) ax=axc
      xC1 = -dasin(ax)/Pi-1.0d0/6.d0
      xC2 = dasin(ax)/Pi+1.d0/1.2d0
      xNumPoints2 = 3*xNumPoints/4
      xNumPoints3 = xNumPoints/6
      xNumPoints1 = xNumPoints-xNumPoints2-xNumPoints3
      xDelt1 = (xC1-xMin)/dfloat(xNumPoints1-1)
      do i1 = 1,xNumPoints1
         xPoints(i1) = ((i1-1)*xDelt1+xMin)*Pi
      enddo
      xDelt2 = (xC2-xC1)/dfloat(xNumPoints2)
      do i2 = 1,xNumPoints2
         xPoints(xNumPoints1+i2) = (i2*xDelt2+xC1)*Pi
      enddo
      xDelt3 = (xMax-xC2)/dfloat(xNumPoints3)
      do i3 = 1,xNumPoints3
         xPoints(xNumPoints1+xNumPoints2+i3) = (i3*xDelt3+xC2)*Pi
      enddo
c     Smooth Grid (jpdincao 10-06-03)
      do k=1,xNumPoints/5
         do ix=2,xNumPoints-1
            xPoints(ix)=(xPoints(ix-1)+2.d0*xPoints(ix)+xPoints(ix+1))/4.d0
         enddo
      enddo

      
      r0y = r0yfactor*Rscale

      ay = 1.d0-(dsqrt(2.d0)*r0y/(d31*R))**2
      ayc = 1.d0/dsqrt(2.d0)
      if(ay.le.ayc) ay=ayc
      yC = dasin(ay)/Pi
      yNumPoints2 = yNumPoints/2
      yNumPoints1 = yNumPoints-yNumPoints2
      yDelt1 = (yC-yMin)/dfloat(yNumPoints1-1)
      do j1 = 1,yNumPoints1
         yPoints(j1) = ((j1-1)*yDelt1+yMin)*Pi
      enddo
      yDelt2 = (yMax-yC)/dfloat(yNumPoints2)
      do j2 = 1,yNumPoints2
         yPoints(yNumPoints1+j2) = (j2*yDelt2+yC)*Pi
      enddo
c     Smooth Grid (jpdincao 10-06-03)
      do k=1,yNumPoints/5
         do iy=2,yNumPoints-1
            yPoints(iy)=(yPoints(iy-1)+2.d0*yPoints(iy)+yPoints(iy+1))/4.d0
         enddo
      enddo

      return
      end


      subroutine GridMakerNEW(R,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints,d12,d31,d23,
     .     phi12,phi31,phi23,RscaleAB,RscaleAA,r0xfactor,r0yfactor) 

      integer xNumPoints,yNumPoints,ncount
      double precision R,xMin,xMax,yMin,yMax,RscaleAB,RscaleAA,r0xfactor,r0yfactor
      double precision d12,d31,d23
      double precision phi12,phi31,phi23
      double precision xPoints(xNumPoints),yPoints(yNumPoints)

      integer xNumPoints1,xNumPoints2,xNumPoints3,xNumPoints4,xNumPoints5,xNumPoints6,xNumPoints7,yNumPoints1,yNumPoints2
      integer i1,i2,i3,i4,i5,i6,i7,j1,j2,k,ix,iy,i
      double precision Pi
      double precision xC1,xC2,xC3,xC4,xC5,xC6,xC7,xDelt1,xDelt2,xDelt3,xDelt4,xDelt5,xDelt6,xDelt7
      double precision yC,yDelt1,yDelt2
      double precision ax,axc,ay,ayc,r0,r0xAA,r0xAB,r0y 
      double precision phi12c,phi31c,phi23c

      Pi = 3.1415926535897932385d0

      phi12c = MOD(Pi-phi12,2.d0*Pi)
      phi23c = Pi

      r0xAA = r0xfactor*RscaleAA
      r0xAB = r0xfactor*RscaleAB
      
      ax = (dsqrt(2.d0)*r0xAB/(d12*R))**2-1.d0
      axc = dsin(phi12c/2.d0-Pi/2.d0) !-dsqrt(3.d0)/2.d0
      if (ax.ge.axc) ax = axc
      xC1 = -dasin(ax)/Pi-1.d0/2.d0+phi12c/Pi
      xC2 =  dasin(ax)/Pi+1.d0/2.d0+phi12c/Pi
      ax = (dsqrt(2.d0)*r0xAA/(d23*R))**2-1.d0
      axc = -dsqrt(3.d0)/2.d0
      if (ax.ge.axc) ax = axc
      xC3 = -dasin(ax)/Pi-1.d0/2.d0+phi23c/Pi
      xC4 =  xMax
      xNumPoints1 = xNumPoints/8
      xNumPoints2 = 3*xNumPoints/8
      xNumPoints3 = xNumPoints/8
      xNumPoints4 = xNumPoints-xNumPoints1-xNumPoints2-xNumPoints3
      xDelt1 = (xC1-xMin)/(xNumPoints1-1)
      do i1 = 1,xNumPoints1
         xPoints(i1) = ((i1-1)*xDelt1+xMin)*Pi
      enddo
      xDelt2 = (xC2-xC1)/(xNumPoints2)
      do i2 = 1,xNumPoints2 
         xPoints(xNumPoints1+i2) = (i2*xDelt2+xC1)*Pi
      enddo
      xDelt3 = (xC3-xC2)/(xNumPoints3)
      do i3 = 1,xNumPoints3
         xPoints(xNumPoints1+xNumPoints2+i3) = (i3*xDelt3+xC2)*Pi
      enddo
      xDelt4 = (xC4-xC3)/(xNumPoints4)
      do i4 = 1,xNumPoints4
         xPoints(xNumPoints1+xNumPoints2+xNumPoints3+i4) = (i4*xDelt4+xC3)*Pi
      enddo
c     Smooth Grid (jpdincao 10-06-03)
      do k=1,xNumPoints/5
         do ix=2,xNumPoints-1 
            xPoints(ix) = (xPoints(ix-1)+2.d0*xPoints(ix)+xPoints(ix+1))/4.d0
         enddo
      enddo

!      do k=1,xNumPoints
!         write(777,*)k,xPoints(k)
!      enddo

      r0y = r0yfactor*(RscaleAA+RscaleAB)/2.d0

      dij = min(d12,d31,d23)
      ay = 1.d0-(r0y*dsqrt(2.d0)/dij/R)**2
      ayc = dsqrt(2.d0)/2.d0
      if (ay.le.ayc) ay = ayc
      yC = dasin(ay)/Pi
      yNumPoints2 = yNumPoints/2
      yNumPoints1 = yNumPoints-yNumPoints2
      yDelt1 = (yC-yMin)/(yNumPoints1-1)
      do j1=1,yNumPoints1
         yPoints(j1) = ((j1-1)*yDelt1+yMin)*Pi
      enddo
      yDelt2 = (yMax-yC)/(yNumPoints2)
      do j2=1,yNumPoints2
         yPoints(yNumPoints1+j2) = (j2*yDelt2+yC)*Pi
      enddo
c     Smooth Grid (jpdincao 10-06-03)
      do k=1,yNumPoints/5
         do iy=2,yNumPoints-1
            yPoints(iy) = (yPoints(iy-1)+2.d0*yPoints(iy)+yPoints(iy+1))/4.d0
         enddo
      enddo

      return
      end




      subroutine Telescope(NumEntries,MatchTol,Energies,Psi,PsiLeadDim)

      integer NumEntries,PsiLeadDim
      double precision MatchTol,Energies(NumEntries),Psi(PsiLeadDim,NumEntries)

      integer i,j,k,l,m

      l = NumEntries
      i = 1
      do while (i .le. l)
         j = 0
         do while ((dabs((Energies(i)-Energies(i+j+1))/Energies(i)) .gt. MatchTol) .AND. (i+j+1 .le. l))
            j = j + 1
         enddo
         if (j .ne. l-i) then
            do k = i+j+1,l
               Energies(k-j-1) = Energies(k)
               do m = 1,PsiLeadDim
                  Psi(m,k-j-1) = Psi(m,k)
               enddo
            enddo
            l = l - (j+1)
         endif
         i = i + 1
      enddo

      NumEntries = l

      return
      end
      
      subroutine MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,Psi,Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)

      integer NumStateInc,MatrixDim,LeadDim,HalfBandWidth,NumStates,NumFirst
      double precision Energies(NumStates,2),Psi(MatrixDim,NumStates)
      double precision H(HalfBandWidth+1,MatrixDim),S(HalfBandWidth+1,MatrixDim),LUFac(LeadDim,MatrixDim)
      double precision Shift

      integer i,j,TempNumStates,iLoop
      integer CurrentNumStates,NumStatesAsk
      integer ncv,iparam(11),lworkl,info
      integer, allocatable :: iwork(:)
      logical, allocatable :: Select(:)
      double precision Tol,TempShift,OldTempShift,MatchTol,Shift2
      double precision, allocatable :: workd(:),workl(:),Residuals(:)
      double precision, allocatable :: TempPsi(:,:),TempEnergies(:,:)
      double precision, allocatable :: AllPsi(:,:),AllEnergies(:)
      common /Rvalue/ Rvalue

      MatchTol = 1.0d-8

      H = (Rvalue**2.d0)*H
c     S = (Rvalue**2.d0)*S
      Shift = (Rvalue**2.d0)*Shift
      Shift2 = (Rvalue**2.d0)*Shift2

c     ncv = 2*NumStateInc
      ncv = 2*NumStates
      lworkl = ncv*ncv+8*ncv
      allocate(Select(ncv),iwork(MatrixDim),workd(3*MatrixDim),workl(lworkl),Residuals(MatrixDim))
      allocate(TempPsi(MatrixDim,ncv),TempEnergies(ncv,2))
      allocate(AllPsi(MatrixDim,NumStates+NumStateInc),AllEnergies(NumStates+NumStateInc))

      TempShift = Shift
      info = 0
      if (NumStates .lt. NumStateInc) then
         NumStatesAsk = NumStates
      else
         NumStatesAsk = NumFirst
      endif
      call MyDsband(Select,TempEnergies,TempPsi,MatrixDim,TempShift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,
     >     HalfBandWidth,NumStatesAsk,Tol,Residuals,ncv,TempPsi,MatrixDim,iparam,workd,workl,lworkl,iwork,info)
      TempNumStates = min(NumStatesAsk,iparam(5))
      CurrentNumStates = TempNumStates
c     OldTempShift = TempShift
c     TempShift = TempEnergies(TempNumStates,1)
c     TempShift = 2.0d0*TempEnergies(TempNumStates,1)-TempEnergies(TempNumStates-5,1)
c     if (dabs((TempShift-OldTempShift)/OldTempShift) .le. 1.0d-10) then
c     TempShift = TempShift + TempEnergies(TempNumStates,1) - TempEnergies(TempNumStates-2,1)
c     endif
      do i = 1,TempNumStates
         AllEnergies(i) = TempEnergies(i,1)
         do j = 1,MatrixDim
            AllPsi(j,i) = TempPsi(j,i)
         enddo
      enddo
c     TempShift = AllEnergies(CurrentNumStates-1)

      TempShift = Shift2
      NumStatesAsk = NumStateInc

      iLoop = 0
      do while ((CurrentNumStates .lt. NumStates) .AND. (iLoop .lt. 5))

         iLoop = iLoop+1
c     write(*,*)iLoop

!write(6,*) 'Loop:',iLoop,CurrentNumStates
c     !write(6,*) '          ',CurrentNumStates,TempShift
         info = 0
!write(6,*)NumStatesAsk,CurrentNumStates,TempShift 
c     if (NumStates-CurrentNumStates+2 .lt. NumStateInc) then
c     NumStatesAsk = NumStates-CurrentNumStates+2
c     else
c     NumStatesAsk = NumStateInc
c     endif
         call MyDsband(Select,TempEnergies,TempPsi,MatrixDim,TempShift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,
     >        HalfBandWidth,NumStatesAsk,Tol,Residuals,ncv,TempPsi,MatrixDim,iparam,workd,workl,lworkl,iwork,info)
!write(6,*)'Current',iparam(5)

c     if (TempEnergies(1,1) .le. AllEnergies(CurrentNumStates)+MatchTol*dabs(AllEnergies(CurrentNumStates))) then

         TempNumStates = min(NumStatesAsk,iparam(5))
c     OldTempShift = TempShift
c     TempShift = TempEnergies(TempNumStates,1)
c     TempShift = 2.0d0*TempEnergies(TempNumStates,1)-TempEnergies(TempNumStates-5,1)
c     if (dabs((TempShift-OldTempShift)/OldTempShift) .le. 1.0d-10) then
c     TempShift = TempShift + TempEnergies(TempNumStates,1) - TempEnergies(TempNumStates-2,1)
c     endif

         do i = 1,TempNumStates
            AllEnergies(CurrentNumStates+i) = TempEnergies(i,1)
c     !write(6,*) i,TempEnergies(i,1)
            do j = 1,MatrixDim
               AllPsi(j,CurrentNumStates+i) = TempPsi(j,i)
            enddo
         enddo
c     !write(6,*)

         CurrentNumStates = CurrentNumStates+TempNumStates
         
         call Telescope(CurrentNumStates,MatchTol,AllEnergies,AllPsi,MatrixDim)

         TempShift = AllEnergies(CurrentNumStates-1)

c     else

c     OldTempShift = TempShift
c     TempShift = AllEnergies(CurrentNumStates-1)
c     if (dabs((TempShift-OldTempShift)/OldTempShift) .le. 1.0d-10) then
c     TempShift = TempEnergies(TempNumStates-1,1)
c     c         TempShift = TempShift - (TempEnergies(TempNumStates,1) - TempEnergies(TempNumStates-1,1))
c     endif
c     !write(6,*) TempShift
c     TempShift = AllEnergies(CurrentNumStates-5)

c     endif

      enddo

      do i = 1,NumStates
         Energies(i,1) = AllEnergies(i)
         do j = 1,MatrixDim
            Psi(j,i) = AllPsi(j,i)
         enddo
      enddo

      call CalcEigenErrors(info,MatrixDim,H,S,HalfBandWidth,NumStates,Psi,Energies,NumStates)

      deallocate(Select,iwork,workd,workl,Residuals)
      deallocate(TempPsi,TempEnergies)
      deallocate(AllPsi,AllEnergies)

      H = H/(Rvalue**2.d0)
c     S = S/(Rvalue**2.d0)
      Shift = Shift/(Rvalue**2.d0)
      Shift2 = Shift2/(Rvalue**2.d0)
      Energies = Energies/(Rvalue**2.d0)


      return
      end

      subroutine CalcEigenErrors(info,MatrixDim,H,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

      integer info,MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
      double precision H(HalfBandWidth+1,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

      integer j
      double precision dnrm2
      double precision, allocatable :: HPsi(:),SPsi(:)

      if ( info .eq. 0) then

c     Compute the residual norm: ||  A*x - lambda*x ||

         allocate(HPsi(MatrixDim))
         allocate(SPsi(MatrixDim))
         do j = 1,NumStates
            call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,HalfBandWidth+1,Psi(1,j),1,0.0d0,HPsi,1)
            call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
            call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
            Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
            Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
         enddo
         deallocate(HPsi)
         deallocate(SPsi)
      else
!write(6,*) ' '
!write(6,*) ' Error with _sband, info= ', info
!write(6,*) ' Check the documentation of _sband '
!write(6,*) ' '
      end if

      return
      end

      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      double precision arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      double precision a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         do 12 j=l+1,ir
            a=arr(j)
            do 11 i=j-1,1,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
 11         continue
            i=0
 2          arr(i+1)=a
 12      continue
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
         endif
         i=l+1
         j=ir
         a=arr(l)
 3       continue
         i=i+1
         if(arr(i).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         goto 3
 5       arr(l)=arr(j)
         arr(j)=a
         jstack=jstack+2
         if(jstack.gt.NSTACK) then
            write(6,*) "press enter to continue, NSTACK too small to sort"
            read(*,*) 
         endif
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END


      subroutine CalcRadialCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q)

      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision P(NumStates,NumStates),Q(NumStates,NumStates)

      integer i,j,k
      double precision aP,aQ,ddot
      double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:)

      allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim))

      aP = 0.5d0/RDelt
      aQ = aP*aP

      do j = 1,NumStates
         do k = 1,MatrixDim
            rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
         do i = 1,NumStates
            P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
            do k = 1,MatrixDim
               lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
            enddo
            Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
         enddo
      enddo

      deallocate(lDiffPsi,rDiffPsi,TempPsi)

      return
      end

      subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)

      integer NumStates,HalfBandWidth,MatrixDim,ncv
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)

      integer i,j
      double precision Phase,ddot
      double precision, allocatable :: TempPsi(:)

      allocate(TempPsi(MatrixDim))

      do i = 1,NumStates
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
         Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
         if (Phase .lt. 0.0d0) then
            do j = 1,MatrixDim
               rPsi(j,i) = -rPsi(j,i)
            enddo
         endif
      enddo

      deallocate(TempPsi)

      return
      end

c$$$  
c$$$  subroutine CalcPlotBasisFuncs(Left,Right,Order,xPoints,x,xNumPlotPoints,MatrixDim,xBounds,xNumPoints,Deriv,u)
c$$$  
c$$$  integer Left,Right,Order,LegPoints,xNumPlotPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
c$$$  double precision xPoints(*),x(*)
c$$$  double precision u(MatrixDim,xNumPlotPoints)
c$$$  
c$$$  integer i,k,l,Count
c$$$  integer, allocatable :: t(:)
c$$$  double precision BSpline
c$$$  double precision, allocatable :: tx(:),b(:)
c$$$  
c$$$  allocate(t(xNumPoints+2*Order),tx(xNumPoints+2*Order),b(xNumPoints+Order))
c$$$  
c$$$  do i = 1,Order
c$$$  t(i) = 1
c$$$  tx(i) = xPoints(1)
c$$$  enddo
c$$$  do i = 1,xNumPoints
c$$$  t(i+Order) = i
c$$$  tx(i+Order) = xPoints(i)
c$$$  enddo
c$$$  do i = 1,Order
c$$$  t(i+Order+xNumPoints) = xNumPoints
c$$$  tx(i+Order+xNumPoints) = xPoints(xNumPoints)
c$$$  enddo
c$$$  
c$$$  do i = 1,xNumPoints+2*Order
c$$$  xBounds(i) = 0.0d0
c$$$  enddo
c$$$  
c$$$  select case (Left)
c$$$  case (0:1)
c$$$  select case (Right)
c$$$  case (0:1)
c$$$  do i = 2,xNumPoints+2*Order-1
c$$$  xBounds(i-1) = t(i)
c$$$  enddo
c$$$  case (2)
c$$$  do i = 2,xNumPoints+2*Order
c$$$  xBounds(i-1) = t(i)
c$$$  enddo
c$$$  end select
c$$$  case (2)
c$$$  select case (Right)
c$$$  case (0:1)
c$$$  do i = 1,xNumPoints+2*Order-1
c$$$  xBounds(i) = t(i)
c$$$  enddo
c$$$  case (2)
c$$$  do i = 1,xNumPoints+2*Order
c$$$  xBounds(i) = t(i)
c$$$  enddo
c$$$  end select
c$$$  end select
c$$$  
c$$$  deallocate(t)
c$$$  
c$$$  Count = 1
c$$$  select case (Left)
c$$$  case (0)
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(k))
c$$$  enddo
c$$$  Count = Count + 1
c$$$  case (1)
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,1,tx,b,x(k))+BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(k))
c$$$  enddo
c$$$  Count = Count + 1
c$$$  case(2)
c$$$  do i = 1,2
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
c$$$  enddo
c$$$  Count = Count + 1
c$$$  enddo
c$$$  end select
c$$$  
c$$$  do i = 3,xNumPoints+Order-3
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
c$$$  enddo
c$$$  Count = Count + 1
c$$$  enddo
c$$$  
c$$$  select case (Right)
c$$$  case (0)
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(k))
c$$$  enddo
c$$$  case (1)
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(k))+
c$$$  >                BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,tx,b,x(k))
c$$$  enddo
c$$$  case(2)
c$$$  do i = xNumPoints+Order-2,xNumPoints+Order-1
c$$$  do k = 1,xNumPlotPoints
c$$$  u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
c$$$  enddo
c$$$  Count = Count + 1
c$$$  enddo
c$$$  end select
c$$$  
c$$$  deallocate(tx,b)
c$$$  
c$$$  return
c$$$  end
