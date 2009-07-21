!>     This program calculates the evolution of an ice sheet on a
!!     bed with constant slope.  The surface mass balance is taken
!!     constant.
!!
!!     Boundary conditions are an ice divide at I = 1.  At the other
!!     end, the flotation criterion is used to determine the position
!!     of the grounding line.  The profile is extended two gridpoints
!!     using the equilibrium profile of a free-floating ice shelf.
!!
!!     For a description of the model, and discussion of the results,
!!     see Applied Glaciology, Section 8.3.
!!     This version was used for the EISMINT model intercomparison
!!     experiments (Sept. 1997).

program eism

  implicit none
  integer,parameter :: dp = kind(1.0d0) !< double precision

  integer, parameter :: domainSize = 51  !< domain size
  real(kind=dp), dimension(domainSize) :: h    !< ice thickness
  real(kind=dp), dimension(domainSize) :: diff !< diffusivity
  real(kind=dp), dimension(domainSize) :: flux !< ice flux
  real(kind=dp), dimension(domainSize) :: dhdt !< change of ice thickness with time
  real(kind=dp), dimension(domainSize) :: hb   !< ice bed elevation
  real(kind=dp), dimension(domainSize) :: hs   !< ice surface elevation
  real(kind=dp), dimension(domainSize) :: ux   !< 
  real(kind=dp), dimension(domainSize) :: us   !< 

  !---  Input parameters
  real(kind=dp), parameter ::  SMALL=1.d-20    !< a small value
  real(kind=dp), parameter ::  AFL=1.0d-18     !< flow parameter
  real(kind=dp), parameter ::  ROI=910.0       !< density of ice
  real(kind=dp), parameter ::  ROW=1028.0      !< density of water
  real(kind=dp), parameter ::  G=9.81          !< acceleration due to gravity
  real(kind=dp), parameter ::  CONST = ((2.0/5.0)*AFL)*(ROI*G**3) !< a constant

  real(kind=dp), parameter ::  ALPHA=5.0*0.001
  real(kind=dp), parameter ::  HBNUL=-250.0
  real(kind=dp), parameter ::  ACC=0.3         !< surface mass balance
  real(kind=dp), parameter ::  CFLOAT=1.0-(ROI/ROW)
  real(kind=dp), parameter ::  C1=(ROW-ROI)/(4.0*ROW)
  real(kind=dp), parameter ::  CSHELF=(2.6*CONST)*(C1**3)

  real(kind=dp), parameter ::  DX=2000.0  !< horizontal grid spacing
  real(kind=dp), parameter ::  DT=0.50    !< time step

  real(kind=dp), parameter ::  TMAX=15000.0 !< final time
  real(kind=dp), parameter ::  CDIFF=CONST/((2.0*DX)**2)
  real(kind=dp), parameter ::  DTPR=250.0

  real(kind=dp) :: time=0.0 !< the current time
  real(kind=dp) :: tpr=0.0

  real(kind=dp) :: dist,fl,exx,dhsdx,depth,hfloat,vol
  real(kind=dp) :: hbase, ugr
  real(kind=dp) :: x1,x2,x3,x4a,x4b
  integer i,igr

  open(2,file='sshapenb.cvv',status='new')
  open(4,file='tshapenb.cvv',status='new')
  open(3,file='sshastnb.cvv',status='old')
  
  !---  Initial geometry.
  do i=1,domainSize
     read(3,1000) dist,h(i),hs(i),ux(i),us(i),fl,exx
     dist=real(i-1,kind=dp)*DX
     hb(i)=HBNUL-(ALPHA*dist)
!     h(I)=570.0
!     hs(i)=hb(i)+h(i)
  end do

  !---  start time integration
  time=0.0
  do while(time<TMAX)
     write(*,*) time,igr

     !---  Calculate diffusivity.
     diff(1)=0.0
     do i=2,domainSize
        if (abs(h(i))<SMALL) then
           diff(i) = 0.
        else
           dhsdx=hs(i+1)-hs(i-1)
           diff(I)=(CDIFF*dhsdx*dhsdx)*(h(i)**5.0)
        end if
     end do

     !---  Calculate ice flux.
     do i=1,domainSize-1
        flux(i)=((diff(i)+diff(i+1))/(2.0*dx))*(hs(i+1)-hs(i))
     end do

     !---  Calculate change in ice thickness.
     dhdt(1)=((2.0*flux(1))/DX)+ACC
     do i=2,domainSize-1
        dhdt(i)=((flux(i)-flux(i-1))/DX)+ACC
     end do
     dhdt(domainSize)=0.0

     !---  Calculate new ice thickness.
     do i=1,domainSize-1
        h(i)=h(i)+(dhdt(i)*DT)
        hfloat=(roi/row)*h(i)
        depth=-hb(i)
        if (hfloat<depth .and. depth>0.) then
           igr=i-1
           exit
        end if
        h(i) = max(h(i),0.0)
        hs(i)=hb(i)+h(i)
        igr=i
     end do

     !---  Extend ice-shelf profile one gridpoint.
     do i=igr,domainSize
        !265   CONTINUE
        dhsdx=(hs(i)-hs(i-1))/DX
        ugr=((dhsdx*dhsdx*dhsdx)*(h(i)**4.0))*CONST
        ugr=-ugr

        x1=CSHELF/ACC
        x2=((h(i)**4.0)*x1)-1.0
        x3=(ugr*ugr*ugr*ugr)*x2
        x4a=(ACC*DX)+(h(i)*ugr)
        x4a=x4a**4.0
        x4b=(ACC*2.0*DX)+(h(i)*ugr)
        x4b=x4b**4.0


        h(i+1)=(x1-(x3/x4a))**(-0.25)
        h(i+1) = min(h(i+1),0.577*h(i))
        hs(i+1)=CFLOAT*h(i)


        hfloat=(ROI/ROW)*h(i+1)
        depth=-hb(i+1)
        if (hfloat>depth) then
           hs(igr+1)=h(igr+1)+hb(igr+1)
           igr=igr+1
        else
           exit
        end if
     end do

     do i=igr+2,domainSize
        h(i)=0.0
        hs(i)=0.0
     end do

     !---  Write results to file.

     time=time+DT
     if (time>=tpr) then
        vol=sum(h(:igr))*DX
        !?????
        ux(21)=-(flux(21)+flux(20))*(0.5/h(21))
        ux(31)=0.0

        write(4,2000) time,vol,igr,h(1),ux(21),ux(31)
        tpr=tpr+DTPR
      end if
   end do


   do  i=2,domainSize
      if(abs(h(i))<SMALL) then
         ux(i)=0.0
         us(i)=0.0
      else
         ux(i)=-(flux(i)+flux(i-1))*(0.5/h(i))
         us(i)=0.0
      end if
   end do


   do i=1,domainSize
      dist=(real(I-1))*(DX/1000.0)
      hbase=hs(i)-h(i)
      fl=(ux(i)+us(i))*h(i)
      if(i==1) then
         exx=(ux(2)+us(2))/DX
      else
         if(i.eq.domainSize) then
            exx=0.0
         else
            exx=((ux(i+1)+us(i+1))-(ux(i-1)+us(i-1)))/(2.0*DX)
         end if
      end if
      write(2,1000) dist,h(i),hs(i),ux(i),us(i),fl,exx
   end do

1000 format(5(1x,f8.2),2(1xe14.6))
2000 format(1x,f8.0,1x,e14.6,1x,i2,3(1x,f8.2))

 end program eism
