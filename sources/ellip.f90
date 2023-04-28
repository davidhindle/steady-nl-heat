! program ellipt-1
! 1d, non-linear elliptic, differential equation (diffusion) in steady state
! (D(u)u')' = -q  ; u(0)=a, u(L)=b or u'(L)=bb , x E (0,..,L)
! boundary conditions - Dirichlet, surface; Dirichlet or Neumann at base
! Neumann defined as +ve for heat flow into domain
! -------------------------------------------------------------------------------------------
! METHOD
! 
! equation reduces to a finite difference approximation.
! D(i+1)
! leads to linear system
! a u = q
! solve with L.U decompostion of a - this is the method used to "invert" the matrix a.
!
! -------------------------------------------------------------------------------------------
! INPUT VARIABLES
! dx - grid spacing in metres
! Nuemann true of false - whether to use Neumann bc at x=L
! bc1, bc2, boundary conditions at x=0, x=L
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
! the following are derived from heat-in.dat 
! depth to top of lower crust (for Chapman)
! number of layers
! Layer depths, (to base of layer), layer standard thermal conductivity (ko), layer radiogenic heat production (ao)
!--------------------------------------------------------------------------------------------
! INTERMEDIATE VARIABLES
! tr,tr1 = temperature function in iteration
! q(i) - source term (radiogenic heat production by node)
! D(i) - variable diffusion coefficient in space (x dimension)
! ko(i) - initial STP values of thermal conductivity of layers
!--------------------------------------------------------------------------------------------
! OUTPUT VARIABLES 
! tro = initial solution for linear problem
! trj = Jaupart non-linear variant solution (Jaupart and Mareschal papers)
! trc = Chapman non-linear solution (Chapman 1986 paper)
! hfo, hfj, hfc - heat flows as above
! D0, dj, dc - thermal conductivity values as above
!_____________________________________________________________________________________________

use TDMA
implicit none

double precision, parameter :: zero=0.0d0
integer, parameter :: dp = kind(1.d0)
integer i,j,l,ct,ierr,nsub,n,kon
integer, dimension(:), allocatable :: kdn
double precision, dimension(:), allocatable :: u,q,inis,tr,tr1,tro,trj,trc,p1,p2,p3,D,diff,hf
double precision, dimension(:), allocatable :: hfo,hfc,hfj,dc,dj,D0
double precision :: qv,D1,D2,depth,iter,b,c,lc
double precision, dimension(:), allocatable :: ko,k1,ao
double precision :: dx,dx2
double precision ::  bc1, bc2, tol, kr
logical :: bckind, radiative


character out1*11, out2*11, cnla*3


! read in model parameters ko data - 
open (12,file='heat-in.dat')
read (12,*) radiative
read (12,*) bckind
read (12,*) tol 
read (12,*) kr
read (12,*) dx
dx2=dx**2
read (12,*) bc1
read (12,*) bc2
read (12,*) lc
read (12,*) kon
allocate (k1(kon),kdn(kon+1),ao(kon))
kdn(1)=1
i=1
do while (i.le.kon)
15 read (12,*) k1(i),ao(i),depth
kdn(i+1)=int(depth/dx) 
print *, 'node number of base of layer ', i, ' = ',  kdn(i+1)
i=i+1
end do

50 continue
nsub=kdn(i)+1
allocate (u(nsub),tr(nsub),tr1(nsub),q(nsub),inis(nsub),p1(nsub),p2(nsub),p3(nsub),D(nsub),ko(nsub),diff(nsub),hf(nsub))
allocate (tro(nsub),trj(nsub),trc(nsub),hfo(nsub),hfc(nsub),hfj(nsub),D0(nsub),dc(nsub),dj(nsub))
! set normal Ko, Ao, values

do j=1,kon
print *, 'layer top node = ', kdn(j), 'layer ko = ', k1(j)
 do i=kdn(j),kdn(j+1)
 ko(i)=k1(j)
 q(i)=-ao(j)
 D(i)=ko(i)
 end do
end do

! initialise vectors
p1 = zero ; p2= zero ; p3=zero 
!q(1:nsub)=0.d0

! call subroutine for linear, steady state solution

call matrix(D,nsub,bckind,dx,bc1,bc2,ao,q,p1,p2,p3)

!   solve linear problem using tdma  
!   solution stored as "tro"
call tri(nsub,p1,p2,p3,q,tr1) !n,a,b,c,q,u

do i=2,nsub-1
    hf(i)=((D(i+1)+D(i)+D(i-1))/3.d0) * (tr1(i+1) - tr1(i-1))/(2.d0*dx)
end do
hf(1)=hf(2)
hf(nsub)=hf(nsub-1)
d0=D
hfo=hf
hf=zero
D=zero
tr=tr1
tro=tr1


! define D(T,z) for different, temperature dependent models and solve
!________________________________________________________________________________
! ---- this is the Jaupart algorithm --------------------------------------------
! ---- includes radiative heat contribution to thermal conductivity--------------
! _______________________________________________________________________________
iter = 10.d0
do while (iter.gt.1.d-3)
    do i=1,nsub
     D(i) =  2.26d0 - 618.251d0/tr(i) + ko(i)*(355.576/tr(i) - 0.30247)
     if (radiative .eqv. .true.) then
      if (tr(i).gt.kr) then
       D(i) = D(i) + (tr(i)**3.d0) * 0.37d-9
      end if
     end if 
    end do

  D(nsub)=D(nsub-1)

!  q is the source terms and boundary conditions

  call matrix(D,nsub,bckind,dx,bc1,bc2,ao,q,p1,p2,p3)

!   initial solve of system 
!   solve system using tridiagonal matrix algorithm
    call tri(nsub,p1,p2,p3,q,tr1) !n,a,b,c,q,u
    diff=abs(tr1-tr)
    tr(1:nsub)=tr1(1:nsub)
    iter=sum(diff)
    print *, 'jaupart convergence ', iter
    
end do    

!-------------------------------------------------------------------------------
! store Jaupart solution as trj

do i=2,nsub-1
  hf(i)=((D(i+1)+D(i)+D(i-1))/3.d0) * (tr1(i+1) - tr1(i-1))/(2.d0*dx)  
end do
hf(1)=hf(2)
hf(nsub)=hf(nsub-1)

dj=D
hfj=hf
hf=zero
D=zero
trj=tr1
!________________________________________________________________________________

!________________________________________________________________________________
! ------ this is the Chapman algorithm - depth and pressure dependent ---------
! coefficients c - for pressure, is constant over all depth
! b for temperature has a jump at the depth of the top of the lower crust
! ------------------------------------------------------------------------------

iter = 10.d0
do while (iter.gt.tol)
    c=1.5d-6
    b=1.5d-3
    do i=1,nsub
     depth=dx*(i-1)  
      if (depth.gt.lc) then
       b=1.d-4
      end if
     D(i) = ko(i)*(1+(c*dx*(i-1)))/(1+b*(tr(i)-273.d0))
    end do

  D(nsub)=D(nsub-1)

!  q is the source terms and boundary conditions

  call matrix(D,nsub,bckind,dx,bc1,bc2,ao,q,p1,p2,p3)

!   initial solve of system 
!   solve system using tridiagonal matrix algorithm
    call tri(nsub,p1,p2,p3,q,tr1) !n,a,b,c,q,u
    diff=abs(tr1-tr)
    tr(1:nsub)=tr1(1:nsub)
    iter=sum(diff)
    print *, 'chapman convergence ', iter
    
end do    

! store solution as trc
do i=2,nsub-1
  hf(i)=((D(i+1)+D(i)+D(i-1))/3.d0) * (tr1(i+1) - tr1(i-1))/(2.d0*dx)
end do
hf(1)=hf(2)
hf(nsub)=hf(nsub-1)
dc=D
D=zero
hfc=hf
hf=zero
trc = tr1

! write results
! -----------------------------------------------------------------------------
    open (14, file='out1.dat') ! geotherms
    open (15, file='out2.dat') ! heat flow curves
    open (13, file='out3.dat') ! thermal conductivity profiles

    i=1

 200 continue
    write (15, '(4(1pe16.7))') dx*(i-1)/1000, hfo(i), hfj(i), hfc(i)    
    write (14,'(4(1pe16.7))') dx*(i-1)/1000, tro(i) -273.d0, trj(i) -273.d0, trc(i)-273.d0
    write (13, '(4(1pe16.7))') dx*(i-1)/1000, D0(i), dj(i), dc(i)       
    i=i+1
    if (i.lt.nsub+1) then
      goto 200
    else
      continue
      close (14)
      close (13)
      close (15)
    end if


end program

! *********************************************************************************************
! SUBROUTINES
! *********************************************************************************************
! subtroutine matrix
! set the matrix diagonals and boundary conditions for the tdma matrix solver

subroutine matrix(D,nsub,bckind,dx,bc1,bc2,ao,q,p1,p2,p3)

implicit none
integer :: nsub,nv,i,j
double precision, parameter :: zero=0.0d0
double precision :: D(nsub),ao(nsub),q(nsub),p1(nsub),p2(nsub),p3(nsub),rho(nsub),cp(nsub),v(nsub)
double precision :: dx, dx2, bc1, bc2
logical :: bckind
dx2 = dx**2.d0
p1 = zero ; p2= zero ; p3=zero 

! matrix diagonals
! p1 = principal (a), p2 = upper (b), p3 = lower (c)
! a,b,c correspond to form of tridiagonal solver in module_tdma.f90
!-----------------------------------------------------------------------------------------------
 do i=2,nsub-1
  p2(i)=(0.5d0/dx2)*(D(i)+D(i+1)) ! upper = b, b(1) = part of bc, or zero
 end do

! p2(1) - part of boundary condition

 do i=2,nsub-1 
  p1(i)=(-0.5d0/dx2)*(D(i+1)+2.d0*D(i)+D(i-1))  ! principal = a, a(1) = bc, a(N) = bc 
 end do

! p1(1) part of boundary condition
! p1(nsub) part of boundary condition

 do i=2,nsub
  p3(i)=(0.5d0/dx2)*(D(i)+D(i-1))  ! lower = c, c(1)= 0, c(N) = part of bc or zero
 end do

p1(nsub)=p1(nsub-1)
p2(nsub)=p2(nsub-1)

! ---------------------------------------------------------------------------------------------------
! boundary conditions: choice between top dirichlet and base either dirichlet or neumann
! neumann at base = +ve value = heat flow in at base
! bc temps in °C 
! ---------------------------------------------------------------------------------------------------
if (bckind .eqv. .false.) then

!! this sets dirichlet boundary conditions on both sides of model (bckind = false)
!!---------------------------------------------------------------------------------------------------

p1(1) = 1.d0 !a(1)
p1(nsub) = 1.d0 !a(n)
p2(1) = 0.d0 !b(1)
p3(nsub) = 0.d0 !c(n)

!! and
!! ***


q(1)=bc1 + 273.d0
q(nsub)=bc2 + 273.d0

print *, 'dirichlet selected'

elseif (bckind .eqv. .true.) then

!! this sets basal neumann boundary condition (bckind = true)
!!-------------------------------------------------------------------------------------------------------

p1(1) = 1.d0 !a(1) dirichlet, z=0
p2(1) = 0.d0 !b(1) dirichlet, z=0
p1(nsub) = 1.d0/dx !a(n) neumann, z=B, u(n) 
p3(nsub) = -1.d0/dx !c(n) numeann, z=B, u(n-1)


!! and
!! ***

q(1)=bc1 + 273.d0
q(nsub)=bc2/D(nsub)

print *, 'neumann selected'

end if


end subroutine



