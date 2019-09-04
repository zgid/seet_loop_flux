subroutine inigrd(nts)
  use gf2, only : nuniform
  use gf2, only : power_mesh, power_mesh_small, power_mesh_small_int,max_points_mesh
  implicit none
  integer :: nts

  call initialize_omegas
  call initialize_tau_grid
  nts=-10
  !write(*,*)nuniform  
 call intitialize_small_tau_grid(nts)
  write(*,*)"Number of points in the small grid", nts,allocated(power_mesh_small),allocated(power_mesh_small_int)
  call flush(6)
  if (.not. allocated(power_mesh_small)) then
     allocate(power_mesh_small(nts))
  end if
  if (.not. allocated(power_mesh_small_int)) then
     allocate(power_mesh_small_int(nts))
  end if
  write(*,*)"Number of points in the small grid", nts,allocated(power_mesh_small),allocated(power_mesh_small_int),max_points_mesh,allocated(power_mesh),power_mesh(max_points_mesh)
  call flush(6)
  call intitialize_small_tau_grid(nts)
end subroutine inigrd

subroutine initialize_omegas
  use gf2, only : iwmax, beta, omegas
  implicit none
  integer :: i
  double precision, save :: pi = 4.*atan(1.0d0)

  if(.not.allocated(omegas)) allocate(omegas(iwmax))

  do i=1,iwmax
     omegas(i)=((2.0d0*(i-1)+1)*pi* dcmplx(0.,1.))/beta
  end do

end subroutine initialize_omegas

subroutine initialize_tau_grid
  use gf2, only: beta, power_tau,uniform_tau,max_points_mesh,power_mesh
  implicit none
  integer:: max_power_points
  max_power_points=2*(power_tau+1)+1
  max_points_mesh=(uniform_tau)*(max_power_points-1)+1
  write(*,*)"max_points_tau_mesh",max_points_mesh
  call flush(6)
  if (.not. allocated(power_mesh)) then
     allocate(power_mesh(max_points_mesh))
  end if

  call generate_mesh(power_mesh,max_points_mesh,power_tau, uniform_tau, beta)
end subroutine initialize_tau_grid

subroutine intitialize_small_tau_grid(nts)
  use gf2, only : max_points_mesh, power_mesh, power_mesh_small, power_mesh_small_int, uniform_tau, power_tau,nuniform
  implicit none
  integer:: bbd,pstep,ntaul,index1
  integer:: nts
  integer:: p,n,s


  bbd=nuniform
  pstep=uniform_tau/nuniform
  write(*,*)"all",allocated(power_mesh_small),allocated(power_mesh)

  if (nts>0) then
    write(*,*)"sizeof power mesh: ",sizeof(power_mesh_small_int)/4
    call flush(6)
    power_mesh_small_int(1)=8
    power_mesh_small_int(353)=5
    write(*,*)"value at 1: ",power_mesh_small_int(1)
    write(*,*)"value last: ",power_mesh_small_int(353)
     n=0
     do p=1,2*power_tau+2
        ntaul=(p-1)*uniform_tau


        if (p==(2*power_tau+2)) bbd=nuniform+1
               
        do s=1,bbd
           index1=ntaul+(s-1)*pstep+1
           !write(*,*)n,index1,power_mesh(index1) 
           !write(*,*)n
           !call flush(6) 
           n=n+1
           !write(*,*)n
           !if(n > 353) write(*,*)"problem with n",n
          power_mesh_small_int(n)=index1
          power_mesh_small(n)=power_mesh(index1)  
        end do
     end do
  else
     nts=0
     do p=1,2*power_tau+2
        ntaul=(p-1)*uniform_tau

        
        if (p==(2*power_tau+2)) bbd=nuniform+1

        do s=1,bbd
           index1=ntaul+(s-1)*pstep+1
           nts=nts+1
        end do
     end do
  end if


end subroutine intitialize_small_tau_grid

subroutine generate_mesh(power_mesh, max_points_mesh,power, uniform, beta)     
!generates a power mesh with 2*(power+1)+1 intervals that are further divided
!into uniform intervals
!total number of points in the grid is
!(2*(power+1)+1)*uniform=(2*power+3)*uniform
implicit none
integer:: max_points_mesh
double precision:: power_mesh(max_points_mesh)
integer:: power, uniform
double precision:: beta
double precision, dimension(:), allocatable:: power_points,weights
integer:: max_power_points
integer:: i,j,k,l
double precision:: dtau

max_power_points=2*(power+1)+1

if (.not.allocated(power_points)) allocate(power_points(max_power_points))


power_points(1)=0.
k=2
do i=power,0,-1
   power_points(1+power-i+1)=beta*0.5d0*2**(-i*1.)
   k=k+1
end do
do j=power,1,-1
   power_points(k)=beta*(1.0d0-0.5d0*2**(-j*1.))
   k=k+1
end do


power_points(max_power_points)=beta


call Bubble(power_points,max_power_points)

if (mod(uniform,2)/=0)then
   write(*,*)"Simpson weights in power grid: please choose even uniform spacing."
   call flush(6)
   stop
else
   k=0
   do i=0,max_power_points-2
         dtau=(power_points(i+2)-power_points(i+1))/(uniform*1.)
      do j=0,uniform-1
         power_mesh(k+1)=power_points(i+1)+dtau*j
         k=k+1
      end do
   end do
end if
power_mesh(max_points_mesh)=power_points(max_power_points)
deallocate(power_points)
end subroutine generate_mesh

subroutine Order(p,q)
implicit none
double precision:: p,q,temp
  if (p>q) then
    temp=p
    p=q
    q=temp
  end if
  return
end subroutine Order

subroutine Bubble(A, n)
  implicit none
  integer n,i,j
  double precision :: A(1:n),p,q
  do i=1, n
     do j=n, i+1, -1
        ! call order(A(j-1),A(j))
        p = A(j-1)
        q = A(j)
        call Order(p,q)
        A(j-1) = p
        A(j)   = q
     end do
  end do
end subroutine Bubble

subroutine Bubble_int(A,A_int,n)
  implicit none
  integer n,i,j
  integer A_int(1:n)
  double precision :: A(1:n),p,q

  do i=1, n
     do j=n, i+1, -1
        p = A(j-1)
        q = A(j)
        A_int(j)=j
        A_int(j-1)=j-1
        call Order(p,q)
        A(j-1) = p
        A(j)   = q
        A_int(j-1)=j
        A_int(j)=j-1
     end do
  end do
end subroutine Bubble_int

subroutine initialize_Gr_tau(nts)
  use gf2, only: ncell,nao,Gr_tau,Sigma_tau
  integer nts

  if (.not. allocated(Gr_tau)) then
     allocate(Gr_tau(nao,nao,ncell,nts))
  endif

  if (.not. allocated(Sigma_tau)) then
     allocate(Sigma_tau(nao,nao,ncell,nts))
  endif
  Gr_tau(:,:,:,:) = 0.0d0
  Sigma_tau(:,:,:,:) = 0.0d0

end subroutine initialize_Gr_tau
