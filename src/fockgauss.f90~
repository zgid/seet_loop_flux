subroutine fock_gauss(dm_r,E_S0,iterc)
  use gf2, only: nao,ncell,nk
  implicit none
  complex*16, dimension(nao,nao,ncell)::dm_r
  double precision, dimension(:,:,:),allocatable::Pr_0g

  integer::i,j,k,l,iterc
  integer::nf,ic,jj
  logical::exist
  double precision val, E_S0, E_HF,nel
!  INTEGER OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!  integer::n_threads

  allocate(Pr_0g(nao,nao,ncell))


  do ic=1,ncell
      do i=1,nao
         do j=1,nao            
            Pr_0g(i,j,ic) = 0.50d0*dble(dm_r(i,j,ic))
!            write(*,*)"D2",ic,i,j,0.50d0*dble(dm_r(i,j,ic)) 
!            call flush(6) 
         end do
      end do
   end do

  inquire(file="P_real_long1.txt",exist=exist)
  If(exist) then
  open(1,file="P_real_long1.txt",status="old")
  else
  open(1,file="P_real_long1.txt",status="new")
  Endif

  rewind(1)
  do ic = 1, ncell
    do i = 1, nao
      do j = 1, i
        write(1,"(F18.14)") Pr_0g(i,j,ic)
      enddo
    enddo
  enddo
  close(1)

 
  call system("cd gaussian_tmp; cp ../P_real_long1.txt P_real_long.txt; g09 recontract.com; cp Fock-real.txt ../; cp Energy.txt ../; cp Fock-k.txt ../")
  open(1,file='Energy.txt',status='old')
  rewind(1)
  read(1,*) E_S0
  close(1)
end subroutine fock_gauss 

