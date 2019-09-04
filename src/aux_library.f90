subroutine invert_complex_matrix(ndim,matrix)
  implicit none
  integer:: ndim
  integer, dimension(:), allocatable:: ipiv
  complex*16,dimension(ndim,ndim)::matrix
  integer::lwork
  complex*16, dimension(:),allocatable::work
  integer:: io
  lwork=ndim*ndim
  allocate(ipiv(ndim))
  allocate(work(lwork))

  ipiv(:)=0

  call ZGETRF(ndim,ndim,matrix,ndim,ipiv(:),io)
  call ZGETRI(ndim,matrix,ndim,ipiv(:),work(1),lwork,io)

  if (io .ne. 0) then
     write(*,*)"something went wrong with the complex matrix inversion(invert_complex_matrix containing ZGETRF and ZGETRI)"
     call flush(6)
     stop
  end if

  deallocate(ipiv)
  deallocate(work)

end subroutine invert_complex_matrix

subroutine char4(itau,cht)
  implicit none
  integer itau
  character(len=4)cht
  character(len=4)chaux
  if (itau<10) then
     write(chaux,'(I1)')itau
     cht='000'//trim(chaux)
  end if
  if (itau>=10 .and. itau<100) then
     write(chaux,'(I2)')itau
     cht='00'//trim(chaux)
  end if
  if (itau>=100 .and. itau<1000) then
     write(chaux,'(I3)')itau
     cht='0'//trim(chaux)
  end if
  if (itau>=1000 .and. itau<10000) then
     write(chaux,'(I4)')itau
     cht=trim(chaux)
  end if
end subroutine char4

subroutine compute_ortho_basis_complex(overlap,xmat,ndim,eigval)
  implicit none
  integer:: ndim
  complex*16,dimension(ndim,ndim)::overlap,xmat
  complex*16,dimension(:,:),allocatable::aux
  complex*16,dimension(:,:),allocatable::aux1
  integer::lwork,LRWORK,LIWORK
  complex*16, dimension(:),allocatable::work
  double precision, dimension(:),allocatable::rwork
  double precision, dimension(ndim)::eigval
  integer, dimension(:),allocatable::iwork
  integer:: io
  complex*16::alpha,betaz
  integer::i,j
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)
  allocate(aux(ndim,ndim),aux1(ndim,ndim))
  lwork=2*ndim + ndim**2
  allocate(work(lwork))
  LRWORK=1 + 5*ndim + 2*ndim**2
  LIWORK=3 + 5*ndim
  allocate(rwork(LRWORK))
  allocate(iwork(LIWORK))
  aux(:,:)=overlap(:,:)

  eigval(:)=0.0
  call ZHEEVD('V','U',ndim,aux,ndim,eigval,work(1),lwork,rwork(1),LRWORK, IWORK, LIWORK,io)
!ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
!     $                   LRWORK, IWORK, LIWORK, INFO )
! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
! complex Hermitian matrix A.

  if (io .ne. 0) then
     write(*,*)"something went wrong with compute_ortho_basis_complex (ZHEEVD)"
     call flush(6)
     stop
  end if


  aux1(:,:)=0.0d0
  do i=1,ndim
!     write(*,*)"eee",i,eigval(i)
    aux1(i,i)=1./sqrt(eigval(i))
    do j=1,ndim
      xmat(j,i)=aux(j,i)*aux1(i,i)
    end do
  end do
  deallocate(aux,aux1)
  deallocate(work)
  deallocate(rwork,iwork)

end subroutine compute_ortho_basis_complex

subroutine compute_diag_basis_complex(overlap,matrix,xmat,ndim)
  implicit none
  integer:: ndim
  complex*16,dimension(ndim,ndim)::overlap,xmat,matrix
  complex*16,dimension(:,:),allocatable::aux
  complex*16,dimension(:,:),allocatable::aux1
  integer::lwork,LRWORK,LIWORK
  complex*16, dimension(:),allocatable::work
  double precision, dimension(:),allocatable::eigval,rwork
  integer, dimension(:),allocatable::iwork
  integer:: io
  complex*16::alpha,betaz
  integer::i,j
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)

  allocate(aux(ndim,ndim),aux1(ndim,ndim))
  allocate(eigval(ndim))
  lwork=2*ndim + ndim**2
  allocate(work(lwork))
  LRWORK=1 + 5*ndim + 2*ndim**2
  LIWORK=3 + 5*ndim  
  allocate(rwork(LRWORK))
  allocate(iwork(LIWORK))
  aux(:,:)=matrix(:,:) ! Fock matrix

  call compute_ortho_basis_complex(overlap,xmat,ndim) ! Lowdin: xmat
  aux1(:,:)=dcmplx(0.0d0,0.0d0)
 call ZGEMM('t','n',ndim,ndim,ndim,alpha,conjg(xmat),ndim,aux(:,:),ndim,betaz,aux1,ndim)
 call ZGEMM('n','n',ndim,ndim,ndim,alpha,aux1,ndim,xmat,ndim,betaz,aux(:,:),ndim)
 ! Fock transformed to Lowdin basis
 aux1(:,:)=xmat(:,:)

  eigval(:)=0.0
  call ZHEEVD('V','U',ndim,aux,ndim,eigval,work(1),lwork,rwork(1),LRWORK, IWORK,LIWORK,io)
!ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
!     $                   LRWORK, IWORK, LIWORK, INFO )
! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
! complex Hermitian matrix A.
  if (io .ne. 0) then
     write(*,*)"something went wrong with compute_diag_basis_complex (ZHEEVD)"
     call flush(6)
     stop
  end if

 call ZGEMM('n','n',ndim,ndim,ndim,alpha,aux1,ndim,aux,ndim,betaz,xmat(:,:),ndim) ! xmat diagonalizes Fock

  deallocate(aux,aux1)
  deallocate(eigval)
  deallocate(work)
  deallocate(rwork,iwork)

end subroutine compute_diag_basis_complex


subroutine ft_w_tau_inner(matrix_in,matrix_out,ndim,w,tau,beta)
  implicit none
  integer::ndim
  integer::w
  double precision::tau,beta
  complex*16,dimension(ndim,ndim)::matrix_in,matrix_out
  double precision::cs,sn
  
  cs=cos(-1.0d0*(2.0d0*(w-1)+1.0d0)*tau)
  sn=sin(-1.0d0*(2.0d0*(w-1)+1.0d0)*tau)
  
  matrix_out(:,:)=matrix_out(:,:)+2.0/beta*(cs*dble(matrix_in(:,:))-sn*dimag(matrix_in(:,:))) 

end subroutine ft_w_tau_inner


subroutine open_file_ip_write(k,ip)
  implicit none
  integer k,ip
  character(len=4)chk
!  write(*,*)"IP open write",ip,k
  call flush(6)
  call char4(k,chk)   ! writing k into character*4
  open(ip,file='grtbin.'//chk,status='replace',form='unformatted')
!  open(ip,file='grtbin.'//chk,status='new',form='unformatted')
  rewind(ip)
end subroutine open_file_ip_write

subroutine close_file_ip_write(k,ip)
  implicit none
  integer k,ip
  character(len=4)chk
!  write(*,*)"IP open write",ip,k
  call flush(6)
  call char4(k,chk)   ! writing k into character*4
  close(ip,status='delete')
end subroutine close_file_ip_write

subroutine check_sigma_cells(overlap)
  use gf2, only: nao,ncell,nk
  implicit none
  complex*16, dimension(nao,nao,2*nk) ::overlap
  complex*16, dimension(:,:,:),allocatable::aux
  complex*16, dimension(:,:,:),allocatable::G1
  double precision:: thrf
  integer nzc

  allocate(aux(nao,nao,2*nk))
  allocate(G1(nao,nao,ncell))
  call calc_G1(overlap,aux,G1)
  thrf=1.d-10     
  call screen_any_matrix(thrf,dble(G1),nzc)  
  write(*,*)"number of non-zero cell in slef-energy matrix is",nzc
  deallocate(aux)
  deallocate(G1)

end subroutine check_sigma_cells

subroutine screen_any_matrix(thrf,matrix,nzc)
  use gf2, only: nao,ncell
  implicit none
  integer nzc
  integer ndimf
  integer ip,i,j,k,kk,ic
  double precision:: thrf
  double precision, dimension(nao,nao,ncell)::matrix
  double precision, dimension(:), allocatable::flat_matrix
  integer, dimension(:), allocatable::flat_cell
  double precision::tmp

  ndimf=(nao*(nao+1)/2)*ncell
  allocate(flat_matrix(ndimf))
  allocate(flat_cell(ndimf))

  k=1
  do ic = 1, ncell
    do i = 1, nao
      do j = 1, i
         flat_matrix(k)=matrix(i,j,ic)
         flat_cell(k)=ic
         k=k+1
      enddo
    enddo
  enddo


  do k=ndimf,1,-1
     if (abs(flat_matrix(k))>=thrf) then
        nzc=flat_cell(k)
        exit
     end if
  end do

  write(*,*)"number of non-zero cells in matrix is ",nzc
  deallocate(flat_matrix)
  deallocate(flat_cell)
end subroutine screen_any_matrix

subroutine screen_real_overlap_matrix(thrf)

  use gf2, only: nao,ncell,ncs
  implicit none
  integer nzcf
  integer ndimf
  integer ip,i,k,kk,ic
  double precision:: thrf
  double precision, dimension(:), allocatable::flat_fock
  integer, dimension(:), allocatable::flat_cell
  double precision::tmp

  ndimf=(nao*(nao+1)/2)*ncell
  allocate(flat_fock(ndimf))
  allocate(flat_cell(ndimf))

  ip=3111
  open(ip,file='S-real.txt',status='old')
  rewind(ip)

  k=1
  do i=1,ndimf
     read(ip,*)ic,kk,tmp
     flat_fock(k)=tmp
     flat_cell(k)=ic
     k=k+1
  end do

  close(ip)


  do k=ndimf,1,-1
     if (abs(flat_fock(k))>=thrf) then
        ncs=flat_cell(k)
        exit
     end if
  end do

  write(*,*)"number of non-zero cells in overlap matrix is ",ncs
  deallocate(flat_fock)
  deallocate(flat_cell)
end subroutine screen_real_overlap_matrix


subroutine symmetric_ortho_basis(overlap,xmat,ndim)
  implicit none
  integer:: ndim
  complex*16,dimension(ndim,ndim)::overlap,xmat
  complex*16,dimension(:,:),allocatable::aux
  complex*16,dimension(:,:),allocatable::aux1,aux2
  integer::lwork,LRWORK,LIWORK
  complex*16, dimension(:),allocatable::work
  double precision, dimension(:),allocatable::eigval,rwork
  integer, dimension(:),allocatable::iwork
  integer:: io
  complex*16::alpha,betaz
  integer::i,j
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)
  allocate(aux(ndim,ndim),aux1(ndim,ndim),aux2(ndim,ndim))
  allocate(eigval(ndim))
  lwork=2*ndim + ndim**2
  allocate(work(lwork))
  LRWORK=1 + 5*ndim + 2*ndim**2
  LIWORK=3 + 5*ndim
  allocate(rwork(LRWORK))
  allocate(iwork(LIWORK))
  aux(:,:)=overlap(:,:)

  eigval(:)=0.0
  call ZHEEVD('V','U',ndim,aux,ndim,eigval,work(1),lwork,rwork(1),LRWORK, IWORK, LIWORK,io)
!ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,                                                                                                                                                                                                                                      
!     $                   LRWORK, IWORK, LIWORK, INFO )                                                                                                                                                                                                                                     
! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a                                                                                                                                                                                                                         
! complex Hermitian matrix A.                                                                                                                                                                                                                                                               

  if (io .ne. 0) then
     write(*,*)"something went wrong with compute_ortho_basis_complex (ZHEEVD)"
     call flush(6)
     stop
  end if


  aux1(:,:)=0.0d0
  do i=1,ndim
    aux1(i,i)=1.0d0/dsqrt(eigval(i))
    write(*,*)"eigval",i,eigval(i)
    do j=1,ndim
      xmat(j,i)=aux(j,i)*aux1(i,i)
      aux2(i,j) = dconjg(aux(j,i))
    end do
  end do

  xmat = matmul(xmat,aux2)
  deallocate(aux,aux1,aux2)
  deallocate(eigval)
  deallocate(work)
  deallocate(rwork,iwork)

end subroutine symmetric_ortho_basis

