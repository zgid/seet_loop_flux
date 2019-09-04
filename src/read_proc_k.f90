subroutine read_process(F_kk,S_kk)
  use gf2, only : nao, nk
  use gf2, only : F_k, S_k
  implicit none
  integer nn
  complex*16, dimension(nao,nao,2*nk) :: S_kk,F_kk
  write(*,*)"before read_process SFK"
  call flush(6)  
  call read_SFk
  write(*,*)"before read_process four coeff gaussian"
  call flush(6)  
  call read_four_coeff_gaussian
  write(*,*)"before read_process read_k_weights_gaussian"
  call flush(6)  
  call read_k_weights_gaussian
  write(*,*)"before read_process matrices"
  call flush(6)  
  do nn = 1, nk
    S_kk(:,:,2*nn-1) = S_k(:,:,nn)
    S_kk(:,:,2*nn) = dconjg(S_k(:,:,nn))

    F_kk(:,:,2*nn-1) = F_k(:,:,nn)
    F_kk(:,:,2*nn) = dconjg(F_k(:,:,nn))
  end do
  write(*,*)"before read_process end"
  call flush(6)  
end subroutine read_process

subroutine update_Fk(F_kk)
  use gf2, only : nao, nk
  use gf2, only : F_k, S_k
  implicit none
  integer nn
  integer k,  i,  j
  integer kd, id, jd
  double precision x,y
  complex*16, dimension(nao,nao,2*nk) :: F_kk
  F_k = dcmplx(0.0d0,0.0d0)
  F_kk = dcmplx(0.0d0,0.0d0)
  open(42,file='Fock-k.txt',status='old')
  rewind(42)
  kfloop : do k=1,nk
     ifloop : do i=1,nao
        jfloop : do j=1,i
           read(42,*)kd, id, jd, x, y
          ! write(*,*)"L",kd, id, jd, x, y
           !call flush(6)
           F_k(i, j, kd) = dcmplx(x, y)
           F_k(j, i, kd) = dcmplx(x,-y)
        enddo jfloop
     enddo ifloop
  enddo kfloop
  close(42)
  do nn = 1, nk
    F_kk(:,:,2*nn-1) = F_k(:,:,nn)
    F_kk(:,:,2*nn) = dconjg(F_k(:,:,nn))
  end do
end subroutine update_Fk

subroutine read_SFk
  use gf2, only : nk,  nao, ncell, F_k, S_k, vertex
  implicit none 
  integer k,  i,  j, ic, k_read
  integer kd, id, jd
  double precision x,y, diff

  write(*,*)"before allocate"
  call flush(6)  
  if (.not. allocated(S_k)) allocate(S_k(nao, nao, nk))
  if (.not. allocated(F_k)) allocate(F_k(nao, nao, nk))
    write(*,*)"after allocate"
  call flush(6)  
  open(41,file='S-k.txt',status='old')
  rewind(41)
  write(*,*)"before overlap read"
  call flush(6)  
  kloop : do k=1,nk
     iloop : do i=1,nao
        jloop : do j=1,i
           read(41,*)kd, id, jd, x, y
           S_k(i, j, kd) = dcmplx(x,  y) 
           S_k(j, i, kd) = dcmplx(x, -y)
        enddo jloop
     enddo iloop
  enddo kloop
  write(*,*)"after overlap read"
  call flush(6)  
  write(*,*)"before fock read"
  call flush(6)  
  close(41)
  open(42,file='Fock-k.txt',status='old')
  rewind(42)
  kfloop : do k=1,nk
     ifloop : do i=1,nao
        jfloop : do j=1,i
           read(42,*)kd, id, jd, x, y
           !write(*,*)"L2",kd, id, jd, x, y
           !call flush(6)
           F_k(i, j, kd) = dcmplx(x, y)
           F_k(j, i, kd) = dcmplx(x,-y)
        enddo jfloop
     enddo ifloop
  enddo kfloop
  close(42)
  write(*,*)"after fock read"
  call flush(6)  
  call check_hermiticity
  write(*,*)"after herm check"
  call flush(6)  
contains
  subroutine check_S_k
    implicit none    
    kvertex : if(vertex(k)) then
       if(abs(y)>1d-10) then
          write(6,*)' problem, S_k at k=',k,' should be real but for i,j ',i,j,' it is ',S_k(i,j,k)
          call flush(6)
          stop
       endif
    end if kvertex
  end subroutine check_S_k

  subroutine check_F_k
    implicit none
    
    kvertex : if(vertex(k)) then
       if(abs(y)>1d-10) then
          write(6,*)' problem, F_k at k=',k,' should be real but for i,j ',i,j,' it is ',F_k(i,j,k)
          call flush(6)
          stop
       endif
    end if kvertex
  end subroutine check_F_k

  subroutine check_hermiticity
    implicit none
    do k=1, nk
       diff = maxval(abs(S_k(:,:,k) - transpose(dconjg(S_k(:,:,k)))))
       if(diff>1d-7) then
          write(6,*) 'problem, for k=',k,' S diff from hermitcity ',diff
       endif
    enddo
    
    do k=1, nk
       diff = maxval(abs(F_k(:,:,k) - transpose(dconjg(F_k(:,:,k)))))
       if(diff>1d-7) then
          write(6,*) 'problem, for k=',k,' S diff from hermitcity ',diff
       endif
    enddo
  end subroutine check_hermiticity
end subroutine read_SFk

 subroutine read_four_coeff_gaussian
   use gf2, only : nk, ncell, ncell, four_coef 
   implicit none
   integer:: k, kd, icel, iceld
   double precision x,y
   integer kk, ic

   if(.not.allocated(four_coef)) then;  allocate(four_coef(nk, ncell))
   end if

   open(43,file='Fourier_factors.txt',status='old')
   rewind(43)

   do k=1,nk
      do icel=1,ncell
         read(43,*)kk,ic,x,y
         four_coef(kk,ic)=cmplx(x,y)
      end do
   end do

 end subroutine read_four_coeff_gaussian


 subroutine read_k_weights_gaussian
   use gf2, only : nk, k_weights 
   use gf2, only : nkpw, w_func
   implicit none
   integer k, kd, icel, iceld,np,i
   double precision x,y,nrmlz,ii
   double precision,dimension(:),allocatable :: k_w


   open(44,file='k_weights.txt ',status='old')

   np=0  
   rewind(44)
   do  
      read(44,*,end=11)ii
      np=np+1
   end do
11  continue
   
   if (np .ne. nk) then
      write(*,*)"k_weights_gaussian np not equal nk"
      stop
   end if

   if(.not.allocated(k_weights)) then;  allocate(k_weights(nk))
   end if
   if(.not.allocated(w_func)) then; allocate(w_func(nk)) 
   end if
   if(.not.allocated(k_w)) then;  allocate(k_w(nk))
   end if
   
   w_func(:)=0.0d0

   rewind(44)                                                                  
  
   do k=1,nk
      read(44,*)x                                                              
      k_weights(k)=x
      k_w(k)=x
   end do
 
   nrmlz=1.0d0/k_weights(1)
   nkpw=0.0d0
   do k=1,nk
     k_w(k)=k_w(k)*nrmlz

     if (abs(k_w(k)-1.0d0)>1.0e-7) then
        w_func(k) = 1.0d0
     end if
     nkpw = nkpw + k_w(k)
   end do

   nkpw = 1.0d0/nkpw
 
 deallocate(k_w)
 end subroutine read_k_weights_gaussian

subroutine read_tnl
  use gf2, only:nleg,iwmax,tnl
  implicit none
  integer::i,n,l
  double precision:: valr,vali

  open(47,file='tnl.dat',status='old')
  rewind(47)
  do i=1,nleg*iwmax
     read(47,*)l,n,valr,vali
     tnl(n+1,l+1)=dcmplx(valr,vali)
  end do
  close(47)

end subroutine read_tnl

