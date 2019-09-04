subroutine chkmott(nts,uniform,iterc,mott)
  use gf2, only:nao,nk,ncell,iwmax,Sigma_tau,omegas,nc,nleg
  use gf2, only:sig11
  implicit none
  integer :: nts, uniform, iterc
  double precision :: sdiff, sdiff1 
  complex*16, dimension(:,:,:),allocatable :: rsig,isig,sigr,sigma_w,sigma_tau_k
  complex*16, dimension(:,:,:,:),allocatable :: sigl
  character*256 file_name
  integer:: i,j,k,ic,w,kk,t,ip,ii
  logical mott

  allocate(sigma_w(nao,nao,2*nk))
  allocate(rsig(nao,nao,0:nleg),isig(nao,nao,0:nleg),sigl(nao,nao,0:nleg,2*nk))
  allocate(sigr(nao,nao,ncell))
  allocate(sigma_tau_k(nao,nao,nts))

  do k=1,2*nk

     do t = 1, nts
       call FT_r_to_single_k(sigma_tau_k(:,:,t),Sigma_tau(:,:,:,t),k)
     enddo
     call calc_sigl_from_sigtau(nts,uniform,rsig,dble(sigma_tau_k(:,:,:)))
     call calc_sigl_from_sigtau(nts,uniform,isig,dimag(sigma_tau_k(:,:,:)))


     sigl(:,:,:,k)=rsig(:,:,:)+dcmplx(0.0d0,1.0d0)*isig(:,:,:)
  end do

  write(file_name,fmt='(A6,I0,A4)') 'sig11_',iterc,'.txt'
  open(950,file=file_name,status='replace')

  do w = 1, iwmax
    call print_se11_inner(w,nts,uniform,sigl,sig11)
  end do
  call flush(6)
  close(950)

  deallocate(sigma_w)
  deallocate(sigr)
  deallocate(rsig,isig,sigl)

  mott = .true.
  do w = 2, 100
    if(sig11(w-1).ge.sig11(w)) then
      mott = .false.
      write(*,*)'Mott behavior violated, 1st criterion, at iter ', iterc
      exit
    endif
  enddo


  if(mott) then
  sdiff = abs(sig11(2) - sig11(1))
  do w = 3, 101
    sdiff1 = abs(sig11(w) - sig11(w-1))
    If(sdiff1.ge.sdiff) then
      mott = .false.
      write(*,*)'Mott behavior violated, 2nd criterion, at iter ', iterc
      exit
    Else
      sdiff = sdiff1
    Endif
  enddo
  endif  

end subroutine chkmott

subroutine print_sigma_11(nts,uniform,iterc,iters,Sigma_tau_prev,E_diff,damp,prtsig,mott)
  use gf2, only:nao,nk,ncell,iwmax,Sigma_tau,omegas,nc,nleg
  use gf2, only:sig11,Sigma_tau_tmp
  implicit none
  integer::nts,uniform,iterc,iters
  double precision :: E_diff,damp,sdiff,sdiff1
  complex*16, dimension(:,:,:),allocatable :: rsig,isig,sigr,sigma_w,sigma_tau_k
  complex*16, dimension(:,:,:,:),allocatable :: sigl
  double precision,dimension(nao,nao,ncell,nts) :: Sigma_tau_prev

  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer:: i,j,k,ic,w,kk,t,ip,ii
  character*256 file_name
  logical mott, prtsig, update
  allocate(sigma_w(nao,nao,2*nk))
  allocate(rsig(nao,nao,0:nleg),isig(nao,nao,0:nleg),sigl(nao,nao,0:nleg,2*nk))

  allocate(sigr(nao,nao,ncell))
  allocate(sigma_tau_k(nao,nao,nts))

  update = .false.
  !If(prtsig.and.mott) update = .true.
  If(prtsig) update = .true.


  do k=1,2*nk

     do t = 1, nts
       call FT_r_to_single_k(sigma_tau_k(:,:,t),Sigma_tau(:,:,:,t),k)
     enddo
     call calc_sigl_from_sigtau(nts,uniform,rsig,dble(sigma_tau_k(:,:,:)))
     call calc_sigl_from_sigtau(nts,uniform,isig,dimag(sigma_tau_k(:,:,:)))


     sigl(:,:,:,k)=rsig(:,:,:)+dcmplx(0.0d0,1.0d0)*isig(:,:,:)
  end do

  If(prtsig) then
  write(file_name,fmt='(A6,I0,A4)') 'sig11_',iterc,'.txt'
  open(950,file=file_name,status='replace')
  
  do w = 1, iwmax
    call print_se11_inner(w,nts,uniform,sigl,sig11)
  end do
  call flush(6)
  close(950)
  Endif
  
  call flush(6)
!  close(950) 
  deallocate(sigma_w)
  deallocate(sigr)
  deallocate(rsig,isig,sigl)
 
  mott = .true. 
  do w = 2, 100
    if(sig11(w-1).gt.sig11(w)) then
      mott = .false.
      write(*,*)'Mott behavior violated, 1st criterion, at iter ', iterc, prtsig
      exit
    endif
  enddo
  

  if(mott.and.iterc.ge.0.and.prtsig) then
  sdiff = abs(sig11(2) - sig11(1))
  do w = 3, 101
    sdiff1 = abs(sig11(w) - sig11(w-1))
    If(sdiff1.ge.sdiff) then
      mott = .false.
      write(*,*)'Mott behavior violated, 2nd criterion, at iter ', iterc
      exit
    Else
      sdiff = sdiff1
    Endif
  enddo  
  endif

  if(update.and.mott.and.iterc.ge.10) then 
    write(*,*)'update Sigma_tau_prev at iteration ', iterc
    Sigma_tau_prev(:,:,:,:) = Sigma_tau(:,:,:,:)
    damp = 7.5d-01
  endif

end subroutine print_sigma_11 


subroutine print_se11_inner(w,nts,uniform,sigl,sig11)
  use gf2, only:nao,nk,ncell,iwmax,nleg,omegas
  implicit none
  integer:: w,nts,uniform
  complex*16, dimension(:,:),allocatable :: sigrw
  complex*16, dimension(:,:,:),allocatable :: sigma_w,sigma
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  double precision, dimension(iwmax) :: sig11
  integer:: i,j,k,ic,kk,t

  allocate(sigrw(nao,nao))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(sigma(nao,nao,nts))
!
  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))
  end do
  call FT_k_to_r_ic(sigma_w(:,:,:), sigrw(:,:),1) 
  write(950,*) w, dimag(sigrw(1,1))
  sig11(w) = dimag(sigrw(1,1))
  deallocate(sigrw)
  deallocate(sigma_w)
  deallocate(sigma)
!

end subroutine print_se11_inner 
