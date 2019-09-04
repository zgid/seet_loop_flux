subroutine pbc_nao(S_kk,P_0g,wannier,xmatk)
  use gf2, only: nao, ncell, diffcell, nk
  use gf2, only: four_coef, k_weights, nkpw, w_func
  implicit none
  double precision nel
  double precision, dimension(nao,nao,ncell) :: P_0g
  complex*16, allocatable, dimension(:,:,:) :: Sr_0g
  complex*16, dimension(nao,nao,2*nk) :: S_kk, xmatk

  complex*16, allocatable, dimension(:,:) :: P_tmp
  
  double precision, allocatable, dimension(:,:,:) :: ABr, ABr1, U_SAOR, U_SAORI
  

  complex*16, allocatable, dimension(:,:,:) :: U_SAO, U_SAOI, BB, UV, UVI, Pk_ao, Pk_sao, IUV, IUVI, BB1
  complex*16, dimension(:),allocatable :: work
  double precision, dimension(:),allocatable :: rwork
  double precision, dimension(:,:),allocatable :: eigval
  integer, dimension(:),allocatable :: iwork
  integer::lwork, LRWORK, LIWORK, io

  integer :: i, j, k, l, p, q, r, g, h, gg, hh, g1, nf, nfk, k1, kk, ic, gg1, ncell1

  double precision :: alpha, beta
  
  double precision :: vnorm
  double precision, dimension(ncell*nao,nao) :: wannier 


  alpha = 1.0d0
  beta =  0.0d0

  allocate(ABr(nao,nao,ncell),ABr1(nao,nao,ncell))
  allocate(Sr_0g(nao,nao,ncell))
  nf = nao*ncell

  allocate(BB(nao,nao,ncell),BB1(nao,nao,ncell))
  allocate(U_SAO(nao,nao,2*nk),UV(nao,nao,2*nk),U_SAOI(nao,nao,2*nk),UVI(nao,nao,2*nk))
  allocate(IUV(nao,nao,2*nk),IUVI(nao,nao,2*nk))
  allocate(Pk_ao(nao,nao,2*nk),Pk_sao(nao,nao,2*nk))
  allocate(U_SAOR(nao,nao,ncell),U_SAORI(nao,nao,ncell))
  allocate(eigval(nao,2*nk)) 
  allocate(P_tmp(nao,nao))
  lwork=2*nao + nao**2 
  allocate(work(lwork))
  LRWORK=1 + 5*nao + 2*nao**2
  LIWORK=3 + 5*nao
  allocate(rwork(LRWORK))
  allocate(iwork(LIWORK))

  Sr_0g(:,:,:) = dcmplx(0.0d0,0.0d0)

  do ic = 1, ncell
    call FT_k_to_r_ic(S_kk,Sr_0g(:,:,ic),ic)
  enddo
  U_SAO(:,:,:) = dcmplx(0.0d0,0.0d0) !S_kk(:,:,:)


  write(*,*)'Sk eigenvectors'
  do k = 1, 2*nk
    call symmetric_ortho_basis(S_kk(:,:,k),U_SAO(:,:,k),nao)
    U_SAOI(:,:,k) = U_SAO(:,:,k)
    call invert_complex_matrix(nao,U_SAOI(:,:,k)) 
    call FT_r_to_single_k(Pk_ao(:,:,k),P_0g,k)
    call ZGEMM('N','N',nao,nao,nao,alpha,U_SAOI(:,:,k),nao,Pk_ao(:,:,k),nao,beta,P_tmp,nao)
    call ZGEMM('N','N',nao,nao,nao,alpha,P_tmp,nao,U_SAOI(:,:,k),nao,beta,Pk_sao(:,:,k),nao) 
    call ZHEEVD('V','U',nao,Pk_sao(:,:,k),nao,eigval(:,k),work(1),lwork,rwork(1),LRWORK,IWORK,LIWORK,io)
    Pk_ao(:,:,k) = dcmplx(0.0d0,0.0d0)
    do i = 1, nao
      Pk_ao(i,i,k) = eigval(i,k)
    enddo
!    write(*,*)'k = ', k
!    write(*,*) eigval(:,k) 
    do i = 1, nao
      Pk_sao(:,i,k) = Pk_sao(:,i,k)/Pk_sao(1,i,k)
      vnorm = 0.0d0
      do j = 1, nao
        vnorm = vnorm + Pk_sao(j,i,k)*dconjg(Pk_sao(j,i,k))
      enddo      
      vnorm = sqrt(vnorm)
      Pk_sao(:,i,k) = Pk_sao(:,i,k)/vnorm
    enddo
    call ZGEMM('N','N',nao,nao,nao,alpha,U_SAO(:,:,k),nao,Pk_sao(:,:,k),nao,beta,UV(:,:,k),nao)
    call ZGEMM('C','C',nao,nao,nao,alpha,Pk_sao(:,:,k),nao,U_SAO(:,:,k),nao,beta,UVI(:,:,k),nao) 
    xmatk(:,:,k) = UV(:,:,k)
    IUV(:,:,k) = UV(:,:,k)
    IUVI(:,:,k) = UVI(:,:,k)
    call invert_complex_matrix(nao,IUV(:,:,k))
    call invert_complex_matrix(nao,IUVI(:,:,k))
  enddo

  
  BB1=dcmplx(0.0d0,0.0d0)
  write(*,*)'P real via Fourier'
  call FT_k_to_r(Pk_ao,BB1)
  do gg = 1, ncell
    write(*,*)'gg = ', gg
    write(*,*) real(BB1(:,:,gg))    
  enddo


  BB=dcmplx(0.0d0,0.0d0) 
  call FT_k_to_r(UVI,BB)
  do gg = 1, ncell
    U_SAORI(:,:,gg) = real(BB(:,:,gg))
  enddo

!  do gg = 1, ncell
!    k = 0
!    do i = 1, nao
!      do j = 1, nao
!        if(abs(U_SAORI(i,j,gg)).ge.1.0d-6) k = k + 1
!      enddo
!    enddo
!    if(k.eq.0) U_SAORI(:,:,gg) = 0.0d0
!  enddo 

  do i = 1, nao
    do j = 1, nao
      U_SAOR(j,i,1) = U_SAORI(i,j,1)
    enddo
  enddo 
  do gg = 2, ncell
    gg1 = gg
    if(mod(gg,2).gt.0) then
      gg1 = gg - 1
    else
      gg1 = gg + 1
    endif
    do i = 1, nao
      do j = 1, nao
        U_SAOR(i,j,gg) = U_SAORI(j,i,gg1) 
      enddo
    enddo
  enddo 
  
  ABr = 0.0d0
  ABr1 = 0.0d0

  write(*,*)'Check the Wannier orbitals orthogonality'
  do g = 1, ncell
    do h = 1, ncell
      gg = diffcell(h,g)
      if(gg.gt.0) then
        call DGEMM('N','N',nao,nao,nao,alpha,Real(Sr_0g(:,:,h)),nao,U_SAOR(:,:,gg),nao,alpha,ABr(:,:,g),nao)
      endif
    enddo
  enddo

  do g = 1, ncell
    do h = 1, ncell
      gg = diffcell(h,g)
      if(gg.gt.0) then
        call DGEMM('N','N',nao,nao,nao,alpha,U_SAORI(:,:,h),nao,ABr(:,:,gg),nao,alpha,ABr1(:,:,g),nao)
      endif
    enddo
  enddo
  !stop
  !write(*,*) 'Print ABr1'

  do g = 1, ncell
    write(*,*) 'g = ', g
    write(*,*) real(ABr1(:,:,g))
  enddo
  !stop

!  write(*,*)'Print Wannier natural orbitals (NAO)'
!  do i = 1, nao
!    write(*,*)'i = ', i
!    do g = 1, ncell
!      do j = 1, nao
!        write(*,*) U_SAORI(i,j,g)
!      enddo
!    enddo
!  enddo

  do i = 1, nao
    do g = 1, ncell
      do j = 1, nao
        wannier((g-1)*nao+j,i) = U_SAORI(i,j,g)
      enddo
    enddo
  enddo

  BB=dcmplx(0.0d0,0.0d0)
  call FT_k_to_r(IUV,BB)
  do gg = 1, ncell
    U_SAOR(:,:,gg) = real(BB(:,:,gg))
  enddo

  do i = 1, nao
    do j = 1, nao
      U_SAORI(j,i,1) = U_SAOR(i,j,1)
    enddo
  enddo
  do gg = 2, ncell
    gg1 = gg
    if(mod(gg,2).gt.0) then
      gg1 = gg - 1
    else
      gg1 = gg + 1
    endif
    do i = 1, nao
      do j = 1, nao
        U_SAORI(i,j,gg) = U_SAOR(j,i,gg1)
      enddo
    enddo
  enddo

  ABr = 0.0d0
  ABr1 = 0.0d0

  do g = 1, ncell
    do h = 1, ncell
      gg = diffcell(h,g)
      if(gg.gt.0) then
        call DGEMM('N','N',nao,nao,nao,alpha,P_0g(:,:,h),nao,U_SAORI(:,:,gg),nao,alpha,ABr(:,:,g),nao)
      endif
    enddo
  enddo

  do g = 1, ncell
    do h = 1, ncell
      gg = diffcell(h,g)
      if(gg.gt.0) then
        call DGEMM('N','N',nao,nao,nao,alpha,U_SAOR(:,:,h),nao,ABr(:,:,gg),nao,alpha,ABr1(:,:,g),nao)
      endif
    enddo
  enddo

  write(*,*) 'Check P(r) via FT against that via the real-space transformation'

  nel = 0.0d0
  do g = 1, ncell
    write(*,*) 'g = ', g
    write(*,*) real(ABr1(:,:,g)) - real(BB1(:,:,g))
    do i = 1, nao
      nel = nel + real(ABr1(i,i,g))
    enddo
  enddo 
 
  write(*,*) "nel = ", nel

  deallocate(U_SAO,U_SAOI,BB,UV,UVI,Pk_ao,Pk_sao,IUV,IUVI,P_tmp)
  deallocate(ABr,ABr1,U_SAOR,U_SAORI,Sr_0g,work,iwork,rwork,eigval)

end subroutine pbc_nao
