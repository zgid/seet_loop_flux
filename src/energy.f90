subroutine calc_energy_freq_depend_old(mu,Fk,Sk,X_k,sig_sol_loc_ndc,sigma_inf_loc,E_corr,nts,uniform,iterc)
  use gf2, only:nao,nk,ncell,iwmax,Sigma_tau,beta,omegas,nc,logm,nleg,four_coef,n_threads,ncellgf2
  implicit none
  integer::nts,uniform,iterc
  double precision,save::pi=4.D0*DATAN(1.D0)
  double precision:: mu
  double precision:: E_corr
  complex*16, dimension(nao,nao,2*nk) :: Sk,Fk,X_k
  complex*16, dimension(nao,nao)::sigma_inf_loc
  complex*16, dimension(nao,nao,iwmax)::sig_sol_loc_ndc
  complex*16, dimension(:,:),allocatable:: Xr
  complex*16, dimension(:,:,:),allocatable:: Sr_inv
  complex*16, dimension(:,:,:),allocatable :: Gk,auxx,rsig,isig,sigr,sigma_w,sigma_tau_k
  complex*16, dimension(:,:,:,:),allocatable :: sigl

  double precision, dimension(:),allocatable::aux,erg
  double precision, dimension(:,:,:), allocatable :: aux1
  double precision, dimension(:,:),allocatable::sigma1
  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer:: i,j,k,l,ic,w,kk,t,ip,ii
  double precision:: energy,e1,tmpr
  complex*16::om
  integer:: iwmax_max
  integer::tn,nth!,n_threads
  character*256 file_name
  allocate(sigma_w(nao,nao,2*nk))
  allocate(rsig(nao,nao,0:nleg),isig(nao,nao,0:nleg),sigl(1:nao,1:nao,0:nleg,1:2*nk))
  nth=n_threads
  call omp_set_num_threads(nth) 

  iwmax_max=100000  
  write(*,*)'enter energy calculation'
  call flush(6)
  allocate(Sr_inv(nao,nao,ncell))
  allocate(sigr(nao,nao,ncell))
  allocate(auxx(nao,nao,2*nk))
  allocate(Xr(nao,nao))
  allocate(Gk(nao,nao,2*nk))
  allocate(erg(nth))
!  allocate(aux(iwmax_max))
!  allocate(aux1(nao,nao,iwmax_max))
  allocate(sigma1(nao,nao))
  allocate(sigma_tau_k(nao,nao,nts))

  do k=1,2*nk
     auxx(:,:,k)=Sk(:,:,k)
     call invert_complex_matrix(nao,auxx(:,:,k))
  end do

  Sr_inv=dcmplx(1.0d0,0.0d0)
  call FT_k_to_r(auxx,Sr_inv)

  do k=1,2*nk

     do t = 1, nts
       call FT_r_to_single_k(sigma_tau_k(:,:,t),Sigma_tau(:,:,:,t),k)
     enddo
     call calc_sigl_from_sigtau(nts,uniform,rsig,dble(sigma_tau_k(:,:,:)))
     call calc_sigl_from_sigtau(nts,uniform,isig,dimag(sigma_tau_k(:,:,:)))


     sigl(:,:,:,k)=rsig(:,:,:)+dcmplx(0.0d0,1.0d0)*isig(:,:,:)

           if (k<=2) then
              do i=0,nleg,2
                 write(*,*)i,k,"slk d",sigl(1,1,i,k)
              end do

              do i=0,nleg,2
                 write(*,*)i,k,"slk od",sigl(1,2,i,k)
              end do
           end if

        end do



  ! open(950,file="sigr.txt",status="replace")
  ! open(951,file="sigl.txt",status="replace")
  ! do l = 0,nleg
  !   do k = 1, 2*nk
  !     do i = 1, nao
  !       do j = 1, nao
  !         write(unit=950,fmt='(4I4XX,F15.10X)') l, k, i, j, dble(sigl(i,j,l,k))
  !         write(unit=951,fmt='(4I4XX,F15.10X)') l, k, i, j, dimag(sigl(i,j,l,k))
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(950)
  ! close(951)


  auxx(1:nao,1:nao,1:2*nk)=X_k(1:nao,1:nao,1:2*nk)
  do k=1,2*nk
     call invert_complex_matrix(nao,auxx(:,:,k))
  end do


  energy=0.0d0
!  aux(:)=0.0d0
  erg(:)=0.0d0

!  open(950,file='sigma_11.txt',status='replace')
!  If(mod(iterc,2).eq.0) then 
!    write(file_name,fmt='(A7,I0,A4)') 'sig11e_',iterc,'.txt'
!    open(950,file=file_name,status='replace')
!  else
!    write(file_name,fmt='(A7,I0,A4)') 'sig11o_',iterc,'.txt'
!    open(950,file=file_name,status='replace')
!  Endif

!1$OMP PARALLEL DO DEFAULT(none) &                                                                                 
!1$OMP SHARED(Fk,Sk,auxx,mu,iwmax,erg,nts,uniform,sigl,sig_sol_loc_ndc,sigma_inf_loc) &        
!1$OMP PRIVATE(w,tn)
  
  do w = 1, iwmax
        tn=omp_get_thread_num()  
        call calc_energy_inner_new1(w,mu,Sk,Fk,auxx,nts,uniform,erg(tn+1),sigl,sig_sol_loc_ndc(:,:,w),sigma_inf_loc)
  end do

  energy=0.0d0  
  do i=1,nth
     energy=energy+erg(i)
  end do

  e1=energy/beta


  w=iwmax
  
  write(*,*)"e1",e1, iwmax
  call flush(6)
  
 ! do k=1,2*nk
 !    call calc_sigw_using_tnl_inner(iwmax,sigma_w(:,:,k),sigl(:,:,:,k))
 ! end do

 !  sigr=dcmplx(0.0d0,0.0d0)
  !  do ic=1,ncellgf2
  !     call FT_k_to_r_ic(sigma_w(:,:,:), sigr(:,:,ic),ic) 
  !     if (ic==1) then
  !        sigr(:,:,ic)=sigr(:,:,ic)+sig_sol_loc_ndc(:,:,iwmax)
  !     end if
  !     aux(:)=0.0d0
  !     aux1(:,:,:)=0.0d0
  !     sigma1(:,:)=-1*dimag(sigr(:,:,ic))*dimag(omegas(iwmax))
  !     do w = iwmax+1,iwmax_max
  !        om=dcmplx(0.0d0,(2.0d0*(w-1)+1)*pi/beta)
  !        do i=1,nao
  !           do j=1,nao
  !              do k = 1, nao
  !              aux1(i,j,w)=aux1(i,j,w)-dimag(Sr_inv(i,k,ic)/om)*dimag(sigma1(j,k)/om)
  !              enddo
  !           end do
  !        end do
  !        do i = 1, nao
  !        energy=energy+2.0d0*aux1(i,i,w)
  !     enddo
  !  end do
  !  end do

 !  energy=energy/beta

  energy=0.0d0
  call calc_energy_inner_hf_tail(w,mu,Sk,Fk,auxx,Sr_inv,energy,sigl,sig_sol_loc_ndc(:,:,w),sigma_inf_loc,iwmax_max)

  write(*,*)"HF tail energy", energy/beta
  write(*,*)"energy",energy/beta+e1

  E_corr=energy/beta+e1

  call flush(6)
  
!  close(950) 
  deallocate(sigma_w)
!write(*,*)"1"
call flush(6) 
deallocate(rsig,isig,sigl)
!write(*,*)"10"
call flush(6) 

 deallocate(Sr_inv)
!write(*,*)"2"
call flush(6) 
  deallocate(sigr)
!write(*,*)"3"
call flush(6) 
  deallocate(auxx)
!write(*,*)"4"
call flush(6) 
  deallocate(Xr)
!write(*,*)"5"
call flush(6) 
  deallocate(Gk)
!write(*,*)"6"
call flush(6) 
!if (allocated(erg)) then
!   write(*,*) "erg allocated" 
!   call flush(6)
!else
!      write(*,*) "erg NOT allocated" 
!   call flush(6)
!end if
  deallocate(erg)
!write(*,*)"11"
call flush(6) 
!  deallocate(aux)
!write(*,*)"7"
call flush(6) 
!  deallocate(aux1)
!write(*,*)"8"
call flush(6) 
  deallocate(sigma1)
!  write(*,*)"9"
call flush(6)
deallocate(sigma_tau_k)
!write(*,*)"12"
call flush(6) 

end subroutine calc_energy_freq_depend_old


subroutine calc_energy_inner_new1(w,mu,Sk,Fk,Xk,nts,uniform,e1,sigl,sig_sol_loc_ndc,sigma_inf_loc)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,omegas,ncellgf2
  implicit none
  integer:: w,nts,uniform
  double precision:: mu
  complex*16, dimension(nao,nao,2*nk) :: Sk,Fk,Xk
  complex*16, dimension(nao,nao) :: sig_sol_loc_ndc,sigma_inf_loc 
  complex*16, dimension(:,:),allocatable:: Xr,sigrw
  complex*16, dimension(:,:,:),allocatable :: Gk,auxx,sigma_w
  complex*16, dimension(1:nao,1:nao,0:nleg,1:2*nk):: sigl
  double precision, dimension(1)::e1
  complex*16, dimension(:,:),allocatable::sigma1,aux,aux1,aux2
  integer:: i,j,k,ic,kk,t
  double precision:: energy
  complex*16::om,muomega
  integer:: iwmax_max
  complex*16::alpha,betaz
  alpha=dcmplx(1.0d0,0.0d0)  
  betaz=dcmplx(0.0d0,0.0d0) 

  allocate(sigrw(nao,nao))
  allocate(Xr(nao,nao))
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(aux(nao,nao),aux1(nao,nao),aux2(nao,nao))

  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))

     aux2=dcmplx(0.0d0,0.0d0)
     aux1(:,:)=sig_sol_loc_ndc(:,:)

     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux1(:,:),nao,betaz,aux2,nao)
     aux1=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux2,nao,Xk(:,:,k),nao,betaz,aux1(:,:),nao)

     sigma_w(:,:,k)=sigma_w(:,:,k)+aux1(:,:)

     aux1(:,:)=sigma_inf_loc(:,:)

     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux1(:,:),nao,betaz,aux2,nao)
     aux1=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux2,nao,Xk(:,:,k),nao,betaz,aux1(:,:),nao)

     muomega =omegas(w)+mu
     Gk(:,:,k) = (muomega*Sk(:,:,k)-Fk(:,:,k)-aux1(:,:)-sigma_w(:,:,k))
     call invert_complex_matrix(nao,Gk(:,:,k))
     
  end do
  energy=0.0d0
  aux(:,:)=dcmplx(0.0d0,0.0d0)

  do ic=1, ncellgf2
     Xr(:,:)=dcmplx(0.0d0,0.0d0)
     sigrw(:,:)=dcmplx(0.0d0,0.0d0)  
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:),ic) 
     call FT_k_to_r_ic(sigma_w(:,:,:), sigrw(:,:),ic) 
!     if(ic .eq. 2) then
!        if(w < 201) write(*,*) w, Xr(1,1), sigrw(1,1)
!     end if
     do i=1,nao
        do j=1,nao
           do k = 1, nao 
             aux(i,j)=aux(i,j)+dble(Xr(i,k))*dble(sigrw(j,k))
             aux(i,j)=aux(i,j)-dimag(Xr(i,k))*dimag(sigrw(j,k))
           enddo
        end do
     end do
  end do
  do i = 1, nao 
     e1(1)=e1(1)+2.0d0*dble(aux(i,i))
  enddo

  deallocate(sigrw)
  deallocate(Xr)
  deallocate(Gk)
  deallocate(sigma_w)
  deallocate(aux,aux1,aux2)
!

end subroutine calc_energy_inner_new1








subroutine calc_energy_inner_hf_tail(w,mu,Sk,Fk,Xk,Sr_inv,energy,sigl,sig_sol_loc_ndc,sigma_inf_loc,iwmax_max)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,omegas,ncellgf2,beta
  implicit none
  integer:: w,iwmax_max
  double precision:: mu
  complex*16, dimension(nao,nao,2*nk) :: Sk,Fk,Xk
  complex*16, dimension(nao,nao,ncell)::Sr_inv
  complex*16, dimension(nao,nao) :: sig_sol_loc_ndc,sigma_inf_loc 
  complex*16, dimension(:,:),allocatable:: Xr,sigrw
  complex*16, dimension(:,:,:),allocatable :: Gk,sigma_w
  complex*16, dimension(1:nao,1:nao,0:nleg,1:2*nk):: sigl
  double precision, dimension(:,:),allocatable::sigma1
  complex*16,dimension(:,:),allocatable::aux,aux1,aux2
  double precision, dimension(:,:,:), allocatable :: aux11
  integer:: i,j,k,ic,kk,t
  double precision:: energy
  complex*16::om,muomega
  double precision, save::pi=4.D0*DATAN(1.D0) 
  complex*16::alpha,betaz
  alpha=dcmplx(1.0d0,0.0d0)  
  betaz=dcmplx(0.0d0,0.0d0) 

  allocate(sigrw(nao,nao))
  allocate(Xr(nao,nao))
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(aux(nao,nao),aux1(nao,nao),aux2(nao,nao),aux11(nao,nao,iwmax_max))

  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))

     aux2=dcmplx(0.0d0,0.0d0)
     aux1(:,:)=sig_sol_loc_ndc(:,:)

     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux1(:,:),nao,betaz,aux2,nao)
     aux1=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux2,nao,Xk(:,:,k),nao,betaz,aux1(:,:),nao)

     sigma_w(:,:,k)=sigma_w(:,:,k)+aux1(:,:)

     aux1(:,:)=sigma_inf_loc(:,:)

     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux1(:,:),nao,betaz,aux2,nao)
     aux1=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux2,nao,Xk(:,:,k),nao,betaz,aux1(:,:),nao)

     muomega =omegas(w)+mu
     Gk(:,:,k) = (muomega*Sk(:,:,k)-Fk(:,:,k)-aux1(:,:)-sigma_w(:,:,k))
     call invert_complex_matrix(nao,Gk(:,:,k))
     
  end do

  energy=0.0d0
  aux11=0.0d0
   aux1=dcmplx(0.0d0,0.0d0)
   do ic=1, ncellgf2
      Xr(:,:)=dcmplx(0.0d0,0.0d0)
      sigrw(:,:)=dcmplx(0.0d0,0.0d0)  
      call FT_k_to_r_ic(sigma_w(:,:,:), sigrw(:,:),ic) 
      aux11(:,:,:)=0.0d0
      aux1(:,:)=-1*dimag(sigrw(:,:))*dimag(omegas(iwmax))  !sigma1
      do w = iwmax+1,iwmax+iwmax_max
         om=dcmplx(0.0d0,(2.0d0*(w-1)+1)*pi/beta)
         do i=1,nao
            do j=1,nao
               do k = 1, nao
                  aux11(i,j,w-iwmax)=aux11(i,j,w-iwmax)-dimag(Sr_inv(i,k,ic)/om)*dimag(aux1(j,k)/om)
!                  write(*,*)w,Sr_inv(i,k,ic),om,aux1(j,k),aux11(i,j,w-iwmax)
               enddo
            end do
         end do
         do i = 1, nao
            energy=energy+2.0d0*aux11(i,i,w-iwmax)
!            write(*,*)"ENG",energy,aux11(i,i,w-iwmax)
         enddo
      end do
   enddo

  deallocate(sigrw)
  deallocate(Xr)
  deallocate(Gk)
  deallocate(sigma_w)
  deallocate(aux,aux1,aux2,aux11)
!

end subroutine calc_energy_inner_hf_tail
