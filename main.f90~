program gf2pbc
  use seet, only : cells_seet
  use gf2, only : nk, nao, ncell, nel_cell, ncellgf2, seet_gf2_iter
  use gf2, only : uniform_tau, power_tau, max_points_mesh, power_mesh_small_int 
  use gf2, only : nuniform
  use gf2, only : nleg
  use gf2, only : iwmax

  use gf2, only : beta, mu
 
  use gf2, only : indices_cell

  use gf2, only : k_weights

  use gf2, only : power_mesh, power_mesh_small 

  use gf2, only : S_k, F_k
  use gf2, only : four_coef, tnl
  use gf2, only : Sigma, Gr_tau, Sigma_tau, Sigma_tau_tmp
  
  use gf2, only : E_thr
  use gf2, only : itermax
  use gf2, only : damp, dampd
  use gf2, only : hf_temp
  use gf2, only : debug_p
  use gf2, only : rst,rst_old_sig

  use gf2, only : sig11
 

  implicit none
  integer iterc, nts, iters
  integer i, j, ic, w, ww,it,k
  integer it1, ic1, i1, j1
  logical exist, mott, mott1, haveprev,symmetric,print_sigl,nodm
  complex*16, allocatable, dimension(:,:,:) :: S_kk, F_kk,unityk
  complex*16, allocatable, dimension(:,:,:) :: dm_r, dm_r_prev, dm_r_prev_e, dm_r_prev_o,sig_sol_loc_ndc
  complex*16, allocatable, dimension(:,:) :: sigma_inf_loc 
  double precision, allocatable, dimension(:,:,:,:) :: Sigma_tau_prev, Sigma_tau_prev_e, Sigma_tau_prev_o, Gr_tau_prev
  double precision E_corr, E_S0, Energy, Energy_prev, E_diff
  double precision Energy_e, Energy_prev_e, E_diff_e, Energy_o, Energy_prev_o, E_diff_o,thrf

  character(2)::pt
  character(3)::ut
  character(2)::na
  character(10)::bb
  character(10)::iw

  double precision beta0, beta1, betainc,tmpi,tmpr
  integer itt, icc, ii, jj, its


  haveprev = .false.
  cells_seet =1


  call timest('Top of periodic GF2')
! Read input parameters from gf2.inp
  call read_params
  call timest('gf2.inp input file read successfully')
  call flush(6)
!
! Initialize omega and tau grids, Gr_tau
  write(*,*)"before inigrids"
  call flush(6)
  call inigrd(nts)
  call flush(6)
  write(*,*)"after inigrids"
  call flush(6)
  call timest('Omega and tau grids initialized')
  call flush(6)
  allocate(sig_sol_loc_ndc(nao,nao,iwmax))
  allocate(sigma_inf_loc(nao,nao))
  sig_sol_loc_ndc=dcmplx(0.0d0,0.0d0)
  sigma_inf_loc=dcmplx(0.0d0,0.0d0)
! Read Gaussian files: F_k, S_k, k_weights, Fourier_coeff, expand F_k, S_k to
! F_kk, S_kk
!

  write(*,*)"before allocation",power_mesh_small(353),power_mesh_small(1)
  call flush(6)
  
  allocate(S_kk(nao,nao,2*nk))
  allocate(F_kk(nao,nao,2*nk))
  allocate(unityk(nao,nao,2*nk))
  unityk=dcmplx(0.0d0,0.0d0)  
  do k = 1, 2*nk
     do i = 1, nao 
        unityk(i,i,k)=dcmplx(1.0d0,0.0d0)
     end do
  end do
  write(*,*)"after allocation"
  call flush(6)
  write(*,*)"before read_process"
  call flush(6)  
  call read_process(F_kk,S_kk)
  write(*,*)"after read_process"
  call flush(6)
  inquire(file="S_kk.txt", exist=exist)
  write(*,*)'exist = ', exist
  If(exist) then
     open(1,file="S_kk.txt",status="old")
  else
     open(1,file="S_kk.txt",status="new")
  endif
  do k = 1, nk
     do i = 1, nao
        do j = 1, nao
           write(1,*)k, i, j, dble(S_kk(i,j,2*k-1)), dimag(S_kk(i,j,2*k-1))
        enddo
     enddo
  enddo
  close(1)
  call flush(6)

  call timest('Gaussian files read successfully')
  call flush(6)
  thrf=1.d-16
  call screen_real_overlap_matrix(thrf)
  call flush(6)
  call calc_S_eigval(S_kk)
  call check_sigma_cells(S_kk)
  
  !call read_gl      
  !call flush(6)
  
  
  allocate(tnl(iwmax,nleg))
  allocate(sig11(iwmax))
  tnl(:,:)=dcmplx(0.0d0,0.0d0)
  call read_tnl
  call timest('Legendre coefficients read successfully')
  call flush(6)

!
! Here the iterations should start
  iterc = 0
  iters = 0
  mu = -7.587542976269955E-002 ! -2.681950099423187E-002 !0.0d0
  allocate(Sigma_tau_prev(nao,nao,ncell,nts))
  allocate(Gr_tau_prev(nao,nao,ncell,nts))
  Sigma_tau_prev(:,:,:,:)=0.0d0

  write(*,*)"check allocated",allocated(Sigma_tau),allocated(Gr_tau_prev)
  call flush(6)
  Gr_tau_prev(:,:,:,:)=0.0d0

  allocate(dm_r(nao,nao,ncell))
  allocate(dm_r_prev(nao,nao,ncell))

  dm_r = dcmplx(0.0d0,0.0d0)
  dm_r_prev =  dcmplx(0.0d0,0.0d0)
!  dm_r_prev_e = dcmplx(0.0d0,0.0d0)
!  dm_r_prev_o = dcmplx(0.0d0,0.0d0)
  call initialize_Gr_tau(nts)
  Sigma_tau(:,:,:,:) = 0.0d0

  write(pt,'(I2)')power_tau
  write(ut,'(I3)')uniform_tau
  write(na,'(I2)')nao
  write(iw,'(I5)')iwmax
  write(bb,'(F10.1)')beta
!  call comp_sig_fock(F_kk,nts)
!  call flush(6)
!  stop
! Make PBC specific arrangements: read indices, create array of differences
  call dffcll
  call addcll
  call timest('Done with PBC arrangements')
  call flush(6)
  if (hf_temp==0) then
!     call read_in_v_integrals_short_n_1
!     write(*,*)"int read"
!     call flush(6)
!     call integrals_ce_n_emb
!     write(*,*)"int read2"
!     call flush(6)
  end if
  call timest('Done with integrals') 
  call flush(6)
  mu=-5.411822589464305E-002 !-4.185678761392139E-002 
  Energy = 0.0d0
  Energy_prev = 0.0d0 
  E_diff = 0.0d0
  Energy_e = 0.0d0
  Energy_prev_e = 0.0d0
  E_diff_e = 0.0d0
  Energy_o = 0.0d0
  Energy_prev_o = 0.0d0
  E_diff_o = 0.0d0

  E_S0 = 0.0d0
  call update_Fk(F_kk)

  Sigma_tau = 0.0d0
  call inigrd(nts)
  If(rst.gt.0) then
     open(1,file="SE_nonsym.txt",status="old")
     do it = 1, nts
        do ic = 1, ncellgf2
           do i = 1, nao
              do j = 1, nao
                read(1,*) it1, ic1, i1, j1, Sigma_tau(i1,j1,ic1,it1)
                ! write(unit=1,fmt='(4I4XX,F15.10X)') it, ic, i, j,
                ! Sigma_tau(i,j,ic,it)
              enddo
           enddo
        enddo
     enddo
     close(1)
  Endif
  if (rst_old_sig==1) then
     do w=1,iwmax  
        do i=1,nao 
           do j=1,nao  
              read(6000,*)ww,ii,jj,tmpr,tmpi
              sig_sol_loc_ndc(i,j,w)=dcmplx(tmpr,tmpi)
           end do
        end do
     end do

     do i=1,nao 
        do j=1,nao  
           read(7000,*)ii,jj,tmpr,tmpi
           sigma_inf_loc(i,j)=dcmplx(tmpr,tmpi)
        end do
     end do

  end if
  print_sigl=.false.
  nodm=.false.
  call mu_dens(F_kk,S_kk,unityk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,print_sigl,nodm) 
  write(*,*)'initial chemical potential search: mu = ', mu 
  call flush(6)

        write(*,*)"------------old_dm"
        write(*,*) dm_r(:,:,1)
        write(*,*)"------------"
        call flush(6)
   
!  Energy_prev = Energy

!  call calc_ft_gwk_to_gtk2(S_kk,F_kk,mu,nts,nuniform,dm_r) 
  call calc_ft_gwk_to_gtk2(S_kk,F_kk,sig_sol_loc_ndc,sigma_inf_loc,mu,nts,nuniform,dm_r)
 
!!BGN check SEET on HF
!  call seet_main_loop(F_kk,S_kk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,E_S0,nts)
!  write(*,*)"after seet_main_loop"
!  call flush(6)
!  stop
!!END
  

 call timest('Initial GF built') 
  call flush(6)

!!  write(*,*)'Read integrals for GF2'
!!  call read_in_v_integrals_short_n_1
!!  write(*,*)"int read"
!!  call flush(6)
!!  call integrals_ce_n_emb
!!  write(*,*)"int read2"
!!  call flush(6)
!!  write(*,*)'Integrals have been successfully read'
  write(*,*) "seet_gf2_iter", seet_gf2_iter

  do its=1,seet_gf2_iter  

100  write(*,*)' +++++++++++++ top of GF2-SEET iteration ', iterc , '++++++++++++++++'
     call flush(6) 
  

       if(rst /= 1) then
          if (hf_temp==0) then
 
           call sigma_outer(Gr_tau,nts)
           call timest('Self-energy calculated')
           call flush(6)

        !!---symmetryzation---!
        !  do t=1,nts
        !     call check_symmetry(Sigma_tau(:,:,:,t),symmetric)
        !
        !     call restore_symmetry(Sigma_tau(:,:,:,t),symmetric)
        !  end do
        
!        call write_self_energy
      
        write(*,*)'damping factor at iter ', iterc, ' is ', damp
        call flush(6)
        If(iterc.gt.0) then 
           Sigma_tau(:,:,:,:) = (1.0d0-damp)*Sigma_tau(:,:,:,:) + damp*Sigma_tau_prev(:,:,:,:)
        Endif
      
        Sigma_tau_prev(:,:,:,:) = Sigma_tau(:,:,:,:)
      
        !mu=-2.00936576889430
        print_sigl=.false.
        nodm=.true.
        dm_r(:,:,:)=dcmplx(0.0d0,0.0d0)   
        call mu_dens(F_kk,S_kk,unityk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,print_sigl,nodm)
        print_sigl=.false.
        nodm=.false.
        call mu_dens(F_kk,S_kk,unityk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,print_sigl,nodm)
     end if

     If(iterc.gt.0) dm_r(:,:,:) = dampd*dm_r_prev(:,:,:) + (1.0d0-dampd)*dm_r(:,:,:)

        write(*,*)"------------old_dm"
        write(*,*) dm_r(:,:,1)
        write(*,*)"------------"
        call flush(6)
       
     dm_r_prev(:,:,:) = dm_r(:,:,:)  
   
     call fock_gauss(dm_r,E_S0,iterc)
     call update_Fk(F_kk)
   !    mu=-7.062712295716458E-002
     print_sigl=.true.
     nodm=.true.
     !mu=-2.30012130321793 
     dm_r(:,:,:)=dcmplx(0.0d0,0.0d0)
     call mu_dens(F_kk,S_kk,unityk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,print_sigl,nodm)
     dm_r(:,:,:)=dcmplx(0.0d0,0.0d0)  
     call calc_ft_gwk_to_gtk2(S_kk,F_kk,sig_sol_loc_ndc,sigma_inf_loc,mu,nts,nuniform,dm_r)
   
  end if
   
     if (debug_p==1) then
        open(950,file="SE_tau.txt",status="replace")
        open(951,file="GF_tau.txt",status="replace")
        do it = 1, nts
           do ic = 1, ncell
              do i = 1, nao
                 do j = 1, nao
                    write(unit=950,fmt='(4I4XX,F15.10X)') it, ic, i, j, Sigma_tau(i,j,ic,it)
                    write(unit=951,fmt='(4I4XX,F15.10X)') it, ic, i, j, Gr_tau(i,j,ic,it)
                enddo
             enddo
          enddo
       enddo
       close(950)
       close(951)
    end if
    
    E_corr = 0.0d0
    if (hf_temp==0) then
       call calc_energy_freq_depend_old(mu,F_kk,S_kk,unityk,sig_sol_loc_ndc,sigma_inf_loc,E_corr,nts,nuniform,iterc)
    end if
    write(*,*)'E_corr = ', E_corr
    call flush(6)
    write(*,*)'E_S0 = ', E_S0
    call flush(6)
    Energy_prev = Energy
    Energy = E_corr + E_S0
    write(*,*)'Total energy: ', Energy
    call flush(6)

    if (rst == 1) then
       write(*,*) "restarting from converged GF2 =>> perform only 1 iteration"
       goto 101
    end if
    
    E_diff = Abs(Energy-Energy_prev)
    write(*,*)'Energy difference = ', E_diff  
    call flush(6)
!    stop
    !    write(*,*)'printing GF to file'
!    call flush(6)
    !    call print_g_se(S_kk,F_kk,mu,nts,iterc)
    
    iterc = iterc + 1
    If(E_diff.le.E_thr) then !.and.E_diff_o.le.E_thr) then
       write(*,*)'GF2 has converged in ', iterc, ' iteations'
       !    call update_Fk(F_kk)
       !    call print_g_se(S_kk,F_kk,mu,nts,iterc) 
       goto 101
    Else
       If(iterc.le.itermax) then
          goto 100
       Else
          write(*,*)'GF2 failed to converge in ', iterc, ' iterations'
          !    call system("rm fort.* grtbin.* ")
!          call print_g_se(S_kk,F_kk,mu,nts,iterc)
          goto 101
       Endif
    Endif

101 call flush(6)    

    write(*,*) "=========="
    write(*,*) dm_r(:,:,1)
    write(*,*) "=========="
    call flush(6)
    
    call seet_main_loop(F_kk,S_kk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,E_S0,nts)
    write(*,*)"after seet_main_loop"
    itermax=1
    call flush(6)
   
 end do


  !call system("rm fort.* grtbin.* ")
    deallocate(sig_sol_loc_ndc)
    deallocate(sigma_inf_loc)
    deallocate(Sigma_tau_prev)
    !  deallocate(Sigma_tau_prev_e)
    !  deallocate(Sigma_tau_prev_o)
    deallocate(Gr_tau_prev)
    deallocate(Gr_tau)
    deallocate(S_kk)
    deallocate(F_kk) 
    deallocate(unityk)
    deallocate(dm_r)
    deallocate(dm_r_prev)
    !  deallocate(dm_r_prev_e)
    !  deallocate(dm_r_prev_o)
    deallocate(Sigma_tau)
    deallocate(sig11)
    !  deallocate(Sigma_tau_tmp)
    call timest('End of GF2')
  end program gf2pbc


!--------------------------------------------------


  subroutine write_self_energy
    use gf2, only : Sigma_tau

        open(1,file="SE.txt",status="replace")
        do it = 1, nts
           do ic = 1, ncellgf2/2
              do i = 1, nao
                 do j = 1, nao
                    if(dabs(Sigma_tau(i,j,ic*2,it)-Sigma_tau(j,i,ic*2+1,it)).gt.1.0d-6) then
                       write(1,*) it, ic*2, i, j, (Sigma_tau(i,j,ic*2,it)-Sigma_tau(j,i,ic*2+1,it)), Sigma_tau(i,j,ic*2,it)
                    endif
                 enddo
              enddo
           enddo
        enddo
        close(1)
        open(1,file="SE_1.txt",status="replace")
        do it = 1, nts
         do ic = 1, 1
            do i = 1, nao
               do j = 1, nao
                  if(dabs(Sigma_tau(i,j,ic,it)-Sigma_tau(j,i,ic,it)).gt.1.0d-6) then
                     write(1,*) it, ic, i, j, (Sigma_tau(i,j,ic,it)-Sigma_tau(j,i,ic,it)), Sigma_tau(i,j,ic,it)
                  endif
               enddo
            enddo
         enddo
      enddo
      close(1)
      
      !     open(1,file="SE_nonsym.txt")
      !open(1,file="SE_nonsym.txt",status="replace")
      !do it = 1, nts
      !   do ic = 1, ncellgf2
      !      do i = 1, nao
      !         do j = 1, nao
      !            !                read(1,*) itt, icc, ii, jj, Sigma_tau(i,j,ic,it)
      !            write(unit=1,fmt='(4I4XX,F15.10X)') it, ic, i, j, Sigma_tau(i,j,ic,it)
      !         enddo
      !      enddo
      !   enddo
      !enddo
      !close(1)

     !If(symse.gt.0) then
     !  symmetrize the central cell, then other cells
      write(*,*)'symmetrizing the self-energy'
      call flush(6)
      do it = 1, nts
         do i = 1, nao
            do j = 1, i
               Sigma_tau(i,j,1,it) = 5.0d-01*(Sigma_tau(i,j,1,it)+Sigma_tau(j,i,1,it))
               Sigma_tau(j,i,1,it) = Sigma_tau(i,j,1,it)
            enddo
         enddo
      enddo
      do it = 1, nts
         do ic = 1, ncellgf2/2
            do i = 1, nao
               do j = 1, i
                  Sigma_tau(i,j,ic*2,it) = 5.0d-01*(Sigma_tau(i,j,ic*2,it)+Sigma_tau(j,i,ic*2+1,it))
                  Sigma_tau(j,i,ic*2+1,it) = Sigma_tau(i,j,ic*2,it)
                  
                  Sigma_tau(j,i,ic*2,it) = 5.0d-01*(Sigma_tau(j,i,ic*2,it)+Sigma_tau(i,j,ic*2+1,it))
                  Sigma_tau(i,j,ic*2+1,it) = Sigma_tau(j,i,ic*2,it)
               enddo
            enddo
         enddo
      enddo
    end subroutine write_self_energy
