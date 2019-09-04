subroutine transform_ints(wannier,vint,norb,imp)
  use gf2, only: nao, ncell
  use seet
  implicit none
  integer imp
  integer norb
  double precision,dimension(nao*ncell,nao*ncell,nao*ncell,nao*ncell)::vintf
  double precision,dimension(norb,norb,norb,norb)::vint
  double precision,dimension(nao*ncell,nao)::wannier
  double precision,dimension(:,:,:,:),allocatable::vita,vitb,vitc

  integer i,j,k,l,l1,k1,j1,i1
  integer ii,jj,kk,ll
  integer ic,jc,kc,lc,nf,inc,mi

  double precision, allocatable, dimension(:) :: vi_c1i
  integer, allocatable, dimension(:,:) :: i_vi_c1i
  integer, allocatable, dimension(:) :: n_vi_c
  integer inmax_c
  CHARACTER(LEN=16) :: filename
  INTEGER :: file_id
  CHARACTER(LEN=4) :: dsetname_ind1
  CHARACTER(LEN=4) :: dsetname_ind2
  CHARACTER(LEN=4) :: dsetname_ind3



  filename = "integrals.h5"
  dsetname_ind1 = "ind1"

  nf=nao*ncell
  allocate(n_vi_c(nf))


  call  size_ints_init(nf,inmax_c,filename,file_id,n_vi_c)


  allocate(vi_c1i(inmax_c))
  allocate(i_vi_c1i(4,inmax_c))

  allocate(vita(nf,nf,nf,norb))  
  vita=0.0d0 


  do ic=1,ncell
     do i=1,nao
 
        inc=(ic-1)*nao+i
        call read_next(filename,file_id,inc,inmax_c,i_vi_c1i,vi_c1i)

        do mi = 1, n_vi_c(inc) 
           j1 = i_vi_c1i(2,mi)
           jc=(j1-1)/nao+1
           j=j1-(jc-1)*nao
 
           k1 = i_vi_c1i(3,mi)
           kc=(k1-1)/nao+1
           k=k1-(kc-1)*nao

           l1 = i_vi_c1i(4,mi)
           lc=(l1-1)/nao+1
           l=l1-(lc-1)*nao        


           do ll=1,norb
              l1=imp_orb_numbers(imp,ll)
           
              vita((ic-1)*nao+i,(jc-1)*nao+j,(kc-1)*nao+k,ll)=vita((ic-1)*nao+i,(jc-1)*nao+j,(kc-1)*nao+k,ll)+vi_c1i(mi)*wannier((lc-1)*nao+l,l1)

           end do

        end do
     end do

  end do

 allocate(vitb(nf,nf,norb,norb))  
 vitb=0.0d0

 do kk=1,norb
    k1=imp_orb_numbers(imp,kk)

    do ll=1,norb

       do ic=1,ncell
          do i=1,nao

             do jc=1,ncell
                do j=1,nao
                
                   do kc=1,ncell
                      do k=1,nao
                      
                         vitb((ic-1)*nao+i,(jc-1)*nao+j,kk,ll)=vitb((ic-1)*nao+i,(jc-1)*nao+j,kk,ll)+vita((ic-1)*nao+i,(jc-1)*nao+j,(kc-1)*nao+k,ll)*wannier((kc-1)*nao+k,k1) 
 
                      end do
                   end do

                end do
             end do

          end do
       end do

    end do

 end do

 deallocate(vita)

 allocate(vitc(nf,norb,norb,norb))
 vitc=0.0d0

 do jj=1,norb
    j1=imp_orb_numbers(imp,jj)

    do ll=1,norb

       do kk=1,norb 

          do ic=1,ncell
             do i=1,nao

                do jc=1,ncell
                   do j=1,nao
                
                      vitc((ic-1)*nao+i,jj,kk,ll)=vitc((ic-1)*nao+i,jj,kk,ll)+vitb((ic-1)*nao+i,(jc-1)*nao+j,kk,ll)*wannier((jc-1)*nao+j,j1) 
                            
                   end do
                end do

             end do
          end do
       
       end do

    end do

 end do

 deallocate(vitb)
 vint=0.0d0

 do ii=1,norb
    i1=imp_orb_numbers(imp,ii)

    do ll=1,norb

       do kk=1,norb 

          do jj=1,norb  

             do ic=1,ncell
                do i=1,nao

                   vint(ii,jj,kk,ll)=vint(ii,jj,kk,ll)+vitc((ic-1)*nao+i,jj,kk,ll)*wannier((ic-1)*nao+i,i1) 
                            
                end do
             end do

          end do
       
       end do

    end do

 end do

 deallocate(vitc)






  deallocate(vi_c1i)
  deallocate(i_vi_c1i)
  deallocate(n_vi_c)


end subroutine transform_ints
  




subroutine read_next(filename,file_id,nc,inmax_c,i_vi_c1i,vi_c1i)
  integer, intent(in) :: nc,inmax_c
  character*6 :: dsetnum
  character*9 :: dsetnumi, dsetnumj
  character*9 :: dsetnum1, dsetnum2, dsetnum3
  double precision, dimension(inmax_c) :: vi_c1i
  integer, dimension(4,inmax_c) :: i_vi_c1i
  INTEGER :: error
  CHARACTER(LEN=16) :: filename
  INTEGER :: file_id
  INTEGER, DIMENSION(1) :: data_dims=(/1/)
  INTEGER, DIMENSION(1) :: data_dim1
  INTEGER, DIMENSION(2) :: data_dim2
  INTEGER, DIMENSION(3) :: data_dim3
  integer::rank
      
  ! Read next Coulomb part
  write(dsetnum,"(I6)") nc
  dsetnum = adjustl(dsetnum)
  dsetnumi = 'fi_'//trim(dsetnum)
  dsetnum1 = 'f1_'//trim(dsetnum)
  rank = 2
  data_dim2(1) = 4
  data_dim2(2) = inmax_c
  call h5ltread_dataset_int_f(file_id,dsetnum1,i_vi_c1i(:,:),data_dim2,error)
  rank = 1
  data_dim1(1) = inmax_c
  call h5ltread_dataset_double_f(file_id,dsetnumi,vi_c1i(:),data_dim1,error)
end subroutine read_next

subroutine size_ints_init(dim_size,inmax_c,filename,file_id,n_vi_c)
  use gf2, only: nao, ncell, hf_temp, ncellgf2  
  integer, intent(IN) :: dim_size
  integer, intent(OUT) :: inmax_c
  integer, dimension(dim_size)::n_vi_c
  INTEGER, DIMENSION(1) :: data_dims=(/1/)
  INTEGER, DIMENSION(1) :: data_dim1
  INTEGER, DIMENSION(2) :: data_dim2
  INTEGER, DIMENSION(3) :: data_dim3

  CHARACTER(LEN=16) :: filename
  CHARACTER(LEN=4) :: dsetname_ind1
  CHARACTER(LEN=4) :: dsetname_ind2
  CHARACTER(LEN=4) :: dsetname_ind3
  INTEGER :: file_id





  CALL h5open_f(error)
  CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
  data_dim1(1) =dim_size
  call h5ltread_dataset_int_f(file_id,dsetname_ind1,n_vi_c,data_dim1,error)  

  ! get the maximum possible index for coulomb integral
  inmax_c = maxval(n_vi_c)
  write(*,*)'inmax_c = ', inmax_c

  ! perform memory allocation

  CALL h5fclose_f(file_id,error)
  CALL h5close_f(error)
end subroutine size_ints_init
