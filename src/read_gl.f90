subroutine read_gl
  use gf2, only: nk,nao,nleg,ncell
  implicit none
  integer i,j,k,l,ii,jj,kk,ll,ic
  complex*16,dimension(:,:,:,:), allocatable::sigl
  complex*16,dimension(:,:,:), allocatable::siglr
  double precision:: tmpr,tmpi

  allocate(sigl(nao,nao,2*nk,0:nleg))
  allocate(siglr(nao,nao,0:nleg))

     open(950,file='GF.bin',action='read',form='unformatted')
     do k=1,2*nk 
        do l=0,nleg
           do i=1,nao
              do j=1,nao
                read(950)kk,ll,ii,jj,tmpr,tmpi 
                sigl(ii,jj,kk,ll)=dcmplx(tmpr,tmpi)
              end do
           end do
        end do
     end do
     close(950)
     write(*,*)"I AM HERE"
     CALL FLUSH(6)
     siglr=dcmplx(0.0d0,0.0d0)
     ic=1
     do l=0,nleg-2  
     write(*,*)"I AM HERE2",l
     CALL FLUSH(6)
        call FT_k_to_r_ic(sigl(1:nao,1:nao,1:2*nk,l), siglr(1:nao,1:nao,l),ic)
     end do

     !even polynomials

     do l=2,nleg-2,2
        write(11,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l)),dble(siglr(8,8,l)),dble(siglr(9,10,l))
        write(14,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l)),dble(siglr(9,12,l))  
     end do


     !odd polynomials
     do l=1,nleg-2,2
        write(111,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l))
        write(141,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l))
     end do

     siglr=dcmplx(0.0d0,0.0d0)
     ic=9
     do l=0,nleg-2
     write(*,*)"I AM HERE3",l
     CALL FLUSH(6)  
        call FT_k_to_r_ic(sigl(:,:,:,l), siglr(:,:,l),ic)
     end do

do l=2,nleg-2,2
       write(91,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l))
        write(94,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l))
     end do

     !odd polynomials
     do l=1,nleg-2,2
        write(911,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l))
        write(941,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l))
     end do



     siglr=dcmplx(0.0d0,0.0d0)
     ic=20
     do l=0,nleg-2
     write(*,*)"I AM HERE3",l
     CALL FLUSH(6)  
        call FT_k_to_r_ic(sigl(:,:,:,l), siglr(:,:,l),ic)
     end do

do l=2,nleg-2,2
       write(21,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l))
        write(24,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l))
     end do

     !odd polynomials
     do l=1,nleg-2,2
        write(211,*)l, dble(siglr(1,1,l)),dble(siglr(1,3,l)),dble(siglr(1,8,l))
        write(241,*)l, dble(siglr(4,1,l)),dble(siglr(4,3,l)),dble(siglr(4,8,l))
     end do

     deallocate(sigl)
     deallocate(siglr)
   end subroutine read_gl
