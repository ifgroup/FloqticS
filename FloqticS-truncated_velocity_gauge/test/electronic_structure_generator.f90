  PROGRAM   commutatorgen
    implicit none

    integer,parameter :: dp=kind(0.0d0)

    real(kind=dp):: pi,tpi,hte,bta,unit_cell_len,fermi,omega,absorption,total_volume
    real(kind=dp)::aa,bb,cc,dd,ee,unit_cell_volume
    real(kind=dp)::lattvec(3,3),reclattvec(3,3)
    real(kind=dp),allocatable :: hopping_matrix(:,:,:),dipole_wannier(:,:), dipole_wannier2(:,:,:,:)
    real(kind=dp),allocatable :: bandeng(:),kpoint(:,:),directionindex(:)

    complex(kind=dp) :: iota
    complex(kind=dp),allocatable :: ham(:,:),commutator_matrix(:,:,:,:),vn(:,:,:)
    complex(kind=dp),allocatable :: derham(:,:,:),derxi(:,:,:),v3(:,:)
    complex(kind=dp),allocatable :: v1(:,:),v2(:,:)
    complex(kind=dp),allocatable :: com1(:,:,:),com2(:,:,:),delh(:,:,:)
    complex(kind=dp),allocatable :: eigvec(:,:),expval(:)
    
    integer :: i,j,k,l,m,n,o,p,io,neib,ka,direction
    integer :: hdim,total_term,nkpt
    integer,allocatable ::  rindex(:,:),trackrind(:,:)

    iota = (0.0_dp,1.0_dp)
    hte  = 27.211386245_dp
    bta  = 0.5291772109_dp
    pi   = 3.1415926535_dp
    tpi  = 2.0_dp*3.1415926535_dp
    
 

! ---------------------------------------------------------------------
! enter a, lattice vectors, and reciprocal lattice vectors from QE here
! ---------------------------------------------------------------------

    unit_cell_volume = 4.7168_dp
    lattvec(:,1) = [ 1.0_dp ,  0.000000_dp ,  0.000000_dp ]
    lattvec(:,2) = [ 0.0_dp ,  4.006410_dp ,  0.000000_dp ]
    lattvec(:,3) = [ 0.0_dp ,  0.000000_dp ,  4.006410_dp ]
    lattvec = lattvec*unit_cell_volume

    reclattvec(:,1) = [ 1.0_dp, 0.000000_dp, 0.000000_dp ]
    reclattvec(:,2) = [ 0.0_dp, 0.249600_dp, 0.000000_dp ]
    reclattvec(:,3) = [ 0.0_dp, 0.000000_dp, 0.249600_dp ]
    reclattvec = reclattvec*tpi/unit_cell_volume

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------



    open(60,file='input.txt') 
      read(60,*)hdim
      read(60,*)direction
      read(60,*)i
      read(60,*)aa
      read(60,*)aa
      read(60,*)i
      read(60,*)fermi
      read(60,*)nkpt
      read(60,*)total_term
      read(60,*)unit_cell_volume
    close(60)

    if (direction /= 1) then
        if (direction /= 2 ) then
           if (direction /= 3 ) then
              write(*,*)'direction should be 1,2,3'
              call exit(0)
           endif
        endif
    endif


    open(10,file='parameter.dat')
        n=0
        do
          read(10,*,iostat=io)
          if(io == 0) then
             n=n+1
          else
             exit
          endif
        enddo
        rewind(10)
    
        neib=int(n/hdim/hdim)
        allocate(rindex(neib,3))
        allocate(hopping_matrix(neib,hdim,hdim))
        allocate(dipole_wannier2(3,neib,hdim,hdim))
         

        total_volume=unit_cell_volume*nkpt
        do i=1,neib
           do j=1,hdim
              do k=1,hdim
                 read(10,*)rindex(i,1:3),ee,ee,hopping_matrix(i,k,j),dipole_wannier2(1:3,i,k,j)
              enddo
           enddo
        enddo
    close(10)


   allocate(dipole_wannier(hdim,hdim))
   dipole_wannier = 0.0_dp
   dipole_wannier(:,:) = dipole_wannier2(direction,(neib+1)/2,:,:)


   deallocate(dipole_wannier2)

    hopping_matrix=hopping_matrix/hte
    dipole_wannier=dipole_wannier/bta

   allocate(directionindex(neib))
   do i=1,neib
      directionindex(i) = dot_product(lattvec(direction,:),rindex(i,:))
   enddo

    open(20,file='kpoint.txt')
       allocate(kpoint(nkpt,3))
       do i=1,nkpt
          read(20,*)kpoint(i,1:3)
       enddo
    close(20)
   

    open(40,file='bandstructure.txt')
    open(50,file='fieldfreeabsspectra.txt')

    allocate(commutator_matrix(total_term,nkpt,hdim,hdim))
    commutator_matrix=(0.0_dp,0.0_dp)

    do ka=1,nkpt

       allocate(expval(neib))
       do i=1,neib
          expval(i)=cdexp(iota*dot_product(rindex(i,:),kpoint(ka,:))*tpi)
       enddo

       allocate(ham(hdim,hdim))
       ham=(0.0_dp,0.0_dp)
       do i=1,neib
          ham(:,:)=ham(:,:) + expval(i)*hopping_matrix(i,:,:)
       enddo

       deallocate(expval)

       allocate(eigvec(hdim,hdim))
       allocate(bandeng(hdim))

          eigvec=(0.0_dp,0.0_dp)
          bandeng=0.0_dp
          eigvec=ham
          call  CALL_ZHEEV(eigvec, hdim,bandeng)

       write(40,'(*(f15.8,1X))')kpoint(ka,1:3),bandeng(1:hdim)

       deallocate(ham)

       call multiplymatrixleft(dcmplx(hopping_matrix),dipole_wannier*(-iota),com1)
       call multiplymatrixright(dcmplx(hopping_matrix),(-iota)*dipole_wannier,com2)

       call derivmatrix(dcmplx(hopping_matrix),directionindex,delh)
          
       allocate(expval(neib))
       do i=1,neib
          expval(i)=cdexp(iota*dot_product(rindex(i,:),kpoint(ka,:))*tpi)
       enddo

       allocate(vn(neib,hdim,hdim))
       vn=delh+com1-com2

       do i=1,neib 
         commutator_matrix(1,ka,:,:)=commutator_matrix(1,ka,:,:) + expval(i)*vn(i,:,:)
       enddo

        deallocate(expval)
        deallocate(com1)
        deallocate(com2)
        deallocate(delh)
     
       do j=2,total_term
          
          call multiplymatrixleft(vn,dipole_wannier*(-iota), com1)
          call multiplymatrixright(vn, (-iota)*dipole_wannier, com2)
          
          call derivmatrix(vn,directionindex, delh)

          allocate(expval(neib))
          do i=1,neib
             expval(i)=cdexp(iota*dot_product(rindex(i,:),kpoint(ka,:))*tpi)
          enddo

          deallocate(vn)
          allocate(vn(neib,hdim,hdim))

          vn=delh+com1-com2
          vn=vn/dble(j)
          do i=1,neib 
             commutator_matrix(j,ka,:,:)=commutator_matrix(j,ka,:,:) + expval(i)*vn(i,:,:)
          enddo

             
           deallocate(expval)
           deallocate(com1)
           deallocate(com2)
           deallocate(delh)

       enddo
           deallocate(vn)
          
        do i=1,total_term
           commutator_matrix(i,ka,:,:)=matmul(transpose(conjg(eigvec)),matmul(commutator_matrix(i,ka,:,:),eigvec))
        enddo


        


         do i=1,hdim
              do j=1,hdim 
                 if (bandeng(i) <= fermi .and. bandeng(j) > fermi) then
                     omega=bandeng(j)-bandeng(i)
               absorption=abs(commutator_matrix(1,ka,i,j)*commutator_matrix(1,ka,j,i))
                    absorption=absorption/omega*54440750.148254760_dp/total_volume 
                      if (abs(absorption) > 1e-2_dp) then
                   write(50,'(5(f20.7,1X))')kpoint(ka,1:3), omega*hte,absorption
                      endif
                endif
              enddo
           enddo

       deallocate(bandeng)
       deallocate(eigvec)
         
  enddo



    open(30,file='commutator_matrix.txt')
   do i=1,nkpt
      do j=1,total_term
         do k=1,hdim
            do l=1,hdim
         write(30,'(2f50.15)') commutator_matrix(j,i,l,k)
            enddo
         enddo
      enddo
   enddo

   close(30)


    deallocate(commutator_matrix)
    deallocate(rindex)
    deallocate(hopping_matrix,dipole_wannier)
    deallocate(kpoint)
    deallocate(directionindex)
 
  contains 

      subroutine  multiplymatrixright(mat1,mat2,resmat)
      implicit none

      integer,parameter:: dp=kind(1.0d0)
      integer :: dim1,dim2,dim3,i,j
      complex(kind=dp),dimension(:,:,:),intent(in) :: mat1
      complex(kind=dp),dimension(:,:),intent(in) :: mat2
      complex(kind=dp),dimension(:,:,:),intent(out),allocatable :: resmat(:,:,:)
  
        dim1 = size(mat1, 1)
        dim2 = size(mat1, 2)
        dim3 = size(mat1, 3)

        if (dim2 /= size(mat2, 1) .or. dim3 /= size(mat2, 2)) then
            write(*,*) "mat1,mat2 dimension dont match"
            call exit(0)
        endif

        allocate(resmat(dim1,dim2,dim3))
        do i=1,dim1
           resmat(i,:,:) = matmul(mat1(i,:,:),mat2(:,:))
        enddo
      end subroutine
  

      subroutine  multiplymatrixleft(mat1,mat2,resmat)
      implicit none

      integer,parameter:: dp=kind(1.0d0)
      integer :: dim1,dim2,dim3,i,j
      complex(kind=dp),dimension(:,:,:),intent(in) :: mat1
      complex(kind=dp),dimension(:,:),intent(in) :: mat2
      complex(kind=dp),dimension(:,:,:),intent(out),allocatable :: resmat(:,:,:)
  
        dim1 = size(mat1, 1)
        dim2 = size(mat1, 2)
        dim3 = size(mat1, 3)

        if (dim2 /= size(mat2, 1) .or. dim3 /= size(mat2, 2)) then
            write(*,*) "mat1,mat2 dimension dont match"
            call exit(0)
        endif

        allocate(resmat(dim1,dim2,dim3))
        do i=1,dim1
           resmat(i,:,:) = matmul(mat2(:,:),mat1(i,:,:))
        enddo
      end subroutine


      subroutine  derivmatrix(mat1,dindex,resmat)
        implicit none

        integer,parameter:: dp=kind(1.0d0)
        integer :: dim1,dim2,dim3,i,j,k,neigh
        complex(kind=dp),dimension(:,:,:),intent(in) :: mat1
        real(kind=dp),dimension(:),intent(in) :: dindex
        complex(kind=dp),dimension(:,:,:),intent(out),allocatable:: resmat
        complex(kind=dp) :: iota
    
        iota=(0.0_dp,1.0_dp)

        dim1 = size(mat1, 1)
        dim2 = size(mat1, 2)
        dim3 = size(mat1, 3)

        neigh=size(dindex,1)

        if (dim1 /= neigh) then
            write(*,*) "mat1,dindex dimension dont match"
            call exit(0)
        else
           allocate(resmat(dim1,dim2,dim3))
           
           do i=1,neigh
              resmat(i,:,:)=iota*dindex(i)*mat1(i,:,:)
           enddo
        endif

      end subroutine


      SUBROUTINE CALL_ZHEEV(H, NDIM,EIGVALS)
        IMPLICIT NONE
     
        CHARACTER (LEN = 1), PARAMETER  :: JOBZ="V"
        CHARACTER (LEN = 1), PARAMETER  :: UPLO="U"
        INTEGER          :: NDIM,INFO,LDA,LWORK
        COMPLEX(kind=dp)         :: H(NDIM,NDIM)
        REAL(kind=dp)            :: EIGVALS(NDIM)
        COMPLEX(kind=dp), ALLOCATABLE :: WORK(:)
        REAL(kind=dp),    ALLOCATABLE :: RWORK(:)
     
             LDA = NDIM
             LWORK=2*NDIM+10
             ALLOCATE(WORK(LWORK))
             ALLOCATE(RWORK(3*NDIM-2))
        CALL ZHEEV (JOBZ,UPLO,NDIM,H,LDA,EIGVALS,WORK, LWORK,RWORK,INFO)
             DEALLOCATE(WORK)
             DEALLOCATE(RWORK)
     
            return
  
          END subroutine

  end program

