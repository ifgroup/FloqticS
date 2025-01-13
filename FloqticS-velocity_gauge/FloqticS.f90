      PROGRAM elpatestmpidist                                           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Part of the FloqticS code from
!! Floquet theory and computational method for the optical absorption of laser-dressed solids
!!
!! Code author: Vishal Tiwari
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use elpa

      implicit none
      
      include "mpif.h"

      integer                          :: nblk
      integer                          :: np_rows, np_cols,myrow,mycol
      integer                          :: na_rows, na_cols
      integer                          :: ddarray
      integer                          :: infile, icontxt
      integer                          :: myid, nprocs,my_prow,my_pcol
      integer                          :: mpi_comm_rows, mpi_comm_cols
      integer                          :: mpierr, my_blacs_ctxt
      integer                          :: sc_desc(9), info,nprow,npcol
      integer                          :: jl,jpc,il,ipr
      integer*8                        :: locsize
      integer(kind=MPI_ADDRESS_KIND)   :: lb,locextent
      integer(kind=MPI_OFFSET_KIND)    :: disp
      integer, external                :: numroc,INDXG2L,INDXG2P,INDXL2G
      integer                   :: dims(2),distribs(2),dargs(2),pdim(2)
      integer                          :: STATUS
      integer,dimension(MPI_STATUS_SIZE) :: mpistatus
      integer                          :: success, incontxt
      character(len=8)                 :: task_suffix
      integer, parameter               :: error_unit = 0
      class(elpa_t), pointer           :: elpap

      integer*4                        :: na, nev, hdim, flm, zero
      integer*4                        :: i,j,k,l,n,m,o,p,ka
      integer*4                        :: io,di,kpt
      integer*4, allocatable           :: bzindex(:),occupation(:)
      complex*16                       :: iota
      complex*16, allocatable          :: drivedipole(:,:)
      complex*16, allocatable          :: probedipole(:,:)
      complex*16, allocatable          :: ham(:,:), z(:,:)
      complex*16, allocatable          :: a2(:,:),eigvec(:,:)
      complex*16, allocatable          :: coff2(:,:),fouriermme(:,:,:)
      complex*16, allocatable          :: testc(:)
      real*8                           :: awm,omega,pi,totvol
      real*8                           :: EO, OMGA
      real*8, allocatable              :: syseng(:), ev(:), popfact(:,:)
      real*8, allocatable              :: kval(:)
      real*8, allocatable              :: testr(:)
      character(len=20)                :: filekval,filebandeng
      character(len=20)                :: filedrivemme,fileprobemme
      character(len=20)                :: fileoccup





      call mpi_init(mpierr)
      call mpi_comm_rank(mpi_comm_world,myid,mpierr)
      call mpi_comm_size(mpi_comm_world,nprocs,mpierr)


        pi=3.1415926535d0
        iota=(0d0,1d0)
 
        filekval     =  'kvector.txt'
        filebandeng  =  'bandeng.txt'
        filedrivemme =  'dipoledrive.txt'
        fileprobemme =  'dipoleprobe.txt'
        fileoccup    =  'initialoccup.txt'
      
        if (myid .eq. 0) then
            open(10,file='input.txt')
                read(10,*)hdim
                read(10,*)flm
                read(10,*)EO
                read(10,*)OMGA
                read(10,*)kpt
                read(10,*)totvol
             close(10)
        endif

        nblk=6

        call mpi_bcast(hdim,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
        call mpi_bcast(flm,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
        call mpi_bcast(EO,1,MPI_REAL8,0,mpi_comm_world,mpierr)
        call mpi_bcast(OMGA,1,MPI_REAL8,0,mpi_comm_world,mpierr)
        call mpi_bcast(nblk,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
        call mpi_bcast(kpt,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
        call mpi_bcast(totvol,1,MPI_REAL8,0,mpi_comm_world,mpierr)

        call mpi_barrier(mpi_comm_world, mpierr)


        if (myid .eq. 0) then

!!!!!!!!!!!!!!! quick input file checks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(testr(hdim))
        allocate(testc(hdim))

        open(1,file=filekval)
            do i = 1, kpt
                read(1, *, iostat=io)(testr(j),j=1,3)
            enddo
            if (io /= 0.0) then
               write(*,*) "error reading k-vectors file"
               call exit(0)
            endif
        close(1)

        open(2,file=filebandeng)
            do i = 1, kpt
               read(2, *, iostat=io)(testr(j),j=1,hdim)
            enddo
            if (io /= 0) then
                write(*,*) "error reading bandeng file"
                call exit(0)
            endif
        close(2)

        open(3,file=filedrivemme)
            do i = 1, kpt*hdim
               read(3, *, iostat=io)(testc(j),j=1,hdim)
            enddo
            if (io /= 0) then
               write(*,*) "error reading dipoledrive file"
               call exit(0)
            end if
        close(3)

        open(4,file=fileprobemme)
            do i = 1, kpt*hdim
               read(4, *, iostat=io)(testc(j),j=1,hdim)
            enddo
            if (io /= 0) then
               write(*,*) "error reading dipoleprobe file"
               call exit(0)
            end if
        close(4)

        open(5,file=fileoccup)
            do i = 1, kpt
               read(5, *, iostat=io)(testr(j),j=1,hdim)
            enddo
            if (io /= 0) then
               write(*,*) "error reading occupation file"
               call exit(0)
            end if
        close(5)
        

!!!!!!!!!!!!!!!!!!!!  openeing input and output files !!!!!!!!!!!!!!!!!

        open(1,file=filekval)
        open(2,file=filebandeng)
        open(3,file=filedrivemme)
        open(4,file=fileprobemme)
        open(5,file=fileoccup)

        open(10,file='band.txt')
        open(20,file='quasienergy.txt')

        endif



!!!!!!!!!!!!!!!!!! reading the input files one k-vector at a time !!!!!!!!!

!!!!!!!!!!!!!!!!!!!!! k loop starts !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ka=1,kpt


        allocate(kval(3))
        allocate(syseng(hdim))
        allocate(drivedipole(hdim,hdim))
        allocate(probedipole(hdim,hdim))
        allocate(occupation(hdim))

        if (myid == 0) then
               
          read(1,*)(kval(j),j=1,3)

          read(2,*)(syseng(j),j=1,hdim)

          do i=1,hdim
          read(3,*)(drivedipole(i,j),j=1,hdim)
          enddo

          do i=1,hdim
          read(4,*)(probedipole(i,j),j=1,hdim)
          enddo
        
          read(5,*)(occupation(j),j=1,hdim)

        write(10,'(*(f15.5,1X))')kval(:)/0.529177d0, &
              syseng(:)*27.211386245988d0

        endif


!!!!!!!!!!!!!!!!!!!!! sending them to all nodes !!!!!!!!!!!!!!!!!!!!!!!!!!!

        call mpi_barrier(mpi_comm_world, mpierr)

        call mpi_bcast(syseng,hdim,MPI_REAL8,0,mpi_comm_world,mpierr)
        call mpi_bcast(drivedipole,hdim*hdim,    &
                  MPI_COMPLEX16,0,mpi_comm_world,mpierr)
        call mpi_bcast(probedipole,hdim*hdim,    &
                  MPI_COMPLEX16,0,mpi_comm_world,mpierr)

        call mpi_barrier(mpi_comm_world, mpierr)

        na=hdim*(2*flm+1)




!!!!!!!!!!!!!!! set up Floquet-Hamiltonian for ELPA !!!!!!!!!!!!!!!!!!!

        pdim=0
        call MPI_Dims_create(nprocs, 2, pdim, mpierr)


        np_rows=pdim(1)
        np_cols=pdim(2)

        if (nblk > (na/np_rows)) nblk = na/np_rows
        if (nblk > (na/np_cols)) nblk = na/np_cols

        call mpi_barrier(mpi_comm_world, mpierr)
 


        call blacs_get(-1, 0, icontxt)
        call blacs_gridinit(icontxt, 'R', np_rows, np_cols)
        call blacs_gridinfo(icontxt, nprow, npcol, my_prow, my_pcol)


        na_rows = numroc(na, nblk, my_prow, 0, np_rows)
        na_cols = numroc(na, nblk, my_pcol, 0, np_cols)



         call mpi_barrier(mpi_comm_world, mpierr)

         call descinit(sc_desc, na, na, nblk, nblk, 0, 0,  &
                   icontxt, na_rows, info)

      if (info .ne. 0) then
        write(error_unit,*) 'Error in BLACS descinit! info=',info
        write(error_unit,*) 'more MPI tasks than are possible for your'
        write(error_unit,*) 'problem size (matrix size and blocksize)!'
        write(error_unit,*) 'Try reducing the number of MPI tasks...'
        call MPI_ABORT(mpi_comm_world, 1, mpierr)
      endif



      allocate(ham(na_rows,na_cols))
      ham=(0d0,0d0)


!!!!!! diganoal part of Floquet-Bloch Hamiltonian !!!!!!!!

      do k=1,2*flm+1
         do l=1,hdim
            j = (k-1)*hdim+l
            il = indxg2l(j,nblk,0,0,np_rows)
            jl = indxg2l(j,nblk,0,0,np_cols)
            ipr = indxg2p(j,nblk,0,0,np_rows)
            jpc = indxg2p(j,nblk,0,0,np_cols)
            if (ipr==my_prow .and. jpc==my_pcol) then
               HAM(il,jl)=syseng(l)+dble(k-flm-1)*OMGA
            endif
         enddo
      enddo



!!!!!!!!!!!!! off-diagonal part of Floquet-Bloch Hamiltonian !!!!!!!!!!!!!!!!!!

      do k=1,2*flm
         do l=1,hdim
            do m=1,hdim

               i = (k-1)*hdim+l
               j = k*hdim+m

               il = indxg2l(i,nblk,0,0,np_rows)
               jl = indxg2l(j,nblk,0,0,np_cols)

               ipr = indxg2p(i,nblk,0,0,np_rows)
               jpc = indxg2p(j,nblk,0,0,np_cols)

            if (ipr==my_prow .and. jpc==my_pcol) then
               HAM(il,jl)=-drivedipole(l,m)*EO/2d0/OMGA/iota
            endif


            enddo
         enddo
      enddo

      do k=1,2*flm
         do l=1,hdim
            do m=1,hdim

               i = k*hdim+l
               j = (k-1)*hdim+m

               il = indxg2l(i,nblk,0,0,np_rows)
               jl = indxg2l(j,nblk,0,0,np_cols)

               ipr = indxg2p(i,nblk,0,0,np_rows)
               jpc = indxg2p(j,nblk,0,0,np_cols)

            if (ipr==my_prow .and. jpc==my_pcol) then
              HAM(il,jl)=conjg(-drivedipole(m,l)*EO/2d0/OMGA/iota)
            endif

            enddo
         enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!! initializing ELPA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpi_barrier(mpi_comm_world,mpierr)
 

        nev=na
        allocate(ev(na))
        allocate(z(na_rows,na_cols))
        z=(0d0,0d0)
        

!-------------------------------------------------------------------------------

      if (elpa_init(20170403) /= elpa_ok) then
        write(*,*)
        print *, "ELPA API version not supported"
        write(*,*)
        stop
      endif
      elpap => elpa_allocate()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! set parameters decribing the matrix and it's MPI distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call elpap%set("na", na, success)
      call elpap%set("nev", nev, success)
      call elpap%set("local_nrows", na_rows, success)
      call elpap%set("local_ncols", na_cols, success)
      call elpap%set("nblk", nblk, success)
      call elpap%set("mpi_comm_parent", mpi_comm_world, success)
      call elpap%set("process_row", my_prow, success)
      call elpap%set("process_col", my_pcol, success)

      success = elpap%setup()

      call elpap%set("solver", elpa_solver_1stage, success)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate eigenvalues/eigenvectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           
      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only

      call elpap%eigenvectors(ham, ev, z, success)
        
      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only

      deallocate(ham)

!!!!!!!!!!!!!!! Floquet-Bloch Hamiltonian diagonalization done !!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!! collecting eigenvalues and eigenvectors in fundamental FBZ  !!!!!!!!!!!!!!!!!

      allocate(eigvec(na,hdim))
      allocate(a2(na,hdim))
      allocate(bzindex(hdim))

       

      if (myid .eq. 0) then


        j=0
         do i=1,na
            if ( abs(ev(i)) .le.  OMGA/2d0 ) then
                j=j+1
                bzindex(j)=i
            endif
         enddo

         if ( j .ne. hdim) then
            write(*,*)'error no. of quasienergy is not equal to no. of &
             bands. Increase the time-periodic functions'
            call MPI_ABORT(mpi_comm_world, 1, mpierr)
         endif
 
         write(20,'(*(f9.5,1X))')kval(:),ev(bzindex(:))*27.211d0


        endif

      call mpi_bcast(bzindex,hdim,MPI_INTEGER,0,mpi_comm_world,mpierr)


        call mpi_barrier(mpi_comm_world, mpierr)

       a2=(0d0,0d0)
       eigvec=(0d0,0d0)
!
      do k=1,hdim
        do jl=1,na_cols
         j=INDXL2G(jl,nblk,my_pcol,0,np_cols)
         
            if (j .eq. bzindex(k) ) then
       
               do il=1,na_rows
    
                  i=INDXL2G(il,nblk,my_prow,0,np_rows)
                  a2(i,k)=z(il,jl)
    
               enddo
            endif

        enddo
      enddo

        call mpi_reduce(a2,EIGVEC,na*hdim,MPI_COMPLEX16,MPI_SUM,0,   & 
           mpi_comm_world,mpierr)



        call mpi_barrier(mpi_comm_world,mpierr)

       deallocate(a2)
       deallocate(z)




!!!!!!!!!!!!!!!!! computing population factor and Fourier component of MME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (myid==0) then



          allocate(coff2(hdim,hdim))
          coff2=(0d0,0d0)


          do i=1,hdim
             do j=1,hdim
                do k=1,2*flm+1
                      coff2(j,i)=coff2(j,i)+eigvec(j+(k-1)*hdim,i)
                enddo
                 coff2(j,i)=coff2(j,i)*conjg(coff2(j,i))
             enddo
          enddo




          allocate(popfact(hdim,hdim))
          popfact=0d0
          do i=1,hdim
             do j=1,hdim
                do k=1,hdim
                   do l=1,hdim
                      popfact(i,j)=popfact(i,j)+dble(coff2(k,j)        &
                *(coff2(l,i)))*dble(occupation(l)*(1-occupation(k)))
                   enddo
                enddo
             enddo
          enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! only compute fourier components 0 - 2*flm+1 due to property of  MME
!!!!!!!!!!!!!! \mathcal{P}_{\alpha\beta \mathbf{k}}^{(n)}=\mathcal{P}_{\beta\alpha\mathbf{k}}^{(-n)*}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          allocate(fouriermme(hdim,hdim,4*flm+1))
          fouriermme=(0d0,0d0)
          do i=1,hdim
             do j=1,hdim
                do k=0,2*flm
                   do l=k,2*flm
                      do m=1,hdim
                         do n=1,hdim
                fouriermme(i,j,k+2*flm+1)=fouriermme(i,j,k+2*flm+1)+   &
                 probedipole(m,n)*conjg(eigvec((l-k)*hdim+m,i))*      &
                 eigvec(l*hdim+n,j)
                         enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo





          fouriermme=conjg(fouriermme)*fouriermme


!!!!!!!!!!!!!!!!!!!!!!!!! inter-FBZ transitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! convert into realunits 54440750. = e^2*hbar*(mom au-si)^2/4/pi/m^2/(permvacu)/c/(length au-si)**3/(enrgy au-si)**2

          do i=1,hdim
            do j=1,hdim
               do k=1,2*flm
                  omega=ev(bzindex(i))-ev(bzindex(j))+k*OMGA
                  awm=fouriermme(i,j,k+2*flm+1)      &
                    *(popfact(i,j)-popfact(j,i))
                  awm=awm/omega*54440750.148254760d0/totvol
         if (abs(awm) > 1e-6 ) then
             write(*,'(3(f9.5,1X),2(f15.5,1X),3(I3,1X),5(f15.5,1X))')          &
          kval(:)/0.529177d0,omega*27.211d0,awm,                &  
              i,j,k,                                               &  
             ev(bzindex(i))*27.211d0, ev(bzindex(j))*27.211d0,     &
             abs(fouriermme(i,j,k+2*flm+1))*(1.2438408262263d0)**2, &
                      popfact(i,j),popfact(j,i)
                 endif
               enddo
            enddo
         enddo


!!!!!!!!!!!!!!!!!!!!!!!!! intra-FBZ transitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          do i=1,hdim
             do j=1,hdim
                if (i > j) then
                  omega=ev(bzindex(i))-ev(bzindex(j))
                  awm=real(fouriermme(i,j,2*flm+1))      &
                    *(popfact(i,j)-popfact(j,i))
                  awm=awm/omega*54440750.148254760d0/totvol    
              if (abs(awm) > 1e-6 ) then
             write(*,'(3(f9.5,1X),2(f15.5,1X),3(I3,1X),5(f15.5,1X))')  &
          kval(:)/0.529177d0,omega*27.211d0,awm,         & 
              i,j,0,                                        & 
             ev(bzindex(i))*27.211d0, ev(bzindex(j))*27.211d0, &
                 abs(fouriermme(i,j,2*flm+1))*(1.2438408262263d0)**2, &
                      popfact(i,j),popfact(j,i)
                  endif
                endif
            enddo
         enddo

          deallocate(coff2)
          deallocate(popfact)
          deallocate(fouriermme)



       endif
 
       call mpi_barrier(mpi_comm_world,mpierr)
        


!!!!!!!!!!!!!!!!!!!!!! deallocations and closing files !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       deallocate(eigvec)
       deallocate(bzindex)
       deallocate(ev)


       call elpa_deallocate(elpap)
       call elpa_uninit()
       
       deallocate(kval)
       deallocate(syseng)
       deallocate(drivedipole)
       deallocate(probedipole)
       deallocate(occupation)

      call mpi_barrier(mpi_comm_world,mpierr)


        


!!!!!!!!!!!!!!!!!!!!! k loop ends !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo

        close(10)
        close(20)
    
        if (myid == 0) then
           close(1)
           close(2)
           close(3)
           close(4)
           close(5)
       endif


      call mpi_finalize(mpierr)


      end   

        
