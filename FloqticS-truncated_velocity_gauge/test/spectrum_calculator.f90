   PROGRAM elpatestmpidist
   use elpa
   implicit none
   include "mpif.h"
 
   integer, parameter   :: dp = kind(0.0d0)
   integer              :: nblk,flmpow
   integer              :: na, nev, hdim, flm, zero
 
   integer              :: np_rows, np_cols,myrow,mycol
   integer              :: na_rows, na_cols
   integer              :: ddarray
   integer              :: infile, icontxt
 
   integer              :: myid, nprocs,my_prow,my_pcol
   integer              :: mpi_comm_rows, mpi_comm_cols
   integer              :: mpierr, my_blacs_ctxt,me
   integer              :: sc_desc(9), info,nprow,npcol
 
   integer              :: jl,jpc,il,ipr
   integer              :: totsize,locn,locdisp
 
   integer*8           :: locsize
   integer(kind=MPI_ADDRESS_KIND) :: lb,locextent
 
   integer(kind=MPI_OFFSET_KIND)  :: disp
   integer, external              :: numroc,INDXG2L,INDXG2P,INDXL2G
 
   real(kind=dp), allocatable     :: syseng(:), ev(:),popfact(:,:)
   real(kind=dp), allocatable     :: val(:),meane(:)
   complex(kind=dp), allocatable  :: dipole(:,:,:),temp(:,:)
   complex(kind=dp), allocatable  :: ham(:,:), z(:,:)
   complex(kind=dp), allocatable  :: a2(:,:),eigvec(:,:)
   complex(kind=dp), allocatable  :: coff2(:,:),muflq(:,:,:)
   complex(kind=dp), allocatable  :: locmuflq(:)
   integer, allocatable           :: bzindex(:)
   integer                        :: dims(2),distribs(2),dargs(2),pdim(2)
   real(kind=dp)                  ::   sums,awm,omega
   real(kind=dp)                  ::   totvol,nck
 
   integer                        :: STATUS
   integer,dimension(MPI_STATUS_SIZE) :: mpistatus
   integer                        :: success, incontxt
   character(len=8)               :: task_suffix
   character(len=20)              :: filekpt
   integer                        :: i,j,k,l,n,m,o,p,ka
   integer                          :: ii,jj,kk,ll,mm,nn
 
   real(kind=dp)            :: EO, OMGA,tstart,ctime,fermi
   real(kind=dp)            :: aa,bb,unit_cell_volume
 
 
   class(elpa_t), pointer   :: elpap
 
   INTEGER                  :: io,di,kpt,totpow,power,term
   complex(kind=dp)         :: tr,iota
   real(kind=dp)            :: br,pi,hte,bta
 
   INTEGER,allocatable :: neindex(:,:),reccount(:),discount(:)
   REAL(kind=8),allocatable :: intercoup(:,:,:),dmat(:,:,:),kval(:,:)
 
   call mpi_init(mpierr)
   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)
 
   iota=(0.0_dp,1.0_dp)
   hte=27.211386245_dp
   bta =0.5291772109_dp
   pi  =3.1415926535_dp
 
   tstart=mpi_wtime()
 
     if (myid == 0) then
        open(10,file='input.txt')
            read(10,*)hdim
            read(10,*)i
            read(10,*)flm
            read(10,*)EO
            read(10,*)OMGA
            read(10,*)nblk
            read(10,*)fermi
            read(10,*)kpt
            read(10,*)totpow
            read(10,*)unit_cell_volume
         close(10)
    endif


    call mpi_bcast(hdim,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
    call mpi_bcast(flm,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
    call mpi_bcast(EO,1,MPI_REAL8,0,mpi_comm_world,mpierr)
    call mpi_bcast(OMGA,1,MPI_REAL8,0,mpi_comm_world,mpierr)
    call mpi_bcast(nblk,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
    call mpi_bcast(fermi,1,MPI_REAL8,0,mpi_comm_world,mpierr)
    call mpi_bcast(kpt,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
    call mpi_bcast(totpow,1,MPI_INTEGER,0,mpi_comm_world,mpierr)
    call mpi_bcast(unit_cell_volume,1,MPI_REAL8,0,mpi_comm_world,mpierr)

    totvol=unit_cell_volume*kpt

    allocate(kval(kpt,3))

    if (myid == 0) then

    open(12,file="kpoint.txt")
    do i=1,kpt
      read(12,*)kval(i,1),kval(i,2),kval(i,3)
    enddo
    close(12)

    open(30,file="commutator_matrix.txt",action='read')
    open(40,file="bandstructure.txt",action='read')
    open(20,file="bandeng.txt")

    endif

    call mpi_bcast(kval,kpt*3,MPI_REAL8,0,mpi_comm_world,mpierr)

!        ############loopforkval########################################################
  do ka=1,kpt

    allocate(syseng(hdim))
    allocate(dipole(totpow,hdim,hdim))

    if (myid == 0) then
               
        do j=1,totpow
          do k=1,hdim
             do l=1,hdim
                read(30,*)aa,bb
                dipole(j,l,k)=aa+iota*bb
             enddo
           enddo
        enddo
    
        read(40,*)aa,aa,aa,(syseng(j),j=1,hdim)
    
        write(20,'(*(f9.5,1X))')kval(ka,:),syseng(:)*hte

        

    endif


    call mpi_bcast(syseng,hdim,MPI_REAL8,0,mpi_comm_world,mpierr)
    call mpi_bcast(dipole,totpow*hdim*hdim,MPI_COMPLEX16,0,mpi_comm_world,mpierr)

    na=hdim*(2*flm+1)

    pdim=0
    call MPI_Dims_create(nprocs, 2, pdim, mpierr)

    np_rows=pdim(1)
    np_cols=pdim(2)

    if (nblk > (na/np_rows)) nblk = na/np_rows
    if (nblk > (na/np_cols)) nblk = na/np_cols


    call blacs_get(-1, 0, icontxt)
    call blacs_gridinit(icontxt, 'R', np_rows, np_cols)
    call blacs_gridinfo(icontxt, nprow, npcol, my_prow, my_pcol)

    na_rows = numroc(na, nblk, my_prow, 0, np_rows)
    na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

    call descinit(sc_desc, na, na, nblk, nblk, 0, 0, icontxt, na_rows, info)

    if (info .ne. 0) then
        write(*,*) 'Error in BLACS descinit! info=',info
        write(*,*) 'more MPI tasks than are possible for your'
        write(*,*) 'problem size (matrix size and blocksize)!'
        write(*,*) 'Try reducing the number of MPI tasks...'
        call MPI_ABORT(mpi_comm_world, 1, mpierr)
    endif

    allocate(ham(na_rows,na_cols))
    ham=(0.0_dp,0.0_dp)

      do k=1,2*flm+1
         do l=1,hdim
            j = (k-1)*hdim+l
            jl = indxg2l(j,nblk,0,0,np_cols)
            il = indxg2l(j,nblk,0,0,np_rows)
            ipr = indxg2p(j,nblk,0,0,np_rows)
            jpc = indxg2p(j,nblk,0,0,np_cols)
            if (ipr==my_prow .and. jpc==my_pcol) then
               HAM(il,jl)=syseng(l)+dble(k-flm-1)*OMGA
            endif
         enddo
      enddo



      do power=1,totpow

      do k=0,power
        do l=0,2*flm
            m=l-power+2*k 
            if (m >= 0 .and. m <=2*flm) then
               do n=1,hdim
               do o=1,hdim
                  i = l*hdim+n
                  j = m*hdim+o

                  il = indxg2l(i,nblk,0,0,np_rows)
                  jl = indxg2l(j,nblk,0,0,np_cols)

                  ipr = indxg2p(i,nblk,0,0,np_rows)
                  jpc = indxg2p(j,nblk,0,0,np_cols)

                  if (ipr==my_prow .and. jpc==my_pcol) then
                      HAM(il,jl) = Ham(il,jl)+(iota*eo/omga/2.0_dp)**power &
                              *nck(power,k)*(-1.0_dp)**k*dipole(power,n,o)
                  endif

               enddo
               enddo
            endif
         enddo
      enddo

      enddo


    call mpi_barrier(mpi_comm_world,mpierr)

    nev=na
    allocate(ev(na))
    allocate(z(na_rows,na_cols))
        

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

    ctime=mpi_wtime()

!!!!!!!!!checkeigvals separated by omega!!!!!!!!!!!!!!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    allocate(bzindex(hdim))

    bzindex=1

    if (myid .eq. 0) then

       j=0
       do i=1,na
          if ( abs(ev(i)) <= OMGA/2.0_dp .and. j <= hdim ) then
              j=j+1
              bzindex(j)=i
          endif
       enddo
       if ( j /= hdim) then
          write(*,*)'error fbz not equal to hdim',j,hdim
          call MPI_ABORT(mpi_comm_world, 1, mpierr)
       endif

    endif

    call mpi_bcast(bzindex,hdim,MPI_INTEGER,0,mpi_comm_world,mpierr)

    allocate(a2(na,hdim))
    a2=(0.0_dp,0.0_dp)

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

    deallocate(z)

    allocate(eigvec(na,hdim))
    eigvec=(0.0_dp,0.0_dp)

    call  mpi_allreduce(a2,eigvec,na*hdim,MPI_COMPLEX16,MPI_SUM,mpi_comm_world,mpierr)

    deallocate(a2)

    if (myid==0) then
          
      do j=1,hdim

          sums=dot_product(eigvec(:,j),eigvec(:,j))

          if ( abs(sums-1.0_dp) >= 1e-5) then
             write(*,*)'error in norm',sums,'for col no.',j
             call MPI_ABORT(mpi_comm_world, 1, mpierr)
          endif

      enddo

      allocate(coff2(hdim,hdim))
      coff2=(0.0_dp,0.0_dp)

      do i=1,hdim
         do j=1,hdim
            do k=1,2*flm+1
                  coff2(j,i)=coff2(j,i)+eigvec(j+(k-1)*hdim,i)
            enddo
             coff2(j,i)=coff2(j,i)*conjg(coff2(j,i))
         enddo
      enddo


    allocate(popfact(hdim,hdim))
    popfact=0.0_dp
    do i=1,hdim
       do j=1,hdim
          do k=1,hdim
             do l=1,hdim
      if(syseng(k) <= fermi .and. syseng(l) > fermi) then
                popfact(i,j)=popfact(i,j)+dble(coff2(k,j) * (coff2(l,i)))
               endif
             enddo
          enddo
       enddo
    enddo


   endif




     flmpow=4*flm+1+2*(totpow-1)

        totsize=hdim*hdim*flmpow

        locn = totsize / nprocs
        locdisp = myid*(totsize / nprocs) + min(myid, mod(totsize, nprocs)) + 1

        if (myid < mod(totsize, nprocs)) locn = locn + 1


      allocate(locmuflq(locn))
      locmuflq=(0.0_dp,0.0_dp)

       do m=1,locn
           n=locdisp+m-1
           i=int(n-1)/int(hdim*flmpow)+1
           j=mod(n-1,int(hdim*flmpow))/flmpow+1
           o=mod(mod(n-1,int(hdim*flmpow)),flmpow)+1

           if (o >= 2*flm+1+(totpow-1))  then

              do power=1,totpow
              do kk=0,power-1
                 do ii=1,2*flm+1
                    jj=o-2*flm-totpow+ii-power+1+2*kk
                    if (jj >= 1 .and. jj <= 2*flm+1) then
             locmuflq(m)=locmuflq(m)                                       &
             + dot_product(eigvec((ii-1)*hdim+1:ii*hdim,i),                &
               matmul(dipole(power,:,:),eigvec((jj-1)*hdim+1:jj*hdim,j)))  &
              *power*(iota*eo/2.0_dp/omga)**(power-1)*nck(power-1,kk)    &
                *(-1.0_dp)**kk

                    endif
                 enddo
               enddo
               enddo
            endif

        enddo


      locmuflq=conjg(locmuflq)*locmuflq

      allocate(reccount(nprocs))
      allocate(discount(nprocs))
      do i = 1, nprocs
         reccount(i) = totsize / nprocs
         if (i <= mod(totsize, nprocs)) reccount(i) = reccount(i) + 1
            discount(i) = (i - 1) * (totsize / nprocs) + min(i - 1, mod(totsize, nprocs))
      enddo

      allocate(muflq(hdim,hdim,flmpow))
      muflq=(0.0_dp,0.0_dp)

      call MPI_Gatherv(locmuflq, locn, MPI_COMPLEX16, muflq,   &
            reccount ,discount, MPI_COMPLEX16, 0,  mpi_comm_world, mpierr)
 
         deallocate(reccount)
         deallocate(discount)

     deallocate(locmuflq)

    if (myid == 0) then 
        muflq=reshape(muflq,(/hdim,hdim,flmpow/),order=(/3,2,1/))


          do i=1,hdim
            do j=1,hdim
               do k=1,2*flm+totpow-1
                  omega=ev(bzindex(i))-ev(bzindex(j))+k*OMGA
                  awm=dble(muflq(i,j,k+2*flm+totpow))             &
                         *(popfact(i,j)-popfact(j,i))
                  awm=awm/omega*54440750.148254760d0/totvol
                  if (abs(awm) > 1e-2) then
               write(*,'(5(f20.7,1X),3(I3,1X),5(f15.7,1X))')         &
                     kval(ka,1:3)/bta,                                   &
                     omega*hte,awm,i,j,k,                        &
               ev(bzindex(i))*hte,ev(bzindex(j))*hte,               &
           abs(muflq(i,j,k+2*flm+totpow))*(1.2438408262263d0)**2,  &
                 popfact(i,j),popfact(j,i)
                 endif
               enddo
            enddo
         enddo

          do i=1,hdim
             do j=1,hdim
                if (i > j) then
                  omega=ev(bzindex(i))-ev(bzindex(j))
                  awm=dble(muflq(i,j,2*flm+totpow))         &
                      *(popfact(i,j)-popfact(j,i))
                    awm=awm/omega*54440750.148254760d0/totvol
                  if (abs(awm) > 1e-2) then
              write(*,'(5(f20.7,1X),3(I3,1X),5(f15.7,1X))')            &
             kval(ka,1:3)/bta,omega*hte,awm,                             &
              i,j,0,                                                   &
             ev(bzindex(i))*hte, ev(bzindex(j))*hte,                   &
         abs(muflq(i,j,2*flm+totpow))*(1.2438408262263d0)**2,        &
           popfact(i,j),popfact(j,i)
                  endif
                endif
            enddo
         enddo

       deallocate(popfact)
       deallocate(coff2)


       endif
 
       deallocate(muflq)

       deallocate(eigvec)
       deallocate(bzindex)
       deallocate(ev)
        
       deallocate(syseng)
       deallocate(dipole)

       call mpi_barrier(mpi_comm_world,mpierr)
 
       call elpa_deallocate(elpap)
       call elpa_uninit()

!!!!!!!!!!!!!!!!!!!!!k loopends!!!!!!!!!!!!!!!!
        enddo

      if ( myid == 0) then
        close(20)
        close(30)
        close(40)
       endif

!      call blacs_gridexit(my_blacs_ctxt)
      call mpi_finalize(mpierr)


      end   

        

        real(kind=8) function nck(k,l)

        integer,intent(in)::k,l
        integer(kind=8) :: fact
 
        if (l >=0 .and. l<=k) then
           nck=dble(fact(k)/fact(l)/fact(k-l))
        else
           nck=0.0d0
        endif
        
        return
        end function


        integer*8 function fact(n)

        integer,intent(in)::n
        integer::j
 
        fact=1
        if (n>=0) then
           if (n == 0) then
              fact=1
           else
              do j=1,n
                 fact=fact*j
              enddo
           endif
        endif

        return
        end function



        
        
 

