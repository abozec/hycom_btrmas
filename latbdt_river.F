      subroutine latbdt_river(n,lll)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_nc         ! Netcdf output (debug)
      implicit none
c
      integer n,lll,il,jl,imid,jmid
c
c --- Apply lateral boundary conditions to   barotropic  flow field 
c
c --- River version:
c --- Uses the 'Browning and Kreiss' MICOM open boundary condition.
c
c --- Note that 'speed' is in fact thref*SQRT(gH) which represents
c --- c1*thref/g in the notation of (Bleck and Sun, Open boundary conditions
c --- for MICOM).  The thref/g allows for the use of pressure fields.
c
c --- The first call is made during initialization.
c
c --- Alan J. Wallcraft,  NRL,  July, 2001.
c --- Alexandra Bozec, COAPS, based on Remy Baraille, SHOM, Apr, 2016 
c
      integer, parameter :: npas=18   
      integer, parameter :: nchar=120
c
      character*10 cfall
      logical      lfatal,lfatalp,rarwave, boundary_mpi,debitf
      integer      i,j,k,isec,ifrst,ilast,l,np,npf,npi,npl
      integer      ite,ios,jname,nrec,margin,nts_day,urivr
      real         aline(nchar)
      real         crs,fin,fatal,xf,xfder,dtime
      character*3  char3
      character*100 cline
c
      integer  ll1, ll2
      logical ifirst

c --- Ports properties and locations      
      integer, save :: nports
      integer, save, allocatable :: kdport(:),
     &                              ifport(:),ilport(:),
     &                              jfport(:),jlport(:),lnport(:)
      real   , save, allocatable :: xpport(:)
      character*10, save, allocatable ::  nameriver(:)

c --- River properties 
      real :: xtemp,xsaln,xdeb
      real,    save, allocatable  :: xtemp1(:),xsaln1(:),
     &                               xtemp2(:),xsaln2(:),
     &                               xtime1(:),xdeb1(:),
     &                               xtime2(:),xdeb2(:)
      real,    save, allocatable  :: pline(:),uline(:,:),xline(:)

c --- Boundary conditions properties 
      real, save, allocatable :: speedw(:,:),rspedw(:,:),speede(:,:),
     &                           rspede(:,:),speedn(:,:),rspedn(:,:),
     &                           speeds(:,:),rspeds(:,:)
c
c --- Loop indices of ports      
      integer, save, allocatable :: ibeg(:), iend(:),jbeg(:), jend(:)
      integer, save, allocatable :: indport(:)
      integer, save, allocatable :: ifport_loc(:), ilport_loc(:),
     &                              lnport_loc(:),
     &                              jfport_loc(:), jlport_loc(:)

c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
      integer lcount
      save    lcount
      data    lcount / 0 /

      if (n.eq.0.and.lcount.ne.0) return

c --- Conservation of transport
c     rarwave = true for invariants, false for characteristics
c     debitf = true for transport conservation, false for SSH conservation
      rarwave = .true.
      boundary_mpi = .false.
      debitf = .true.
c
      lcount = lcount + 1
      urivr  = 1000+19
      margin = 5

c --- Get time 
      nts_day = nint(86400.0d0/baclin)     !no. time steps per day
      dtime=(nstep/nts_day)+mod(nstep,nts_day)*(baclin/86400.0d0)

c
c --- The first call just initializes data structures.
c
      if     (lcount.eq.1) then
c
        open(unit=uoff+99,file=trim(flnminp)//'rivers.list')
c
c --- 'nports' = number of boundary port sections.
        call blkini(nports,'nports')

c --- Array allocations       
        allocate(
     &      nameriver(nports),
     &      kdport(nports),
     &      ifport(nports),
     &      ilport(nports),
     &      jfport(nports),
     &      jlport(nports),
     &      xpport(nports),
     &      lnport(nports) )
        
c --- Read the rivers.list
        do l=1,nports
          call blkinc(nameriver(l),'namriv','(a6," = ",a10)')
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          call blkinr(xpport(l),'xpport','(a6," =",f10.4," %")')
        enddo
c
c --- Sanity check.
c
        do l= 1,nports
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) ' port #', l
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) ' port #', l
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
            write(lp,*) 
              write(lp,*) ' port #', l
            write(lp,*) 'error in latbdt - port',
     &                   ' location is not consistent'
            write(lp,*) 
            call flush(lp)
            endif
            call xcstop('(latbdt)')
                   stop '(latbdt)'
          endif
        enddo
c
        close(unit=99)
c
c --- Check ports against masks,
c --- Mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c --- Western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &                j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &                iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c --- Eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &                j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &                iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c --- Northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &                j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &                iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c --- Southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &                j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &                iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          endif
c
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdt - port ',l,' mislocated',
     &                  '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp
        enddo  !l=1,nports
c
c --- Local lfatal to global lfatal
c
        if     (lfatal) then
          fatal = 1.0
        else
          fatal = 0.0
        endif
        call xcmaxr(fatal)
        lfatal = fatal.gt.0.5
c
c --- Write out  -iu-  and -iv- arrays, if they are not too big
c --- Data are written in strips nchar points wide
c
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
c
          util1(1:ii,1:jj) = iv(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif  ! small region
c
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdt - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdt)')
                 stop '(latbdt)'
        endif
c
c ---  Restore iu and iv, and zero iuopn and ivopn.
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo
c
c --- Modification of the storage of loop indices in case of complex boundaries
c

c --- Array allocations        
        allocate(
     &         jbeg(nports),ibeg(nports),
     &         jend(nports),iend(nports),
     &      indport(nports))
        
c --- Western port
        jbeg(:) = 0
        jend(:) = 0
        ll1 = 0
        do l = 1, nports
           if (kdport(l).eq.4) then
             i = ifport(l)
             j = jfport(l)
             do ll2 = 1, ll1
                if ((j.ge.jbeg(ll1).or.j+lnport(l).ge.jbeg(ll1)).and.
     &              (j.le.jend(ll1).or.j+lnport(l).le.jend(ll1))) then
                    ll1 = ll1 + 1
                endif
             enddo
             ll1 = max(1,ll1)
             jbeg(ll1)  = j
             jend(ll1)  = j+lnport(l)
             indport(l) = ll1
           endif
        enddo
c --- Eastern port
        jbeg(:) = 0
        jend(:) = 0
        ll1 = 0
        do l = 1, nports
           if (kdport(l).eq.3) then
             i = ifport(l)
             j = jfport(l)
             do ll2 = 1, ll1
                if ((j.ge.jbeg(ll1).or.j+lnport(l).ge.jbeg(ll1)).and.
     &              (j.le.jend(ll1).or.j+lnport(l).le.jend(ll1))) then
                    ll1 = ll1 + 1
                endif
             enddo
             ll1 = max(1,ll1)
             jbeg(ll1)  = j
             jend(ll1)  = j+lnport(l)
             indport(l) = ll1
           endif
        enddo
c --- Southern port
        ibeg(:) = 0
        iend(:) = 0
        ll1 = 0
        do l = 1, nports
           if (kdport(l).eq.2) then
             i = ifport(l)
             j = jfport(l)
             do ll2 = 1, ll1
                if ((i.ge.ibeg(ll1).or.i+lnport(l).ge.ibeg(ll1)).and.
     &              (i.le.iend(ll1).or.i+lnport(l).le.iend(ll1))) then
                    ll1 = ll1 + 1
                endif
             enddo
             ll1 = max(1,ll1)
             ibeg(ll1)  = i
             iend(ll1)  = i+lnport(l)
             indport(l) = ll1
           endif
        enddo
c --- Northern port
        ibeg(:) = 0
        iend(:) = 0
        ll1 = 0
        do l = 1, nports
           if (kdport(l).eq.1) then
             i = ifport(l)
             j = jfport(l)
             do ll2 = 1, ll1
                if ((i.ge.ibeg(ll1).or.i+lnport(l).ge.ibeg(ll1)).and.
     &              (i.le.iend(ll1).or.i+lnport(l).le.iend(ll1))) then
                    ll1 = ll1 + 1
                endif
             enddo
             ll1 = max(1,ll1)
             ibeg(ll1)  = i
             iend(ll1)  = i+lnport(l)
             indport(l) = ll1
           endif
        enddo
c
        i = -10
        do l = 1, nports
           i = max(i,indport(l))
        enddo

c
c ---  Initialize the ports
c

c --- Array allocations        
        allocate (speedw(jtdm,i) , rspedw(jtdm,i) , speede(jtdm,i) , 
     &            rspede(jtdm,i) , speedn(itdm,i) , rspedn(itdm,i) , 
     &            speeds(itdm,i) , rspeds(itdm,i))
        allocate(
     &           pline(itdm+jtdm),
     &           uline(itdm+jtdm,2),
     &           xline(itdm+jtdm))
c
        do l= 1,nports
          ll1 = indport(l)
          if     (kdport(l).eq.4) then
c
c --- Western port
c
              xf= 0.
              i = ifport(l)
              j = jfport(l)
              if (boundary_mpi) then
                  call xclget(pline(j)     ,lnport(l),depths,i,j,0,1,0)
                  call xclget(rspedw(j,ll1),lnport(l),onetai,i,j,0,1,0)
                  call xclget(xline(j)     ,lnport(l),  scuy,i,j,0,1,0) !!Alex
              else
                  if (i.ge.i0+1-nbdy.and.i.le.i0+ii+nbdy) then
                  do k = max(1-nbdy,j-j0), min(jj+nbdy,j-j0+lnport(l)-1)
                      pline(k+j0)      = depths(i-i0,k)  
                      rspedw(k+j0,ll1) = onetai(i-i0,k)  
                      xline(k+j0)      = scuy(i-i0,k) !!Alex
                  enddo
                  endif
              endif
              do j= jfport(l),jlport(l)
                xf = xf + xline(j)
              enddo
              do j= jfport(l),jlport(l)
c
                if     (i.ge.i0+ 1-nbdy .and.
     &                  i.le.i0+ii+nbdy .and.
     &                  j.ge.j0+ 1-nbdy .and.
     &                  j.le.j0+jj+nbdy      ) then
                    speedw(j,ll1) = rspedw(j,ll1)*pline(j)*onem
                    rspedw(j,ll1) = sqrt(thref/(onem*pline(j)))
                    iuopn(i-i0,j-j0) = 2 ! river
                else
                    speedw(j,ll1) = onem
                    rspedw(j,ll1) = sqrt(thref/onem)
                endif
              enddo
#if defined(MPI)
              call gdsum(xf,xfder) !! send xf and xfder to all processors
#else
              xfder=xf
#endif
              xpport(l) =  0.01*onem*xpport(l)/max(1.,xf)
c
          elseif (kdport(l).eq.3) then
c
c --- Eastern port
c
              xf= 0.
              i = ifport(l)-1
              j = jfport(l)
              if (boundary_mpi) then
                  call xclget(pline(j)     ,lnport(l),depths,i ,j,0,1,0)
                  call xclget(rspede(j,ll1),lnport(l),onetai,i ,j,0,1,0)
                  call xclget(xline(j)     ,lnport(l), scuy,i+1,j,0,1,0) !!Alex
!RB                  call xclget(xline(j)     ,lnport(l), scuy,i,j,0,1,0)
              else
                  if (i.ge.i0+1-nbdy.and.i.le.i0+ii+nbdy) then
                  do k = max(1-nbdy,j-j0), min(jj+nbdy,j-j0+lnport(l)-1)
                    pline(k+j0)      = depths(i-i0,k)
                    rspede(k+j0,ll1) = onetai(i-i0,k)
                    xline(k+j0)      = scuy(i+1-i0,k) !!Alex
!RB                    xline(k+j0) = scuy(i-i0,k)
                  enddo
                  endif
              endif
              do j= jfport(l),jlport(l)
                xf = xf + xline(j)
              enddo
              do j= jfport(l),jlport(l)
c
                if     (i+1.ge.i0+ 1-nbdy .and.
     &                  i+1.le.i0+ii+nbdy .and.
     &                  j  .ge.j0+ 1-nbdy .and.
     &                  j  .le.j0+jj+nbdy      ) then
                    speede(j,ll1) = rspede(j,ll1)*pline(j)*onem
                    rspede(j,ll1) = sqrt(thref/(onem*pline(j)))
                    iuopn(i-i0+1,j-j0) = 2 ! river
                else
                    speede(j,ll1) = onem
                    rspede(j,ll1) = sqrt(thref/onem)
                endif
              enddo
#if defined(MPI)
              call gdsum(xf,xfder) !! send xf and xfder to all processors
#else
              xfder=xf
#endif
              xpport(l) = -0.01*onem*xpport(l)/max(1.,xf)
c
          elseif (kdport(l).eq.1) then
c
c --- Northern port
c
              xf= 0.
              j = jfport(l)-1
              i = ifport(l)
              if (boundary_mpi) then
                  call xclget(pline(i)     ,lnport(l),depths,i,j ,1,0,0)
                  call xclget(rspedn(i,ll1),lnport(l),onetai,i,j ,1,0,0)
                  call xclget(xline(i)     ,lnport(l), scvx,i,j+1,1,0,0) !!Alex
!RB                  call xclget(xline(i)     ,lnport(l), scvx,i,j,1,0,0)
              else
                  if (j.ge.j0+1-nbdy.and.j.le.j0+jj+nbdy) then
                  do k = max(1-nbdy,i-i0), min(ii+nbdy,i-i0+lnport(l)-1)
                        pline(k+i0)      = depths(k,j-j0)
                        rspedn(k+i0,ll1) = onetai(k,j-j0)
                        xline(k+i0)      = scvx(k,j+1-j0)  !!Alex
!RB                        xline(k+i0) = scvx(k,j-j0)  
                  enddo
                  endif
              endif
              do i= ifport(l),ilport(l)
                xf = xf + xline(i)
              enddo
              do i= ifport(l),ilport(l)
c
                if     (i  .ge.i0+ 1-nbdy .and.
     &                  i  .le.i0+ii+nbdy .and.
     &                  j+1.ge.j0+ 1-nbdy .and.
     &                  j+1.le.j0+jj+nbdy      ) then
                    speedn(i,ll1) = rspedn(i,ll1)*pline(i)*onem
                    rspedn(i,ll1) = sqrt(thref/(onem*pline(i)))
                    ivopn(i-i0,j-j0+1) = 2 ! river
                else
                    speedn(i,ll1) = onem
                    rspedn(i,ll1) = sqrt(thref/onem)
                endif
              enddo
#if defined(MPI)
              call gdsum(xf,xfder) !! send xf and xfder to all processors
#else
              xfder=xf
#endif
              xpport(l) = -0.01*onem*xpport(l)/max(1.,xf)
c
          elseif (kdport(l).eq.2) then
c
c --- Southern port
c
              xf= 0.
              j = jfport(l)
              i = ifport(l)
              if (boundary_mpi) then
                  call xclget(pline(i)     ,lnport(l),depths,i,j,1,0,0)
                  call xclget(rspeds(i,ll1),lnport(l),onetai,i,j,1,0,0)
                  call xclget(xline(i)     ,lnport(l),  scvx,i,j,1,0,0) !!Alex
              else
                  if (j.ge.j0+1-nbdy.and.j.le.j0+jj+nbdy) then
                  do k = max(1-nbdy,i-i0), min(ii+nbdy,i-i0+lnport(l)-1)
                        pline(k+i0)      = depths(k,j-j0) 
                        rspeds(k+i0,ll1) = onetai(k,j-j0) 
                        xline(k+i0)      = scvx(k,j-j0)  !!Alex
                  enddo
                  endif
              endif
              do i= ifport(l),ilport(l)
                xf = xf + xline(i)
              enddo
              do i= ifport(l),ilport(l)
c
                if     (i.ge.i0+ 1-nbdy .and.
     &                  i.le.i0+ii+nbdy .and.
     &                  j.ge.j0+ 1-nbdy .and.
     &                  j.le.j0+jj+nbdy      ) then
                    speeds(i,ll1) = rspeds(i,ll1)*pline(i)*onem 
                    rspeds(i,ll1) = sqrt(thref/(onem*pline(i)))
                    ivopn(i-i0,j-j0) = 2 ! river
                else
                    speeds(i,ll1) = onem
                    rspeds(i,ll1) = sqrt(thref/onem)
                endif
              enddo
c
#if defined(MPI)
              call gdsum(xf,xfder) !! send xf and xfder to all processors
#else
              xfder=xf
#endif
              xpport(l) =  0.01*onem*xpport(l)/max(1.,xf)
          endif
        enddo                   !l=1,nports
c
c -- Update halo of iuopn and ivopn
        util1(:,:) = iuopn(:,:)
        util2(:,:) = ivopn(:,:)
        call xctilr(util1,1,1, nbdy,nbdy, halo_us)
        call xctilr(util2,1,1, nbdy,nbdy, halo_vs)
        iuopn(:,:) = util1(:,:)
        ivopn(:,:) = util2(:,:)
c
c --- Modification of loop indices for ports
c
c --- Arrays allocation 
        allocate(
     &      ifport_loc(nports),ilport_loc(nports),
     &      jfport_loc(nports),jlport_loc(nports),
     &      lnport_loc(nports))

        allocate(
     &          xtemp1(nports),xtemp2(nports),
     &          xsaln1(nports),xsaln2(nports),
     &           xdeb1(nports), xdeb2(nports),
     &          xtime1(nports),xtime2(nports))
        
        do l = 1,nports

c --- Western/Eastern port

           if (kdport(l).eq.4.or.kdport(l).eq.3) then
              i = ifport(l)
              ifirst = .false.
              if (i.ge.1+i0-nbdy.and.i.le.ii+i0+nbdy) then
                 do j = jfport(l),jlport(l)
                    if (j.ge.1+j0-nbdy.and.j.le.jj+j0+nbdy) then
                        if (.not.ifirst) then
                           ifrst = j
                           ifirst = .true.
                        endif
                        ilast = j
                    endif
                 enddo
              endif
              if (ifirst) then
                 jfport_loc(l) = ifrst
                 jlport_loc(l) = ilast
                 lnport_loc(l) = ilast-ifrst+1
                          else
                 lnport_loc(l) = 0
                 jfport_loc(l) = 10000
                 jlport_loc(l) = -10000
              endif
           endif

c --- Southern/Northern port
           
           if (kdport(l).eq.1.or.kdport(l).eq.2) then
              j = jfport(l)
              ifirst = .false.
              if (j.ge.1+j0-nbdy.and.j.le.jj+j0+nbdy) then
                 do i = ifport(l),ilport(l)
                    if (i.ge.1+i0-nbdy.and.i.le.ii+i0+nbdy) then
                        if (.not.ifirst) then
                           ifrst = i
                           ifirst = .true.
                        endif
                        ilast = i
                    endif
                 enddo
              endif
              if (ifirst) then
                 ifport_loc(l) = ifrst
                 ilport_loc(l) = ilast
                 lnport_loc(l) = ilast-ifrst+1
                          else
                 lnport_loc(l) = 0
                 ifport_loc(l) = 10000
                 ilport_loc(l) = -10000
              endif
           endif
        enddo
c
c --- End of initialization
        call xcsync(flush_lp)

        xtime1(:) = 0.0
        xtime2(:) = 0.0
        return
        

      endif ! lcount.eq.1
      
c
c --- 'wellposed' treatment of pressure and velocity fields.
c --- Alternate order of ports in case corners are open.
c
      if     (mod(lll,2).eq.1) then
        npf =  1
        npl =  nports
        npi =  1
      else
        npf =  nports
        npl =  1
        npi = -1
      endif
c
      
      do l = npf, npl, npi
c
        ll1 = indport(l)
c
        if (lll.eq.1.and. (xtime1(l).gt.dtime .or.
     &                     dtime    .gt.xtime2(l)) ) then
            
            open(unit=urivr,file=trim(nameriver(l))//'.r')
c
            j = 0
            do nrec = 1, 9999999
              read(urivr,'(a)',iostat=ios) cline
              if (j.eq.0) then
                  if (cline(1:3).eq.'---') j=1
              else
                  if (ios.ne.0) then
                      if (mnproc.eq.1) write(lp,*)'missing data '//
     &                    'for river:',nameriver(l)
                      call xcstop('(latbdt)')
                      stop '(latbdt)'
                  endif
                  
                  read (cline(1:),*) xtime2(l), xdeb2(l),
     &                               xtemp2(l),xsaln2(l)                  

                  if (j.eq.1.and.
     &                dtime-0.5*baclin/86400..le.xtime2(l)) then
                      if (mnproc.eq.1) write(lp,*)'missing data '//
     &                    'for river:',nameriver(l)
                      call xcstop('(latbdt)')
                      stop '(latbdt)'
                  endif
                  if (dtime-0.5*baclin/86400..le.xtime2(l)) then
                      close(urivr)
                      goto 200
                  endif
                  xtime1(l) = xtime2(l)
                   xdeb1(l) =  xdeb2(l)
                  xtemp1(l) = xtemp2(l)
                  xsaln1(l) = xsaln2(l)
                  j = j + 1
              endif
            enddo
 200        continue
        endif ! lll.eq.1 .and. xtime
            
c --- Calculate the transport of the river for each point of the port
        xdeb  = ( xdeb2(l)+(    dtime-xtime2(l))  *(xdeb1(l)-xdeb2(l))/
     &                     (xtime1(l)-xtime2(l)) )* xpport(l)


        if (lll.eq.1) then

c --- Calculate the T and S of the river for each point of the port
            xtemp = xtemp2(l) + (dtime-xtime2(l))*(xtemp1(l)-xtemp2(l))
     &                                           /(xtime1(l)-xtime2(l))

            xsaln = xsaln2(l) + (dtime-xtime2(l))*(xsaln1(l)-xsaln2(l))
     &                                           /(xtime1(l)-xtime2(l))
c
c --- Western/Eastern port
c --- Fill triver and sriver arrays for tsadvc (advem_fct2c)
c
            if (kdport(l).eq.4.or.kdport(l).eq.3) then
                i = ifport(l)
                do j = jfport(l),jlport(l)
                  pline(j) = xtemp
                enddo
                j = jfport(l)
                do k = 1, kk
                  call xclput(pline(j),lnport(l),
     &                        triver(1-nbdy,1-nbdy,k), i, j, 0, 1)
                enddo
                do j = jfport(l),jlport(l)
                  pline(j) = xsaln
                enddo
                j = jfport(l)
                do k = 1, kk
                  call xclput(pline(j),lnport(l),
     &                        sriver(1-nbdy,1-nbdy,k), i, j, 0, 1)
                enddo
            endif
c
c --- Southern/Northern port
c           
            if (kdport(l).eq.1.or.kdport(l).eq.2) then
                j = jfport(l)
                do i = ifport(l),ilport(l)
                  pline(i) = xtemp
                enddo
                i = ifport(l)
                do k = 1, kk
                  call xclput(pline(i),lnport(l),
     &                        triver(1-nbdy,1-nbdy,k), i, j, 1, 0)
                enddo
                do i = ifport(l),ilport(l)
                  pline(i) = xsaln
                enddo
                i = ifport(l)
                do k = 1, kk
                  call xclput(pline(i),lnport(l),
     &                        sriver(1-nbdy,1-nbdy,k), i, j, 1, 0)
                enddo
            endif
        endif ! lll.eq.1
        

c
        if     (kdport(l).eq.4) then
c
c --- Western port
c
            i = ifport(l)
            j = jfport(l)
            if (boundary_mpi) then
                call xclget(uline(j,1),lnport(l),
     &              ubavg(1-nbdy,1-nbdy,n), i+1,j,  0,1, 0)
                call xclget(uline(j,2),lnport(l),
     &              ubavg(1-nbdy,1-nbdy,n), i+2,j,  0,1, 0)
                call xclget(pline(j),  lnport(l),
     &              pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
            else
                if (i.ge.i0+1-margin.and.i.le.i0+ii+margin) then
                    do k=max(1 -margin,j-j0),
     &                   min(jj+margin,j-j0+lnport(l)-1)
                      pline(k+j0) = pbavg(i-i0,k,n)
                enddo
                endif
                if (i.ge.i0-margin.and.i+1.le.i0+ii+margin) then
                    do k=max(1 -margin,j-j0),
     &                   min(jj+margin,j-j0+lnport(l)-1)
                      uline(k+j0,1) = ubavg(i+1-i0,k,n)
                   enddo
                endif
                if (i+1.ge.i0-margin.and.i+2.le.i0+ii+margin) then
                    do k=max(1 -margin,j-j0),
     &                   min(jj+margin,j-j0+lnport(l)-1)
                      uline(k+j0,2) = ubavg(i+2-i0,k,n)
                    enddo
                endif
            endif ! boundary_mpi
            if (rarwave) then
                do j= max(j0+1 -margin,jfport_loc(l)),
     &                min(jj+j0+margin,jlport_loc(l))
                  fin  = uline(j,1)
     &              - 2.*sqrt(max(1.e-6,(speedw(j,ll1)+pline(j))*thref))
                  do ite = 1, npas
                    xf = 0.25*uline(j,1)*(fin-uline(j,1))*
     &                                   (fin-uline(j,1))/thref - xdeb
                    xfder = 0.25*(fin-   uline(j,1))*
     &                           (fin-3.*uline(j,1))/thref
                    uline(j,1) = uline(j,1) - xf/xfder
                  enddo
                  pline(j) = 0.25*(fin-uline(j,1))
     &                           *(fin-uline(j,1))/thref - speedw(j,ll1)

c
c --- Solver convergence problem 
c 
                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(j)   = -speedw(j,ll1)
                      uline(j,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                        'xfder=',xfder
                      endif
                  endif
                enddo
            else                
                do j= max(j0+1 -margin,jfport_loc(l)),
     &                min(j0+jj+margin,jlport_loc(l))
                  fin  = 1.5*uline(j,1)-0.5*uline(j,2)
     &                                 -rspedw(j,ll1)*pline( j)
                  do ite = 1, npas
                    xf = (speedw(j,ll1)-(fin-uline(j,1))/rspedw(j,ll1))
     &                                      *uline(j,1) - xdeb
                    xfder = speedw(j,ll1)-(fin-2.*uline(j,1))
     &                                                  /rspedw(j,ll1)
                    uline(j,1) = uline(j,1) - xf/xfder
                  enddo
                  pline(j) = -(fin - uline(j,1))/rspedw(j,ll1)

c
c --- Solver convergence problem 
c 
                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(j)   = -speedw(j,ll1)
                      uline(j,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                              'xfder=',xfder
                      endif
                  endif
                enddo
            endif ! rarwave
            j = jfport(l)
            if (boundary_mpi) then
                call xclput(pline(j),  lnport(l),
     &              pbsavu(1-nbdy,1-nbdy  ), i,  j,  0,1)
                call xclput(uline(j,1),lnport(l),
     &              ubavg (1-nbdy,1-nbdy,n), i,  j,  0,1) ! normal
            else
                if (i.ge.i0+1-margin.and.i.le.i0+ii+margin) then
                 do k=max(1-margin,j-j0),min(jj+margin,j-j0+lnport(l)-1)
                      pbsavu(i-i0,k  ) = pline(k+j0)  
                      ubavg (i-i0,k,n) = uline(k+j0,1)
                 enddo
                endif
            endif ! boundary_mpi
c
c
        elseif (kdport(l).eq.3) then
c
c --- Eastern port
c
            i = ifport(l)-1
            j = jfport(l)
            if (boundary_mpi) then
                call xclget(uline(j,1),lnport(l),
     &              ubavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
                call xclget(uline(j,2),lnport(l),
     &              ubavg(1-nbdy,1-nbdy,n), i-1,j,  0,1, 0)
                call xclget(pline(j),  lnport(l),
     &              pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
            else
                if (i.ge.i0+1-margin.and.i.le.i0+ii+margin) then
                    do k = max(1 -margin,j-j0),
     &                     min(jj+margin,j-j0+lnport(l)-1)
                      pline(k+j0)   = pbavg(i-i0,k,n)
                      uline(k+j0,1) = ubavg(i-i0,k,n)
                    enddo
                endif
                if (i.ge.i0+2-margin.and.i-1.le.i0+ii+margin) then
                    do k = max(1 -margin,j-j0),
     &                     min(jj+margin,j-j0+lnport(l)-1)
                      uline(k+j0,2) = ubavg(i-1-i0,k,n)
                    enddo
                endif
            endif ! boundary_mpi
            if (rarwave) then
                do j= max(j0+1 -margin,jfport_loc(l)),
     &                min(j0+jj+margin,jlport_loc(l))
                  
                  fin  = uline(j,1)
     &              + 2.*sqrt(max(1.e-6,(speede(j,ll1)+pline(j))*thref))
                  do ite = 1, npas
                    xf = 0.25*uline(j,1)*(fin-uline(j,1))*
     &                                   (fin-uline(j,1))/thref - xdeb
                    xfder = 0.25*(fin-   uline(j,1))*
     &                           (fin-3.*uline(j,1))/thref
                    uline(j,1) = uline(j,1) - xf/xfder
                  enddo
                  pline(j) = 0.25*(fin-uline(j,1))
     &                           *(fin-uline(j,1))/thref - speede(j,ll1)

c
c --- Solver convergence problem 
c
                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                        'xfder=',xfder
                      endif
                      pline(j) = -speede(j,ll1)
                      uline(j,1) = 0.
                  endif
                enddo
            else
                do j= max(j0+1 -margin,jfport_loc(l)),
     &                min(j0+jj+margin,jlport_loc(l))
                  fin  = 1.5*uline(j,1)-0.5*uline(j,2)
     &                      +rspede(j,ll1)*pline( j)
                  do ite = 1, npas
                    xf = (speede(j,ll1)+(fin-uline(j,1))/rspede(j,ll1))
     &                                      *uline(j,1) - xdeb
                    xfder = speede(j,ll1)+(fin-2.*uline(j,1))
     &                                                  /rspede(j,ll1)
                    uline(j,1) = uline(j,1) - xf/xfder
                  enddo
                  pline(j) = (fin - uline(j,1))/rspede(j,ll1)
c
c --- Solver convergence problem 
c 
                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(j) = -speede(j,ll1)
                      uline(j,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder,
     &                              'xf =',xf,
     &                              'xfder =',xfder
                      endif
                  endif
                enddo
            endif ! rarwave
            j = jfport(l)
            if (boundary_mpi) then
               call xclput(pline(j),  lnport(l),
     &              pbsavu(1-nbdy,1-nbdy  ), i+1,j,0,1)
                call xclput(uline(j,1),lnport(l),
     &              ubavg (1-nbdy,1-nbdy,n), i+1,j,0,1) ! normal
            else
                if (i.ge.i0-margin.and.i+1.le.i0+ii+margin) then
                    do k=max(1 -margin,j-j0),
     &                   min(jj+margin,j-j0+lnport(l)-1)
                      pbsavu(i+1-i0,k  ) = pline(k+j0)
                      ubavg (i+1-i0,k,n) = uline(k+j0,1)
                    enddo
                endif
            endif ! boundary_mpi
c
        elseif (kdport(l).eq.1) then
c
c --- Northern port
c
            j = jfport(l)-1
            i = ifport(l)
            if (boundary_mpi) then
                call xclget(uline(i,1),lnport(l),
     &              vbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
                call xclget(uline(i,2),lnport(l),
     &              vbavg(1-nbdy,1-nbdy,n), i,  j-1,1,0, 0)
                call xclget(pline(i),  lnport(l),
     &              pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
            else
                if (j.ge.j0+1-margin.and.j.le.j0+jj+margin) then
                    do k = max(1 -margin,i-i0),
     &                     min(ii+margin,i-i0+lnport(l)-1)
                      pline(k+i0)   = pbavg(k,j-j0,n)
                      uline(k+i0,1) = vbavg(k,j-j0,n)
                    enddo
                endif
                if (j.ge.j0+2-margin.and.j-1.le.j0+jj+margin) then
                    do k = max(1 -margin,i-i0),
     &                     min(ii+margin,i-i0+lnport(l)-1)
                      uline(k+i0,2) = vbavg(k,j-1-j0,n)
                    enddo
                endif
            endif ! boundary_mpi
            if (rarwave) then
                do i= max(i0+1 -margin,ifport_loc(l)),
     &                min(i0+ii+margin,ilport_loc(l))
                  fin  = uline(i,1)
     &              + 2.*sqrt(max(1.e-6,(speedn(i,ll1)+pline(i))*thref))
                  do ite = 1, npas
                    xf = 0.25*uline(i,1)*(fin-uline(i,1))*
     &                                   (fin-uline(i,1))/thref -xdeb
                    xfder = 0.25*(fin-uline(i,1))*
     &                           (fin-3.*uline(i,1))/thref
                    uline(i,1) = uline(i,1) - xf/xfder
                  enddo
                  pline(i) = 0.25*(fin-uline(i,1))*
     &                            (fin-uline(i,1))/thref - speedn(i,ll1)
c
c --- Solver convergence problem 
c 

                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(i) = -speedn(i,ll1)
                      uline(i,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*) 'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                        'xfder=',xfder
                      endif
                  endif
                enddo
            else
                do i= max(i0+1 -margin,ifport_loc(l)),
     &                min(i0+ii+margin,ilport_loc(l))
                  fin  = 1.5*uline(i,1)-0.5*uline(i,2)
     &                       +rspedn(i,ll1)*pline( i)
                  do ite = 1, npas
                    xf = (speedn(i,ll1)+(fin-uline(i,1))/rspedn(i,ll1))
     &                                      *uline(i,1) - xdeb
                    xfder = speedn(i,ll1)+(fin-2.*uline(i,1))
     &                                                  /rspedn(i,ll1)
                    uline(i,1) = uline(i,1) - xf/xfder
                  enddo
                  pline(i) = (fin - uline(i,1))/rspedn(i,ll1)
c
c --- Solver convergence problem 
c  

                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(i) = -speedn(i,ll1)
                      uline(i,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                              'xfder=',xfder
                      endif
                  endif
                enddo
            endif ! rarwave
            i = ifport(l)
            if (boundary_mpi) then
                call xclput(pline(i),  lnport(l),
     &              pbsavv(1-nbdy,1-nbdy  ), i, j+1,1,0)
                call xclput(uline(i,1),lnport(l),
     &              vbavg (1-nbdy,1-nbdy,n), i, j+1,1,0) ! normal
            else
                if (j.ge.j0-margin.and.j+1.le.j0+jj+margin) then
                    do k=max(1 -margin,i-i0),
     &                   min(ii+margin,i-i0+lnport(l)-1)
                      pbsavv(k,j+1-j0  ) = pline(k+i0)
                      vbavg (k,j+1-j0,n) = uline(k+i0,1)
                    enddo
                endif
            endif ! boundary_mpi
        elseif (kdport(l).eq.2) then
c
c --- Southern port
c
            j = jfport(l)
            i = ifport(l)
            if (boundary_mpi) then
                call xclget(uline(i,1),lnport(l),
     &              vbavg(1-nbdy,1-nbdy,n), i,  j+1,1,0, 0)
                call xclget(uline(i,2),lnport(l),
     &              vbavg(1-nbdy,1-nbdy,n), i,  j+2,1,0, 0)
                call xclget(pline(i),  lnport(l),
     &              pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
            else
                if (j.ge.j0+1-margin.and.j.le.j0+jj+margin) then
                    do k = max(1 -margin,i-i0),
     &                     min(ii+margin,i-i0+lnport(l)-1)
                      pline(k+i0) = pbavg(k,j-j0,n)
                    enddo
                endif
                if (j.ge.j0-margin.and.j+1.le.j0+jj+margin) then
                    do k = max(1 -margin,i-i0),
     &                     min(ii+margin,i-i0+lnport(l)-1)
                      uline(k+i0,1) = vbavg(k,j+1-j0,n)
                    enddo
                endif
                if (j+1.ge.j0-margin.and.j+2.le.j0+jj+margin) then
                    do k = max(1 -margin,i-i0),
     &                     min(ii+margin,i-i0+lnport(l)-1)
                      uline(k+i0,2) = vbavg(k,j+2-j0,n)
                    enddo
                endif
            endif ! boundary_mpi
            if (rarwave) then
                do i= max(i0+1 -margin,ifport_loc(l)),
     &                min(i0+ii+margin,ilport_loc(l))
                  fin  = uline(i,1)
     &              - 2.*sqrt(max(1.e-6,(speeds(i,ll1)+pline(i))*thref))
                  do ite = 1, npas
                    xf = 0.25*uline(i,1)*(fin-uline(i,1))*
     &                                   (fin-uline(i,1))/thref - xdeb
                    xfder = 0.25*(fin-uline(i,1))*
     &                           (fin-3.*uline(i,1))/thref
                    uline(i,1) = uline(i,1) - xf/xfder
                  enddo
                  pline(i) = 0.25*(fin-uline(i,1))*
     &                            (fin-uline(i,1))/thref - speeds(i,ll1)

c
c --- Solver convergence problem 
c 

                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(i) = -speeds(i,ll1)
                      uline(i,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                              'xfder=',xfder
                      endif
                  endif
                enddo
            else
                do i= max(i0+1 -margin,ifport_loc(l)),
     &                min(i0+ii+margin,ilport_loc(l))
                  fin  = 1.5*uline(i,1)-0.5*uline(i,2)
     &                                 -rspeds(i,ll1)*pline( i)
                  do ite = 1, npas
                    xf = (speeds(i,ll1)-(fin-uline(i,1))/rspeds(i,ll1))
     &                                      *uline(i,1) - xdeb
                    xfder = speeds(i,ll1)-(fin-2.*uline(i,1))
     &                                                  /rspeds(i,ll1)
                    uline(i,1) = uline(i,1) - xf/xfder
                  enddo
                  pline(i) = -(fin - uline(i,1))/rspeds(i,ll1)

c
c --- Solver convergence problem 
c

                  if (debitf.and.abs(xf/xfder).ge.1e-3) then
                      pline(i) = -speeds(i,ll1)
                      uline(i,1) = 0.
                      if (mnproc.eq.1) then
                          write(6,*) 'Solver convergence pb in ',i,',',j
                          write(6,*)'xf/xfder = ', xf/xfder, 'xf=',xf,
     &                              'xfder=',xfder
                      endif
                  endif
                enddo
            endif ! rarwave
            i = ifport(l)
            if (boundary_mpi) then
                call xclput(pline(i),  lnport(l),
     &              pbsavv(1-nbdy,1-nbdy), i,  j, 1, 0)
                call xclput(uline(i,1),lnport(l),
     &              vbavg (1-nbdy,1-nbdy,n), i,  j, 1, 0) ! normal
            else
                if (j.ge.j0+1-margin.and.j.le.j0+jj+margin) then
                    do k=max(1 -margin,i-i0),
     &                   min(ii+margin,i-i0+lnport(l)-1)
                      pbsavv(k,j-j0  ) = pline(k+i0)   
                      vbavg (k,j-j0,n) = uline(k+i0,1)
                    enddo
                endif
            endif ! boundary_mpi
c
        endif ! kdports     
c
      enddo !l=1,nports
c
      
      return
      end
