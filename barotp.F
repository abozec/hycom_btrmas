#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif
      subroutine barotp(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
      use mod_tides      ! HYCOM tides
#if defined(STOKES)
      use mod_stokes     !    HYCOM Stokes Drift
#endif
      use mod_nc         ! Netcdf output (debug)
      implicit none
c
      integer m,n
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations.
c ---   on entry: -n- is time t-dt, -m- is time t
c ---   on exit:                    -m- is time t, -n- is time t+dt
c ---   time level 3 is only used internally (n and m are always 1 or 2).
c
c --- LeapFrog version based on:
c ---   Y. Morel, Baraille, R., Pichon A. (2008) "Time splitting and
c ---   linear stability of the slow part of the barotropic component", 
c ---   Ocean Modeling, 23, pp 73-81.
c --- ------------------------------------------------------------------------
c
      logical    lpipe_barotp
      parameter (lpipe_barotp=.false.)
      logical    ldebug_barotp
      parameter (ldebug_barotp=.false.)
c
      real    q,pbudel,pbvdel,utndcy,vtndcy,wblpf,z1
      real    d11,d12,d21,d22,ubp,vbp
      real*8  sump
      integer i,j,l,k,lll,ml,nl,mn,lstep1,margin,mbdy,icof
c	 & ,iffstep
      logical ldrag
c	  data iffstep/0/
c	  save iffstep
c
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::
     &        pbavo,ubavo,vbavo,displd,gslpr,gtide,
     &        flxloc,flyloc,uflxba,vflxba
c
      if     (.not.allocated(pbavo)) then
        allocate(
     &          pbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          ubavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         displd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          gslpr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), 
     &          gtide(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         flxloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         flyloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), 
     &         uflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         vflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 6*(idm+2*nbdy)*(jdm+2*nbdy) )
                pbavo = r_init
                ubavo = r_init
                vbavo = r_init
               displd = r_init
                gslpr = r_init
                gtide = r_init
               flxloc = r_init
               flyloc = r_init
               uflxba = r_init
               vflxba = r_init
      endif
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        pbavo,ubavo,vbavo,displd,gslpr,gtide,
     &        flxloc,flyloc,uflxba,vflxba
#endif
c
      mbdy = 6
c
      margin = mbdy
c --- atmospheric pressure forcing
c
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if     (mslprf) then
            if     (natm.eq.2) then
              gslpr(i,j) =( mslprs(i,j,l0)*w0+
     &                      mslprs(i,j,l1)*w1 )*thref
            else
              gslpr(i,j) =( mslprs(i,j,l0)*w0+
     &                      mslprs(i,j,l1)*w1+
     &                      mslprs(i,j,l2)*w2+
     &                      mslprs(i,j,l3)*w3 )*thref
            endif !natm
          else
            gslpr(i,j) = 0.0
          endif !mslprf
c
          if (SEA_P) then
c
c ---       tidal body forcing, including Scalar SAL,
c ---       SAL should only be applied to the mass anomally, i.e. to
c ---       non-steric SSH, so use sshflg>0 except in tides-only cases.
            if     (tidflg.gt.0 .and. sshflg.eq.0) then !tides
              gtide(i,j)=-g*etide(i,j)
     &                   -salfac(i,j)* srfhgt(i,j)
            elseif (tidflg.gt.0 .and. sshflg.ne.0) then !tides
              gtide(i,j)=-g*etide(i,j)
     &                   -salfac(i,j)*(srfhgt(i,j)-steric(i,j))
            else
              gtide(i,j)=0.0
            endif !tides (sshflg)
          endif !ip
        enddo !i
      enddo !j
c
c --- utotn,vtotn from momtum is time step t-1 to t+1 barotropic tendency
      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_vv)

cCL/RB MPI bug correction 2011-01 !!Alex
      call xctilr(oneta(  1-nbdy,1-nbdy,  n),1, 1, 6,6, halo_ps)
      call xctilr(oneta(  1-nbdy,1-nbdy,  m),1, 1, 6,6, halo_ps)
c
      if     (lpipe .and. lpipe_barotp) then
c ---   compare two model runs.
        call pipe_compare_sym2(utotn, iu,'barotp:utotn',
     &                         vtotn, iv,'barotp:vtotn')
        call pipe_compare_sym1(pvtrop,iq,'barotp:pvtrp')
      endif
c
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
c
      if (btrmas) then
          uflxba(:,:)   = 0.0
          vflxba(:,:)   = 0.0
          flxloc(:,:)   = 0.0
          flyloc(:,:)   = 0.0
          onetao(:,:,n) = oneta(:,:,n)
      endif

      if     (btrlfr) then !btrmas implies btrlfr
        if     (delt1.ne.baclin) then  !not on very 1st time step
c ---     start at time level t-dt and go to t+dt.
          lstep1 = lstep + lstep  !more stable, but also more expensive
          icof   = 2
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = pbavg(i,j,n)  !save t-1 for RA filter
              ubavo(i,j)   = ubavg(i,j,n)  !save t-1 for RA filter
              vbavo(i,j)   = vbavg(i,j,n)  !save t-1 for RA filter
c
              pbavg(i,j,3) = pbavg(i,j,n)
              ubavg(i,j,3) = ubavg(i,j,n)
              vbavg(i,j,3) = vbavg(i,j,n)
            enddo !i
          enddo !j
        else !1st time step
c ---     start at time level t and go to t+dt.
          lstep1 = lstep
          icof   = 1
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = 0.0 !makes correct mean height safe
              pbavg(i,j,n) = pbavg(i,j,m)
              ubavg(i,j,n) = ubavg(i,j,m)
              vbavg(i,j,n) = vbavg(i,j,m)
              pbavg(i,j,3) = pbavg(i,j,m)
              ubavg(i,j,3) = ubavg(i,j,m)
              vbavg(i,j,3) = vbavg(i,j,m)
            enddo !i
          enddo !j
        endif !usual:1st time step
      else
c ---   start at time level t    and go to t+dt.
        lstep1 = lstep          !original, less stable, method
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            pbavo(i,j)   = 0.0 !makes correct mean height safe
            pbavg(i,j,n) = pbavg(i,j,m)
            ubavg(i,j,n) = ubavg(i,j,m)
            vbavg(i,j,n) = vbavg(i,j,m)
          enddo !i
        enddo !j
      endif !btrlfr
c
      ldrag = tidflg.gt.0 .and. drgscl.ne.0.0 .and. thkdrg.eq.0.0
c
      if     (ldrag) then
        displd(:,:) = 0.0
      endif
c
c --- time step loop
c
      if     (btrlfr) then
        wblpf = 0.0   !1st minor time step, lll=1, only
      else
        wblpf = wbaro
      endif
c
      do 840 lll=1,lstep1,2
c
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- odd minor time step.
c
      ml=n
      nl=3
c
c --- continuity equation, and tidal drag on p-grid
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      if     (btrmas) then
c
        margin = mbdy
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
              uflxba(i,j) = uflxba(i,j) +
     &            coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
            enddo !i
            i=ifu(j,l)-1
            if (i.ge.1-margin) then
                if (iuopn(i,j).ne.0) then
                    flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
                    uflxba(i,j) = uflxba(i,j) +
     &                  coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
                endif !iuopn
            endif !i
            i=ilu(j,l)+1
            if (i.le.ii+margin) then
                if (iuopn(i,j).ne.0) then
                    flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
                    uflxba(i,j) = uflxba(i,j) +
     &                  coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
                endif !iuopn
            endif !i
          enddo !l
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
              vflxba(i,j) = vflxba(i,j) +
     &            coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
            enddo !i  
          enddo !l
        enddo !j
c
        do i=1-margin,ii+margin
          do l=1,jsv(i)
            j=jfv(i,l)-1
            if (j.ge.1-margin) then
                if (ivopn(i,j).ne.0) then
                    flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
                    vflxba(i,j) = vflxba(i,j) +
     &                  coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
                endif !ivopn
            endif !j
            j=jlv(i,l)+1
            if (j.le.jj+margin) then
                if (ivopn(i,j).ne.0) then
                    flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
                    vflxba(i,j) = vflxba(i,j) +
     &                  coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
                endif !ivopn
            endif !j
          enddo !l
        enddo !i
c
        margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ubp,vbp,d11,d12,d21,d22,q)
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
!           do l=1,isp(j)
!             do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           do i=1-margin,ii+margin
             if (SEA_P) then
                pbavg(i,j,nl)=
     &            ((1.0-wblpf)*pbavg(i,j,ml)+
     &                  wblpf *pbavg(i,j,nl) )-
     &             (1.0+wblpf)*dlt*(flxloc(i+1,j)-flxloc(i,j) +
     &                              flyloc(i,j+1)-flyloc(i,j)  )*
     &                             scp2i(i,j)
c
                if     (ldrag) then
c
c ---             tidal drag tensor on p-grid:
c ---               ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---               vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---             solve implicitly by inverting the matrix:
c ---                1+(dlt/H)*t.11    (dlt/H)*t.12
c ---                  (dlt/H)*t.21  1+(dlt/H)*t.22
c ---             use depths (H) rather than onem*pbavg (h) for stability.
c
                  ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
                  vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
                  d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
                  d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
                  d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
                  d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
                  q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---             set util5,util6 to the ubavg,vbavg drag increment
                  util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
                  util6(i,j) = q*(ubp*d21+vbp*(1.0-d12)) - vbp
c ---             add an explicit antidrag correction
*                 util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                       d12*vntide(i,j) )
*                 util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                       d22*vntide(i,j) )
c
*                 if (ldebug_barotp .and.
*    &                i.eq.itest.and.j.eq.jtest) then
*                   write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &                nstep,i+i0,j+j0,lll,
*    &                'ubp,new,vbp,new =',
*    &              ubp,ubp+util5(i,j),
*    &              vbp,vbp+util6(i,j)
*                 endif !debug
                else
                  util5(i,j) = 0.0
                  util6(i,j) = 0.0
                endif !ldrag
            endif 
          enddo !i
        enddo !j
      else !.not.btrmas

          margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel,
!$OMP&                    ubp,vbp,d11,d12,d21,d22,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
c
c   Barotropic Stokes flow included here
c
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))*
     &                    (depthu(i+1,j)*scuy(i+1,j))
     &              -(ubavg(i,  j,ml)+usdbavg(i,  j))*
     &                    (depthu(i,  j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))*
     &                    (depthv(i,j+1)*scvx(i,j+1))
     &              -(vbavg(i,j,  ml)+vsdbavg(i,j  ))*
     &                    (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))
#endif
            pbavg(i,j,nl)=
     &        ((1.-wblpf)*pbavg(i,j,ml)+
     &             wblpf *pbavg(i,j,nl) )-
     &         (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
c
            if     (ldrag) then
c
c ---         tidal drag tensor on p-grid:
c ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---         solve implicitly by inverting the matrix:
c ---            1+(dlt/H)*t.11    (dlt/H)*t.12
c ---              (dlt/H)*t.21  1+(dlt/H)*t.22
c ---         use depths (H) rather than onem*pbavg (h) for stability.
c
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
c ---         add an explicit antidrag correction
*             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                   d12*vntide(i,j) )
*             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                   d22*vntide(i,j) )
c ---         dissipation per m^2
              displd(i,j) = displd(i,j) +
     &                      (ubp*util5(i,j) + vbp*util6(i,j))*
     &                      depths(i,j)*qthref/dlt
c
*             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &            nstep,i+i0,j+j0,lll,
*    &            'ubp,new,vbp,new =',
*    &          ubp,ubp+util5(i,j),
*    &          vbp,vbp+util6(i,j)
*             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
      endif !btrmas:else
c
      mn=ml
c
c --- u momentum equation, 1st
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)-
     &                    (gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)-
     &                    (gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wblpf)*ubavg(i,j,ml)+
     &             wblpf *ubavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(utndcy+utotn(i,j))+
     &                0.5*(util5(i,j)+util5(i-1,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'u_old,u_new,p_grad,t_g,m_g,corio,u_star,drag =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &           -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &                 -(gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)*dlt,
*    &                 -(gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt,
*    &          0.5*(util5(i,j)+util5(i-1,j))
*           endif !debug
          endif !iu
        enddo !i
      enddo !j
c
      mn = nl
c
c --- v momentum equation, 2nd
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &                    (gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)-
     &                    (gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wblpf)*vbavg(i,j,ml)+
     &             wblpf *vbavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+
     &                0.5*(util6(i,j)+util6(i,j-1))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'v_old,v_new,p_grad,t_g,m_g,corio,v_star,drag =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &                -(gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)*dlt,
*    &                -(gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt,
*    &          0.5*(util6(i,j)+util6(i,j-1))
*           endif !debug
          endif !iv
        enddo !i
      enddo !j
c
*     if     (ldebug_barotp) then
*       call xcsync(flush_lp)
*     endif
c
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll)
      endif
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- even minor time step.
c
      ml=3
      nl=n
      wblpf = wbaro  !used for all subsequent time steps: lll=2,lstep1
c
c --- continuity equation
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      if     (btrmas) then
         margin = mbdy
c
!$OMP    PARALLEL DO PRIVATE(j,l,i)
!$OMP&            SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
!           do l=1,isu(j)
!             do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
           do i=1-margin,ii+margin
             if (SEA_U) then
                 flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
                 uflxba(i,j) = uflxba(i,j) +
     &               coeflx(lll+2,icof)*(1.0+wblpf)*flxloc(i,j)
             endif 
           enddo !i
         enddo !j
c
!$OMP    PARALLEL DO PRIVATE(j,l,i)
!$OMP&            SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
!           do l=1,isv(j)
!             do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
           do i=1-margin,ii+margin
             if (SEA_V) then
                flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
                vflxba(i,j) = vflxba(i,j) +
     &                        coeflx(lll+2,icof)*(1.0+wblpf)*flyloc(i,j)
            endif 
          enddo !i
        enddo !j
c
         margin = mbdy - 1
c
!$OMP    PARALLEL DO PRIVATE(j,l,i,ubp,vbp,d11,d12,d21,d22,q)
!$OMP&            SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
!           do l=1,isp(j)
!             do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           do i=1-margin,ii+margin
             if (SEA_P) then
                pbavg(i,j,nl)=
     &            ((1.0-wblpf)*pbavg(i,j,ml)+
     &                  wblpf *pbavg(i,j,nl) )-
     &             (1.0+wblpf)*dlt*(flxloc(i+1,j)-flxloc(i,j) +
     &                              flyloc(i,j+1)-flyloc(i,j)  )*
     &                             scp2i(i,j)
c
                if     (ldrag) then
c
c ---             tidal drag tensor on p-grid:
c ---               ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---               vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---             solve implicitly by inverting the matrix:
c ---                1+(dlt/H)*t.11    (dlt/H)*t.12
c ---                  (dlt/H)*t.21  1+(dlt/H)*t.22
c ---             use depths (H) rather than onem*pbavg (h) for stability.
c
                  ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
                  vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
                  d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
                  d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
                  d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
                  d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
                  q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---             set util5,util6 to the ubavg,vbavg drag increment
                  util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
                  util6(i,j) = q*(ubp*d21+vbp*(1.0-d12)) - vbp
c ---             add an explicit antidrag correction
*                 util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                       d12*vntide(i,j) )
*                 util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                       d22*vntide(i,j) )
c
*                 if (ldebug_barotp .and.
*    &                i.eq.itest.and.j.eq.jtest) then
*                   write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &                nstep,i+i0,j+j0,lll,
*    &                'ubp,new,vbp,new =',
*    &              ubp,ubp+util5(i,j),
*    &              vbp,vbp+util6(i,j)
*                 endif !debug
                else
                  util5(i,j) = 0.0
                  util6(i,j) = 0.0
                endif !ldrag
            endif 
          enddo !i
        enddo !j
      else !.not.btrmas
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel,
!$OMP&                    ubp,vbp,d11,d12,d21,d22,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
c
c           Barotropic Stokes flow included here
c
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))*
     &                    (depthu(i+1,j)*scuy(i+1,j))
     &              -(ubavg(i,  j,ml)+usdbavg(i,  j))*
     &                     (depthu(i ,j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))*
     &                    (depthv(i,j+1)*scvx(i,j+1))
     &              -(vbavg(i,j,  ml)+vsdbavg(i,j  ))*
     &                    (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))

#endif
            pbavg(i,j,nl)=
     &        ((1.-wblpf)*pbavg(i,j,ml)+
     &             wblpf *pbavg(i,j,nl) )-
     &         (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
c
            if     (ldrag) then
c ---         tidal drag tensor on p-grid:
c ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---         solve implicitly by inverting the matrix:
c ---            1+(dlt/H)*t.11    (dlt/H)*t.12
c ---              (dlt/H)*t.21  1+(dlt/H)*t.22
c ---         use depths (H) rather than onem*pbavg (h) for stability.
c
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
c ---         add an explicit antidrag correction
*             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                   d12*vntide(i,j) )
*             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                   d22*vntide(i,j) )
c ---         dissipation per m^2
              displd(i,j) = displd(i,j) +
     &                      (ubp*util5(i,j) + vbp*util6(i,j))*
     &                      depths(i,j)*qthref/dlt
c
*             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &            nstep,i+i0,j+j0,lll+1,
*    &            'ubp,new,vbp,new =',
*    &          ubp,ubp+util5(i,j),
*    &          vbp,vbp+util6(i,j)
*             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
      endif !btrmas:else
c
      mn=ml
c
c --- v momentum equation, 1st
c
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &                    (gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)-
     &                    (gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wblpf)*vbavg(i,j,ml)+
     &             wblpf *vbavg(i,j,nl))+
     &         (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+
     &                0.5*(util6(i,j)+util6(i,j-1))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'v_old,v_new,p_grad,t_g,m_g,corio,v_star,drag =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &                -(gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)*dlt,
*    &                -(gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt,
*    &          0.5*(util6(i,j)+util6(i,j-1))
*           endif !debug
          endif !iv
        enddo !i
      enddo !j
c
      mn=nl
c
c --- u momentum equation, 2nd
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)-
     &                    (gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)-
     &                    (gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wblpf)*ubavg(i,j,ml)+
     &             wblpf *ubavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(utndcy+utotn(i,j))+
     &                0.5*(util5(i,j)+util5(i-1,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,7g15.6)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'u_old,u_new,p_grad,t_g,m_g,corio,u_star,drag =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &                -(gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)*dlt,
*    &                -(gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt,
*    &          0.5*(util5(i,j)+util5(i-1,j))
*           endif !debug
          endif !iu
        enddo !i
      enddo !j
c
*     if     (ldebug_barotp) then
*       call xcsync(flush_lp)
*     endif
c
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll+1)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll+1)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll+1)
      endif
c
 840  continue  ! lll=1,lstep1,2
c
      if     (ldrag) then  !disp_count updated in momtum
        displd_mn(:,:) = displd_mn(:,:) + displd(:,:)/real(lstep1)
      endif
c
      if     (lbflag.eq.1) then
c
c ---   correct mean height.
c ---   this should not be required - so there may be a bug in the bc.
c
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1-margin,ii+margin
            if (SEA_P) then
              util1(i,j) = pbavg(i,j,n)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
        call xcsum(sump, util1,ipa)
        q = sump/area
c
c ---   rhs: pbavg
c ---   lhs: pbavg
c
        margin = 0
c
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pbavo(i,j)   = pbavo(i,j)   - q
              pbavg(i,j,1) = pbavg(i,j,1) - q
              pbavg(i,j,2) = pbavg(i,j,2) - q
              pbavg(i,j,3) = pbavg(i,j,3) - q
            endif !ip
          enddo !i
        enddo !j
      endif ! lbflag.eq.1
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barotp:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barotp:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barotp:pbav3')
      endif
c
      if     (btrlfr .and. delt1.ne.baclin) then  !not on very 1st time step
c ---   Robert-Asselin time filter 
!$OMP   PARALLEL DO PRIVATE(j,i,q)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            q = 0.5*ra2fac*(    pbavo(i,j)   +   !t-1
     &                          pbavg(i,j,n) -   !t+1
     &                      2.0*pbavg(i,j,m)  )  !t
            pbavg(i,j,m) = pbavg(i,j,m) + q
            q = 0.5*ra2fac*(    ubavo(i,j)   +   !t-1
     &                          ubavg(i,j,n) -   !t+1
     &                      2.0*ubavg(i,j,m)  )  !t
            ubavg(i,j,m) = ubavg(i,j,m) + q
            q = 0.5*ra2fac*(    vbavo(i,j)   +   !t-1
     &                          vbavg(i,j,n) -   !t+1
     &                      2.0*vbavg(i,j,m)  )  !t
            vbavg(i,j,m) = vbavg(i,j,m) + q
          enddo !i
        enddo !j
      endif !btrlfr & not on very 1st time step
c
      if (btrmas) then
c
c ---   update 1 + heta
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
!          do l=1,isp(j)
!            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          do i=1-margin,ii+margin
            if (SEA_P) then
                oneta(i,j,n)  = onetai(i,j) + pbavg(i,j,n)/pbot(i,j)
            endif 
          enddo !i
        enddo !j
c
c ---   Corrections imposed on dp in order to conserve mass
        do j=1-margin,jj+margin
!          do l=1,isp(j)
!            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          do i=1-margin,ii+margin
            if (SEA_P) then
                onetao(i,j,m) = onetai(i,j) + pbavg(i,j,m)/pbot(i,j)
            endif 
          enddo !i
        enddo !j
c
        uflxba(:,:) = dlt*uflxba(:,:)/delt1
        vflxba(:,:) = dlt*vflxba(:,:)/delt1
c
        call xctilr(uflxba(  1-nbdy,1-nbdy),  1, 1, 2,2, halo_uv)
        call xctilr(vflxba(  1-nbdy,1-nbdy),  1, 1, 2,2, halo_vv)
        call xctilr(onetacnt(1-nbdy,1-nbdy),  1, 1, 2,2, halo_ps)
        call xctilr(oneta(   1-nbdy,1-nbdy,n),1, 1, 6,6, halo_ps)
        call xctilr(onetao(  1-nbdy,1-nbdy,m),1, 1, 6,6, halo_ps)
c
c                Be sure that the sum over the vertical of the baroclinic
c                mass fluxes is exactly the same as the barotropic fluxes
c                --------------------------------------------------------
c
        margin = 1
        util1(:,:) = 0.0
        util2(:,:) = 0.0
        do k =1,kk
!$OMP     PARALLEL DO PRIVATE(j,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
!            do l=1,isp(j)
!              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do i=1-margin,ii+margin
              if (SEA_P) then
                  dp(i,j,k,n)=onetacnt(i,j)*dp(i,j,k,n)
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              endif  
            enddo !i
            do l=1,isu(j)
              do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
                  util1(i,j) = util1(i,j) + uflx(i,j,k)
                enddo !i
                i=ifu(j,l)-1
                if (i.ge.1-margin) then
                    if (iuopn(i,j).ne.0) then
                        util1(i,j) = util1(i,j)+uflx(i,j,k)
                    endif
                endif
                i=ilu(j,l)+1
                if (i.le.ii+margin) then
                    if (iuopn(i,j).ne.0) then
                        util1(i,j) = util1(i,j)+uflx(i,j,k)
                    endif
                endif
              enddo !l
              do l=1,isv(j)
                do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
                  util2(i,j) = util2(i,j) + vflx(i,j,k)
                enddo !i
              enddo !l
            enddo !j
!$OMP     END PARALLEL DO
c
            do i=1-margin,ii+margin
              do l=1,jsv(i)
                j=jfv(i,l)-1
                if (j.ge.1-margin) then
                    if (ivopn(i,j).ne.0) then
                        util2(i,j) = util2(i,j)+vflx(i,j,k)
                    endif
                endif
                j=jlv(i,l)+1
                if (j.le.jj+margin) then
                    if (ivopn(i,j).ne.0) then
                        util2(i,j) = util2(i,j)+vflx(i,j,k)
                    endif
                endif
              enddo !l
            enddo !i
c
          enddo  !k
c
        do k= 1,kk
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,z1)
!$OMP&           SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
c
            do l=1,isu(j)
              do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
                if (uflxba(i,j)-util1(i,j).ge.0.) then
                    z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)
                else
                    z1=dp(i  ,j,k,n)/p(i  ,j,kk+1)
                endif
                uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
              enddo             !i
              i=ifu(j,l)-1
              if (i.ge.1-margin  ) then
                  if (iuopn(i,j).ne.0) then
                      z1=dp(i,j,k,n)/p(i,j,kk+1)
                      uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                      uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
                  endif
              endif
              i=ilu(j,l)+1
              if (i.le.ii+margin ) then
                  if (iuopn(i,j).ne.0) then
                      z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)
                      uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                      uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
                  endif
              endif
            enddo 
c
            do l=1,isv(j)
              do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
                if (vflxba(i,j)-util2(i,j).ge.0.) then
                    z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)
                else
                    z1=dp(i,j  ,k,n)/p(i,j  ,kk+1)
                endif
                vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
              enddo !i
            enddo !l
c
          enddo !j
!$OMP     END PARALLEL DO
c
          do i=1-margin,ii+margin
            do l=1,jsv(i)
              j=jfv(i,l)-1
              if (j.ge.1-margin  ) then
                  if (ivopn(i,j).ne.0) then
                      z1=dp(i,j,k,n)/p(i,j,kk+1)
                      vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                      vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
                  endif
              endif
              j=jlv(i,l)+1
              if (j.le.jj+margin  ) then
                  if (ivopn(i,j).ne.0) then
                      z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)
                      vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                      vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
                  endif
              endif
            enddo !l
          enddo !i
c
!$OMP     PARALLEL DO PRIVATE(j,l,i,dpmin)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
!            do l=1,isp(j)
!              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do i=1-margin,ii+margin
              if (SEA_P) then
                  dp(i,j,k,n)=dp(i,j,k,n)-
     &                   ((uflux(i+1,j)-uflux(i,j))+
     &                    (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              endif 
            enddo !i
          enddo !j
c
        enddo !k
c
!$OMP   PARALLEL DO PRIVATE(j,l,k,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
!          do l=1,isp(j)
            do k=1,kk
!              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              do i=1-margin,ii+margin
                if (SEA_P) then
                    if (p(i,j,kk+1) .gt. 0.0) then
                        dp(i,j,k,n)=dp(i,j,k,n)*pbot(i,j)/p(i,j,kk+1)
                    else
                        dp(i,j,k,n)=0.
                        dp(i,j,1,n)=pbot(i,j)
                    endif
                    p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                    if (isopyc .and. k.eq.1) then
                        dpmixl(i,j,n)=dp(i,j,k,n)
                    endif
                endif 
              enddo !i
            enddo !k
          enddo !j
!$OMP   END PARALLEL DO
c
        call xctilr(dp(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
c
      endif  !btrmas

      return
      end subroutine barotp
c
      subroutine barotp_dry_and_nodry(m,n)
      use mod_xc        ! HYCOM communication interface
      use mod_pipe      ! HYCOM debugging interface
      use mod_cb_arrays ! HYCOM saved arrays
      use mod_nc ! Netcdf output debug
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations from baroclinic time level -m- to level -n
c --- ------------------------------------------------------------------------
c --- ------------------------------------------------------------------------
c --- advance barotropic equations.
c ---   on entry: -n- is time t-dt, -m- is time t
c ---   on exit:                    -m- is time t, -n- is time t+dt
c ---   time level 3 is only used internally (n and m are always 1 or 2).
c
c --- LeapFrog version based on:
c ---   Y. Morel, Baraille, R., Pichon A. (2008) "Time splitting and
c ---   linear stability of the slow part of the barotropic component",
c ---   Ocean Modeling, 23, pp 73-81.
c
c ---  Wetting and Drying only at needed points (where depths < depth_wet)
c
c===========================================================--============
c
      implicit none

      real  z1, wblpf
c
      integer m, n, ml, nl, mn
      integer lstep1, icof
      integer i, j, l, k, lll, ia, ib, ja, jb,margin,mbdy
c
c     Mass conservation
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::
     &        pbavo,ubavo,vbavo,tmp,
     &        flxloc,flyloc,uflxba,vflxba
c
      if     (.not.allocated(pbavo)) then
          allocate(
     &          tmp  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          pbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          ubavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         flxloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         flyloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         uflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         vflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 6*(idm+2*nbdy)*(jdm+2*nbdy) )
                pbavo = r_init
                ubavo = 0.0
                vbavo = r_init
               flxloc = 0.0
               flyloc = 0.0
               uflxba = r_init
               vflxba = r_init
               tmp=0.0
      endif
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        pbavo,ubavo,vbavo,
     &        flxloc,flyloc,uflxba,vflxba
#endif
c
c=======================================================================
c
      if (btrlfr .and. nstep .gt. 1) then
         lstep1 = lstep+lstep
         icof = 2
      else
         lstep1 = lstep
         icof = 1
      endif
c
!      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 2,2, halo_uv)
!      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 2,2, halo_vv)

      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_vv)
      call xctilr(oneta(  1-nbdy,1-nbdy,1  ),1, 2, 6,6, halo_ps)
      call xctilr(onetai(  1-nbdy,1-nbdy   ),1, 1, 6,6, halo_ps)
c
      margin = nbdy
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            onetao(i,j,n) = oneta(i,j,n)
      enddo
      enddo
      enddo
c
      uflxba(:,:) = 0.
      vflxba(:,:) = 0.


c                                                 Save t-1 for RA filter
c                                                 ----------------------
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubavo(i,j) = (1.-ra2fac)*ubavg(i,j,m)+.5*ra2fac*ubavg(i,j,n)
         enddo
         enddo
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbavo(i,j) = (1.-ra2fac)*vbavg(i,j,m)+.5*ra2fac*vbavg(i,j,n)
         enddo
         enddo
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbavo(i,j) = (1.-ra2fac)*pbavg(i,j,m)+.5*ra2fac*pbavg(i,j,n)
         enddo
         enddo
      enddo
c
      if (btrlfr) then
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
            do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               ubavg(i,j,m) = ubavg(i,j,n)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               vbavg(i,j,m) = vbavg(i,j,n)
            enddo
            enddo
            do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               pbavg(i,j,m) = pbavg(i,j,n)
            enddo
            enddo
         enddo
c                           wbaro is set to zero for the first time step
c                           --------------------------------------------
         wblpf = 0.
                  else
         wblpf = wbaro
      endif
c                               ml = n : t-1, mn = 3 : t, nl = m : t+1
c                               - - - - - - - - - - - - - - - - - - - 
      ml = n
      mn = 3
      nl = m

c
      do 840 lll = 1, lstep1
c
         i  = ml
         ml = mn
         mn = nl
         nl = i
c
      call xctilr(pbavg(  1-nbdy,1-nbdy,mn ),1, 1, 4,4, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,mn ),1, 1, 3,3, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,mn ),1, 1, 3,3, halo_vv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,ml ),1, 1, 1,1, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,ml ),1, 1, 1,1, halo_vv)

c
c                         Continuity
c                         ==========
c                                                   Boundary conditions
c                                                   --------------------
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll)
      endif
c
      margin = 3
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               util3(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,mn))
            enddo
         enddo
      enddo
c                              first order u-flux for the wet./dry. area
c                              -----------------------------------------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
               if(ubavg(i,j,mn).ge.0.) then
                  uflux(i,j) = util3(i-1,j)*ubavg(i,j,mn)*scuy(i,j)
               else
                  uflux(i,j) = util3(i,j)*ubavg(i,j,mn)*scuy(i,j)
               endif
               uflux2(i,j) = 0.5*ubavg(i,j,mn)*scuy(i,j)*
     .                       (util3(i,j)+util3(i-1,j)) - uflux(i,j)
               flxloc(i,j) = uflux(i,j)
            enddo
         enddo
      enddo
c                          Second order u-flux for the non wet./dry. area
c                          ----------------------------------------------
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_nw(j)
            do i=max(1-margin,ifu_nw(j,l)),min(ii+margin,ilu_nw(j,l))
               uflux(i,j) = 0.5*ubavg(i,j,mn)*scuy(i,j)*
     .                       (util3(i,j)+util3(i-1,j)) 
               flxloc(i,j) = uflux(i,j)
            enddo
         enddo
      enddo
c                                                    Boundary conditions
c                                                    -------------------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu(j)
            i=max(1-nbdy,ifu(j,l)-1)
            if (iuopn(i,j).ne.0) then
c In the SHOM version, boundary values are stored in pbsavu, pbsavv
c Probably stored in pbavg(.,mn) in the standard version
c For SHOM  uflux(i,j)=max(0.,pbot(i,j)*onetai(i,j)+pbsavu(i,j))
                uflux(i,j)=max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,mn))
     .                     *ubavg(i,j,mn)*scuy(i,j)
                                                      else
                uflux(i,j)=0.
            endif
            flxloc(i,j) = uflux(i,j)
            uflux2(i,j) = 0.
            i=min(ii+nbdy,ilu(j,l)+1)
            if (iuopn(i,j).ne.0) then
c For SHOM  uflux(i,j)=max(0.,pbot(i-1,j)*onetai(i-1,j)+pbsavu(i,j))
                uflux(i,j)=max(0.,pbot(i-1,j)*onetai(i-1,j)+
     .                     pbavg(i,j,mn))
     .                     *ubavg(i,j,mn)*scuy(i,j)
                                                      else
                uflux(i,j)=0.
            endif
            flxloc(i,j) = uflux(i,j)
            uflux2(i,j) = 0.
         enddo
      enddo
c                              first order v-flux for the wet./dry. area
c                              -----------------------------------------
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
               if(vbavg(i,j,mn).ge.0.) then
                  vflux(i,j)=util3(i,j-1)*vbavg(i,j,mn)*scvx(i,j)
               else
                  vflux(i,j)=util3(i,j)*vbavg(i,j,mn)*scvx(i,j)
               endif
               vflux2(i,j) = 0.5*scvx(i,j)*vbavg(i,j,mn)*
     .                       (util3(i,j)+util3(i,j-1)) - vflux(i,j)
               flyloc(i,j) = vflux(i,j)
            enddo
         enddo
      enddo
c                          Second order v-flux for the non wet./dry. area
c                          ----------------------------------------------
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_nw(j)
            do i=max(1-margin,ifv_nw(j,l)),min(ii+margin,ilv_nw(j,l))
               vflux(i,j) = 0.5*scvx(i,j)*vbavg(i,j,mn)*
     .                       (util3(i,j)+util3(i,j-1))
               flyloc(i,j) = vflux(i,j)
            enddo
         enddo
      enddo
c                                                    Boundary conditions
c                                                    -------------------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
      do i=1-margin,ii+margin
         do l=1,jsv(i)
            j=max(1-nbdy,jfv(i,l)-1)
            if (ivopn(i,j).ne.0) then
ccc For SHOM   vflux(i,j)=max(0.,pbot(i,j)*onetai(i,j)+pbsavv(i,j))
               vflux(i,j)=max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,mn))
     .                     *vbavg(i,j,mn)*scvx(i,j)
                                                      else
                vflux(i,j)=0.
            endif
            flyloc(i,j) = vflux(i,j)
            vflux2(i,j)=0.
            j=min(jj+nbdy,jlv(i,l)+1)
            if (ivopn(i,j).ne.0) then
ccc For SHOM   vflux(i,j)=max(0.,pbot(i,j-1)*onetai(i,j-1)+pbsavv(i,j))
               vflux(i,j)=max(0.,pbot(i,j-1)*onetai(i,j-1) +
     .                      pbavg(i,j,mn))
     .                     *vbavg(i,j,mn)*scvx(i,j)
                                                      else
               vflux(i,j)=0.
            endif
            flyloc(i,j) = vflux(i,j)
            vflux2(i,j)=0.
         enddo
      enddo
c                                                           Integration 
c                                                           -----------
      margin = 0
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              pbavg(i,j,nl)=
     .        ((1.-wblpf)*pbavg(i,j,mn)+
     .             wblpf *pbavg(i,j,ml) ) 
     .                    - (uflux(i+1,j)-uflux(i,j) 
     .                    +  vflux(i,j+1)-vflux(i,j)
     .                    )*dlt*scp2i(i,j)*(1.+wblpf)
            enddo
         enddo
      enddo
!      call xctilr(pbavg(  1-nbdy,1-nbdy,nl ),1, 1, 3,3, halo_ps)
      call xctilr(pbavg(  1-nbdy,1-nbdy,nl ),1, 1, 6,6, halo_ps)
c
c			   End of the FCT scheme for the wet./dry. area
c                          --------------------------------------------
      margin = 3
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp_w(j)
            do i=max(1-margin,ifp_w(j,l)),min(ii+margin,ilp_w(j,l))
               util4(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
            enddo
            i = max(1-nbdy,ifp_w(j,l)-1)
            if (ip(i,j).eq.1)
     .         util4(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
            i = min(ii+nbdy,ilp_w(j,l)+1)
            if (ip(i,j).eq.1)
     .         util4(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
         enddo
      enddo
      do i = 1-margin, ii+margin
         do l = 1, jsp_w(i)
            j = max(1-nbdy,jfp_w(i,l)-1)
            if (ip(i,j).eq.1)
     .         util4(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
            j = min(jj+nbdy,jlp_w(i,l)+1)
            if (ip(i,j).eq.1)
     .         util4(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
         enddo
      enddo
c                                                     R+,R- calculations
c                                                     ------------------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib,ja,jb)
!$OMP&           SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp_w(j)
            do i=max(1-margin,ifp_w(j,l)),min(ii+margin,ilp_w(j,l))
              ia=i-1
              if (ip(ia,j).eq.0) ia=i
              ib=i+1
              if (ip(ib,j).eq.0) ib=i
              ja=j-1
              if (ip(i,ja).eq.0) ja=j
              jb=j+1
              if (ip(i,jb).eq.0) jb=j
              util1(i,j)=max(util4(i,j),util4(ia,j),util4(ib,j),
     .                                   util4(i,ja),util4(i,jb))
              util2(i,j)=max(0.,
     .                   min(util4(i,j),util4(ia,j),util4(ib,j),
     .                                   util4(i,ja),util4(i,jb)))
c
              util1(i,j)=(util1(i,j)-util4(i,j))
     .        /(((max(0.,uflux2(i,j))-min(0.,uflux2(i+1,j)))
     .          +(max(0.,vflux2(i,j))-min(0.,vflux2(i,j+1)))+epsil)
     .        *dlt*(1.+wblpf)*scp2i(i,j))
c
              util2(i,j)=(util2(i,j)-util4(i,j))
     .        /(((min(0.,uflux2(i,j))-max(0.,uflux2(i+1,j)))
     .          +(min(0.,vflux2(i,j))-max(0.,vflux2(i,j+1)))-epsil)
     .        *dlt*(1.+wblpf)*scp2i(i,j))
            enddo
         enddo
      enddo
c                                                           Final fluxes
c                                                           ------------
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
              if (uflux2(i,j).ge.0.) then
                uflux(i,j) = uflux2(i,j)
     .                       *min(1.,util1(i,j),util2(i-1,j))
              else
                uflux(i,j) = uflux2(i,j)
     .                       *min(1.,util2(i,j),util1(i-1,j))
              endif
              flxloc(i,j) = flxloc(i,j) + uflux(i,j)
           enddo
        enddo
c
        do l=1,isv_w(j)
           do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
              if (vflux2(i,j).ge.0.) then
                vflux(i,j) = vflux2(i,j)
     .                       *min(1.,util1(i,j),util2(i,j-1))
              else
                vflux(i,j) = vflux2(i,j)
     .                       *min(1.,util2(i,j),util1(i,j-1))
              endif
              flyloc(i,j) = flyloc(i,j) + vflux(i,j)
            enddo
         enddo
      enddo 
c                                       set fluxes to 0 at the interface
c                                       --------------------------------
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            i=max(1-nbdy,ifu_w(j,l)-1)
            uflux(i,j)=0.
            i=min(ii+nbdy,ilu_w(j,l)+1)
            uflux(i,j)=0.
         enddo
      enddo
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
      do i=1-margin,ii+margin
         do l=1,jsv_w(i)
            j=max(1-nbdy,jfv_w(i,l)-1)
            vflux(i,j)=0.
            j=min(jj+nbdy,jlv_w(i,l)+1)
            vflux(i,j)=0.
         enddo
      enddo
c                         Second order integration for the wet./dry.area
c                         ----------------------------------------------
      margin = 0
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp_w(j)
            do i=max(1-margin,ifp_w(j,l)),min(ii+margin,ilp_w(j,l))
               pbavg(i,j,nl)=pbavg(i,j,nl)
     .             - ((iu_w(i+1,j)*uflux(i+1,j)-iu_w(i,j)*uflux(i,j))+
     .                (iv_w(i,j+1)*vflux(i,j+1)-iv_w(i,j)*vflux(i,j)))
     .                *dlt*(1.+wblpf)*scp2i(i,j)
            enddo
         enddo
      enddo
c                                 Barotropic fluxes for mass conservation
c                                 ---------------------------------------
      do j=1-margin,jj+margin
         do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               uflxba(i,j) = uflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                          * flxloc(i,j)
             enddo
         enddo
      enddo
      do j=1-margin,jj+margin
         do l=1,isu(j)
            i=max(1-nbdy,ifu(j,l)-1)
            uflxba(i,j) = uflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                       * flxloc(i,j)
            i=min(ii+nbdy,ilu(j,l)+1)
            uflxba(i,j) = uflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                       * flxloc(i,j)
         enddo
      enddo
      do j=1-margin,jj+margin
         do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              vflxba(i,j) = vflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                          * flyloc(i,j)
            enddo
         enddo
      enddo
      do i=1-margin,ii+margin
         do l=1,jsv(i)
            j=max(1-nbdy,jfv(i,l)-1)
            vflxba(i,j) = vflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                        * flyloc(i,j)
            j=min(jj+nbdy,jlv(i,l)+1)
            vflxba(i,j) = vflxba(i,j) + coeflx(lll+1,icof)*(1.+wblpf)
     .                        * flyloc(i,j)
         enddo
      enddo
c
!      call xctilr(pbavg(  1-nbdy,1-nbdy,nl ),1, 1, 3,3, halo_ps)
      call xctilr(pbavg(  1-nbdy,1-nbdy,nl ),1, 1, 6,6, halo_ps)
c
c                        Momentum equations
c                        ==================
c                            u-equation
c                            - - - - - -
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               util2(i,j) = (-utotn(i,j)
     .               +thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)
     .               -((vbavg(i,j,mn)*depthv(i,j)
     .                 +vbavg(i,j+1,mn)*depthv(i,j+1))+
     .                 (vbavg(i-1,j,mn)*depthv(i-1,j)
     .                 +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     .                 (0.125*(pvtrop(i,j)+pvtrop(i,j+1))))

           enddo
         enddo
      enddo
c                                      Integration for non wet./wry.area
c                                      ---------------------------------
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_nw(j)
            do i=max(1-margin,ifu_nw(j,l)),min(ii+margin,ilu_nw(j,l))
               ubavg(i,j,nl) = ((1.-wblpf)*ubavg(i,j,mn)+
     .                wblpf*ubavg(i,j,ml))-(1.+wblpf)*dlt*util2(i,j)
           enddo
         enddo
      enddo
      call xctilr(ubavg(  1-nbdy,1-nbdy,nl ),1, 1, 6,6, halo_uv) !!Alex
c
      margin = 3
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isp_w(j)
            do i=max(1-margin,ifp_w(j,l)-1),min(ii+margin,ilp_w(j,l)+1)
               util3(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
            enddo
         enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-margin,ii+margin
         do l=1,jsp_w(i)
            j=max(1-margin,jfp_w(i,l)-1)
            util3(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
            j=min(jj+margin,jlp_w(i,l)+1)
            util3(i,j) = max(0.,pbot(i,j)*onetai(i,j)+pbavg(i,j,nl))
         enddo
      enddo
c                                                                Weights
c                                                                -------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)-1),min(ii+margin,ilu_w(j,l)+1)
               util1(i,j)=max(0.,
     .                        min(util3(i,j),util3(i-1,j),hcpu(i,j)))
               util2(i,j)=util1(i,j)*util2(i,j)
           enddo
         enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-margin,ii+margin
         do l=1,jsu_w(i)
            j = max(1-margin,jfu_w(i,l)-1)
            if (iu_w(i,j).ne.1) then !!Alex add iu_w from Remy
            util1(i,j)=max(0.,min(util3(i,j),util3(i-1,j),hcpu(i,j)))
            util2(i,j)=util1(i,j)*util2(i,j)
            endif
            j = min(jj+margin,jlu_w(i,l)+1)
            if (iu_w(i,j).ne.1) then !!Alex add iu_w from Remy
            util1(i,j)=max(0.,min(util3(i,j),util3(i-1,j),hcpu(i,j)))
            util2(i,j)=util1(i,j)*util2(i,j)
            endif
         enddo
      enddo
c
c
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
               ia=i-1
               if (iu(ia,j).ne.1) ia=i
               ib=i+1
               if (iu(ib,j).ne.1) ib=i
               ja=j-1
               if (iu(i,ja).ne.1) ja=j
               jb=j+1
               if (iu(i,jb).ne.1) jb=j
             util4(i,j)=(util2(i,j)+(hcpu(i,j)-util1(i,j))*
     .        (util2 (ia,j)+util2 (ib,j)+util2 (i,ja)+util2 (i,jb))/
     .        (util1(ia,j)+util1(ib,j)+util1(i,ja)+util1(i,jb)+
     .          epspu(i,j)))*invhcpu(i,j)
           enddo
         enddo
      enddo
c                                          Integration for wet./wry.area
c                                          -----------------------------
      util5(:,:) = 0. ! Probably not necessary - to be checked later
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
               util5(i,j)= ((1.-wblpf)*ubavg(i,j,mn)+
     .                wblpf*ubavg(i,j,ml))-(1.+wblpf)*dlt*util4(i,j)
           enddo
         enddo
      enddo
c                            v-equation
c                            - - - - - -
c
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               util2(i,j) = (-vtotn(i,j)
     .               +thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)
     .               +((ubavg(i,j,mn)*depthu(i,j)
     .                 +ubavg(i+1,j,mn)*depthu(i+1,j))+
     .                 (ubavg(i,j-1,mn)*depthu(i,j-1)
     .                 +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     .                 (0.125*(pvtrop(i,j)+pvtrop(i+1,j))))

           enddo
         enddo
      enddo
c                                      Integration for non wet./wry.area
c                                      ---------------------------------
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_nw(j)
            do i=max(1-margin,ifv_nw(j,l)),min(ii+margin,ilv_nw(j,l))
               vbavg(i,j,nl) = ((1.-wblpf)*vbavg(i,j,mn)+
     .                wblpf*vbavg(i,j,ml))-(1.+wblpf)*dlt*util2(i,j)
           enddo
         enddo
      enddo
      call xctilr(vbavg(  1-nbdy,1-nbdy,nl ),1, 1, 6,6, halo_vv) !!Alex
c                                                                Weights
c                                                                -------
      margin = 2
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)-1),min(ii+margin,ilv_w(j,l)+1)
               util1(i,j)=max(0.,
     .                        min(util3(i,j-1),util3(i,j),hcpv(i,j)))
               util2(i,j)=util1(i,j)*util2(i,j)
           enddo
         enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-margin,ii+margin !!Alex bug corrected from jj to ii 
         do l=1,jsv_w(i)
            j=max(1-margin,jfv_w(i,l)-1)
            if (iv_w(i,j).ne.1) then
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcpv(i,j)))
            util2(i,j)=util1(i,j)*util2(i,j)
            endif
            j=min(jj+margin,jlv_w(i,l)+1)
            if (iv_w(i,j).ne.1) then
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcpv(i,j)))
            util2(i,j)=util1(i,j)*util2(i,j)
            endif
         enddo
      enddo
c
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
               ia=i-1
               if (iv(ia,j).ne.1) ia=i
               ib=i+1
               if (iv(ib,j).ne.1) ib=i
               ja=j-1
               if (iv(i,ja).ne.1) ja=j
               jb=j+1
               if (iv(i,jb).ne.1) jb=j
             util4(i,j)=(util2(i,j)+(hcpv(i,j)-util1(i,j))*
     .        (util2 (ia,j)+util2 (ib,j)+util2 (i,ja)+util2 (i,jb))/
     .        (util1(ia,j)+util1(ib,j)+util1(i,ja)+util1(i,jb)+
     .          epspv(i,j)))*invhcpv(i,j)
           enddo
         enddo
      enddo
c                                          Integration for wet./wry.area
c                                          -----------------------------
      util6(:,:) = 0. ! Probably not necessary - to be checked later
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
               util6(i,j)= ((1.-wblpf)*vbavg(i,j,mn)+
     .                wblpf*vbavg(i,j,ml))-(1.+wblpf)*dlt*util4(i,j)
           enddo
         enddo
      enddo
c                     Velocity relaxation when dp is low
c                     ==================================
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
               util1(i,j)=min(hcu(i,j),util3(i,j),util3(i-1,j))
            enddo
            i = max(1-margin,ifu_w(j,l)-1)
            util1(i,j)=min(hcu(i,j),util3(i,j),util3(i-1,j))
            if (iu_nw(i,j).eq.1) util5(i,j)=ubavg(i,j,nl)
            i = min(ii+margin,ilu_w(j,l)+1)
            util1(i,j)=min(hcu(i,j),util3(i,j),util3(i-1,j))
            if (iu_nw(i,j).eq.1) util5(i,j)=ubavg(i,j,nl)
         enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-margin,ii+margin
         do l=1,jsu_w(i)
            j = max(1-margin,jfu_w(i,l)-1)
            util1(i,j)=max(0.,min(util3(i,j),util3(i-1,j),hcu(i,j)))
            if (iu_nw(i,j).eq.1) util5(i,j)=ubavg(i,j,nl)
            j = min(jj+margin,jlu_w(i,l)+1)
            util1(i,j)=max(0.,min(util3(i,j),util3(i-1,j),hcu(i,j)))
            if (iu_nw(i,j).eq.1) util5(i,j)=ubavg(i,j,nl)
         enddo
      enddo
c                                                              Nudging u
c                                                              ---------
      margin = 0
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu_w(j)
            do i=max(1-margin,ifu_w(j,l)),min(ii+margin,ilu_w(j,l))
               ia=i-1
               if (iu(ia,j).ne.1) ia=i
               ib=i+1
               if (iu(ib,j).ne.1) ib=i
               ja=j-1
               if (iu(i,ja).ne.1) ja=j
               jb=j+1
               if (iu(i,jb).ne.1) jb=j
             ubavg(i,j,nl) = ( util1(i,j)*util5(i,j)+
     .         (hcu(i,j)-util1(i,j))*(util1(ia,j)*util5(ia,j)+
     .                          util1(ib,j)*util5(ib,j)+
     .                          util1(i,ja)*util5(i,ja)+
     .                          util1(i,jb)*util5(i,jb))/
     .         (epsu(i,j)+util1(ia,j)+util1(ib,j)+
     .                util1(i,ja)+util1(i,jb)) )*invhcu(i,j)
           enddo
         enddo
      enddo


c
      margin = 1
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
               util1(i,j)=max(0.,
     .                        min(util3(i,j-1),util3(i,j),hcv(i,j)))
            enddo
            i = max(1-margin,ifv_w(j,l)-1)
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcv(i,j)))
            if (iv_nw(i,j).eq.1) util6(i,j)=vbavg(i,j,nl)
            i = min(ii+margin,ilv_w(j,l)+1)
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcv(i,j)))
            if (iv_nw(i,j).eq.1) util6(i,j)=vbavg(i,j,nl)
         enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-margin,ii+margin !!Alex bug corrected from jj to ii
         do l=1,jsv_w(i)
            j=max(1-margin,jfv_w(i,l)-1)
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcv(i,j)))
            if (iv_nw(i,j).eq.1) util6(i,j)=vbavg(i,j,nl)
            j=min(jj+margin,jlv_w(i,l)+1)
            util1(i,j)=max(0.,min(util3(i,j-1),util3(i,j),hcv(i,j)))
            if (iv_nw(i,j).eq.1) util6(i,j)=vbavg(i,j,nl)
         enddo
      enddo
c                                                              Nudging v
c                                                              ---------
      margin = 0
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isv_w(j)
            do i=max(1-margin,ifv_w(j,l)),min(ii+margin,ilv_w(j,l))
               ia=i-1
               if (iv(ia,j).ne.1) ia=i
               ib=i+1
               if (iv(ib,j).ne.1) ib=i
               ja=j-1
               if (iv(i,ja).ne.1) ja=j
               jb=j+1
               if (iv(i,jb).ne.1) jb=j
             vbavg(i,j,nl) = ( util1(i,j)*util6(i,j)+
     .         (hcv(i,j)-util1(i,j))*(util1(ia,j)*util6(ia,j)+
     .                          util1(ib,j)*util6(ib,j)+
     .                          util1(i,ja)*util6(i,ja)+
     .                          util1(i,jb)*util6(i,jb))/
     .         (epsv(i,j)+util1(ia,j)+util1(ib,j)+
     .                util1(i,ja)+util1(i,jb)) )*invhcv(i,j)
           enddo
         enddo
      enddo
c
      wblpf = wbaro
c
 840  continue  ! lll=1,lstep1
c
      do j=1-margin,jj+margin
         do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               uflxba(i,j) = dlt*uflxba(i,j)/delt1
             enddo
         enddo
      enddo
      do j=1-margin,jj+margin
         do l=1,isu(j)
            i=max(1-nbdy,ifu(j,l)-1)
            uflxba(i,j) = dlt*uflxba(i,j)/delt1
            i=min(ii+nbdy,ilu(j,l)+1)
            uflxba(i,j) = dlt*uflxba(i,j)/delt1
         enddo
      enddo
      do j=1-margin,jj+margin
         do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              vflxba(i,j) = dlt*vflxba(i,j)/delt1
            enddo
         enddo
      enddo
      do i=1-margin,ii+margin
         do l=1,jsv(i)
            j=max(1-nbdy,jfv(i,l)-1)
            vflxba(i,j) = dlt*vflxba(i,j)/delt1
            j=min(jj+nbdy,jlv(i,l)+1)
            vflxba(i,j) = dlt*vflxba(i,j)/delt1
         enddo
      enddo
c
c                                 Boundary conditions for 
c                                 ubavg, oneta in cnuity
c                                 ------------------------
c
                                
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubavg(i,j,n) = ubavg(i,j,nl)
         enddo
         enddo
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbavg(i,j,n) = vbavg(i,j,nl)
         enddo
         enddo
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbavg(i,j,n) = pbavg(i,j,nl)
         enddo
         enddo
      enddo
c                                                             new 1 + eta
c                                                             -----------
      call xctilr(pbavg(1-nbdy,1-nbdy,n),1, 1, nbdy,nbdy, halo_ps)

!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            oneta(i,j,n)  = onetai(i,j) + pbavg(i,j,n) /pbot(i,j)
      enddo
      enddo
      enddo
c                                                          Asselin filter
c                                                          --------------
      do j=1-margin,jj+margin
c
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubavg(i,j,m)      = ubavo(i,j)     +0.5*ra2fac*ubavg(i,j,n)
         enddo
         enddo
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbavg(i,j,m)      = vbavo(i,j)     +0.5*ra2fac*vbavg(i,j,n)
         enddo
         enddo
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbavg(i,j,m)      = pbavo(i,j)     +0.5*ra2fac*pbavg(i,j,n)
         enddo
         enddo
c
      enddo
c
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            onetao(i,j,m) = onetai(i,j) + pbavg(i,j,m) /pbot(i,j)
      enddo
      enddo
      enddo
      call xctilr(oneta(   1-nbdy,1-nbdy,n),1, 1, 2,2, halo_ps)
      call xctilr(onetao(  1-nbdy,1-nbdy,m),1, 1, 2,2, halo_ps)
c
!!Alex      if (btrmas) then
      if (btrmas .and. .not.btrmas) then
ccc      if (btrmas) then !Not Validated cf. Remy Baraille
c                                                 Final integration of dp
c                                                 -----------------------
c
      call xctilr(uflxba(  1-nbdy,1-nbdy),1, 1, 2,2, halo_uv)
      call xctilr(vflxba(  1-nbdy,1-nbdy),1, 1, 2,2, halo_vv)
      call xctilr(onetacnt(1-nbdy,1-nbdy),1,1, 2,2, halo_ps)
c
      margin = 1
      util1(:,:) = 0.
      util2(:,:) = 0.
      do k = 1, kk
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            dp(i,j,k,n)=onetacnt(i,j)*dp(i,j,k,n)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
         enddo
         enddo
         enddo
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               util1(i,j) = util1(i,j) + uflx(i,j,k)
            enddo
            i=ifu(j,l)-1
            if (i.ge.1-margin  ) then
               if (iuopn(i,j).ne.0) util1(i,j)=util1(i,j) + uflx(i,j,k)
            endif
            i=ilu(j,l)+1
            if (i.le.ii+margin ) then
               if (iuopn(i,j).ne.0) util1(i,j)=util1(i,j) + uflx(i,j,k)
            endif
         enddo
         enddo
!$OMP   END PARALLEL DO
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            util2(i,j) = util2(i,j) + vflx(i,j,k)
         enddo
         enddo
         enddo
!$OMP   END PARALLEL DO
c
         do i=1-margin,ii+margin
         do l=1,jsv(i)
            j=jfv(i,l)-1
            if (j.ge.1-margin  ) then
               if (ivopn(i,j).ne.0) util2(i,j)=util2(i,j) + vflx(i,j,k)
            endif
            j=jlv(i,l)+1
            if (j.le.jj+margin  ) then
               if (ivopn(i,j).ne.0) util2(i,j)=util2(i,j) + vflx(i,j,k)
            endif
         enddo
         enddo
      enddo
c      
      util3(:,:) = 0.
      util4(:,:) = 0.
      do k=1,kk
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
c
            do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               if (uflxba(i,j)-util1(i,j).ge.0.) then
                  if (dp(i,j,k,n).gt.0.) then
                     util3(i,j) = util3(i,j) + dp(i-1,j,k,n)
     .                                        /p(i-1,j,kk+1)
                  endif
               else
                  if (dp(i-1,j,k,n).gt.0.) then
                     util3(i,j) = util3(i,j) + dp(i,j,k,n)/p(i,j,kk+1)
                  endif
               endif
            enddo
            enddo
c
            do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               if (vflxba(i,j)-util2(i,j).ge.0.) then
                  if (dp(i,j,k,n).gt.0.) then
                     util4(i,j) = util4(i,j) + dp(i,j-1,k,n)
     .                                        /p(i,j-1,kk+1)
                  endif
               else
                  if (dp(i,j-1,k,n).gt.0.) then
                     util4(i,j) = util4(i,j) + dp(i,j,k,n)/p(i,j,kk+1)
                  endif
               endif
            enddo
            enddo
         enddo
!$OMP   END PARALLEL DO
c
      enddo
c      
      do k=1,kk
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,z1)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
c
            do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               if (uflxba(i,j)-util1(i,j).ge.0.) then
                  if (dp(i,j,k,n).gt.0.) then
                     z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)/util3(i,j)
                                         else
                     z1=0.
                  endif
               else
                  if (dp(i-1,j,k,n).gt.0.) then
                     z1=dp(i  ,j,k,n)/p(i  ,j,kk+1)/util3(i,j)
                                           else
                     z1=0.
                  endif
               endif
               uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
               uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
            enddo
            i=ifu(j,l)-1
            if (i.ge.1-margin  ) then
               if (iuopn(i,j).ne.0) then
                  z1=dp(i,j,k,n)/p(i,j,kk+1)
                  uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                  uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
               endif
            endif
            i=ilu(j,l)+1
            if (i.le.ii+margin ) then
               if (iuopn(i,j).ne.0) then
                  z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)
                  uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                  uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
               endif
            endif
            enddo
c
            do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               if (vflxba(i,j)-util2(i,j).ge.0.) then
                  if (dp(i,j,k,n).gt.0.) then
                     z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)/util4(i,j)
                                         else
                     z1=0.
                  endif
               else
                  if (dp(i,j-1,k,n).gt.0.) then
                     z1=dp(i,j  ,k,n)/p(i,j  ,kk+1)/util4(i,j)
                                         else
                     z1=0.
                  endif
               endif
               vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
               vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
            enddo
            enddo
         enddo
!$OMP   END PARALLEL DO
c
        do i=1-margin,ii+margin
        do l=1,jsv(i)
           j=jfv(i,l)-1
           if (j.ge.1-margin  ) then
              if (ivopn(i,j).ne.0) then
                 z1=dp(i,j,k,n)/p(i,j,kk+1)
                 vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                 vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
              endif
           endif
           j=jlv(i,l)+1
           if (j.le.jj+margin  ) then
              if (ivopn(i,j).ne.0) then
                 z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)
                 vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                 vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
              endif
           endif
         enddo
         enddo
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,dpmin)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            dp(i,j,k,n)=dp(i,j,k,n)-
     &                  ((uflux(i+1,j)-uflux(i,j))+
     &                  (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
ccc            p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            p(i,j,k+1)=p(i,j,k)+max(0.,dp(i,j,k,n))
         enddo
         enddo
         enddo
c
       enddo
c
!$OMP PARALLEL DO PRIVATE(j,l,k,i)
!$OMP&         SCHEDULE(STATIC,jblk)
       do j=1-margin,jj+margin
       do l=1,isp(j)
          do k=1,kk
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              if (p(i,j,kk+1) .gt. 0.) then
ccc              dp(i,j,k,n)=dp(i,j,k,n)*pbot(i,j)/p(i,j,kk+1)
              dp(i,j,k,n)=max(0.,dp(i,j,k,n))*pbot(i,j)/p(i,j,kk+1)
              else
              dp(i,j,k,n)=0.
              dp(i,j,1,n)=pbot(i,j)
              endif
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              if (isopyc .and. k.eq.1) then
                dpmixl(i,j,n)=dp(i,j,k,n)
              endif
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      call xctilr(dp(1-nbdy,1-nbdy,1,n),1,kk, 2,2, halo_ps)
c
      endif  
c
      end subroutine barotp_dry_and_nodry
      
      subroutine check_field_arctic_ok(field,itypegrid,itypefield,car6)
      use mod_xc        ! HYCOM communication interface
      use mod_pipe      ! HYCOM debugging interface
      use mod_cb_arrays ! HYCOM saved arrays
      
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .                field, field1
      real work0(itdm,jtdm)
      integer i,ind, itypegrid, itypefield
      character*6 car6
      real xsign
c     itypegrid  = 1 (grille P),  = 2 (grille U), = 3 (grille V)
c     itypefield = 1 (champ P,T,S),  = 2 (champ U), = 3 (champ V)
c
      if (itypegrid.eq.1) field1(:,:) = ip(:,:) * field(:,:)
      if (itypegrid.eq.2) field1(:,:) = iu(:,:) * field(:,:)
      if (itypegrid.eq.3) field1(:,:) = iv(:,:) * field(:,:)
      xsign = -1.
      if (itypefield.eq.1) xsign = 1.
c
      call xcaget(work0, field1, 1)
      if (mnproc.eq.1) then
         ind = 0
         do i = 1, itdm
            if (itypegrid.eq.1) then
               diff = work0(i,jtdm-1) - xsign*work0(itdm+1-i,jtdm)
               if (abs(diff).ge.1.e-10) then
                  write(6,*)'Point h i,j=',i,jtdm-1,'and i,j=',
     .                        itdm+1-i,jtdm,'not symetric',
     .                       ' val1=',work0(i,jtdm-1),
     .                       ' val2=',work0(itdm+1-i,jtdm),
     .                       ' nstep=',nstep,
     .                       ' field ',car6
                  ind = ind + 1
               endif
            endif
            if (itypegrid.eq.2) then
               if (i.eq.1) then
                   diff = 0.
                           else
                   diff = work0(i,jtdm-1) - xsign*work0(itdm+2-i,jtdm)
               endif
               if (abs(diff).ge.1.e-10) then
                  write(6,*)'Point u i,j=',i,jtdm-1,'and i,j=',
     .                        itdm+2-i,jtdm,'not symetric',
     .                       ' val1=',work0(i,jtdm-1),
     .                       ' val2=',work0(itdm+2-i,jtdm),
     .                       ' nstep=',nstep,
     .                       ' field ',car6
                  ind = ind + 1
               endif
            endif
            if (itypegrid.eq.3) then
               diff = work0(i,jtdm) - xsign*work0(itdm+1-i,jtdm)
               if (abs(diff).ge.1.e-10)  then
                  write(6,*)'Point v i,j=',i,jtdm,'and i,j=',
     .                        itdm+1-i,jtdm,'not symetric',
     .                       ' val1=',work0(i,jtdm),
     .                       ' val2=',work0(itdm+1-i,jtdm),
     .                       ' nstep=',nstep,
     .                       ' field ',car6
                  ind = ind + 1
               endif
            endif
         enddo
         if (ind.eq.0) write(6,*)'field ',car6,
     .                 ' at nstep=',nstep,' OK'
      endif
c
      return
      end
c
c=======================================================================
c
      subroutine make_field_arctic_ok(field,itypegrid,itypefield)
      use mod_xc        ! HYCOM communication interface
      use mod_pipe      ! HYCOM debugging interface
      use mod_cb_arrays ! HYCOM saved arrays
     
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .                field, field1
      real work0(itdm,jtdm)
      integer i,ind, itypegrid, itypefield
      real xsign
c     itypegrid  = 1 (grille P),  = 2 (grille U), = 3 (grille V)
c     itypefield = 1 (champ P,T,S),  = 2 (champ U), = 3 (champ V)
c
      if (itypegrid.eq.1) field1(:,:) = ip(:,:) * field(:,:)
      if (itypegrid.eq.2) field1(:,:) = iu(:,:) * field(:,:)
      if (itypegrid.eq.3) field1(:,:) = iv(:,:) * field(:,:)
      xsign = -1.
      if (itypefield.eq.1) xsign = 1.
c
      call xcaget(work0, field1, 1)
      if (mnproc.eq.1) then
         ind = 0
         do i = 1, itdm
            if (itypegrid.eq.1) then
               work0(i,jtdm-1) = xsign*work0(itdm+1-i,jtdm)
            endif
            if (itypegrid.eq.2) then
               if (i.ne.1) then
                  work0(i,jtdm-1) = xsign*work0(itdm+2-i,jtdm)
               endif
            endif
            if (itypegrid.eq.3) then
               work0(i,jtdm) = xsign*work0(itdm+1-i,jtdm)
            endif
            if (itypegrid.eq.4) then
               if (i.ne.1) then
                  work0(i,jtdm-1) = xsign*work0(itdm+2-i,jtdm-1)
               endif
            endif
         enddo
      endif
c
      call xcaput(work1,field,1)
c
      return
      end
c=======================================================================
c
      subroutine check_field_cstlon_ok(field,itypegrid,car6)
      use mod_xc        ! HYCOM communication interface
      use mod_pipe      ! HYCOM debugging interface
      use mod_cb_arrays ! HYCOM saved arrays

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .                field, field1
      real work0(itdm,jtdm)
      integer i,ind, itypegrid
      character*6 car6
      real x1
c     itypegrid  = 1 (grille P),  = 2 (grille U), = 3 (grille V)
c
      if (itypegrid.eq.1) field1(:,:) = ip(:,:) * field(:,:)
      if (itypegrid.eq.2) field1(:,:) = iu(:,:) * field(:,:)
      if (itypegrid.eq.3) field1(:,:) = iv(:,:) * field(:,:)
c
      call xcaget(work0, field1, 1)
      if (mnproc.eq.1) then
         ind = 0
         do j = 1, jtdm
            x1 = 0.
            do i = 1, itdm-1
               x1 = x1 + abs(work0(i+1,j)-work0(i,j))
               if (x1.ge.1.e-14) then
                  ind = ind + 1
                  write(6,*)'car6=',car6,' j=',j,'i=',i,
     .                      'work1(i+1,j)=',work0(i+1,j),
     .                      'work1(i,j)=',work0(i,j)
               endif
            enddo
         enddo
         if (ind.eq.0) write(6,*)'field ',car6,
     .                 ' at nstep=',nstep,' OK'
      endif
c
      return
      end

      
c> Revision history:
c>
c> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
c>             (loops 33,35)
c> Mar. 1995 - changed order of loop nesting in loop 842
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
c> Jan. 2000 - added latbdp for lateral boundary ports
c> Aug. 2001 - two barotropic time steps per loop, for halo efficiency
c> Nov. 2006 - added lbflag==3 (latbdf) and thref_bt (mod_tides)
c> Nov. 2006 - removed thref_bt (and mod_tides)
c> Apr. 2007 - added btrlfr: leapfrog time step; see also momtum
c> Apr. 2010 - bugfixes for 1st time step and 1st miner time step
c> Apr  2011 - added    Robert-Asselin filtering for btrlfr
c> Aug  2011 - reworked Robert-Asselin filtering for btrlfr
c> Mar. 2012 - added latbdtf for nesting with Flather b.c.'s.
c> Jan. 2013 - added tidal drag tensor
c> June 2013 - added   lbflag==6 for latbdtc
c> Apr. 2014 - replace ip with ipa for mass sums
c> May  2014 - use land/sea masks (e.g. ip) to skip land
c> May  2014 - removed lbflag==6 for latbdtc
c> Apr. 2015 - added atmospheric pressure forcing
c> Apr. 2015 - added tidal body forcing
c> Jan. 2016 - added btrmas option
