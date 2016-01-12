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
      subroutine rasselin(m,n,filt,filu)
      use mod_xc        ! HYCOM communication interface
      use mod_cb_arrays ! HYCOM saved arrays
      implicit none
c
      include 'stmt_fns.h'
c
      real,    parameter :: onemu=9806.e-12 !very small layer thickness
!!Alex      real,    parameter :: onezm=9806.e-20 ! insignificant thickness
      real,    parameter :: onezm=9806.e-18 ! insignificant thickness
c
      integer i,j,k,l,ktr,margin
      integer m,n
      logical filt, filu
      real dpold, dpmid, dpnew, dpmidn, q, qdpmidn
      real dpsold, dpsmid, dpsnew
      real sminny(jdm),smaxxy(jdm),smin,sminn,smaxx
      real xmin(kdm),xmax(kdm)
      logical latemp,lath3d,ldtemp,ldth3d

c
c --- -----------------------------------------------------
c --- Asselin Filter
c --- -----------------------------------------------------
c --- on entry:
c==============
c --- salno,dpo,uo(:,:,:,n) = time step t-1
c --- salno,dpo,uo(:,:,:,m) = time step t without RA
c --- saln,dp,u(:,:,:,n) = time step t+1
c --- saln,dp,u(:,:,:,m) = time step t with possibly intermediate RA
c
c --- onetao(:,:,n) = time step t-1
c --- onetao(:,:,m) = time step t WITH RA!!!!!
c --- oneta(:,:,n) = time step t+1
c --- oneta(:,:,m) = time step t WITHOUT RA!!!!
c
c --- on exit:
c==============
c --- oneta(:,:,m) = time step t with RA
c --- saln,dp,u(:,:,:,m) = time step t with RA
c --- -----------------------------------------------------
c
      if (filt) then !!Alex filter var on T grid
      margin = 0  !after advem
c
      if     ( btrmas ) then 
c
c ---   Robert-Asselin time filter of scalar fields
c ---   Note that this is smoothing dp * oneta *scalar,
c ---    but the filter is not conservative across 3 time levels.
c
c ---   rhs: temp.n, th3d.n, saln.n, dpo, dp.m, dp.n, sold, told
c ---   lhs: temp.n, th3d.n, dp.m, saln.m, temp.m, th3d.m
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,ktr,dpsold,dpsmid,dpsnew,q,       !NOCSD
!$OMP&                      dpold,dpmid,dpnew,dpmidn,qdpmidn,smin)  !NOCSD
!$OMP&           SCHEDULE(STATIC,jblk) !NOCSD
        do k = 1, kk
        latemp =  k.le.nhybrd .and. advflg.eq.0       ! advect temp
        lath3d = (k.le.nhybrd .and. advflg.eq.1) .or.
     &           (k.eq.1      .and. isopyc     ) ! advect th3d
        
        do j=1-margin,jj+margin
        sminny(j)= 999.  !simplifies OpenMP parallelization
        smaxxy(j)=-999.  !simplifies OpenMP parallelization
        do l=1,isp(j)
!DIR$     PREFERVECTOR
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            if     (dp(i,j,k,n).gt.onemm) then
              sminny(j)=min(sminny(j),saln(i,j,k,n))
              smaxxy(j)=max(smaxxy(j),saln(i,j,k,n))
            endif
c
            dpmidn = dp( i,j,k,m)*onetao(i,j,m)
            if     (dpmidn.gt.onezm) then !effectively non-zero
            dpold  = dpo(i,j,k,n)*onetao(i,j,n) !t-1
            dpmid  = dpo(i,j,k,m)*oneta( i,j,m) !t
            dpnew  = dp( i,j,k,n)*oneta( i,j,n) !t+1
c ---       redo the Robert-Asselin thickness filter
            q           = 0.5*ra2fac*(dpold+dpnew-2.0*dpmid)
            dpmidn      = dpmid + q
            dp(i,j,k,m) = dpmidn / onetao(i,j,m)
            qdpmidn     = 1.0/max(dpmidn,onemu)
c
            dpsold = dpold*salno(i,j,k,n)
            dpsmid = dpmid*salno(i,j,k,m)
            dpsnew = dpnew*saln(i,j,k,n)
            q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
            saln(i,j,k,m)= (dpsmid + q) * qdpmidn
            if     (latemp) then
              dpsold = dpold*tempo(i,j,k,n)
              dpsmid = dpmid*tempo(i,j,k,m)
              dpsnew = dpnew*temp(i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              temp(i,j,k,m) = (dpsmid + q) * qdpmidn
c ---         update dependent thermodynamic variable
              th3d(i,j,k,m) = sig(temp(i,j,k,m),saln(i,j,k,m))-thbase
            elseif (lath3d) then
              dpsold = dpold*th3do(i,j,k,n)
              dpsmid = dpmid*th3do(i,j,k,m)
              dpsnew = dpnew*th3d(i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              th3d(i,j,k,m) = (dpsmid + q) * qdpmidn
c ---         update dependent thermodynamic variable
              temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase,
     &                               saln(i,j,k,m))
            else   ! exactly isopycnal layer
              th3d(i,j,k,m) = theta(i,j,k)
c ---         update dependent thermodynamic variable
              temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase,
     &                               saln(i,j,k,m))
            endif
            do ktr= 1,ntracr
c ---         non-negative version that exactly conserves constant tracers
              smin   = min(   tracero(i,j,k,n,ktr),
     &                        tracer (i,j,k,m,ktr),
     &                        tracer (i,j,k,n,ktr) )
              dpsold = dpold*(tracero(i,j,k,n,ktr) - smin)  !>=0
              dpsmid = dpmid*(tracero(i,j,k,m,ktr) - smin)  !>=0
              dpsnew = dpnew*(tracer (i,j,k,n,ktr) - smin)  !>=0
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              tracer(i,j,k,m,ktr) = smin + (dpsmid + q) * qdpmidn
            enddo !ktr
            if (mxlmy) then
              dpsold = dpold*q2o(i,j,k,n)
              dpsmid = dpmid*q2o(i,j,k,m)
              dpsnew = dpnew* q2(i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              q2(i,j,k,m) = (dpsmid + q) * qdpmidn
              dpsold = dpold*q2lo(i,j,k,n)
              dpsmid = dpmid*q2lo(i,j,k,m)
              dpsnew = dpnew* q2l(i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              q2l(i,j,k,m) = (dpsmid + q) * qdpmidn
            endif !mxlmy
            endif !effectively non-zero
          enddo !i
        enddo !l
        enddo !j
        xmin(k) = minval(sminny(1:jj))
        xmax(k) = maxval(smaxxy(1:jj))
        enddo !k
!$OMP   END PARALLEL DO !NOCSD
c
      else  !.not.btrmas
c
c ---   Robert-Asselin time filter of scalar fields
c ---   Note that this is smoothing dp *scalar,
c ---    but the filter is not conservative across 3 time levels.
c
c ---   rhs: temp.n, th3d.n, saln.n, dpo, dp.m, dp.n, sold, told
c ---   lhs: temp.n, th3d.n, dp.m, saln.m, temp.m, th3d.m
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,ktr,dpsold,dpsmid,dpsnew,q,  !NOCSD
!$OMP&                      dpold,dpmid,dpnew,qdpmidn)         !NOCSD
!$OMP&           SCHEDULE(STATIC,jblk) !NOCSD
        do k = 1, kk
        latemp =  k.le.nhybrd .and. advflg.eq.0       ! advect temp
        lath3d = (k.le.nhybrd .and. advflg.eq.1) .or.
     &           (k.eq.1      .and. isopyc     ) ! advect th3d

        do j=1-margin,jj+margin
        sminny(j)= 999.  !simplifies OpenMP parallelization
        smaxxy(j)=-999.  !simplifies OpenMP parallelization
!DIR$     PREFERVECTOR
          do i=1-margin,ii+margin
            if (SEA_P) then                
            if     (dp(i,j,k,n).gt.onemm) then
              sminny(j)=min(sminny(j),saln(i,j,k,n))
              smaxxy(j)=max(smaxxy(j),saln(i,j,k,n))
            endif
c
            if (dp( i,j,k,m).gt.onezm) then !effectively non-zero
              dpold   =     dpo(i,j,k,n)
              dpmid   =     dpo(i,j,k,m)
              dpnew   =     dp( i,j,k,n)
              qdpmidn = 1.0/dp( i,j,k,m)

            dpsold = dpold*salno(i,j,k,n)
            dpsmid = dpmid*saln (i,j,k,m)
            dpsnew = dpnew*saln (i,j,k,n)
            q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
            saln(i,j,k,m)= (dpsmid + q) * qdpmidn
           if     (latemp) then
              dpsold = dpold*tempo(i,j,k,n)
              dpsmid = dpmid*temp (i,j,k,m)
              dpsnew = dpnew*temp (i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              temp(i,j,k,m) = (dpsmid + q) * qdpmidn
c ---         update dependent thermodynamic variable
              th3d(i,j,k,m) = sig(temp(i,j,k,m),saln(i,j,k,m))-thbase
           elseif (lath3d) then
              dpsold = dpold*th3do(i,j,k,n)
              dpsmid = dpmid*th3d(i,j,k,m)
              dpsnew = dpnew*th3d(i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              th3d(i,j,k,m) = (dpsmid + q) * qdpmidn
c ---         update dependent thermodynamic variable
              temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase,
     &                               saln(i,j,k,m))
            else   ! exactly isopycnal layer
              th3d(i,j,k,m) = theta(i,j,k)
c ---         update dependent thermodynamic variable
              temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase,
     &                               saln(i,j,k,m))
            endif
            do ktr= 1,ntracr
              dpsold = dpold*tracero(i,j,k,n,ktr)
              dpsmid = dpmid*tracer(i,j,k,m,ktr)
              dpsnew = dpnew*tracer (i,j,k,n,ktr)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              tracer(i,j,k,m,ktr) = (dpsmid + q) * qdpmidn
            enddo !ktr
            if (mxlmy) then
              dpsold = dpold*q2o(i,j,k,n)
              dpsmid = dpmid*q2o(i,j,k,m)
              dpsnew = dpnew*q2 (i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              q2(i,j,k,m) = (dpsmid + q) * qdpmidn
              dpsold = dpold*q2lo(i,j,k,n)
              dpsmid = dpmid*q2lo(i,j,k,m)
              dpsnew = dpnew*q2l (i,j,k,n)
              q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
              q2l(i,j,k,m) = (dpsmid + q) * qdpmidn
            endif !mxlmy
            endif !effectively non-zero
          endif !ip
        enddo !i
      enddo !j
      xmin(k) = minval(sminny(1:jj))
      xmax(k) = maxval(smaxxy(1:jj))
      enddo !k
!$OMP END PARALLEL DO !NOCSD
c
      endif !btrmas:else

c --- check for negative scalar fields.
c
 101  format (i9,' i,j,k =',2i5,i3,a,2f8.2)
c
      if     (mod(nstep,3).eq.0 .or. diagno) then
        call xcminr(xmin(1:kk))
        call xcmaxr(xmax(1:kk))
c
        do k= 1,kk
          sminn=xmin(k)
          smaxx=xmax(k)
c
          if (sminn.lt.0.0) then
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  if (saln(i,j,k,n).eq.sminn) then
                    write (lp,101) nstep,i+i0,j+j0,k,
     &                ' neg. saln after advem call ',
     &                saln(i,j,k,n)
                  endif
                enddo !i
              enddo !l
            enddo !j
            call xcsync(flush_lp)
          endif
c
          if (diagno) then
            if     (mnproc.eq.1) then
            write (lp,'(i9,i3, a,2f7.3, a,1pe9.2,a)')
     &        nstep,k,
     &        ' min/max of s after advection:',sminn,smaxx,
     &        '   (range:',smaxx-sminn,')'
            call flush(lp)
            endif
          endif
        enddo !k
      endif !every 3 time steps or diagno

      endif                     ! filt
      
            if (filu) then !!Alex filter var on U and V grid
         margin = 0
!$OMP   PARALLEL DO PRIVATE(j,l,i,ktr,dpsold,dpsmid,dpsnew,q,  !NOCSD
!$OMP&                      dpold,dpmid,dpnew,dpmidn,qdpmidn)  !NOCSD
!$OMP&           SCHEDULE(STATIC,jblk) !NOCSD
         do k = 1, kk
         do j = 1-margin, jj+margin
         do l = 1, isu(j)
!DIR$     PREFERVECTOR
         do i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
              dpold = dpuo(i,j,k,n)
              dpmid = dpuo(i,j,k,m)
              dpnew = dpu(i,j,k,n)
              qdpmidn = 1.0/(0.5*ra2fac*(dpold+dpnew)+onemm
     .                       +(1.-ra2fac)*dpmid)
c
              dpsold = dpold*uo(i,j,k,n)
              dpsnew = dpnew*u(i,j,k,n)
              u(i,j,k,m) = qdpmidn*(0.5*ra2fac*(dpsold+dpsnew) +
     .                      ((1.-ra2fac)*dpmid+onemm)*uo(i,j,k,m))
c
           enddo !i
           enddo !l
           enddo !j
           enddo !k
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,ktr,dpsold,dpsmid,dpsnew,q,  !NOCSD
!$OMP&                      dpold,dpmid,dpnew,dpmidn,qdpmidn)  !NOCSD
!$OMP&           SCHEDULE(STATIC,jblk) !NOCSD
         do k = 1, kk
         do j = 1-margin, jj+margin
         do l = 1, isv(j)
!DIR$     PREFERVECTOR
         do i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
              dpold = dpvo(i,j,k,n)
              dpmid = dpvo(i,j,k,m)
              dpnew = dpv(i,j,k,n)
              qdpmidn = 1.0/(0.5*ra2fac*(dpold+dpnew)+onemm
     .                       +(1.-ra2fac)*dpmid)
c
              dpsold = dpold*vo(i,j,k,n)
              dpsnew = dpnew*v(i,j,k,n)
              v(i,j,k,m) = qdpmidn*(0.5*ra2fac*(dpsold+dpsnew) +
     .                      ((1.-ra2fac)*dpmid+onemm)*vo(i,j,k,m))
c
           enddo !i
           enddo !l
           enddo !j
           enddo !k
      endif ! filu
c
      if (filt) then !!Alex filter var on T grid
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            oneta(i,j,m) = onetao(i,j,m)
         enddo
         enddo
         enddo
      endif ! filt
      return
      end
