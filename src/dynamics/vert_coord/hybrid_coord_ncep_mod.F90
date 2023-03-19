module hybrid_coord_ncep_mod

  use flogger
  use const_mod
  use namelist_mod
  use process_mod

  implicit none

  private

  public hybrid_coord_ncep

contains

  subroutine hybrid_coord_ncep(p0, ptop, hyai, hybi)

    real(r8), intent(in) :: p0
    real(r8), intent(in) :: ptop
    real(r8), intent(out) :: hyai(:)
    real(r8), intent(out) :: hybi(:)

    real(r8) ps_min

    ! psig  =99400.0_r8   , &
    ! ppre  =7000.0_r8    , &
    ! dpbot =500.0_r8     , &
    ! dpsig =1200.0_r8    , &
    ! dppre =18000.0_r8   , &
    ! dpupp =60.0_r8      , &
    ! dptop =60.0_r8      , &

    call vcoord_gen(        &
      levs  =size(hyai)-1 , &
      lupp  =size(hyai)-1 , &
      pbot  =p0           , &
      psig  =hybrid_coord_ncep_psig, &
      ppre  =hybrid_coord_ncep_ppre, &
      pupp  =ptop         , &
      ptop  =ptop         , &
      dpbot =hybrid_coord_ncep_dpbot, &
      dpsig =hybrid_coord_ncep_dpsig, &
      dppre =hybrid_coord_ncep_dppre, &
      dpupp =hybrid_coord_ncep_dptop, &
      dptop =hybrid_coord_ncep_dptop, &
      pmin  =ps_min       , &
      ak    =hyai         , &
      bk    =hybi         )

    hyai = hyai / p0

    ! Revert from bottom-top to top-bottom.
    hyai = hyai(size(hyai)::-1)
    hybi = hybi(size(hybi)::-1)

  end subroutine hybrid_coord_ncep

  !> @file
  !! Lower and upper triangular decomposition.
  !! @author Mark Iredell @date 2008-08-01
  
  !> This subprogram decomposes a matrix into a product of
  !! lower and upper triangular matrices.
  !!
  !! @param[inout] a - input: a real(np,np) matrix (will be overwritten) output:
  !!                 - output real(np,np) LU-decomposed matrix
  !!                   (U is upper part of A, including diagonal;
  !!                   L is lower part of A, with 1 as diagonal;
  !!                   L*U equals original A after permuting)
  !! @param[in] n integer order of matrix
  !! @param[in] np integer dimension of matrix
  !! @param[out] indx integer(n) pivot indices
  !!              (original A rows are permuted in order i with indx(i))
  !! @param[out] d real determinant permutation (1 or -1, or 0 if singular)
  !!              (determinant is output diagonal product times d)
  !! @author Mark Iredell @date 2008-08-01

  subroutine ludcmp(a, n, np, indx, d)

    integer, intent(in):: n,np
    real(r8), intent(inout):: a(np,np)
    integer,intent(out):: indx(n)
    real(r8), intent(out):: d

    integer i, j, k, imax
    real(r8) aamax, sum, dum
    real(r8) vv(n)

    d = 1
    do i = 1, n
      aamax = 0
      do j = 1, n
        if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
      end do
      if (aamax == 0) then
        d = 0
        return
      end if
      vv(i) = 1 / aamax
    end do
    do j = 1, n
      do i = 1, j - 1
        sum = a(i,j)
        do k = 1, i - 1
          sum = sum - a(i,k) * a(k,j)
        end do
        a(i,j) = sum
      end do
      aamax = 0.
      do i = j, n
        sum = a(i,j)
        do k = 1, j - 1
          sum = sum - a(i,k) * a(k,j)
        end do
        a(i,j) = sum
        dum = vv(i) * abs(sum)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        end if
      end do
      if (j /= imax) then
        do k = 1, n
          dum = a(imax,k)
          a(imax,k) = a(j,k)
          a(j,k) = dum
        end do
        d = -d
        vv(imax) = vv(j)
      end if
      indx(j) = imax
      if (a(j,j) == 0) then
        d = 0
        return
      end if
      if (j /= n) then
        dum = 1 / a(j,j)
        do i = j + 1, n
          a(i,j) = a(i,j) * dum
        end do
      end if
    end do

  end subroutine ludcmp

  !> Lower and upper triangular back substitution
  !! @author Iredell @date 2008-08-01
  !!
  !! This subprogram back substitutes to solve decomposed
  !! lower and upper triangular matrices as outputted by ludcmp.
  !!
  !! @param[in] a real(np,np) LU-decomposed matrix (from ludcmp)
  !! @param[in] n integer order of matrix
  !! @param[in] np integer dimension of matrix
  !! @param[in] indx integer(n) pivot indices (from ludcmp)
  !! @param[inout] b  - input real(n) rhs vector of linear problem (will be overwritten)
  !!                  - output real(n) solution of linear problem 
  !!                    (original A times output B equals original B)
  !! @author Mark Iredell @date 2008-08-01

  subroutine lubksb(a,n,np,indx,b)

    integer, intent(in) :: n,np
    real(r8), intent(in) :: a(np,np)
    integer, intent(in) :: indx(n)
    real(r8), intent(inout) :: b(n)

    integer i, j, ii, ll
    real(r8) sum
  
    ii = 0
    do i = 1, n
      ll = indx(i)
      sum = b(ll)
      b(ll) = b(i)
      if (ii /= 0) then
        do j = ii, i-1
          sum = sum - a(i,j) * b(j)
        end do
      else if (sum /= 0) then
        ii = i
      end if
      b(i) = sum
    end do
    do i = n, 1, -1
      sum = b(i)
      do j = i + 1, n
        sum = sum - a(i,j) * b(j)
      end do
      b(i) = sum / a(i,i)
    end do

  end subroutine lubksb

  !> @file
  !! @brief Generates hybrid coordinate interface profiles.
  !! @author Mark Iredell @date 2008-08-01
  
  !> This subprogram generates hybrid coordinate interface profiles
  !! from a few given parameters. The hybrid coordinate is intended to start
  !! out at the bottom in pure sigma and end up at the top in pure pressure,
  !! with a smooth transition in between. The pressure thickness is close to
  !! quadratic in pressure, with maximum thicknesses in the middle of the domain.
  !! The coordinate pressure will have continuous second derivatives in level.
  !!
  !! The hybrid coordinate is returned in terms of vectors AK and BK, where
  !! the interface pressure is defined as A+B*ps where ps is surface pressure
  !! Thus A=0 in regions of pure sigma and B=0 in regions of pure sigma.
  !! At the bottom, A(0)=0 and B(0)=1 so that surface pressure is the bottom
  !! boundary condition, while at the top, A(levs)=ptop and B(levs)=0 so that
  !! the constant top pressure (which can be zero) is the top boundary condition.
  !!
  !! The procedure for the calculation is described in the remarks section below.
  !!
  !! @param[in] levs     integer number of levels 
  !! @param[in] lupp     integer number of levels below pupp
  !! @param[in] pbot     real nominal surface pressure (Pa) 
  !! @param[in] psig     real nominal pressure where coordinate changes
  !!              from pure sigma (Pa) 
  !! @param[in] ppre     real nominal pressure where coordinate changes
  !!              to pure pressure (Pa) 
  !! @param[in] pupp     real nominal pressure where coordinate changes
  !!              to upper atmospheric profile (Pa)
  !! @param[in] ptop     real pressure at top (Pa)
  !! @param[in] dpbot    real coordinate thickness at bottom (Pa) 
  !! @param[in] dpsig    real thickness of zone within which coordinate changes
  !!              to pure sigma (Pa) 
  !! @param[in] dppre    real thickness of zone within which coordinate changes
  !!              to pure pressure (Pa) 
  !! @param[in] dpupp    real coordinate thickness at pupp (Pa) 
  !! @param[in] dptop    real coordinate thickness at top (Pa) 
  !! @param[out] pmin     real minimum surface pressure (Pa)
  !! @param[out] ak       real(0:levs) a coordinate values, bottom to top (Pa)
  !! @param[out] bk       real(0:levs) b coordinate values, bottom to top ()
  !!
  !! Subprograms called:
  !! - ludcmp   lower and upper triangular decomposition
  !! - lubksb   lower and upper triangular back substitution
  !!
  !! <pre>
  !!   Example: Create the operational GFS 64-level hybrid coordinate.
  !!     real(8) pmin,ak(0:64),bk(0:64)
  !!     call vcoord_gen(64,64,100000.,99400.,7000.,0.,0.,500.,1200.,18000.,60.,60.,&
  !!                  pmin,ak,bk)
  !!     print '(2i6)',2,64
  !!     print '(f12.3,f12.8)',(ak(k),bk(k),k=0,64)
  !!     end
  !!
  !!   Graphical description of parameters and zones:
  !!     ptop  ---  -----  ----------------------
  !!           ...  dptop
  !!           ---         zone U (upper atmos)
  !!           ...
  !!     pupp  ---  -----  ----------------------
  !!           ...  dpupp
  !!           ---  -----
  !!           ...         zone P (pure pressure)
  !!           ---
  !!           ...
  !!     ppre  ---  -----  ----------------------
  !!           ...
  !!           ---  dppre  zone T1 (transition 1)
  !!           ...
  !!           ---  -----  ----------------------
  !!           ...
  !!           ---
  !!           ...         zone T2 (transition 2)
  !!           ---
  !!           ...
  !!           ---  -----  ----------------------
  !!           ...
  !!           ---  dpsig  zone T3 (transition 3)
  !!           ...
  !!     psig  ---  -----  ----------------------
  !!           ...
  !!           ---  -----  zone S (pure sigma)
  !!           ...  dpbot
  !!     pbot  ---  -----  ----------------------
  !! </pre>
  !!
  !!   Detailed procedure description:
  !!   1 STEP 1.
  !!   The pressure profile is computed with respect to the given reference
  !!   surface pressure pbot. For this surface pressure, the 'sigma' thicknesses
  !!   dsig are assumed to be proportional to a quadratic polynomial in sigma sig
  !!   with zero intercepts sig1 and sig2 somewhere below and above the model
  !!   domain, respectively. That is,
  !!     dsig ~ (sig2-sig)*(sig-sig1)*dk
  !!   Integrating this differential equation gives
  !!     sig = (sig1*exp(c1*k+c2)+sig2)/(exp(c1*k+c2)+1)
  !!   The required boundary conditions sig(0)=1 and sig(levs)=0
  !!   fix the proportionality and integration constants c1 and c2.
  !!   The two crossing parameters (sig1 and sig2) are determined
  !!   by two input sigma thickness conditions dsig/dk at the bottom and top
  !!   which are respectively given as dpbot/(pbot-pupp) and dpupp/(pbot-pupp).
  !!   The crossing parameters are computed using Newton-Raphson iteration.
  !!   This procedure fixes the pressure profile for surface pressure pbot.
  !!   2 STEP 2.
  !!   The pressure profile is computed with respect to a minimum surface pressure.
  !!   This minimum surface pressure pmin is yet to be determined.
  !!   Divide the profile into zones:
  !!     zone U (pure pressure) from pupp to ptop
  !!     zone P (pure pressure) from pupp to ppre
  !!     zone T1 (transition 1) from ppre to ppre+dppre
  !!     zone T2 (transition 2) from ppre+dppre to psig-dpsig
  !!     zone T3 (transition 3) from psig-dpsig to psig
  !!     zone S (pure "sigma") from psig to pmin
  !!     (here sigma=p/ps so that d(ln(p))/dk is horizontally uniform)
  !!   The pressure profile in the pure pressure zone P is set from step 1.
  !!   The pressure thicknesses in zone T1 is set to be quadratic in level k.
  !!   The pressure thicknesses in zone T2 is set to be linear in level k.
  !!   The pressure thicknesses in zone T3 is set to be quadratic in level k.
  !!   The pressure profile in the pure sigma zone S is also set from step 1.
  !!   Thus there are nine unknowns:
  !!     the 3 polynomial coefficients in zone T1
  !!     the 2 polynomial coefficients in zone T2
  !!     the 3 polynomial coefficients in zone T3
  !!     and the 1 minimum surface pressure.
  !!   The nine conditions to determine these unknowns are:
  !!     the thickness and its derivative match at zone P and T1 boundary
  !!     the thickness and its derivative match at zone T1 and T2 boundary
  !!     the thickness and its derivative match at zone T2 and T3 boundary
  !!     the thickness and its derivative match at zone T3 and S boundary
  !!     the sum of the thicknesses of zones T1, T2, T3, and S is pmin-ppre
  !!   The unknowns are computed using standard linear decomposition.
  !!   This procedure fixes the pressure profile for surface pressure pmin.
  !!   3 STEP 3.
  !!   (Step 3 skipped if lupp=levs, in which case pupp=ptop and dpupp=dptop.)
  !!   The pressure in zone U is assumed to be the exponential of a cubic
  !!   polynomial in level k. The function must match the pressure at pupp,
  !!   as well as the thickness and its derivative there, and the pressure
  !!   at ptop+dptop at the second to top level. The latter 3 conditions
  !!   are determined by using standard linear decomposition.
  !!   4 STEP 4.
  !!   For an arbitrary surface pressure, the pressure profile is an linear
  !!   combination of the pressure profiles for surface pressures pbot and pmin
  !! <pre>
  !!     p(psfc)=p(pbot)*(psfc-pmin)/(pbot-pmin)+p(pmin)*(pbot-psfc)/(pbot-pmin)
  !! </pre>
  !!   from which the hybrid coordinate profiles ak and bk are found such that
  !! <pre>
  !!     p(psfc)=ak+bk*psfc
  !! </pre>
  !! @author Mark Iredell @date 2008-08-01

  subroutine vcoord_gen(levs, lupp, pbot, psig, ppre, pupp, ptop, &
                        dpbot, dpsig, dppre, dpupp, dptop, pmin, ak, bk)

    integer, intent(in) :: levs, lupp
    real(r8), intent(in) :: pbot, psig, ppre, pupp, ptop
    real(r8), intent(in) :: dpbot, dpsig, dppre, dpupp, dptop
    real(r8), intent(out) :: pmin, ak(0:levs), bk(0:levs)

    integer, parameter :: lo = 100, li = 10  ! outer and inner N-R iterations
    real(r8) pdif        ! thickness from pbot to pupp
    real(r8) delb        ! delta sig at bot
    real(r8) delt        ! delta sig at top
    real(r8) sig1        ! crossing parameter 1
    real(r8) sig2        ! crossing parameter 2
    real(r8) c1          ! proportionality constant
    real(r8) c2          ! integration constant
    real(r8) sig         ! sig variable
    real(r8) dsig        ! delta sig variable
    real(r8) delbio0     ! initial guess at delta sig at bot
    real(r8) deltio0     ! initial guess at delta sig at top
    real(r8) delbio      ! guess at delta sig at bot
    real(r8) deltio      ! guess at delta sig at top
    real(r8) c1sig1      ! d(c1)/d(sig1)
    real(r8) c1sig2      ! d(c1)/d(sig2)
    real(r8) p(2)        ! rhs in N-R iteration
    real(r8) fjac(2,2)   ! lhs in N-R iteration
    integer indx(2)  ! permutations in N-R iteration
    real(r8) ppred       ! pressure at T1-T2 boundary
    real(r8) spre        ! sig at P-T1 boundary
    real(r8) spred       ! sig at T1-T2 boundary
    real(r8) rkpre       ! level at P-T1 boundary
    real(r8) rkpred      ! level at T1-T2 boundary
    real(r8) pkpre       ! dp/dk at P-T1 boundary
    real(r8) pkkpre      ! d2p/dk2 at P-T1 boundary
    real(r8) psigd       ! pressure at T2-T3 boundary
    real(r8) ssig        ! sig at T3-S boundary
    real(r8) ssigd       ! sig at T2-T3 boundary
    real(r8) rksig       ! level at T3-S boundary
    real(r8) rksigd      ! level at T2-T3 boundary
    real(r8) pksig       ! dp/dk at T3-S boundary
    real(r8) pkksig      ! d2p/dk2 at T3-S boundary
    real(r8) p2sig       ! pressure at T3-S boundary for pmin surface pressure
    real(r8) p2sigd      ! pressure at T2-T3 boundary for pmin surface pressure
    real(r8) p2pred      ! pressure at T1-T2 boundary for pmin surface pressure
    real(r8) x2(9)       ! rhs in linear solver
    real(r8) a2(9,9)     ! lhs in linear solver
    integer indx2(9) ! permutations in linear solver
    real(r8) pkupp       ! dp/dk at U-P boundary
    real(r8) pkkupp      ! d2p/dk2 at U-P boundary
    real(r8) x3(3)       ! rhs in linear solver
    real(r8) a3(3,3)     ! lhs in linear solver
    integer indx3(3) ! permutations in linear solver
    real(r8) p1          ! pressure variable for pbot surface pressure
    real(r8) p2          ! pressure variable for pmin surface pressure
    real(r8) d           ! determinant permutation 
    integer io,ii,k

    !  STEP 1.
    pdif=pbot-pupp
    delb=dpbot/pdif
    delt=dpupp/pdif
    sig1=1+delb
    sig2=-delt
    c1=log(-sig2*(1-sig1)/sig1/(sig2-1))/lupp
    c2=log((sig2-1)/(1-sig1))
    sig=1
    dsig=(sig2-sig)*(sig-sig1)*c1/(sig1-sig2)
    delbio0=-dsig
    sig=0
    dsig=(sig2-sig)*(sig-sig1)*c1/(sig1-sig2)
    deltio0=-dsig
    !  Newton-Raphson iterations
    do io = 1, lo
      delbio=delbio0+(delb-delbio0)*io/lo
      deltio=deltio0+(delt-deltio0)*io/lo
      do ii = 1, li
        c1sig1=-1/(sig1*(1-sig1)*lupp)
        c1sig2=-1/(sig2*(sig2-1)*lupp)
        sig=1
        dsig=(sig2-sig)*(sig-sig1)*c1/(sig1-sig2)
        p(1)=-delbio-dsig
        fjac(1,1)=((-c1*(sig+sig2)+(sig-sig1)*c1sig1*(sig1+sig2)) &
                  *(sig2-sig)/(sig1+sig2)**2)
        fjac(1,2)=((c1*(sig+sig1)+(sig2-sig)*c1sig2*(sig1+sig2)) &
                  *(sig-sig1)/(sig1+sig2)**2)
        sig=0
        dsig=(sig2-sig)*(sig-sig1)*c1/(sig1-sig2)
        p(2)=-deltio-dsig
        fjac(2,1)=((-c1*(sig+sig2)+(sig-sig1)*c1sig1*(sig1+sig2)) &
                  *(sig2-sig)/(sig1+sig2)**2)
        fjac(2,2)=((c1*(sig+sig1)+(sig2-sig)*c1sig2*(sig1+sig2)) &
                  *(sig-sig1)/(sig1+sig2)**2)
        call ludcmp(fjac,2,2,indx,d)
        call lubksb(fjac,2,2,indx,p)
        sig1=sig1+p(1)
        sig2=sig2+p(2)
        c1=log(-sig2*(1-sig1)/sig1/(sig2-1))/lupp
        c2=log((sig2-1)/(1-sig1))
      end do
    end do

    !  STEP 2.
    !  Compute minimum surface pressure 
    ppred=ppre+dppre
    spre=(ppre-pupp)/pdif
    spred=(ppred-pupp)/pdif
    rkpre=(log((spre-sig2)/(sig1-spre))-c2)/c1
    rkpred=(log((spred-sig2)/(sig1-spred))-c2)/c1
    pkpre=pdif*c1*(sig2-spre)*(spre-sig1)/(sig1-sig2)
    pkkpre=pkpre*c1*(sig2+sig1-2*spre)/(sig1-sig2)
    psigd=psig-dpsig
    ssig=(psig-pupp)/pdif
    ssigd=(psigd-pupp)/pdif
    rksig=(log((ssig-sig2)/(sig1-ssig))-c2)/c1
    rksigd=(log((ssigd-sig2)/(sig1-ssigd))-c2)/c1
    pksig=pdif*c1*(sig2-ssig)*(ssig-sig1)/(sig1-sig2)
    pkksig=pksig*c1*(sig2+sig1-2*ssig)/(sig1-sig2)
    x2=0
    a2=0
    x2(1)=pkpre
    a2(1,1)=1
    a2(1,2)=rkpre
    a2(1,3)=rkpre**2
    x2(2)=pkkpre
    a2(2,2)=1
    a2(2,3)=2*rkpre
    a2(3,1)=1
    a2(3,2)=rkpred
    a2(3,3)=rkpred**2
    a2(3,4)=-1
    a2(3,5)=-rkpred
    a2(4,2)=1
    a2(4,3)=2*rkpred
    a2(4,5)=-1
    a2(5,4)=-1
    a2(5,5)=-rksigd
    a2(5,6)=1
    a2(5,7)=rksigd
    a2(5,8)=rksigd**2
    a2(6,5)=-1
    a2(6,7)=1
    a2(6,8)=2*rksigd
    a2(7,6)=1
    a2(7,7)=rksig
    a2(7,8)=rksig**2
    a2(7,9)=-pksig/pbot
    a2(8,7)=1
    a2(8,8)=2*rksig
    a2(8,9)=-pkksig/pbot
    x2(9)=ppre
    a2(9,1)=(rkpre-rkpred)
    a2(9,2)=(rkpre**2-rkpred**2)/2
    a2(9,3)=(rkpre**3-rkpred**3)/3
    a2(9,4)=(rkpred-rksigd)
    a2(9,5)=(rkpred**2-rksigd**2)/2
    a2(9,6)=(rksigd-rksig)
    a2(9,7)=(rksigd**2-rksig**2)/2
    a2(9,8)=(rksigd**3-rksig**3)/3
    a2(9,9)=psig/pbot
    call ludcmp(a2,9,9,indx2,d)
    call lubksb(a2,9,9,indx2,x2)
    pmin=x2(9)
    p2sig=psig/pbot*pmin
    p2sigd=p2sig &
          +x2(6)*(rksigd-rksig) &
          +x2(7)*(rksigd**2-rksig**2)/2 &
          +x2(8)*(rksigd**3-rksig**3)/3
    p2pred=p2sigd &
          +x2(4)*(rkpred-rksigd) &
          +x2(5)*(rkpred**2-rksigd**2)/2

    !  STEP 3.
    if (lupp < levs) then
      pkupp=pdif*c1*(sig2-0)*(0-sig1)/(sig1-sig2)
      pkkupp=pkupp*c1*(sig2+sig1-2*0)/(sig1-sig2)
      x3=0
      a3=0
      x3(1)=pkupp
      a3(1,1)=pupp
      x3(2)=pkkupp*pupp-pkupp**2
      a3(2,2)=pupp**2
      x3(3)=log((ptop+dptop)/pupp)
      a3(3,1)=levs-1-lupp
      a3(3,2)=(levs-1-lupp)**2/2
      a3(3,3)=(levs-1-lupp)**3/3
      call ludcmp(a3,3,3,indx3,d)
      call lubksb(a3,3,3,indx3,x3)
    end if

    !  STEP 4.
    !  Compute hybrid interface values
    ak(0) = 0
    bk(0) = 1
    do k = 1, levs -1
      if (k >= lupp) then
        p1=pupp*exp(x3(1)*(k-lupp)+x3(2)*(k-lupp)**2/2+x3(3)*(k-lupp)**3/3)
      else
        p1=(sig1*exp(c1*k+c2)+sig2)/(exp(c1*k+c2)+1)*pdif+pupp
      end if
      if (k >= rkpre) then
        p2=p1
      else if (k >= rkpred) then
        p2=p2pred+x2(1)*(k-rkpred) &
                 +x2(2)*(k**2-rkpred**2)/2 &
                 +x2(3)*(k**3-rkpred**3)/3
      else if (k >= rksigd) then
        p2=p2sigd+x2(4)*(k-rksigd) &
                 +x2(5)*(k**2-rksigd**2)/2
      else if (k >= rksig) then
        p2=p2sig+x2(6)*(k-rksig) &
                +x2(7)*(k**2-rksig**2)/2 &
                +x2(8)*(k**3-rksig**3)/3
      else
        p2=p1/pbot*pmin
      end if
      ak(k)=(p2*pbot-p1*pmin)/(pbot-pmin)
      bk(k)=(p1-p2)/(pbot-pmin)
    end do
    ak(levs) = ptop
    bk(levs) = 0

  end subroutine vcoord_gen

end module hybrid_coord_ncep_mod
