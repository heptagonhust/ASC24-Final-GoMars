module sigma_coord_mod

  use flogger
  use namelist_mod
  use const_mod
  use latlon_mesh_mod
  use process_mod
  use mars_vert_coord_mod

  implicit none

  private

  public sigma_coord_init
  public sigma_coord_final
  public sigma_coord_calc_mg
  public sigma_coord_calc_mg_lev
  public sigma_coord_calc_dmgdt_lev

  real(r8), allocatable, dimension(:) :: sigi
  real(r8), allocatable, dimension(:) :: sig

  namelist /sigma_coord/ &
    sig, sigi

contains

  subroutine sigma_coord_init(namelist_file, template)

    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: template

    integer k, ierr

    if (allocated(sigi)) deallocate(sigi)
    if (allocated(sig )) deallocate(sig )

    allocate(sigi(nlev+1)); sigi = 0
    allocate(sig (nlev  )); sig  = 0

    if (present(namelist_file)) then
      open(10, file=namelist_file, status='old')
      read(10, nml=sigma_coord, iostat=ierr)
      close(10)
    end if

    if (ierr /= 0 .and. .not. (baroclinic .or. advection) .and. proc%is_root()) then
      call log_notice('Run shallow-water model.')
    end if

    if (present(template)) then
      select case (template)
      case ('nasa24')
        call mars_vert_coord_nasa24(ptop, sigi, sig)
      case default
        call sigma_coord_test_l26(ptop, sigi, sig)
      end select
    else
      call sigma_coord_test_l26(ptop, sigi, sig)
    end if

    do k = 1, nlev
      global_mesh%full_lev(k) = sig(k)
    end do
    do k = 1, nlev + 1
      global_mesh%half_lev(k) = sigi(k)
    end do

  end subroutine sigma_coord_init

  subroutine sigma_coord_final()

    deallocate(sigi)
    deallocate(sig )

  end subroutine sigma_coord_final

  pure real(r8) function sigma_coord_calc_mg(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    res = sig(k) * (mgs - ptop) + ptop

  end function sigma_coord_calc_mg

  pure real(r8) function sigma_coord_calc_mg_lev(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    res = sigi(k) * (mgs - ptop) + ptop

  end function sigma_coord_calc_mg_lev

  pure real(r8) function sigma_coord_calc_dmgdt_lev(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = sigi(k) * dmgsdt

  end function sigma_coord_calc_dmgdt_lev

  subroutine sigma_coord_test_l26(ptop, sigi, sig)

    real(r8), intent(out) :: ptop
    real(r8), intent(out) :: sigi(:)
    real(r8), intent(out) :: sig(:)

    if (nlev /= 26) then
      call log_error('sigma_coord_test_l26: nlev /= 26!')
    end if

    sigi = [               &
      0.00000000000000000, & !  1
      0.00270708138944839, & !  2
      0.00770525729190707, & !  3
      0.0158928127640089 , & !  4
      0.0277039568735249 , & !  5
      0.0425225717851423 , & !  6
      0.0595424438370873 , & !  7
      0.0764861789945283 , & !  8
      0.09037000846462   , & !  9
      0.106703636692459  , & ! 10
      0.125919292600401  , & ! 11
      0.148525517478954  , & ! 12
      0.175120586853872  , & ! 13
      0.206408247111974  , & ! 14
      0.243216666654805  , & ! 15
      0.286519798881535  , & ! 16
      0.33746367336785   , & ! 17
      0.397396487729323  , & ! 18
      0.467904355709719  , & ! 19
      0.550853249115629  , & ! 20
      0.648438367513145  , & ! 21
      0.743820793246579  , & ! 22
      0.830649646952043  , & ! 23
      0.903087677425118  , & ! 24
      0.955900674543326  , & ! 25
      0.985079453568069  , & ! 26
      1.000000000000000    & ! 27
    ]
    sig = [                &
      0.00135354069472419, & !  1
      0.00520616934067773, & !  2
      0.011799035027958  , & !  3
      0.0217983848187669 , & !  4
      0.0351132643293336 , & !  5
      0.0510325078111148 , & !  6
      0.0680143114158078 , & !  7
      0.0834280937295741 , & !  8
      0.0985368225785395 , & !  9
      0.11631146464643   , & ! 10
      0.137222405039678  , & ! 11
      0.161823052166413  , & ! 12
      0.190764416982923  , & ! 13
      0.22481245688339   , & ! 14
      0.26486823276817   , & ! 15
      0.311991736124692  , & ! 16
      0.367430080548587  , & ! 17
      0.432650421719521  , & ! 18
      0.509378802412674  , & ! 19
      0.599645808314387  , & ! 20
      0.696129580379862  , & ! 21
      0.787235220099311  , & ! 22
      0.86686866218858   , & ! 23
      0.929494175984222  , & ! 24
      0.970490064055697  , & ! 25
      0.992539726784035    & ! 26
    ]
    ptop = 219.406699761748

  end subroutine sigma_coord_test_l26

end module sigma_coord_mod
