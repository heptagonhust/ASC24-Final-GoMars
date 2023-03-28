module vert_coord_mod

  use flogger
  use const_mod
  use namelist_mod
  use sigma_coord_mod
  use hybrid_coord_mod
  use mesh_mod
  use process_mod

  implicit none

  private

  public vert_coord_init_stage1
  public vert_coord_init_stage2
  public vert_coord_final
  public vert_coord_calc_mg
  public vert_coord_calc_mg_lev
  public vert_coord_calc_dmgdt_lev
  public hyai, hybi

  interface
    pure real(r8) function vert_coord_calc_mg_interface(k, mgs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: mgs
    end function vert_coord_calc_mg_interface

    pure real(r8) function vert_coord_calc_mg_lev_interface(k, mgs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: mgs
    end function vert_coord_calc_mg_lev_interface

    pure real(r8) function vert_coord_calc_dmgdt_lev_interface(k, dmgsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dmgsdt
    end function vert_coord_calc_dmgdt_lev_interface
  end interface

  procedure(vert_coord_calc_mg_interface), pointer :: vert_coord_calc_mg
  procedure(vert_coord_calc_mg_lev_interface), pointer :: vert_coord_calc_mg_lev
  procedure(vert_coord_calc_dmgdt_lev_interface), pointer :: vert_coord_calc_dmgdt_lev

contains

  subroutine vert_coord_init_stage1(nlev, namelist_file, scheme, template)

    integer, intent(in) :: nlev
    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: scheme
    character(*), intent(in), optional :: template

    if (present(scheme)) then
      if (vert_coord_scheme /= scheme .and. proc%is_root()) then
        call log_notice('Change vert_coord_scheme to ' // trim(scheme) // '.')
      end if
      vert_coord_scheme = scheme
    end if
    if (present(template)) then
      if (vert_coord_template /= template .and. proc%is_root()) then
        call log_notice('Change vert_coord_template to ' // trim(template) // '.')
      end if
      vert_coord_template = template
    else if (baroclinic .and. (vert_coord_scheme == 'hybrid' .and. vert_coord_template == 'N/A')) then
      call log_error('Parameter vert_coord_template is not set!', pid=proc%id)
    end if

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_init_stage1(nlev, namelist_file, vert_coord_template)
      vert_coord_calc_mg => sigma_coord_calc_mg
      vert_coord_calc_mg_lev => sigma_coord_calc_mg_lev
      vert_coord_calc_dmgdt_lev => sigma_coord_calc_dmgdt_lev
    case ('hybrid')
      call hybrid_coord_init_stage1(nlev, namelist_file, vert_coord_template)
      vert_coord_calc_mg => hybrid_coord_calc_mg
      vert_coord_calc_mg_lev => hybrid_coord_calc_mg_lev
      vert_coord_calc_dmgdt_lev => hybrid_coord_calc_dmgdt_lev
    end select

  end subroutine vert_coord_init_stage1

  subroutine vert_coord_init_stage2()

    integer k

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_init_stage2()
    case ('hybrid')
      call hybrid_coord_init_stage2()
    end select

    ! Set vertical level intervals.
    if (global_mesh%full_nlev > 1) then
      do k = global_mesh%full_kds, global_mesh%full_kde
        global_mesh%full_dlev(k) = global_mesh%half_lev(k+1) - global_mesh%half_lev(k)
      end do
    else
      global_mesh%full_dlev(1) = 1
    end if
    do k = global_mesh%half_kds + 1, global_mesh%half_kde - 1
      global_mesh%half_dlev(k) = global_mesh%full_lev(k) - global_mesh%full_lev(k-1)
      global_mesh%half_dlev_upper(k) = global_mesh%half_lev(k) - global_mesh%full_lev(k-1)
      global_mesh%half_dlev_lower(k) = global_mesh%full_lev(k) - global_mesh%half_lev(k)
    end do
    global_mesh%half_dlev(1) = global_mesh%full_lev(1) - global_mesh%half_lev(1)
    global_mesh%half_dlev_lower(1) = global_mesh%half_dlev(1)
    global_mesh%half_dlev(global_mesh%half_kde) = global_mesh%half_lev(global_mesh%half_kde) - global_mesh%full_lev(global_mesh%full_kde)
    global_mesh%half_dlev_upper(global_mesh%half_kde) = global_mesh%half_dlev(global_mesh%half_kde)

  end subroutine vert_coord_init_stage2

  subroutine vert_coord_final()

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_final()
    case ('hybrid')
      call hybrid_coord_final()
    end select

  end subroutine vert_coord_final

end module vert_coord_mod
