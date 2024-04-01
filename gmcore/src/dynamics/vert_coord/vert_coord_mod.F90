! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module vert_coord_mod

  use flogger
  use const_mod
  use namelist_mod
  use sigma_coord_mod
  use hybrid_coord_mod
  use smooth_coord_mod
  use latlon_mesh_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public vert_coord_init
  public vert_coord_final
  public vert_coord_calc_mg
  public vert_coord_calc_mg_lev
  public vert_coord_calc_dmgdt_lev
  public hyai, hybi, hyam, hybm, nlevp

  interface
     real(r8) function calc_mg_interface(k, mgs, ref_ps_perb)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: mgs
      real(r8), intent(in), optional :: ref_ps_perb
    end function calc_mg_interface

     real(r8) function calc_mg_lev_interface(k, mgs, ref_ps_perb)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: mgs
      real(r8), intent(in), optional :: ref_ps_perb
    end function calc_mg_lev_interface

     real(r8) function calc_dmgdt_interface(k, dmgsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dmgsdt
    end function calc_dmgdt_interface

     real(r8) function calc_dmgdt_lev_interface(k, dmgsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dmgsdt
    end function calc_dmgdt_lev_interface

     real(r8) function calc_ddmgdt_interface(k, dmgsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dmgsdt
    end function calc_ddmgdt_interface
  end interface

  procedure(calc_mg_interface       ), pointer :: vert_coord_calc_mg
  procedure(calc_mg_lev_interface   ), pointer :: vert_coord_calc_mg_lev
  procedure(calc_dmgdt_interface    ), pointer :: vert_coord_calc_dmgdt
  procedure(calc_dmgdt_lev_interface), pointer :: vert_coord_calc_dmgdt_lev
  procedure(calc_ddmgdt_interface   ), pointer :: vert_coord_calc_ddmgdt

contains

  subroutine vert_coord_init(namelist_file, scheme, template)

    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: scheme
    character(*), intent(in), optional :: template

    integer k

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
      call sigma_coord_init(namelist_file, vert_coord_template)
      vert_coord_calc_mg        => sigma_coord_calc_mg
      vert_coord_calc_mg_lev    => sigma_coord_calc_mg_lev
      vert_coord_calc_dmgdt     => sigma_coord_calc_dmgdt
      vert_coord_calc_dmgdt_lev => sigma_coord_calc_dmgdt_lev
      vert_coord_calc_ddmgdt    => sigma_coord_calc_ddmgdt
    case ('hybrid')
      call hybrid_coord_init(namelist_file, vert_coord_template)
      vert_coord_calc_mg        => hybrid_coord_calc_mg
      vert_coord_calc_mg_lev    => hybrid_coord_calc_mg_lev
      vert_coord_calc_dmgdt     => hybrid_coord_calc_dmgdt
      vert_coord_calc_dmgdt_lev => hybrid_coord_calc_dmgdt_lev
      vert_coord_calc_ddmgdt    => hybrid_coord_calc_ddmgdt
    case ('smooth')
      call smooth_coord_init()
      vert_coord_calc_mg        => smooth_coord_calc_mg
      vert_coord_calc_mg_lev    => smooth_coord_calc_mg_lev
      vert_coord_calc_dmgdt     => smooth_coord_calc_dmgdt
      vert_coord_calc_dmgdt_lev => smooth_coord_calc_dmgdt_lev
      vert_coord_calc_ddmgdt    => smooth_coord_calc_ddmgdt
    case default
      call log_error('Invalid vert_coord_scheme: ' // trim(vert_coord_scheme) // '!')
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

  end subroutine vert_coord_init

  subroutine vert_coord_final()

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_final()
    case ('hybrid')
      call hybrid_coord_final()
    end select

  end subroutine vert_coord_final

end module vert_coord_mod
