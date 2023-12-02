! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module filter_mod

  use const_mod
  use filter_types_mod
  use latlon_field_types_mod

  implicit none

  private

  public filter_type
  public filter_run

  interface filter_run
    module procedure filter_run_2d
    module procedure filter_run_3d
    module procedure filter_run_4d
  end interface filter_run

contains

  subroutine filter_run_2d(filter, x, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field2d_type), intent(inout) :: x
    type(latlon_field2d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ims:x%mesh%full_ime)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, i, j, n, hn

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do j = js, je
      if (ngrid(j) > 1) then
        n  = ngrid(j)
        hn = (n - 1) / 2
        do i = is, ie
          tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j))
        end do
        if (present(y)) then
          y%d(:,j) = tmp
        else
          x%d(:,j) = tmp
        end if
      else if (present(y)) then
        y%d(:,j) = x%d(:,j)
      end if
    end do

  end subroutine filter_run_2d

  subroutine filter_run_3d(filter, x, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field3d_type), intent(inout) :: x
    type(latlon_field3d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ims:x%mesh%full_ime)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, ks, ke, i, j, k, n, hn

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell', 'lon')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    case ('lat', 'vtx')
      js = x%mesh%half_jds
      je = x%mesh%half_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lat
      ngrid => filter%ngrid_lat
    case ('lev')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%half_kds
      ke = x%mesh%half_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do k = ks, ke
      do j = js, je
        if (ngrid(j) > 1) then
          n  = ngrid(j)
          hn = (n - 1) / 2
          do i = is, ie
            tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y%d(:,j,k) = tmp
          else
            x%d(:,j,k) = tmp
          end if
        else if (present(y)) then
          y%d(:,j,k) = x%d(:,j,k)
        end if
      end do
    end do

  end subroutine filter_run_3d

  subroutine filter_run_4d(filter, x, i4, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field4d_type), intent(inout) :: x
    integer, intent(in) :: i4
    type(latlon_field4d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ims:x%mesh%full_ime)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, ks, ke, i, j, k, n, hn

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell', 'lon')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do k = ks, ke
      do j = js, je
        if (ngrid(j) > 1) then
          n  = ngrid(j)
          hn = (n - 1) / 2
          do i = is, ie
            tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j,k,i4))
          end do
          if (present(y)) then
            y%d(:,j,k,i4) = tmp
          else
            x%d(:,j,k,i4) = tmp
          end if
        else if (present(y)) then
          y%d(:,j,k,i4) = x%d(:,j,k,i4)
        end if
      end do
    end do

  end subroutine filter_run_4d

end module filter_mod
