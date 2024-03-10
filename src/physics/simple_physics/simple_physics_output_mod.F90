! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_output_mod

  use fiona
  use simple_physics_objects_mod

  implicit none

  private

  public simple_physics_add_output
  public simple_physics_output

contains

  subroutine simple_physics_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

  end subroutine simple_physics_add_output

  subroutine simple_physics_output(tag, iblk, start, count)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk
    integer, intent(in) :: start(3)
    integer, intent(in) :: count(3)

  end subroutine simple_physics_output

end module simple_physics_output_mod
