! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_physics_output_mod

  use fiona
  use mars_nasa_objects_mod

  implicit none

  private

  public mars_nasa_physics_add_output
  public mars_nasa_physics_output

contains

  subroutine mars_nasa_physics_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    character(3) :: dims(2) = ['lon', 'lat']

    call fiona_add_var(tag, 'tin' , long_name='Surface thermal inertia'     , units='', dim_names=dims, dtype=dtype)

  end subroutine mars_nasa_physics_add_output

  subroutine mars_nasa_physics_output(tag, iblk, start, count)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk
    integer, intent(in) :: start(3)
    integer, intent(in) :: count(3)

    associate (static => objects(iblk)%static, state => objects(iblk)%state)
    call fiona_output(tag, 'tin' , reshape(static%tin, count(1:2)), start=start(1:2), count=count(1:2))
    end associate

  end subroutine mars_nasa_physics_output

end module mars_nasa_physics_output_mod
