module dycore

  use string_utils, only: to_upper
    
  implicit none

  private

  public dycore_is

contains
    
  logical function dycore_is(name)
    
    ! Determine the dynamical core in use.
    
    character(*), intent(in) :: name
    
    character(len(name)) :: uname
    
    uname = to_upper(name)
    
    if (uname == 'GMCORE') then
      dycore_is = .true.
    else
      dycore_is = .false.
    end if
    
  end function dycore_is
    
end module dycore
    