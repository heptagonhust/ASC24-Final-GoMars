module accum_mod

  use container
  use flogger
  use datetime
  use const_mod
  use time_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public accum_type

  integer, public, parameter :: freq_hourly   = 1
  integer, public, parameter :: freq_daily    = 2
  integer, public, parameter :: freq_monthly  = 3
  integer, public, parameter :: freq_yearly   = 4
  integer, public, parameter :: stat_avg      = 1
  integer, public, parameter :: stat_min      = 2
  integer, public, parameter :: stat_max      = 3

  type accum_type
    character(30) :: name = 'N/A'
    character(30) :: units = 'N/A'
    character(50) :: long_name = 'N/A'
    character(30) :: from = 'N/A'
    character(30) :: var_name = 'N/A'
    logical :: active = .true.
    integer :: freq = 0
    integer :: stat = 0
    integer :: time_steps = 0
    logical :: tick = .false.
    type(timedelta_type) dt
    type(datetime_type) last_time
    real(r8), allocatable, dimension(:,:,:,:) :: array
  contains
    procedure :: init => accum_init
    procedure :: clear => accum_clear
    procedure, private :: run_start => accum_run_start
    procedure, private :: run_end => accum_run_end
    procedure :: accum_run_1d
    procedure :: accum_run_2d
    procedure :: accum_run_3d
    procedure :: accum_run_4d
    generic :: run => accum_run_1d, accum_run_2d, accum_run_3d, accum_run_4d
    final :: accum_final
  end type accum_type

contains

  subroutine accum_init(this, name, units, long_name, from, var_name, freq, stat, array_shape, active)

    class(accum_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: units
    character(*), intent(in) :: long_name
    character(*), intent(in) :: from
    character(*), intent(in) :: var_name
    integer, intent(in) :: freq
    integer, intent(in) :: stat
    integer, intent(in) :: array_shape(:)
    logical, intent(in) :: active

    integer, parameter :: calendar = datetime_noleap_calendar

    call this%clear()

    this%name = name
    this%units = units
    this%long_name = long_name
    this%from = from
    this%var_name = var_name
    this%freq = freq
    this%stat = stat
    this%active = active
    this%last_time = curr_time

    associate (s => array_shape)
    select case (size(array_shape))
    case (1)
      allocate(this%array(s(1),   1,   1,  1 ))
    case (2)
      allocate(this%array(s(1),s(2),   1,  1 ))
    case (3)
      allocate(this%array(s(1),s(2),s(3),  1 ))
    case (4)
      allocate(this%array(s(1),s(2),s(3),s(4)))
    case default
      if (proc%is_root()) then
        call log_error('Unsupportted array dimension!', __FILE__, __LINE__)
      end if
    end select
    end associate

    select case (this%freq)
    case (freq_hourly)
      call this%dt%init(days=0, hours=1)
    case (freq_daily)
      call this%dt%init(days=1)
    case (freq_monthly)
      call this%dt%init(days=days_of_month(this%last_time%year, this%last_time%month, calendar))
    case (freq_yearly)
      call this%dt%init(days=days_of_year(this%last_time%year, calendar))
    end select

    select case (stat)
    case (stat_avg)
      this%array = 0
    case (stat_max)
      this%array = -inf
    case (stat_min)
      this%array = inf
    end select

  end subroutine accum_init

  subroutine accum_clear(this)

    class(accum_type), intent(inout) :: this

    if (allocated(this%array)) deallocate(this%array)

  end subroutine accum_clear

  subroutine accum_run_start(this)

    class(accum_type), intent(inout) :: this

    if (this%tick) then
      select case (this%stat)
      case (stat_avg)
        this%array = 0
      case (stat_max)
        this%array = -inf
      case (stat_min)
        this%array = inf
      end select
      this%tick = .false.
    end if

  end subroutine accum_run_start

  subroutine accum_run_end(this)

    class(accum_type), intent(inout) :: this

    this%time_steps = this%time_steps + 1
    this%tick = curr_time - this%last_time >= this%dt
    if (this%tick) then
      if (this%stat == stat_avg) then
        this%array = this%array / this%time_steps
      end if
      this%last_time = curr_time
      this%time_steps = 0
    end if

  end subroutine accum_run_end

  subroutine accum_run_1d(this, array)

    class(accum_type), intent(inout) :: this
    real(r8), intent(in) :: array(:)

    call this%run_start()
    select case (this%stat)
    case (stat_avg)
      this%array(:,1,1,1) = this%array(:,1,1,1) + array
    case (stat_max)
      this%array(:,1,1,1) = max(this%array(:,1,1,1), array)
    case (stat_min)
      this%array(:,1,1,1) = min(this%array(:,1,1,1), array)
    end select
    call this%run_end()

  end subroutine accum_run_1d

  subroutine accum_run_2d(this, array)

    class(accum_type), intent(inout) :: this
    real(r8), intent(in) :: array(:,:)

    call this%run_start()
    select case (this%stat)
    case (stat_avg)
      this%array(:,:,1,1) = this%array(:,:,1,1) + array
    case (stat_max)
      this%array(:,:,1,1) = max(this%array(:,:,1,1), array)
    case (stat_min)
      this%array(:,:,1,1) = min(this%array(:,:,1,1), array)
    end select
    call this%run_end()

  end subroutine accum_run_2d

  subroutine accum_run_3d(this, array)

    class(accum_type), intent(inout) :: this
    real(r8), intent(in) :: array(:,:,:)

    call this%run_start()
    select case (this%stat)
    case (stat_avg)
      this%array(:,:,:,1) = this%array(:,:,:,1) + array
    case (stat_max)
      this%array(:,:,:,1) = max(this%array(:,:,:,1), array)
    case (stat_min)
      this%array(:,:,:,1) = min(this%array(:,:,:,1), array)
    end select
    call this%run_end()

  end subroutine accum_run_3d

  subroutine accum_run_4d(this, array)

    class(accum_type), intent(inout) :: this
    real(r8), intent(in) :: array(:,:,:,:)

    call this%run_start()
    select case (this%stat)
    case (stat_avg)
      this%array(:,:,:,:) = this%array(:,:,:,:) + array
    case (stat_max)
      this%array(:,:,:,:) = max(this%array(:,:,:,:), array)
    case (stat_min)
      this%array(:,:,:,:) = min(this%array(:,:,:,:), array)
    end select
    call this%run_end()

  end subroutine accum_run_4d

  subroutine accum_final(this)

    type(accum_type), intent(inout) :: this

    call this%clear()

  end subroutine accum_final

end module accum_mod
