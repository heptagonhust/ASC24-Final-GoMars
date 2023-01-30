module gas_mod

  implicit none

  private

  public gas_mixture_init
  public gas_mixture_final
  public gases
  public major_gas
  public minor_gas

  integer, parameter :: r8 = 8

  real(r8), parameter :: kb       = 1.380649e-23_r8     ! Boltzmann constant
  real(r8), parameter :: na       = 6.02214076e23_r8    ! Avogadro constant
  real(r8), parameter :: ru       = kb * na             ! Universal gas constant
  real(r8), parameter :: m_o2     = 32.00e-3_r8         ! O2 molar mass (kg mol-1)
  real(r8), parameter :: m_n2     = 28.01e-3_r8         ! N2 molar mass (kg mol-1)
  real(r8), parameter :: m_ar     = 39.94e-3_r8         ! Ar molar mass (kg mol-1)
  real(r8), parameter :: m_co2    = 44.01e-3_r8         ! CO2 molar mass (kg mol-1)
  real(r8), parameter :: m_h2o    = 18.02e-3_r8         ! H2O molar mass (kg mol-1)

  type gas_type
    character(30) :: name = ''
    integer  :: n_atom  = 0
    real(r8) :: dof     = 0 ! Degree of freedom
    real(r8) :: m_ratio = 0 ! Mass ratio (1)
    real(r8) :: v_ratio = 0 ! Volume ratio or mole ratio (1)
    real(r8) :: m       = 0 ! Molecular mass (kg mol-1)
    real(r8) :: r       = 0 ! Gas constant (J K-1 kg-1)
    real(r8) :: cp      = 0 ! Specific heat capacity under constant pressure (J kg-1 K-1)
    real(r8) :: cv      = 0 ! Specific heat capacity under constant volume (J kg-1 K-1)
    real(r8) :: gamma   = 0 ! cp / cv
    real(r8) :: kappa   = 0 ! r / cp
    real(r8) :: l       = 0 ! Latent heat rate (J kg-1)
  contains
    procedure :: gas_init_1
    procedure :: gas_init_2
    generic :: init => gas_init_1, gas_init_2
  end type gas_type

  type(gas_type), allocatable, target :: gases(:)
  type(gas_type) major_gas, minor_gas

contains

  subroutine gas_mixture_init(planet)

    character(*), intent(in) :: planet

    select case (planet)
    case ('earth')
      allocate(gases(5))
      ! NOTE: DOF is calculated at 293K.
      call gases(1)%init('n2' , n_atom=2, m=m_n2 , v_ratio=0.7808_r8  , dof=5.00607939581368_r8)
      call gases(2)%init('o2' , n_atom=2, m=m_o2 , v_ratio=0.2095_r8  , dof=5.07260684628201_r8)
      call gases(3)%init('ar' , n_atom=1, m=m_ar , v_ratio=0.0093_r8  , dof=3.02464704198518_r8)
      call gases(4)%init('co2', n_atom=3, m=m_co2, v_ratio=0.000385_r8, dof=6.93407411251379_r8)
      call gases(5)%init('h2o', n_atom=3, m=m_h2o, v_ratio=999.0_r8   , dof=6.0_r8)
      call major_gas%init('dry_air', gases(1:4))
      call minor_gas%init('water_vapor', gases(5:5))
    case ('mars')
      allocate(gases(4))
      ! NOTE: DOF is calculated at 200K.
      call gases(1)%init('co2', n_atom=3, m=m_co2, v_ratio=0.9600_r8  , dof=5.78098392778265_r8)
      call gases(2)%init('n2' , n_atom=2, m=m_n2 , v_ratio=0.0189_r8  , dof=5.00042596534376_r8)
      call gases(3)%init('ar' , n_atom=1, m=m_ar , v_ratio=0.0193_r8  , dof=3.02464704198518_r8)
      call gases(4)%init('o2' , n_atom=2, m=m_o2 , v_ratio=0.00145_r8 , dof=5.00466195768837_r8)
      call major_gas%init('dry_air', gases(1:4))
    end select

  end subroutine gas_mixture_init

  subroutine gas_mixture_final()

    if (allocated(gases)) deallocate(gases)

  end subroutine gas_mixture_final

  subroutine gas_init_1(this, name, n_atom, m, v_ratio, dof)

    class(gas_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: n_atom
    real(r8), intent(in) :: m ! kg mol-1
    real(r8), intent(in), optional :: v_ratio ! 1
    real(r8), intent(in), optional :: dof

    real(r8) dof_

    this%name = name
    this%n_atom = n_atom
    this%m = m
    this%r = ru / m
    if (present(dof)) then
      dof_ = dof
    else
      select case (n_atom)
      case (1)
        dof_ = 3
      case (2)
        dof_ = 5
      end select
    end if
    this%cv = dof_ * 0.5_r8 * this%r
    this%cp = this%cv + this%r
    if (present(v_ratio)) this%v_ratio = v_ratio

    this%gamma = this%cp / this%cv
    this%kappa = this%r / this%cp

  end subroutine gas_init_1

  subroutine gas_init_2(this, name, gases)

    class(gas_type), intent(inout) :: this
    character(*), intent(in) :: name
    type(gas_type), intent(inout) :: gases(:)

    real(r8) sum_ratio
    integer i

    this%name = name

    ! Calculate mean molecular mass.
    sum_ratio = sum(gases(:)%v_ratio)
    this%m = 0
    do i = 1, size(gases)
      this%m = this%m + gases(i)%v_ratio * gases(i)%m
    end do
    this%m = this%m / sum_ratio
    this%r = ru / this%m

    ! Calculate volume mixing ratio from mass mixing ratio.
    do i = 1, size(gases)
      gases(i)%m_ratio = gases(i)%v_ratio * gases(i)%m / this%m
    end do

    ! Calculate mean specific heat capacity.
    sum_ratio = sum(gases(:)%m_ratio)
    this%cv = 0
    do i = 1, size(gases)
      this%cv = this%cv + gases(i)%m_ratio * gases(i)%cv
    end do
    this%cv = this%cv / sum_ratio
    this%cp = this%cv + this%r

    this%gamma = this%cp / this%cv
    this%kappa = this%r / this%cp

  end subroutine gas_init_2

end module gas_mod
