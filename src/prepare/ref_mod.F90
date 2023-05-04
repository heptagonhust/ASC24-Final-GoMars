module ref_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use laplace_damp_mod

  implicit none

  private

  public ref_calc_ps

contains

  subroutine ref_calc_ps()

    real(r8) p0, t0, a
    integer i, j, iblk

    if (proc%is_root()) call log_notice('Calculating reference surface pressure and its perturbation.')

    select case (planet)
    case ('earth')
      do iblk = 1, size(blocks)
        associate (block       => blocks(iblk)                   , &
                   mesh        => blocks(iblk)%mesh              , &
                   gzs         => blocks(iblk)%static%gzs        , &
                   ref_ps      => blocks(iblk)%static%ref_ps     , &
                   ref_ps_smth => blocks(iblk)%static%ref_ps_smth, &
                   ref_ps_perb => blocks(iblk)%static%ref_ps_perb)
        p0 = 1.0e5_r8
        t0 = 300.0_r8
        a  = 50.0_r8
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ref_ps(i,j) = p0 * exp(-t0 / a + sqrt((t0 / a)**2 - 2 * gzs(i,j) / a / rd))
          end do
        end do
        call fill_halo(block%halo, ref_ps, full_lon=.true., full_lat=.true.)
        ref_ps_smth = ref_ps
        do i = 1, 50
          call laplace_damp_on_cell(block%mesh, block%halo, 2, ref_ps_smth, coef=1.0e6_r8, lon_coef=decay_from_pole)
        end do
        ref_ps_perb = ref_ps - ref_ps_smth
        end associate
      end do
    end select

  end subroutine ref_calc_ps

end module ref_mod
