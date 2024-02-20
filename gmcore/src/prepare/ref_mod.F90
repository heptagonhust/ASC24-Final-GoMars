module ref_mod

  use const_mod
  use namelist_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod, only: proc
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
            ref_ps%d(i,j) = p0 * exp(-t0 / a + sqrt((t0 / a)**2 - 2 * gzs%d(i,j) / a / rd))
          end do
        end do
        call fill_halo(ref_ps)
        ref_ps_smth%d = ref_ps%d
        do i = 1, 50
          call laplace_damp_run(ref_ps_smth, 2, 1.0e6_r8)
        end do
        ref_ps_perb%d = ref_ps%d - ref_ps_smth%d
        end associate
      end do
    end select

  end subroutine ref_calc_ps

end module ref_mod
