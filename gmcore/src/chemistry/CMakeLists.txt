cmake_policy(SET CMP0076 NEW)

if (USE_CHEM)
  target_sources(gmcore PRIVATE
    mozart/cfc11star.F90
    mozart/chemistry.F90
    mozart/mo_airglow.F90
    mozart/mo_chemini.F90
    mozart/mo_cph.F90
    mozart/mo_heatnirco2.F90
    mozart/mo_neu_wetdep.F90
    mozart/mo_srf_emissions.F90
    mozart/mo_waccm_hrates.F90
    mozart/ocean_emis.F90
    mozart/set_cp.F90
    pp_trop_mam4/chem_mods.F90
    pp_trop_mam4/m_spc_id.F90
    pp_trop_mam4/mo_adjrxt.F90
    pp_trop_mam4/mo_exp_sol.F90
    pp_trop_mam4/mo_exp_sol.F90
    pp_trop_mam4/mo_imp_sol.F90
    pp_trop_mam4/mo_indprd.F90
    pp_trop_mam4/mo_lin_matrix.F90
    pp_trop_mam4/mo_lu_factor.F90
    pp_trop_mam4/mo_lu_solve.F90
    pp_trop_mam4/mo_nln_matrix.F90
    pp_trop_mam4/mo_phtadj.F90
    pp_trop_mam4/mo_prod_loss.F90
    pp_trop_mam4/mo_rxt_rates_conv.F90
    pp_trop_mam4/mo_setrxt.F90
    pp_trop_mam4/mo_sim_dat.F90
  )
else ()
  target_sources(gmcore PRIVATE
    pp_none/chemistry.F90
    pp_none/chem_mods.F90
    pp_none/m_spc_id.F90
    pp_none/mo_adjrxt.F90
    pp_none/mo_exp_sol.F90
    pp_none/mo_imp_sol.F90
    pp_none/mo_indprd.F90
    pp_none/mo_lin_matrix.F90
    pp_none/mo_lu_factor.F90
    pp_none/mo_lu_solve.F90
    pp_none/mo_nln_matrix.F90
    pp_none/mo_phtadj.F90
    pp_none/mo_prod_loss.F90
    pp_none/mo_rxt_rates_conv.F90
    pp_none/mo_setrxt.F90
  )
endif ()

target_sources(gmcore PRIVATE
  aerosol/cldaero_mod.F90
  aerosol/drydep_mod.F90
  aerosol/dust_common.F90
  aerosol/dust_sediment_mod.F90
  aerosol/mo_setsox.F90
  aerosol/soil_erod_mod.F90
  aerosol/sslt_sections.F90
  aerosol/wetdep.F90
  mozart/charge_neutrality.F90
  mozart/chem_prod_loss_diags.F90
  mozart/chlorine_loading_data.F90
  mozart/clybry_fam.F90
  mozart/epp_ionization.F90
  mozart/euvac.F90
  mozart/fire_emissions.F90
  mozart/gas_wetdep_opts.F90
  mozart/gcr_ionization.F90
  mozart/lin_strat_chem.F90
  mozart/linoz_data.F90
  mozart/llnl_O1D_to_2OH_adj.F90
  mozart/m_sad_data.F90
  mozart/mo_aero_settling.F90
  mozart/mo_airmas.F90
  mozart/mo_airplane.F90
  mozart/mo_apex.F90
  mozart/mo_aurora.F90
  mozart/mo_calcoe.F90
  mozart/mo_chem_utls.F90
  mozart/mo_chm_diags.F90
  mozart/mo_drydep.F90
  mozart/mo_extfrc.F90
  mozart/mo_fstrat.F90
  mozart/mo_gas_phase_chemdr.F90
  mozart/mo_ghg_chem.F90
  mozart/mo_inter.F90
  mozart/mo_jeuv.F90
  mozart/mo_jlong.F90
  mozart/mo_jpl.F90
  mozart/mo_jshort.F90
  mozart/mo_lightning.F90
  mozart/mo_lymana.F90
  mozart/mo_mass_xforms.F90
  mozart/mo_mean_mass.F90
  mozart/mo_negtrc.F90
  mozart/mo_params.F90
  mozart/mo_pchem.F90
  mozart/mo_photo.F90
  mozart/mo_photoin.F90
  mozart/mo_ps2str.F90
  mozart/mo_rtlink.F90
  mozart/mo_sad.F90
  mozart/mo_schu.F90
  mozart/mo_setaer.F90
  mozart/mo_setair.F90
  mozart/mo_setcld.F90
  mozart/mo_setext.F90
  mozart/mo_sethet.F90
  mozart/mo_setinv.F90
  mozart/mo_seto2.F90
  mozart/mo_setozo.F90
  mozart/mo_setz.F90
  mozart/mo_snoe.F90
  mozart/mo_sphers.F90
  mozart/mo_strato_rates.F90
  mozart/mo_sulf.F90
  mozart/mo_synoz.F90
  mozart/mo_tgcm_ubc.F90
  mozart/mo_tracname.F90
  mozart/mo_trislv.F90
  mozart/mo_tuv_inti.F90
  mozart/mo_usrrxt.F90
  mozart/mo_waveall.F90
  mozart/mo_wavelab.F90
  mozart/mo_wavelen.F90
  mozart/mo_waveo3.F90
  mozart/mo_xsections.F90
  mozart/mo_zadj.F90
  mozart/noy_ubc.F90
  mozart/photo_bkgrnd.F90
  mozart/rate_diags.F90
  mozart/short_lived_species.F90
  mozart/species_sums_diags.F90
  mozart/spehox.F90
  mozart/sums_utils.F90
  mozart/tracer_cnst.F90
  mozart/tracer_srcs.F90
  mozart/upper_bc.F90
  msis/physics_constants.f90
  msis/physics_msis.f90
  msis/utils_constants.f90
  msis/utils_spline.f90
  utils/aerodep_flx.F90
  utils/aircraft_emit.F90
  utils/apex.F90
  utils/horizontal_interpolate.F90
  utils/input_data_utils.F90
  utils/m_types.F90
  utils/mo_constants.F90
  utils/mo_flbc.F90
  utils/mo_msis_ubc.F90
  utils/mo_util.F90
  utils/modal_aero_calcsize.F90
  utils/modal_aero_deposition.F90
  utils/modal_aero_wateruptake.F90
  utils/prescribed_aero.F90
  utils/prescribed_ghg.F90
  utils/prescribed_ozone.F90
  utils/prescribed_strataero.F90
  utils/prescribed_volcaero.F90
  utils/solar_data.F90
  utils/solar_euv_data.F90
  utils/solar_irrad_data.F90
  utils/solar_parms_data.F90
  utils/solar_wind_data.F90
  utils/time_utils.F90
  utils/tracer_data.F90
)

if (USE_BULK_AERO)
  target_sources(gmcore PRIVATE
    bulk_aero/aero_model.F90
    bulk_aero/aerosol_depvel.F90
    bulk_aero/dust_model.F90
    bulk_aero/mo_aerosols.F90
    bulk_aero/mo_setsoa.F90
    bulk_aero/seasalt_model.F90
    bulk_aero/sox_cldaero_mod.F90
  )
elseif (USE_MODAL_AERO)
  target_sources(gmcore PRIVATE
    modal_aero/aero_model.F90
    modal_aero/dust_model.F90
    modal_aero/modal_aero_coag.F90
    modal_aero/modal_aero_convproc.F90
    modal_aero/modal_aero_data.F90
    modal_aero/modal_aero_gasaerexch.F90
    modal_aero/modal_aero_newnuc.F90
    modal_aero/modal_aero_rename.F90
    modal_aero/seasalt_model.F90
    modal_aero/sox_cldaero_mod.F90
  )
endif ()