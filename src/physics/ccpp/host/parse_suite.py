#!/usr/bin/env python3

'''
Parse CCPP suite definition file in XML format.
'''

import argparse
from xml.etree import ElementTree
import os
import configparser

parser = argparse.ArgumentParser(description='Parse CCPP suite definition file in XML format.')
parser.add_argument('-s', '--suite-xml', dest='suite_xml', required=True)
parser.add_argument('-p', '--physics-root', dest='physics_root', required=True)
parser.add_argument('-o', '--output-code', required=True)
args = parser.parse_args()

default_str_len = 256

refactors = {
	'adjusted_vertical_layer_dimension_for_radiation': { 'name': 'nlev_rad_adj' },
	'aerosol_aware_multiplicative_rain_conversion_parameter_for_deep_convection': { 'name': 'asolfac_deep' },
	'aerosol_aware_multiplicative_rain_conversion_parameter_for_shallow_convection': { 'name': 'asolfac_shallow' },
	'air_pressure_at_interface': { 'name': 'p_lev' },
	'air_pressure_at_interface_for_RRTMGP': { 'name': 'p_lev_rrtmgp' },
	'air_pressure_at_surface_adjacent_layer': { 'name': 'p_bot', 'pointer': True },
	'air_pressure_difference_between_midlayers': { 'name': 'dp' },
	'air_pressure': { 'name': 'p' },
	'air_temperature': { 'name': 't' },
	'air_temperature_at_2m': { 'name': 't_2m' },
	'air_temperature_at_lowest_model_layer_for_diag': { 'name': 't_bot_diag' },
	'air_temperature_at_surface_adjacent_layer_on_radiation_timestep': { 'name': 't_bot_rad_lw' },
	'air_temperature_at_surface_adjacent_layer': { 'name': 't_bot', 'pointer': True },
	'air_temperature_of_new_state_at_surface_adjacent_layer': { 'name': 't_bot_new' },
	'air_temperature_of_new_state': { 'name': 't_new' },
	'air_temperature_save_from_convective_parameterization': { 'name': 't_old_cnv' },
	'air_temperature_save': { 'name': 't_old' },
	'area_type': { 'name': 'area_type' },
	'atmosphere_boundary_layer_thickness': { 'name': 'hpbl' },
	'bounded_surface_roughness_length_for_heat_over_land': { 'name': 'ztmax_lnd' },
	'bulk_richardson_number_at_lowest_model_level': { 'name': 'rib_bot' },
	'bulk_richardson_number_at_lowest_model_level_over_land': { 'name': 'rib_bot_lnd' },
	'ccpp_loop_counter': { 'name': 'ccpp_loop_counter' },
	'cell_area': { 'name': 'cell_area' },
	'convective_cloud_area_fraction': { 'name': 'cldfrac_conv' },
	'convective_cloud_condensate_mixing_ratio': { 'name': 'qc_conv_phy' },
	'convective_cloud_water_mixing_ratio': { 'name': 'qc_conv' },
	'convective_precipitation_rate_on_previous_timestep': { 'name': 'prcp_conv_old' },
	'convective_transportable_tracers': { 'name': 'q_conv' },
	'convexity_of_subgrid_orography_small_scale': { 'name': 'ocss' },
	'control_for_microphysics_scheme': { 'name': 'mp_physics' },
	'cosine_of_solar_zenith_angle_for_daytime_points_on_radiation_timestep': { 'name': 'coszen_day' },
	'cloud_condensed_water_mixing_ratio_save': { 'name': 'qc_old' },
	'cloud_ice_mixing_ratio': { 'name': 'qi' },
	'cloud_ice_mixing_ratio_of_new_state': { 'name': 'qi_new' },
	'cloud_liquid_water_mixing_ratio': { 'name': 'qc' },
	'cloud_liquid_water_mixing_ratio_of_new_state': { 'name': 'qc_new' },
	'cloud_work_function': { 'name': 'cloud_work' },
	'date_and_time_of_forecast_in_united_states_order': { 'name': 'jdate' },
	'daytime_points': { 'name': 'idx_day' },
	'detrainment_conversion_parameter_for_deep_convection': { 'name': 'detrain_deep' },
	'detrainment_conversion_parameter_for_shallow_convection': { 'name': 'detrain_shallow' },
	'dimensionless_exner_function': { 'name': 'pk' },
	'dimensionless_exner_function_at_surface_adjacent_layer': { 'name': 'pk_bot', 'pointer': True },
	'do_merra2_aerosol_awareness': { 'name': 'do_merra2_aer' },
	'effective_radius_of_stratiform_cloud_ice_particle': { 'name': 'reff_cice_strat' },
	'effective_radius_of_stratiform_cloud_liquid_water_particle': { 'name': 'reff_cwat_strat' },
	'effective_radius_of_stratiform_cloud_rain_particle': { 'name': 'reff_crain_strat' },
	'effective_radius_of_stratiform_cloud_snow_particle': { 'name': 'reff_csnow_strat' },
	'entrainment_rate_coefficient_for_deep_convection': { 'name': 'entrain_deep' },
	'entrainment_rate_coefficient_for_shallow_convection': { 'name': 'entrain_shallow' },
	'explicit_precipitation_rate_on_previous_timestep': { 'name': 'prcp_explicit_old' },
	'filename_of_namelist': { 'name': 'namelist_file' },
	'flag_for_aerosol_physics': { 'name': 'flag_aer' },
	'flag_for_calling_longwave_radiation': { 'name': 'flag_rad_lw' },
	'flag_for_calling_shortwave_radiation': { 'name': 'flag_rad_sw' },
	'flag_for_diagnostics': { 'name': 'flag_diag' },
	'flag_for_hydrostatic_solver': { 'name': 'flag_hydrostatic' },
	'flag_for_hydrostatic_solver_for_fast_physics': { 'name': 'flag_hydrostatic_fast' },
	'flag_for_generic_tendency_due_to_planetary_boundary_layer': { 'name': 'flag_generic_tend_pbl' },
	'flag_for_output_of_tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep_assuming_clear_sky': { 'name': 'flag_output_dtdt_rad_lw_clr' },
	'flag_for_output_of_tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep_assuming_clear_sky': { 'name': 'flag_output_dtdt_rad_sw_clr' },
	'flag_for_radar_reflectivity': { 'name': 'flag_radar_ref' },
	'flag_nonzero_land_surface_fraction': { 'name': 'flag_land' },
	'fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height_small_scale': { 'name': 'clxss' },
	'freezing_point_temperature_of_seawater': { 'name': 'const_tfrz' },
	'fresh_liquid_water_density_at_0c': { 'name': 'const_rhow0c' },
	'gas_constant_of_dry_air': { 'name': 'const_rd' },
	'gas_constant_water_vapor': { 'name': 'const_rv' },
	'geopotential_at_interface': { 'name': 'gz_lev' },
	'geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature': { 'name': 'dgz_o_tv' },
	'graupel_mixing_ratio': { 'name': 'qg' },
	'graupel_mixing_ratio_of_new_state': { 'name': 'qg_new' },
	'graupel_precipitation_rate_on_previous_timestep': { 'name': 'graupel_old' },
	'gravitational_acceleration': { 'name': 'const_g' },
	'height_above_ground_at_lowest_model_layer': { 'name': 'z_bot' },
	'height_above_mean_sea_level': { 'name': 'z_sfc' },
	'horizontal_loop_extent': { 'name': 'ncol' },
	'identifier_for_morrison_gettelman_microphysics_scheme': { 'name': 'mp_physics_mg' },
	'index_for_ice_cloud_condensate_vertical_diffusion_tracer': { 'name': 'idx_qi_vdiff' },
	'index_for_liquid_cloud_condensate_vertical_diffusion_tracer': { 'name': 'idx_qc_vdiff' },
	'index_for_rain_water_vertical_diffusion_tracer': { 'name': 'idx_qrain_vdiff' },
	'index_for_turbulent_kinetic_energy_vertical_diffusion_tracer': { 'name': 'idx_tke_vdiff' },
	'index_of_cloud_area_fraction_in_atmosphere_layer_in_tracer_concentration_array': { 'name': 'idx_cld' },
	'index_of_cloud_ice_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qi' },
	'index_of_cloud_liquid_water_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qc' },
	'index_of_graupel_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qg' },
	'index_of_mass_number_concentration_of_cloud_droplets_in_tracer_concentration_array': { 'name': 'idx_nc' },
	'index_of_mass_number_concentration_of_cloud_ice_in_tracer_concentration_array': { 'name': 'idx_ni' },
	'index_of_mass_number_concentration_of_hygroscopic_aerosols_in_tracer_concentration_array': { 'name': 'idx_naer' },
	'index_of_ozone_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qo3' },
	'index_of_rain_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qrain' },
	'index_of_snow_mixing_ratio_in_tracer_concentration_array': { 'name': 'idx_qsnow' },
	'index_of_timestep': { 'name': 'time_step' },
	'ice_water_mixing_ratio_save': { 'name': 'qi_old' },
	'ice_precipitation_rate_on_previous_timestep': { 'name': 'prcp_ice_old' },
	'instantaneous_cosine_of_zenith_angle': { 'name': 'coszen' },
	'instantaneous_surface_upward_latent_heat_flux': { 'name': 'lhflx_sfc' },
	'instantaneous_surface_upward_sensible_heat_flux': { 'name': 'shflx_sfc' },
	'instantaneous_surface_x_momentum_flux': { 'name': 'uflx_sfc' },
	'instantaneous_surface_y_momentum_flux': { 'name': 'vflx_sfc' },
	'instantaneous_tendency_of_specific_humidity_due_to_microphysics': { 'name': 'dqvdt_mp' },
	'joules_per_calorie_constant': { 'name': 'const_jcal' },
	'kinematic_surface_upward_latent_heat_flux_over_ice': { 'name': 'lhflx_ice' },
	'kinematic_surface_upward_latent_heat_flux_over_land': { 'name': 'lhflx_lnd' },
	'kinematic_surface_upward_latent_heat_flux_over_water': { 'name': 'lhflx_wat' },
	'kinematic_surface_upward_sensible_heat_flux_over_ice': { 'name': 'shflx_ice' },
	'kinematic_surface_upward_sensible_heat_flux_over_land': { 'name': 'shflx_lnd' },
	'kinematic_surface_upward_sensible_heat_flux_over_water': { 'name': 'shflx_wat' },
	'kinematic_surface_upward_sensible_heat_flux_reduced_by_surface_roughness_and_vegetation': { 'name': 'shflx' },
	'lagrangian_tendency_of_air_pressure': { 'name': 'omega' },
	'latent_heat_of_fusion_of_water_at_0C': { 'name': 'const_hfus' },
	'latent_heat_of_vaporization_of_water_at_0C': { 'name': 'const_hvap' },
	'latitude': { 'name': 'lat' },
	'level_of_dividing_streamline': { 'name': 'rdxzb' },
	'longitude': { 'name': 'lon' },
	'lw_fluxes_top_atmosphere': { 'name': 'lwflx_toa' },
	'lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep': { 'name': 'prcp_conv_dyn' },
	'lwe_thickness_of_deep_convective_precipitation_amount': { 'name': 'prcp_deep' },
	'lwe_thickness_of_explicit_precipitation_amount': { 'name': 'prcp_explicit' },
	'lwe_thickness_of_explicit_rain_amount': { 'name': 'rain_explicit' },
	'lwe_thickness_of_graupel_amount': { 'name': 'graupel' },
	'lwe_thickness_of_graupel_amount_on_dynamics_timestep': { 'name': 'graupel_dyn' },
	'lwe_thickness_of_ice_amount': { 'name': 'ice' },
	'lwe_thickness_of_ice_amount_on_dynamics_timestep': { 'name': 'ice_dyn' },
	'lwe_thickness_of_shallow_convective_precipitation_amount': { 'name': 'prcp_shallow' },
	'lwe_thickness_of_snow_amount': { 'name': 'snow' },
	'lwe_thickness_of_snow_amount_on_dynamics_timestep': { 'name': 'snow_dyn' },
	'magnitude_of_perturbation_of_vegetation_fraction': { 'name': 'vegfrac_pert_mag' },
	'mass_mixing_ratio_of_aerosol_from_gocart_or_merra2': { 'name': 'qaer' },
	'mass_number_concentration_of_hygroscopic_aerosols': { 'name': 'nwfa' },
	'mass_number_concentration_of_hygroscopic_aerosols_of_new_state': { 'name': 'nwfa_new' },
	'maximum_temperature_at_2m_over_maximum_hourly_time_interval': { 'name': 't_2m_max_hr' },
	'maximum_u_wind_at_10m_over_maximum_hourly_time_interval': { 'name': 'u_10m_max_hr' },
	'maximum_v_wind_at_10m_over_maximum_hourly_time_interval': { 'name': 'v_10m_max_hr' },
	'maximum_wind_at_10m': { 'name': 'gust_10m' },
	'maximum_x_wind_at_10m': { 'name': 'u_10m_max' },
	'maximum_y_wind_at_10m': { 'name': 'v_10m_max' },
	'minimum_temperature_at_2m_over_maximum_hourly_time_interval': { 'name': 't_2m_min_hr' },
	'minimum_value_of_specific_humidity': { 'name': 'const_qv_min' },
	'momentum_transport_reduction_factor_due_to_pressure_gradient_force_for_deep_convection': { 'name': 'pgcon_deep' },
	'momentum_transport_reduction_factor_due_to_pressure_gradient_force_for_shallow_convection': { 'name': 'pgcon_shallow' },
	'mpi_rank': { 'name': 'mpi_rank' },
	'mpi_root': { 'name': 'mpi_root' },
	'Monin_Obukhov_similarity_function_for_heat': { 'name': 'fh' },
	'Monin_Obukhov_similarity_function_for_heat_at_2m_over_land': { 'name': 'fh_2m_lnd' },
	'Monin_Obukhov_similarity_function_for_heat_over_ice': { 'name': 'fh_ice' },
	'Monin_Obukhov_similarity_function_for_heat_over_land': { 'name': 'fh_lnd' },
	'Monin_Obukhov_similarity_function_for_heat_over_water': { 'name': 'fh_wat' },
	'Monin_Obukhov_similarity_function_for_momentum_at_10m_over_land': { 'name': 'fm_10m_lnd' },
	'Monin_Obukhov_similarity_function_for_momentum_at_10m_over_water': { 'name': 'fm_10m_wat' },
	'Monin_Obukhov_similarity_function_for_momentum': { 'name': 'fm' },
	'Monin_Obukhov_similarity_function_for_momentum_over_ice': { 'name': 'fm_ice' },
	'Monin_Obukhov_similarity_function_for_momentum_over_land': { 'name': 'fm_lnd' },
	'Monin_Obukhov_similarity_function_for_momentum_over_water': { 'name': 'fm_wat' },
	'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep': { 'name': 'prcp' },
	'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ice': { 'name': 'prcp_ice' },
	'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_land': { 'name': 'prcp_lnd' },
	'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_water': { 'name': 'prcp_wat' },
	'number_of_chemical_tracers': { 'name': 'ntracers_chem' },
	'number_of_columns_per_RRTMGP_LW_block': { 'name': 'ncol_rrtmgp' },
	'number_of_columns_per_RRTMGP_SW_block': { 'skip': True },
	'number_of_condensate_species': { 'name': 'ntracers_cond' },
	'number_of_days_in_current_year': { 'name': 'jday' },
	'number_of_openmp_threads': { 'name': 'nthreads' },
	'number_of_tracers': { 'name': 'ntracers' },
	'number_of_tracers_for_samf': { 'name': 'ntracers_samf' },
	'number_of_vertical_diffusion_tracers': { 'name': 'ntracers_vdiff' },
	'ocean_mixed_layer_thickness': { 'name': 'zmix_ocn' },
	'period_of_shortwave_radiation_calls': { 'name': 'dt_rad_sw' },
	'perturbation_of_heat_to_momentum_roughness_length_ratio': { 'name': 'z0t_pert' },
	'perturbation_of_leaf_area_index': { 'name': 'lai_pert' },
	'perturbation_of_momentum_roughness_length': { 'name': 'z0m_pert' },
	'perturbation_of_soil_type_b_parameter': { 'name': 'soil_type_b_pert' },
	'perturbation_of_vegetation_fraction': { 'name': 'vegfrac_pert' },
	'pi': { 'name': 'const_pi' },
	'process_split_cumulative_tendency_of_air_temperature': { 'name': 'dtdt' },
	'process_split_cumulative_tendency_of_x_wind': { 'name': 'dudt' },
	'process_split_cumulative_tendency_of_y_wind': { 'name': 'dvdt' },
	'rain_conversion_parameter_for_deep_convection': { 'name': 'c0s_deep' },
	'rain_conversion_parameter_for_shallow_convection': { 'name': 'c0s_shallow' },
	'rain_mixing_ratio_of_new_state': { 'name': 'qrain_new' },
	'random_number': { 'name': 'random_number' },
	'ratio_of_dry_air_to_water_vapor_gas_constants_minus_one': { 'name': 'const_epsm1' },
	'ratio_of_dry_air_to_water_vapor_gas_constants': { 'name': 'const_eps' },
	'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer': { 'name': 'ratio_pk_bot' },
	'ratio_of_gas_constant_dry_air_to_gravitational_acceleration': { 'name': 'const_rd_o_g' },
	'ratio_of_gas_constant_dry_air_to_specific_heat_of_dry_air_at_constant_pressure': { 'name': 'const_rd_o_cpd' },
	'ratio_of_vapor_to_dry_air_gas_constants_minus_one': { 'name': 'const_fvirt' },
	'RRTMGP_lw_flux_profile_downward_allsky_on_radiation_timestep': { 'name': 'lwflxdn_rrtmgp' },
	'RRTMGP_lw_flux_profile_upward_allsky_on_radiation_timestep': { 'name': 'lwflxup_rrtmgp' },
	'sea_ice_area_fraction_of_sea_area_fraction': { 'name': 'seaice_area_frac' },
	'sea_land_ice_mask': { 'name': 'slmsk' },
	'sea_land_ice_mask_cice': { 'name': 'slmsk_cice' },
	'sea_surface_temperature': { 'name': 'sst' },
	'snow_mixing_ratio': { 'name': 'qsnow' },
	'snow_mixing_ratio_of_new_state': { 'name': 'qsnow_new' },
	'snowfall_rate_on_previous_timestep': { 'name': 'snow_old' },
	'soil_type_classification': { 'name': 'soil_type' },
	'specific_heat_of_dry_air_at_constant_pressure': { 'name': 'const_cpd' },
	'specific_humidity_at_surface_adjacent_layer': { 'name': 'qv_bot', 'pointer': True },
	'specific_humidity': { 'name': 'qv' },
	'specific_humidity_of_new_state': { 'name': 'qv_new' },
	'specific_humidity_of_new_state_at_surface_adjacent_layer': { 'name': 'qv_new_bot', 'pointer': True },
	'standard_deviation_of_subgrid_orography_small_scale': { 'name': 'z_sfc_var' },
	'stefan_boltzmann_constant': { 'name': 'const_sbc' },
	'surface_air_pressure_diag': { 'name': 'p_sfc_diag' },
	'surface_air_pressure': { 'name': 'p_sfc' },
	'surface_air_temperature_for_radiation': { 'name': 't_bot_rad' },
	'surface_albedo_diffuse_NIR_over_land': { 'name': 'alb_dif_nir_lnd' },
	'surface_albedo_diffuse_visible_over_land': { 'name': 'alb_dif_vis_lnd' },
	'surface_albedo_direct_NIR_over_land': { 'name': 'alb_dir_nir_lnd' },
	'surface_albedo_direct_visible_over_land': { 'name': 'alb_dir_vis_lnd' },
	'surface_albedo_for_diffused_shortwave_on_radiation_timestep': { 'name': 'alb_dif' },
	'surface_albedo_perturbation': { 'name': 'alb_pert' },
	'surface_dimensionless_exner_function': { 'name': 'pk_sfc' },
	'surface_downwelling_diffuse_nir_shortwave_flux_on_radiation_timestep': { 'name': 'swflxdn_sfc_dif_nir'},
	'surface_downwelling_diffuse_uv_and_vis_shortwave_flux_on_radiation_timestep': { 'name': 'swflxdn_sfc_dif_uv_vis' },
	'surface_downwelling_direct_nir_shortwave_flux_on_radiation_timestep': { 'name': 'swflxdn_sfc_dir_nir' },
	'surface_downwelling_direct_uv_and_vis_shortwave_flux_on_radiation_timestep': { 'name': 'swflxdn_sfc_dir_uv_vis' },
	'surface_downwelling_longwave_flux': { 'name': 'lwflxdn_sfc' },
	'surface_downwelling_longwave_flux_absorbed_by_ground_over_ice': { 'name': 'lwflxdn_sfc_abs_ice' },
	'surface_downwelling_longwave_flux_absorbed_by_ground_over_land': { 'name': 'lwflxdn_sfc_abs_lnd' },
	'surface_downwelling_longwave_flux_absorbed_by_ground_over_water': { 'name': 'lwflxdn_sfc_abs_wat' },
	'surface_downwelling_shortwave_flux': { 'name': 'swflxdn_sfc' },
	'surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice': { 'name': 'ch_ice' },
	'surface_drag_coefficient_for_heat_and_moisture_in_air_over_land': { 'name': 'ch_lnd' },
	'surface_drag_coefficient_for_heat_and_moisture_in_air_over_water': { 'name': 'ch_wat' },
	'surface_drag_coefficient_for_momentum_in_air_over_ice': { 'name': 'cm_ice' },
	'surface_drag_coefficient_for_momentum_in_air_over_land': { 'name': 'cm_lnd' },
	'surface_drag_coefficient_for_momentum_in_air_over_water': { 'name': 'cm_wat' },
	'surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice': { 'name': 'chh_ice' },
	'surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land': { 'name': 'chh_lnd' },
	'surface_drag_mass_flux_for_heat_and_moisture_in_air_over_water': { 'name': 'chh_wat' },
	'surface_drag_wind_speed_for_momentum_in_air_over_ice': { 'name': 'cmm_ice' },
	'surface_drag_wind_speed_for_momentum_in_air_over_land': { 'name': 'cmm_lnd' },
	'surface_drag_wind_speed_for_momentum_in_air_over_water': { 'name': 'cmm_wat' },
	'surface_friction_velocity_over_ice': { 'name': 'ustar_ice' },
	'surface_friction_velocity_over_land': { 'name': 'ustar_lnd' },
	'surface_friction_velocity_over_water': { 'name': 'ustar_wat' },
	'surface_friction_velocity': { 'name': 'ustar' },
	'surface_ground_temperature_for_radiation': { 'name': 't_sfc_rad' },
	'surface_longwave_emissivity_over_ice': { 'name': 'emiss_lw_ice' },
	'surface_longwave_emissivity_over_land': { 'name': 'emiss_lw_lnd' },
	'surface_longwave_emissivity_over_water': { 'name': 'emiss_lw_wat' },
	'surface_lw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep': { 'name': 'lwflx_sfc_clr' },
	'surface_net_downwelling_shortwave_flux': { 'name': 'swflxdn_sfc_net' },
	'surface_roughness_length_from_wave_model': { 'name': 'z0_wav' },
	'surface_roughness_length_over_ice': { 'name': 'z0_ice' },
	'surface_roughness_length_over_land': { 'name': 'z0_lnd' },
	'surface_roughness_length_over_water': { 'name': 'z0_wat' },
	'surface_roughness_length': { 'name': 'z0' },
	'surface_runoff_flux': { 'name': 'rofflx_sfc' },
	'surface_skin_temperature_after_iteration_over_land': { 'name': 't_sfc_iter_lnd' },
	'surface_skin_temperature_after_iteration_over_water': { 'name': 't_sfc_iter_wat' },
	'surface_skin_temperature': { 'name': 't_sfc' },
	'surface_skin_temperature_for_coupling': { 'name': 't_sfc_cpl' },
	'surface_skin_temperature_for_nsst': { 'name': 't_sfc_nsst' },
	'surface_skin_temperature_on_radiation_timestep': { 'name': 't_sfc_radtime' },
	'surface_skin_temperature_over_ice': { 'name': 't_sfc_ice' },
	'surface_skin_temperature_over_land': { 'name': 't_sfc_lnd' },
	'surface_skin_temperature_over_water': { 'name': 't_sfc_wat' },
	'surface_slope_classification': { 'name': 'slope_type' },
	'surface_snow_area_fraction_over_land': { 'name': 'snowfrac_lnd' },
	'surface_snow_thickness_water_equivalent_over_ice': { 'name': 'snowd_ice' },
	'surface_snow_thickness_water_equivalent_over_land': { 'name': 'snowd_lnd' },
	'surface_sw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep': { 'name': 'swflx_sfc_clr' },
	'surface_specific_humidity_over_ice': { 'name': 'qv_sfc_ice' },
	'surface_specific_humidity_over_land': { 'name': 'qv_sfc_lnd' },
	'surface_specific_humidity_over_water': { 'name': 'qv_sfc_wat' },
	'surface_specific_humidity': { 'name': 'qv_sfc' },
	'surface_upward_potential_latent_heat_flux_over_ice': { 'name': 'ep_ice' },
	'surface_upward_potential_latent_heat_flux_over_land': { 'name': 'ep_lnd' },
	'surface_upward_potential_latent_heat_flux_over_water': { 'name': 'ep_wat' },
	'surface_upwelling_diffuse_nir_shortwave_flux_on_radiation_timestep': { 'name': 'swflxup_sfc_dif_nir' },
	'surface_upwelling_diffuse_uv_and_vis_shortwave_flux_on_radiation_timestep': { 'name': 'swflxup_sfc_dif_uv_vis' },
	'surface_upwelling_direct_nir_shortwave_flux_on_radiation_timestep': { 'name': 'swflxup_sfc_dir_nir' },
	'surface_upwelling_direct_uv_and_vis_shortwave_flux_on_radiation_timestep': { 'name': 'swflxup_sfc_dir_uv_vis' },
	'surface_wind_stress_over_water': { 'name': 'stress_wat' },
	'surface_wind_stress_over_land': { 'name': 'stress_lnd' },
	'sw_fluxes_top_atmosphere': { 'name': 'swflx_toa' },
	'temperature_at_2m_for_coupling': { 'name': 't_2m_cpl' },
	'temperature_at_2m_from_noahmp': { 'name': 't_2m_noahmp' },
	'temperature_at_zero_celsius': { 'name': 'const_t0c' },
	'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_timestep': { 'name': 'dtdt_rad_lw_clr' },
	'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep': { 'name': 'dtdt_lw_rad' },
	'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_timestep': { 'name': 'dtdt_rad_sw_clr' },
	'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep': { 'name': 'dtdt_sw_rad' },
	'tendency_of_vertically_diffused_tracer_concentration': { 'name': 'dqdt_vdiff' },
	'tendency_of_x_wind_due_to_blocking_drag': { 'name': 'dudt_blkd' },
	'tendency_of_x_wind_due_to_form_drag': { 'name': 'dudt_formd' },
	'tendency_of_x_wind_due_to_mesoscale_orographic_gravity_wave_drag': { 'name': 'dudt_mgwd' },
	'tendency_of_x_wind_due_to_small_scale_gravity_wave_drag': { 'name': 'dudt_sgwd' },
	'tendency_of_y_wind_due_to_blocking_drag': { 'name': 'dvdt_blkd' },
	'tendency_of_y_wind_due_to_form_drag': { 'name': 'dvdt_formd' },
	'tendency_of_y_wind_due_to_mesoscale_orographic_gravity_wave_drag': { 'name': 'dvdt_mgwd' },
	'tendency_of_y_wind_due_to_small_scale_gravity_wave_drag': { 'name': 'dvdt_sgwd' },
	'timestep_for_dynamics': { 'name': 'dt_dyn' },
	'timestep_for_physics': { 'name': 'dt_phys' },
	'total_runoff': { 'name': 'rof_tot' },
	'tracer_concentration': { 'name': 'q' },
	'tracer_concentration_save': { 'name': 'q_old' },
	'triple_point_temperature_of_water': { 'name': 'const_ttrp' },
	'unfiltered_height_above_mean_sea_level': { 'name': 'z_sfc_uf' },
	'unsmoothed_nonhydrostatic_upward_air_velocity': { 'name': 'w_uf' },
	'updated_tendency_of_air_temperature_due_to_longwave_heating_on_physics_timestep': { 'name': 'dtdt_lw_phys' },
	'upward_heat_flux_in_soil': { 'name': 'hflx_soil' },
	'upward_heat_flux_in_soil_over_ice': { 'name': 'hflx_soil_ice' },
	'upward_heat_flux_in_soil_over_land': { 'name': 'hflx_soil_lnd' },
	'upward_heat_flux_in_soil_over_water': { 'name': 'hflx_soil_wat' },
	'vegetation_area_fraction': { 'name': 'vegfrac' },
	'vegetation_type_classification': { 'name': 'veg_type' },
	'vertically_diffused_tracer_concentration': { 'name': 'q_vdiff' },
	'vertical_dimension_for_fast_physics': { 'name': 'nlev_fast' },
	'vertical_dimension_for_radiation': { 'name': 'nlev_rad' },
	'vertical_dimension_of_soil': { 'name': 'nlev_soil' },
	'vertical_layer_dimension': { 'name': 'nlev' },
	'vertically_integrated_x_momentum_flux_due_to_blocking_drag': { 'name': 'dudt_blkd_col' },
	'vertically_integrated_x_momentum_flux_due_to_form_drag': { 'name': 'dudt_formd_col' },
	'vertically_integrated_x_momentum_flux_due_to_mesoscale_orographic_gravity_wave_drag': { 'name': 'dudt_mgwd_col' },
	'vertically_integrated_x_momentum_flux_due_to_small_scale_gravity_wave_drag': { 'name': 'dudt_sgwd_col' },
	'vertically_integrated_y_momentum_flux_due_to_blocking_drag': { 'name': 'dvdt_blkd_col' },
	'vertically_integrated_y_momentum_flux_due_to_form_drag': { 'name': 'dvdt_formd_col' },
	'vertically_integrated_y_momentum_flux_due_to_mesoscale_orographic_gravity_wave_drag': { 'name': 'dvdt_mgwd_col' },
	'vertically_integrated_y_momentum_flux_due_to_small_scale_gravity_wave_drag': { 'name': 'dvdt_sgwd_col' },
	'water_equivalent_accumulated_snow_depth_over_ice': { 'name': 'weasd_ice' },
	'water_equivalent_accumulated_snow_depth_over_land': { 'name': 'weasd_lnd' },
	'water_vapor_specific_humidity_at_Lagrangian_surface': { 'skip': True },
	'water_vapor_specific_humidity_at_lowest_model_layer_for_diag': { 'name': 'q_bot_diag' },
	'wind_speed_at_lowest_model_layer': { 'name': 'wsp_bot' },
	'x_wind_at_10m_for_coupling': { 'name': 'u_10m_cpl' },
	'x_wind_at_10m': { 'name': 'u_10m' },
	'x_wind_at_lowest_model_layer_for_diag': { 'name': 'u_bot_diag' },
	'x_wind_at_surface_adjacent_layer': { 'name': 'u_bot', 'pointer': True },
	'x_wind_of_new_state': { 'name': 'u_new' },
	'x_wind_of_new_state_at_surface_adjacent_layer': { 'name': 'u_bot_new', 'pointer': True },
	'x_wind_save': { 'name': 'u_old' },
	'x_wind': { 'name': 'u' },
	'y_wind_at_10m_for_coupling': { 'name': 'v_10m_cpl' },
	'y_wind_at_10m': { 'name': 'v_10m' },
	'y_wind_at_lowest_model_layer_for_diag': { 'name': 'v_bot_diag' },
	'y_wind_at_surface_adjacent_layer': { 'name': 'v_bot', 'pointer': True },
	'y_wind_of_new_state': { 'name': 'v_new' },
	'y_wind_of_new_state_at_surface_adjacent_layer': { 'name': 'v_bot_new', 'pointer': True },
	'y_wind_save': { 'name': 'v_old' },
	'y_wind': { 'name': 'v' },
}

ccpp_builtin_vars = ('errmsg', 'errflg', 'blkno', 'ccpp_loop_counter', 'nsteps')

def parse_suite(xml_file, physics_root):
	var_list = {}
	tree = ElementTree.parse(xml_file)
	suite = tree.getroot()
	print(f'[Notice]: Parsing suite {suite.get("name")} version {suite.get("version")}.')
	for group in suite:
		print(f'[Notice]: Parsing group {group.get("name")}.')
		for subcycle in group:
			for scheme in subcycle:
				print(f'[Notice]: Parsing scheme {scheme.text.strip()}.')
				if os.path.isfile(f'{physics_root}/{scheme.text.strip()}.meta'):
					meta = configparser.ConfigParser(strict=False)
					meta.read(f'{physics_root}/{scheme.text.strip()}.meta')
					for section in meta.sections():
						if section == 'ccpp-table-properties':
							print(f'==> {meta[section].get("name")}')
						elif section == 'ccpp-arg-table':
							print(f'*** {meta[section].get("name")}')
						else:
							standard_name = meta[section].get('standard_name')
							if standard_name in refactors and 'skip' in refactors[standard_name] and refactors[standard_name]['skip']: continue
							if section in ccpp_builtin_vars: continue
							units = meta[section].get('units')
							type = meta[section].get('type')
							if '_type' in type:
								if type in ('GFS_control_type', 'GFS_interstitial_type'): continue
								type = f'type({type})'
							kind = meta[section].get('kind')
							if kind != None:
								if 'len=' in kind: kind = kind.split('=')[1]
								if type == 'character' and kind == '*': kind = default_str_len
								type = f'{type}({kind})'
							intent = meta[section].get('intent')
							dimensions = meta[section].get('dimensions')[1:-1].split(',')
							if len(dimensions) == 1 and dimensions[0] == '': dimensions = []
							if standard_name not in var_list:
								var_list[standard_name] = {
									'name': section if standard_name not in refactors else refactors[standard_name]['name'],
									'units': units,
									'type': type,
									'dims': dimensions,
									'intent': intent,
									'pointer': True if standard_name in refactors and 'pointer' in refactors[standard_name] else False,
									'functions': {
										scheme.text.strip(): intent
									}
								}
							else:
								# Check same standard name but different name.
								if var_list[standard_name]['name'] != section and standard_name not in refactors:
									print(f'[Error]: Standard name {standard_name} has different names ({section} vs {var_list[standard_name]["name"]}) in scheme {scheme.text.strip()}.')
									exit(1)
								elif standard_name not in ('errmsg'):
									var_list[standard_name]['functions'][scheme.text.strip()] = intent
	# Check duplicate names.
	for standard_name1, var1 in var_list.items():
		count = 0
		for standard_name2, var2 in var_list.items():
			if var1['name'] == var2['name'] and standard_name1 != standard_name2:
				print(f'[Error]: Variable {var1["name"]} has duplicate names.')
				print(standard_name1, var1)
				print(standard_name2, var2)
				exit(1)
	return var_list	

var_list = parse_suite(args.suite_xml, args.physics_root)

def gen_code(var_list):
	var_statements = []
	for standard_name, var in var_list.items():
		line = var['type']
		if len(var['dims']) == 0:
			line += f' {var["name"]}'
		elif var['pointer']:
			line += f', pointer :: {var["name"]}({",".join([":" for _ in range(len(var["dims"]))])})'
		else:
			line += f', allocatable :: {var["name"]}({",".join([":" for _ in range(len(var["dims"]))])})'
		line += f' ! {standard_name} [{var["units"]}]'
		if len(var['dims']) > 0:
			line += f' ({",".join(var["dims"])})'
		line += f' {{{",".join([f"{key}:{val}" for key, val in var["functions"].items()])}}}'
		var_statements.append('    ' + line)
	var_statements = '\n'.join(sorted(var_statements))
	code = f'''
!!! This code is automatically generated, so do not modify it directly !!!

module gmcore_typedefs_mod

  use machine, only: kind_phys

  implicit none

  private

  public gmcore_ccpp_type

  type gmcore_ccpp_type
{var_statements}
  end type gmcore_ccpp_type

end module gmcore_typedefs_mod'''
	return code

code = gen_code(var_list)
open(args.output_code, 'w').write(code)