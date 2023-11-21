dt_min = 1
dt_max = 1728000 ## 20days. 1 day=86400
starting_time = 0
production_end_time = 9.461e+8 ## 30 years

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -1000
    xmax = 0
    nx = 500
  []
[]

[GlobalParams]
  point = '0 0 0'
  reactor = reactor
[]

[SpatialReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  charge_balance_species = "Cl-"
  constraint_species = "H2O              Ca++              Cl-              SO4--              H+                         HCO3-"
# ASSUME that 1 litre of solution contains:
  constraint_value = "  0.99              0.000663673             0.000868756              0.000234234              2.6136e-06           0.000199944"
  constraint_meaning = "kg_solvent_water bulk_composition  bulk_composition bulk_composition    bulk_composition     bulk_composition"
  constraint_unit = "   kg               moles            moles            moles                moles                   moles"
  initial_temperature = 50
  kinetic_species_name = Calcite
# Per 1 litre (1000cm^3) of aqueous solution (1kg of solvent water), there is 9000cm^3 of Calcite, which means the initial porosity is 0.1.
  kinetic_species_initial_value = 1e-100
  kinetic_species_unit = cm3
  temperature = temperature
  source_species_names = 'H2O                 Ca++                Cl-                SO4--           H+             HCO3-'
  source_species_rates = 'rate_H2O_per_1l     rate_Ca_per_1l    rate_Cl_per_1l rate_SO4_per_1l    rate_H_per_1l  rate_HCO3_per_1l'
  ramp_max_ionic_strength_initial = 0 # max_ionic_strength in such a simple problem does not need ramping
  add_aux_pH = false # there is  H+ in this system
  evaluate_kinetic_rates_always = true # implicit time-marching used for stability
  execute_console_output_on = ''
[]

[UserObjects]
  [rate_Calcite]
    type = GeochemistryKineticRate
    kinetic_species_name = Calcite
    intrinsic_rate_constant = 1.0E-7
    multiply_by_mass = true
    area_quantity = 1
    activation_energy = 72800.0
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = "/projects/moose/modules/geochemistry/database/moose_geochemdb.json"
    basis_species = "H2O SO4-- Ca++ Cl- H+ HCO3-"
    kinetic_minerals = "Calcite"
    kinetic_rate_descriptions = "rate_Calcite"
    remove_all_extrapolated_secondary_species = true
  []
  [nodal_void_volume_uo]
    type = NodalVoidVolume
    porosity = porosity
    execute_on = 'initial timestep_end' # "initial" means this is evaluated properly for the first timestep
  []
[]

[AuxVariables]
  [temperature]
    initial_condition = 50
  []
  [porosity]
    initial_condition = 1
  []
  [nodal_void_volume]
  []
  [free_cm3_calcite]
  []
  [pf_rate_H2O] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_Ca] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_Cl] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_SO4] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_H] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_HCO3] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []

  [rate_H2O_per_1l] # rate per 1 litre of aqueous solution that we consider at each node
  []
  [rate_Ca_per_1l]
  []
  [rate_Cl_per_1l]
  []
  [rate_SO4_per_1l]
  []
  [rate_H_per_1l]
  []
  [rate_HCO3_per_1l]
  []

  [transported_H2O]
  []
  [transported_Ca]
  []
  [transported_Cl]
  []
  [transported_SO4]
  []
  [transported_H]
  []
  [transported_HCO3]
  []
  [transported_mass]
  []

  [massfrac_Ca]
  []
  [massfrac_Cl]
  []
  [massfrac_SO4]
  []
  [massfrac_H2O]
  []
  [massfrac_H]
  []
  [massfrac_HCO3]
  []
[]

[AuxKernels]
  [free_cm3_calcite]
   type = GeochemistryQuantityAux
   variable = free_cm3_calcite
   species = 'Calcite'
   quantity = free_cm3
   execute_on = ' timestep_end'
  []

  [porosity]
    type = ParsedAux
    coupled_variables = free_cm3_calcite
    expression = '1000.0 / (1000.0 + free_cm3_calcite)'
    variable = porosity
    execute_on = 'timestep_end'
  []
  [nodal_void_volume_auxk]
    type = NodalVoidVolumeAux
    variable = nodal_void_volume
    nodal_void_volume_uo = nodal_void_volume_uo
    execute_on = 'initial timestep_end' # "initial" to ensure it is properly evaluated for the first timestep
  []
  [rate_H2O_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_H2O nodal_void_volume'
    variable = rate_H2O_per_1l
# pf_rate = change in kg at every node
# pf_rate * 1000 / molar_mass_in_g_per_mole = change in moles at every node
# pf_rate * 1000 / molar_mass / (nodal_void_volume_in_m^3 * 1000) = change in moles per litre of aqueous solution
    expression = 'pf_rate_H2O / 18.0152 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Ca_per_1l]
    type = ParsedAux
    coupled_variables = 'pf_rate_Ca nodal_void_volume'
    variable = rate_Ca_per_1l
    expression = 'pf_rate_Ca / 40.078 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Cl_per_1l]
    type = ParsedAux
    coupled_variables = 'pf_rate_Cl nodal_void_volume'
    variable = rate_Cl_per_1l
    expression = 'pf_rate_Cl / 35.453 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_SO4_per_1l]
    type = ParsedAux
    coupled_variables = 'pf_rate_SO4 nodal_void_volume'
    variable = rate_SO4_per_1l
    expression = 'pf_rate_SO4 /  96.06 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_H_per_1l]
    type = ParsedAux
    coupled_variables = 'pf_rate_H nodal_void_volume'
    variable = rate_H_per_1l
    expression = 'pf_rate_H / 1.00784 / nodal_void_volume'
    execute_on = ' timestep_begin'
  []
  [rate_HCO3_per_1l]
    type = ParsedAux
    coupled_variables = 'pf_rate_HCO3 nodal_void_volume'
    variable = rate_HCO3_per_1l
    expression = 'pf_rate_HCO3 / 61.0168  / nodal_void_volume'
    execute_on = 'timestep_begin'
  []

  [transported_H2O_auxk]
    type = GeochemistryQuantityAux
    variable = transported_H2O
    species = H2O
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_Ca]
    type = GeochemistryQuantityAux
    variable = transported_Ca
    species = Ca++
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_Cl]
    type = GeochemistryQuantityAux
    variable = transported_Cl
    species = Cl-
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_SO4]
    type = GeochemistryQuantityAux
    variable = transported_SO4
    species = SO4--
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_H]
    type = GeochemistryQuantityAux
    variable = transported_H
    species = H+
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_HCO3]
    type = GeochemistryQuantityAux
    variable = transported_HCO3
    species = HCO3-
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_mass_auxk]
    type = ParsedAux
    coupled_variables = 'transported_H2O transported_Ca transported_Cl transported_SO4 transported_H transported_HCO3'
    variable = transported_mass
    expression = 'transported_H2O * 18.0152 + transported_Ca * 40.078 + transported_Cl * 35.453 + transported_SO4 * 96.06 + transported_H * 1.00784 + transported_HCO3 * 61.0168'
    execute_on = 'timestep_end'
  []

  [massfrac_H2O]
    type = ParsedAux
    coupled_variables = 'transported_H2O transported_mass'
    variable = massfrac_H2O
    expression = 'transported_H2O * 18.0152 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Ca]
    type = ParsedAux
    coupled_variables = 'transported_Ca transported_mass'
    variable = massfrac_Ca
    expression = 'transported_Ca * 40.078 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Cl]
    type = ParsedAux
    coupled_variables = 'transported_Cl transported_mass'
    variable = massfrac_Cl
    expression = 'transported_Cl * 35.453 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_SO4]
    type = ParsedAux
    coupled_variables = 'transported_SO4 transported_mass'
    variable = massfrac_SO4
    expression = 'transported_SO4 * 96.06 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_H]
    type = ParsedAux
    coupled_variables = 'transported_H transported_mass'
    variable = massfrac_H
    expression = 'transported_H * 1.00784  / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_HCO3]
    type = ParsedAux
    coupled_variables = 'transported_HCO3 transported_mass'
    variable = massfrac_HCO3
    expression = 'transported_HCO3 * 61.0168  / transported_mass'
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [cm3_calcite]
    type = PointValue
    variable = free_cm3_calcite
  []
  [porosity]
    type = PointValue
    variable = porosity
  []
  [solution_temperature]
    type = PointValue
    variable = solution_temperature
  []
  [massfrac_H2O]
    type = PointValue
    variable = massfrac_H2O
  []
  [massfrac_Ca]
    type = PointValue
    variable = massfrac_Ca
  []
  [massfrac_Cl]
    type = PointValue
    variable = massfrac_Cl
  []
  [massfrac_SO4]
    type = PointValue
    variable = massfrac_SO4
  []
  [massfrac_H]
    type = PointValue
    variable = massfrac_H
  []
  [massfrac_HCO3]
    type = PointValue
    variable = massfrac_HCO3
  []
[]

[Executioner]
  type = Transient
  start_time = ${starting_time}
  end_time = ${production_end_time}
  dtmax = ${dt_max}
  l_tol = 1e-12
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50
  l_max_its = 10
  solve_type = NEWTON
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = ${dt_min}
    growth_factor = 1.2
  [../]
[]

[Outputs]
  print_linear_residuals = true
  exodus = true
  [csv]
    type = CSV
    execute_on = ' TIMESTEP_END'
  []
[]
