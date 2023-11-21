dt_min = 1
dt_max = 1728000 ## 20days. 1 day=86400
starting_time = 0
production_end_time = 9.461e+8 ## 30 years

injection_rate = -8000 # kg/s/m, negative because injection as a source
production_rate = 8000 # kg/s/m

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -1000
    xmax = 0
    nx = 500
  []
  [injection_node]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = injection_node
    coord = '0 0 0'
  []
  [production_node]
    input = injection_node
    type = ExtraNodesetGenerator
    new_boundary = production_node
    coord = '-1000 0 0'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  multiply_by_density = true
  gravity = '-9.81 0 0'
  porepressure = porepressure
  temperature = temperature
[]

[Variables]
  [f0]
    # initial_condition = 0.002285946
    # scaling = 1e+2
  []
  [f1]
    # initial_condition = 0.0035252
    # scaling = 1e+2
  []
  [f2]
    # initial_condition = 1.3741E-05
    # scaling = 1e+2
  []
  [f3]
    # initial_condition =
    # scaling = 1e+2
  []
  [f4]
    # initial_condition =
    # scaling = 1e+2
  []
  [porepressure]
    scaling = 1E-5
    # initial_condition = 2E6
  []
  [temperature]
    # initial_condition = 50
    scaling = 1e-13 # fluid enthalpy is roughly 1E6
  []
[]

[BCs]
  [injection_temperature]
    type = MatchedValueBC
    variable = temperature
    v = injection_temperature
    boundary = injection_node
  []
  [injection_p]
    type = FunctionDirichletBC
    variable = porepressure
    boundary = right
    function = hydrostatic
  []
[]

[ICs]
  [./pressure_ic]
    type = FunctionIC
    variable = porepressure
    function = hydrostatic
  [../]
  [./temperature_ic]
    type = FunctionIC
    variable = temperature
    function = temperature_ini
  [../]
[]

[DiracKernels]
  [inject_Ca]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_Ca
    point_file = injection.bh
    variable = f0
  []
  [inject_Cl]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_Cl
    point_file = injection.bh
    variable = f1
  []
  [inject_SO4]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_SO4
    point_file = injection.bh
    variable = f2
  []
  [inject_H]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_H
    point_file = injection.bh
    variable = f3
  []
  [inject_HCO3]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_HCO3
    point_file = injection.bh
    variable = f4
  []
  [inject_H2O]
    type = PorousFlowPolyLineSink
    SumQuantityUO = injected_mass
    fluxes = ${injection_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    multiplying_var = injection_rate_massfrac_H2O
    point_file = injection.bh
    variable = porepressure
  []

  [produce_Ca]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_Ca
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 0
    point_file = production.bh
    variable = f0
  []
  [produce_Cl]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_Cl
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 1
    point_file = production.bh
    variable = f1
  []
  [produce_SO4]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_SO4
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 2
    point_file = production.bh
    variable = f2
  []
  [produce_H]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_H
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 3
    point_file = production.bh
    variable = f3
  []
  [produce_HCO3]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_HCO3
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 4
    point_file = production.bh
    variable = f4
  []
  [produce_H2O]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_H2O
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    mass_fraction_component = 5
    point_file = production.bh
    variable = porepressure
  []

  [produce_heat]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_heat
    fluxes = ${production_rate}
    p_or_t_vals = 0.0
    line_length = 1.0
    point_file = production.bh
    variable = temperature
    use_enthalpy = true
  []
[]

[UserObjects]
  [injected_mass]
    type = PorousFlowSumQuantity
  []
  [produced_mass_Ca]
    type = PorousFlowSumQuantity
  []
  [produced_mass_Cl]
    type = PorousFlowSumQuantity
  []
  [produced_mass_SO4]
    type = PorousFlowSumQuantity
  []
  [produced_mass_H]
    type = PorousFlowSumQuantity
  []
  [produced_mass_HCO3]
    type = PorousFlowSumQuantity
  []
  [produced_mass_H2O]
    type = PorousFlowSumQuantity
  []
  [produced_heat]
    type = PorousFlowSumQuantity
  []
[]

[Postprocessors]
  [dt]
    type = TimestepSize
    execute_on = TIMESTEP_BEGIN
  []
  [tot_kg_injected_this_timestep]
    type = PorousFlowPlotQuantity
    uo = injected_mass
  []
  [kg_Ca_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_Ca
  []
  [kg_Cl_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_Cl
  []
  [kg_SO4_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_SO4
  []
  [kg_H_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_H
  []
  [kg_HCO3_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_HCO3
  []
  [kg_H2O_produced_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_mass_H2O
  []

  [mole_rate_Ca_produced]
    type = FunctionValuePostprocessor
    function = moles_Ca
    indirect_dependencies = 'kg_Ca_produced_this_timestep dt'
  []
  [mole_rate_Cl_produced]
    type = FunctionValuePostprocessor
    function = moles_Cl
    indirect_dependencies = 'kg_Cl_produced_this_timestep dt'
  []
  [mole_rate_SO4_produced]
    type = FunctionValuePostprocessor
    function = moles_SO4
    indirect_dependencies = 'kg_SO4_produced_this_timestep dt'
  []
  [mole_rate_H_produced]
    type = FunctionValuePostprocessor
    function = moles_H
    indirect_dependencies = 'kg_H_produced_this_timestep dt'
  []
  [mole_rate_HCO3_produced]
    type = FunctionValuePostprocessor
    function = moles_HCO3
    indirect_dependencies = 'kg_HCO3_produced_this_timestep dt'
  []
  [mole_rate_H2O_produced]
    type = FunctionValuePostprocessor
    function = moles_H2O
    indirect_dependencies = 'kg_H2O_produced_this_timestep dt'
  []

  [heat_joules_extracted_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_heat
  []

  [production_T]
    type = PointValue
    point = '-1000 0 0'
    variable = temperature
  []
  [production_p]
    type = PointValue
    point = '-1000 0 0'
    variable = porepressure
  []
[]

[Functions]
  [moles_Ca]
    type = ParsedFunction
    symbol_names = 'kg_Ca dt'
    symbol_values = 'kg_Ca_produced_this_timestep dt'
    expression = 'kg_Ca * 1000 / 40.078 / dt'
  []
  [moles_Cl]
    type = ParsedFunction
    symbol_names = 'kg_Cl dt'
    symbol_values = 'kg_Cl_produced_this_timestep dt'
    expression = 'kg_Cl * 1000 / 35.453 / dt'
  []
  [moles_SO4]
    type = ParsedFunction
    symbol_names = 'kg_SO4 dt'
    symbol_values = 'kg_SO4_produced_this_timestep dt'
    expression = 'kg_SO4 * 1000 / 96.06 / dt'
  []
  [moles_H]
    type = ParsedFunction
    symbol_names = 'kg_H dt'
    symbol_values = 'kg_H_produced_this_timestep dt'
    expression = 'kg_H * 1000 / 1.00784 / dt'
  []
  [moles_HCO3]
    type = ParsedFunction
    symbol_names = 'kg_HCO3 dt'
    symbol_values = 'kg_HCO3_produced_this_timestep dt'
    expression = 'kg_HCO3 * 1000 / 61.0168 / dt'
  []
  [moles_H2O]
    type = ParsedFunction
    symbol_names = 'kg_H2O dt'
    symbol_values = 'kg_H2O_produced_this_timestep dt'
    expression = 'kg_H2O * 1000 / 18.0152 / dt'
  []

  [./temperature_ini]
    type = ParsedFunction
    expression = '20-0.03*x'
  [../]
  [./hydrostatic]
    type = ParsedFunction
    expression = '1e5-1000*9.81*x' ##density should be close to real value!
  [../]
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0
    bulk_modulus = 2E10
    viscosity = 1E-3
    density0 = 1000
    cv = 4000.0
    cp = 4000.0
  []
[]

[PorousFlowFullySaturated]
  coupling_type = ThermoHydro
  # porepressure = porepressure
  # temperature = temperature
  mass_fraction_vars = 'f0 f1 f2 f3 f4'
  save_component_rate_in = 'rate_Ca rate_Cl rate_SO4 rate_H rate_HCO3 rate_H2O' # change in kg at every node / dt
  fp = the_simple_fluid
  temperature_unit = Celsius
[]

[AuxVariables]
  [injection_temperature]
    initial_condition = 90
  []
  [injection_rate_massfrac_Ca]
    initial_condition = 2.66E-05
  []
  [injection_rate_massfrac_Cl]
    initial_condition = 3.08E-05
  []
  [injection_rate_massfrac_SO4]
    initial_condition = 2.25E-05
  []
  [injection_rate_massfrac_H]
    initial_condition = 8.91251E-9
  []
  [injection_rate_massfrac_HCO3]
    initial_condition = 1.22E-05
  []
  [injection_rate_massfrac_H2O]
    initial_condition = 0.994175112
  []

  [rate_Ca]
  []
  [rate_Cl]
  []
  [rate_SO4]
  []
  [rate_H]
  []
  [rate_HCO3]
  []
  [rate_H2O]
  []

  [./density]
   family = MONOMIAL
   order = FIRST
  [../]
[]

[AuxKernels]
  [./fluid_density]
   type = PorousFlowPropertyAux
   property = density
   variable = density
  [../]
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst # this simulation has no porosity changes from dissolution
    porosity = 1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-1 0 0   0 1E-1 0   0 0 1E-1'
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.5 0 0  0 2.5 0  0 0 2.5'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1000.0
  []
[]

[Preconditioning]
  active = 'typically_efficient'

  [typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
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
  [exodus]
   type = Exodus
   # interval = 5
  []
  [csv]
  type = CSV
  []
   [./console]
    type = Console
    # interval = 500
    # max_rows = 1
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[MultiApps]
  [react]
    type = TransientMultiApp
    input_files = aquifer_geochemistry.i
    clone_parent_mesh = true
    # execute_on = 'timestep_end'
    # catch_up = true
  []
[]

[Transfers]
  [changes_due_to_flow]
    type = MultiAppCopyTransfer
    source_variable = 'rate_H2O rate_Ca rate_Cl rate_SO4 rate_H rate_HCO3 temperature'
    variable = 'pf_rate_H2O pf_rate_Ca pf_rate_Cl pf_rate_SO4 pf_rate_H pf_rate_HCO3 temperature'
    to_multi_app = react
  []
  [massfrac_from_geochem]
    type = MultiAppCopyTransfer
    source_variable = 'massfrac_Ca massfrac_Cl massfrac_SO4 massfrac_H massfrac_HCO3'
    variable = 'f0 f1 f2 f3 f4'
    from_multi_app = react
  []
[]
