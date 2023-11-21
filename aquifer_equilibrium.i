# Finds the equilibrium free molality of SiO2(aq) when in contact with QuartzLike at 50degC
[UserObjects]
  [definition]
    type = GeochemicalModelDefinition
    database_file = "/projects/moose/modules/geochemistry/database/moose_geochemdb.json"
    basis_species = "H2O              Ca++              Cl-              SO4--              H+                         HCO3-"
    equilibrium_minerals = "Calcite"
  []
[]

[TimeIndependentReactionSolver]
  model_definition = definition
  # geochemistry_reactor_name = reactor
  charge_balance_species = "Cl-"
  # swap_out_of_basis = "SiO2(aq)"
  # swap_into_basis = Calcite
  constraint_species = "H2O  Ca++           Cl-           SO4--         H+           HCO3-"
  constraint_value = "  1.0  0.000663673    0.000868756   0.000234234   8.91251E-09  0.000199944" # amount of QuartzLike is unimportant (provided it is positive).  396.685 is used simply because the other geotes_2D input files use this amount
  constraint_meaning = "kg_solvent_water bulk_composition bulk_composition bulk_composition    activity  bulk_composition"
  constraint_unit = "   kg               moles            moles            moles                dimensionless                   moles"
  temperature = 13
  prevent_precipitation = Calcite
  ramp_max_ionic_strength_initial = 0 # max_ionic_strength in such a simple problem does not need ramping
  add_aux_pH = false # there is no H+ in this system
  precision = 5
[]

[Postprocessors]
  [bulk_moles_Calcite]
    type = PointValue
    point = '0 0 0'
    variable = 'bulk_moles_Calcite'
  []
[]

[Outputs]
  csv = true
[]
