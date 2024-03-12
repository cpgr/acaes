# This file is used to model the cushion-bubble formation process,
# and provides initial conditions for subsequent simulations.
#
# Initially, the reservoir is fully saturated.
#
# Air injection rate: 0 - 1 kg/s for the first 30 days,
#				  1 kg/s for the next 30 days,
#                 0 kg/s for another 60 days (Well closure period).
# Air Temp = 150 C to avoid water flashing.
#
# phase 0 : liquid phase
# phase 1 : gas phase
# component 0 : water
# component 1 : air
#
#

rock_density = 2600 # kg/m^3
rock_specific_heat_cap = 920 # J/Kg/K
rock_thermal_cond = 2.54 # W/m/K

biasx = 1.02

[Mesh]
  coord_type = RZ
  rz_coord_axis = Y
  [underburden]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.25
    xmax = 800
    nx = 100
    bias_x = ${biasx}
    ymin = 0
    ymax = 40
    ny = 8
    bias_y = 0.9
  []
  [reservoir]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.25
    xmax = 800
    nx = 100
    bias_x = ${biasx}
    ymin = 40
    ymax = 60
    ny = 20
  []
  [caprock]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.25
    xmax = 800
    nx = 100
    bias_x = ${biasx}
    ymin = 60
    ymax = 100
    ny = 10
    bias_y = 1.1
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'underburden reservoir caprock'
    stitch_boundaries_pairs = 'top bottom; top bottom'
    clear_stitched_boundary_ids = true
    prevent_boundary_ids_overlap = false
  []
  [underburden_id]
    type = ParsedSubdomainMeshGenerator
    input = stitch
    combinatorial_geometry = 'y >= 0 & y <= 40'
    block_id = 1
  []
  [reservoir_id]
    type = ParsedSubdomainMeshGenerator
    input = underburden_id
    combinatorial_geometry = 'y >= 40 & y <= 60'
    block_id = 2
  []
  [overburden_id]
    type = ParsedSubdomainMeshGenerator
    input = reservoir_id
    combinatorial_geometry = 'y >= 60 & y <= 100'
    block_id = 3
  []
  [rename]
    type = RenameBlockGenerator
    input = overburden_id
    old_block = '1 2 3'
    new_block = 'bedrock reservoir caprock'
  []
  [injection_area]
    type = ParsedGenerateSideset
    input = rename
    combinatorial_geometry = 'x=0.25 & y>50 & y<59'
    included_subdomains = reservoir
    new_sideset_name = 'injection_area'
  []
  [skew]
    type = ParsedNodeTransformGenerator
    input = injection_area
    y_function = 'y - 0.03333 * x'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'
  temperature_unit = Celsius
[]

[Variables]
  [pg]
  []
  [zi]
    initial_condition = 0
    scaling = 1e5
  []
  [temperature]
    scaling = 1e-6
  []
[]

[Debug]
  show_var_residual_norms = true
[]

[Kernels]
  [mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pg
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pg
  []
  [mass_air_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = zi
  []
  [flux_air]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zi
  []
  [heat_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [heat_advection]
    type = PorousFlowHeatAdvection
    variable = temperature
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = temperature
  []
[]

[AuxVariables]
  [water_density]
    family = MONOMIAL
    order = CONSTANT
  []
  [gas_density]
    family = MONOMIAL
    order = CONSTANT
  []
  [Yh2o]
    order = CONSTANT
    family = MONOMIAL
  []
  [Xair]
    order = CONSTANT
    family = MONOMIAL
  []
  [pl]
    family = MONOMIAL
    order = CONSTANT
  []
  [sl]
    family = MONOMIAL
    order = CONSTANT
  []
  [sg]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
  [perm]
    family = MONOMIAL
    order = CONSTANT
  []
  [pcap]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [water_density]
    type = PorousFlowPropertyAux
    property = density
    variable = water_density
    execute_on = 'initial timestep_end'
    phase = 0
  []
  [gas_density]
    type = PorousFlowPropertyAux
    property = density
    variable = gas_density
    execute_on = 'initial timestep_end'
    phase = 1
  []
  [Yh2o]
    type = PorousFlowPropertyAux
    property = mass_fraction
    variable = Yh2o
    phase = 1
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  [Xair]
    type = PorousFlowPropertyAux
    property = mass_fraction
    variable = Xair
    phase = 0
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [pl]
    type = PorousFlowPropertyAux
    property = pressure
    variable = pl
    execute_on = 'initial timestep_end'
    phase = 0
  []
  [sl]
    type = PorousFlowPropertyAux
    property = saturation
    variable = sl
    execute_on = 'initial timestep_end'
    phase = 0
  []
  [sg]
    type = PorousFlowPropertyAux
    property = saturation
    variable = sg
    execute_on = 'initial timestep_end'
    phase = 1
  []
  [porosity]
    type = PorousFlowPropertyAux
    property = porosity
    variable = porosity
    execute_on = 'initial'
  []
  [permeability]
    type = PorousFlowPropertyAux
    property = permeability
    variable = perm
    execute_on = 'initial'
  []
  [capillary_pressure]
    type = PorousFlowPropertyAux
    property = capillary_pressure
    variable = pcap
    execute_on = 'initial timestep_end'
  []
[]

[ICs]
  [pg]
    type = FunctionIC
    variable = pg
    function = hydrostatic
    block = 'bedrock reservoir caprock'
  []
  [temperature]
    type = FunctionIC
    variable = temperature
    function = geothermal_gradient
  []
[]

[Functions]
  [hydrostatic]
    type = ParsedFunction
    expression = '1e6 + 9810 * (60-y)'
  []
  [geothermal_gradient]
    type = ParsedFunction
    expression = '30 + 0.025 * (60-y)'
  []
  [injection_rate]
    type = ParsedFunction
    symbol_values = injection_area
    symbol_names = area
    expression = 'if(t>5184000, 0, if(t>2592000, -1/area, -1*t/2592000/area))'
  []
  [inject]
    type = ParsedFunction
    expression = 'if(t>5184000, 0, 1)'
  []
[]

[BCs]
  [bc_for_pressure]
    type = FunctionDirichletBC
    boundary = 'bottom top right'
    variable = pg
    function = hydrostatic
  []
  [bc_for_temperature]
    type = FunctionDirichletBC
    boundary = 'bottom top right'
    variable = temperature
    function = geothermal_gradient
  []
  [air_injection]
    type = PorousFlowSink
    boundary = 'injection_area'
    variable = zi
    fluid_phase = 1
    flux_function = injection_rate
    enable = true
  []
  [air_injection_temp]
    type = DirichletBC
    boundary = injection_area
    variable = temperature
    value = 150
  []
[]

[Controls]
  [inject_heat_on]
    type = ConditionalFunctionEnableControl
    enable_objects = 'BCs::air_injection_temp'
    conditional_function = inject
    implicit = false
    execute_on = 'initial timestep_begin'
  []
[]

[FluidProperties]
  [true_water]
    type = Water97FluidProperties
  []
  [tabulated_water]
    type = TabulatedBicubicFluidProperties
    fp = true_water
  []
  [air]
    type = AirFluidProperties
  []
  [tabulated_air]
    type = TabulatedBicubicFluidProperties
    fp = air
    temperature_min = 275
    pressure_max = 1E7
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = air_coolprop.csv
    num_T = 291
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pg zi temperature'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc_reservoir]
    type = PorousFlowCapillaryPressureVG
    alpha = 7.95e-4 # 1/Pa
    m = 0.457
    sat_lr = 0.15
    pc_max = 1e6
  []
  [pc_caprock]
    type = PorousFlowCapillaryPressureVG
    alpha = 1.5e-4 # 1/Pa
    m = 0.2
    sat_lr = 0.2
    pc_max = 1e6
  []
  [water_air_res]
    type = PorousFlowWaterAir
    capillary_pressure = pc_reservoir
    air_fp = air
    water_fp = tabulated_water
    water97_fp = true_water
  []
  [water_air_cap]
    type = PorousFlowWaterAir
    capillary_pressure = pc_reservoir
    air_fp = air
    water_fp = tabulated_water
    water97_fp = true_water
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability_res]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0 0 1e-12 0 0 0 1e-12'
    block = 'reservoir'
  []
  [permeability_cap]
    type = PorousFlowPermeabilityConst
    permeability = '1e-16 0 0 0 1e-16 0 0 0 1e-16'
    block = 'caprock bedrock'
  []
  [relperm_liq_res]
    type = PorousFlowRelativePermeabilityBC
    lambda = 2
    phase = 0
    s_res = 0.15
    sum_s_res = 0.16
    block = 'reservoir'
  []
  [relperm_gas_res]
    type = PorousFlowRelativePermeabilityBC
    lambda = 2
    phase = 1
    s_res = 0.01
    nw_phase = true
    sum_s_res = 0.16
    block = 'reservoir'
  []
  [relperm_liq_cap]
    type = PorousFlowRelativePermeabilityVG
    m = 0.2
    phase = 0
    s_res = 0.2
    sum_s_res = 0.21
    block = 'caprock bedrock'
  []
  [relperm_gas_cap]
    type = PorousFlowRelativePermeabilityVG
    m = 0.2
    phase = 1
    s_res = 0.01
    sum_s_res = 0.21
    wetting = false
    block = 'caprock bedrock'
  []
  [water_and_air_res]
    type = PorousFlowFluidState
    gas_porepressure = pg
    z = zi
    fluid_state = water_air_res
    capillary_pressure = pc_reservoir
    block = 'reservoir'
  []
  [water_and_air_cap]
    type = PorousFlowFluidState
    gas_porepressure = pg
    z = zi
    fluid_state = water_air_cap
    capillary_pressure = pc_caprock
    block = 'bedrock caprock'
  []
  [internal_energy_res]
    type = PorousFlowMatrixInternalEnergy
    block = 'reservoir'
    density = ${rock_density}
    specific_heat_capacity = ${rock_specific_heat_cap}
  []
  [internal_energy_cap]
    type = PorousFlowMatrixInternalEnergy
    block = 'caprock bedrock'
    density = ${rock_density}
    specific_heat_capacity = ${rock_specific_heat_cap}
  []
  [aq_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    block = 'reservoir'
    dry_thermal_conductivity = '${rock_thermal_cond} 0 0  0 ${rock_thermal_cond} 0  0 0 0'
  []
  [caps_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    block = 'caprock bedrock'
    dry_thermal_conductivity = '${rock_thermal_cond} 0 0  0 ${rock_thermal_cond} 0  0 0 0'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' asm     boomeramg'
    # petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    # petsc_options_value = ' hypre      lu           NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 1.0368e7
  nl_abs_tol = 1e-06
  nl_rel_tol = 1e-05
  nl_max_its = 20
  l_tol = 1e-4
  l_abs_tol = 1e-5
  num_steps = 10
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e3
  []
[]

[Times]
  [days]
    type = TimeIntervalTimes
    start_time = 0
    end_time = 1.0368e7 # 100 days in seconds
    time_interval = 8.64e4
    outputs = 'none'
  []
  [every_1_hours]
    type = TimeIntervalTimes
    start_time = 0
    end_time = 1.0368e7 # 100 days in seconds
    time_interval = 3600
    outputs = 'none'
  []
[]

[Postprocessors]
  [mass_water]
    type = PorousFlowFluidMass
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  [mass_air]
    type = PorousFlowFluidMass
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [injection_area]
    type = AreaPostprocessor
    boundary = injection_area
    execute_on = initial
  []
  [p_bottom]
    type = PointValue
    variable = pg
    point = '0.25 54 0'
    execute_on = 'initial timestep_end'
  []
  [temp_bottom]
    type = PointValue
    variable = temperature
    point = '0.25 54 0'
    execute_on = 'initial timestep_end'
  []
  [sg_bottom]
    type = PointValue
    variable = sg
    point = '0.25 54 0'
    execute_on = 'initial timestep_end'
  []
  [z]
    type = PointValue
    variable = zi
    point = '0.25 54 0'
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = false
  [exodus]
    type = Exodus
    # interval = 10
    # sync_times_object = days
  []
  [csv1]
    type = CSV
    sync_only = true
    sync_times_object = every_1_hours
    time_data = false
  []
[]
