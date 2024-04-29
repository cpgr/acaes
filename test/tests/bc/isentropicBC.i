
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 1
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 10
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp temp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
  []
[]

[Variables]
  [pp]
    initial_condition = 1e6
  []
  [temp]
    initial_condition = 300
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [heat_conduction]
    type = PorousFlowEnergyTimeDerivative
    variable = temp
  []
[]

[FluidProperties]
  [air]
    type = AirFluidProperties
  []
[]

[Materials]
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [air]
    type = PorousFlowSingleComponentFluid
    fp = air
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
  []
  [heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 700
    density = 2500
  []
[]

[BCs]
  [left_p]
    type = PorousFlowSink
    variable = pp
    boundary = left
    flux_function = -1
  []
  [left_T]
    type = ACAESIsentropicCompBC
    variable = temp
    boundary = left
    fp = air
    flux_function = -1
    porepressure_var = pp
  []
[]

[Postprocessors]
  [total_heat_energy]
    type = PorousFlowHeatEnergy
    phase = 0
  []
  [pinj]
    type = SideAverageValue
    variable = pp
    boundary = left
  []
  [Tinj]
    type = SideAverageValue
    variable = temp
    boundary = left
  []
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 20
  end_time = 2e2
[]

[Outputs]
  [csv]
    type = CSV
  []
[]
