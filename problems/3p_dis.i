nx = 100     # number of elements in x
ny = 101     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.00       # ND size of the side in y
a = 0.3     # type A monomer density
b = 0.3     # type B monomer density
chi12 = 2.0   # Flory-Huggins parameter
chi13 = 0.1   # Flory-Huggins parameter
chi23 = 0.1   # Flory-Huggins parameter
N1 = 5       # Degree of polymerisation
N2 = 5       # Degree of polymerisation
N3 = 1       # Degree of polymerisation
M = 1e-0       # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
k = 1e-1    # gradient energy coefficient

R = 8.314  # Universal gas constant
T = 300 # Temperature in Kelvin
beta = 1e-3*R*T
delta = 1e-3

[Mesh]
    [2d]
        # generate a 2D mesh
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
        # uniform_refine = 2
    []
    # Subdomain for ramp
    [c3_domain]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        combinatorial_geometry = 'y > 0.99'
        input = 2d
    []
[]

[Variables]
    # polymer volume fraction
    [c1]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w1]
        order = FIRST
        family = LAGRANGE
    []
    # polymer volume fraction
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w2]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    [c1]
        type = SolutionIC
        from_variable = 'c1_rescale'
        solution_uo = 2phase
        variable = c1
        block = 0
    []
    [c2]
        type = SolutionIC
        from_variable = 'c2_rescale'
        solution_uo = 2phase
        variable = c2
        block = 0
    []
    # [c2]
    #     type = CoupledValueFunctionIC
    #     function = c_2phase
    #     variable = c2
    #     v = c1
    #     block = 0
    # []
    # [top_c1]
    #     type = ConstantIC
    #     value = ${delta}
    #     variable = c1
    #     block = 1
    # []
    # [top_c2]
    #     type = ConstantIC
    #     value = ${delta}
    #     variable = c2
    #     block = 1
    # []
    # [w1]
    #     type = CoupledValueFunctionIC
    #     function = w1_2phase
    #     variable = w1
    #     v = 'c1 c2'
    #     block = 0
    # []
    # [w2]
    #     type = CoupledValueFunctionIC
    #     function = w2_2phase
    #     variable = w2
    #     v = w1
    #     block = 0
    # []
[]

[Functions]
  [c_2phase]
    type = ParsedFunction
    expression = '1 - x - ${delta}'
  []
  [w1_2phase]
    type = ParsedFunction
    expression = '${R}*${T}*(-1 + y*${chi12} - x*${chi13} + (1-x-y)*${chi13} - y*${chi23} + 1/${N1} + log(x)/${N1} - log(1-x-y)/${N3})*${s}'
  []
  [w2_2phase]
    type = ParsedFunction
    expression = '${R}*${T}*(-1 + x*${chi12} - x*${chi13} + (1-x-y)*${chi23} - y*${chi23} + 1/${N2} + log(y)/${N2} - log(1-x-y)/${N3})*${s}'
  []
[]

[UserObjects]
  [2phase]
    type = SolutionUserObject
    mesh = 'output/2phase_copy.e'
    system_variables = 'c1_rescale c2_rescale'
    timestep = LATEST
  []
[]

[Distributions]
    [Normal_a]
        type = Normal
        mean = ${a}
        standard_deviation = 0.02
    []
    [Normal_b]
        type = Normal
        mean = ${b}
        standard_deviation = 0.02
    []
[]

[AuxVariables]
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [c3]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    [w1_dot]
        type = CoupledTimeDerivative
        variable = w1
        v = c1
    []
    [coupled_res1]
        type = SplitCHWRes
        variable = w1
        mob_name = M1
    []
    [coupled_parsed1]
        type = SplitCHParsed
        variable = c1
        coupled_variables = 'c2'
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    [w2_dot]
        type = CoupledTimeDerivative
        variable = w2
        v = c2
    []
    [coupled_res2]
        type = SplitCHWRes
        variable = w2
        mob_name = M2
    []
    [coupled_parsed2]
        type = SplitCHParsed
        variable = c2
        coupled_variables = 'c1'
        f_name = f_mix
        kappa_name = kappa
        w = w2
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa'
        interfacial_vars = 'c1 c2'
    []
    # calculate c3
    [c3]
        type = ParsedAux
        variable = c3
        coupled_variables = 'c1 c2'
        expression = '1 - c1 - c2'
    []
[]

[BCs]
    [top1]
        type = DirichletBC
        variable = c1
        boundary = 2
        value = ${delta}
    []
    [top2]
        type = DirichletBC
        variable = c2
        boundary = 2
        value = ${delta}
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa'
        prop_values = '${fparse k*s}'
    []
    # [mobility]
    #     type = GenericFunctionMaterial
    #     prop_names = 'M1    M2'
    #     prop_values = '${fparse M/s} ${fparse M/s}'
    # []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        expression = '(M*(c1)^2)/s'
        # derivative_order = 2
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c2'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        expression = '(M*(c2)^2)/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c2'
        constant_names = 'R      T       chi12      chi13       chi23     N1        N2      N3       s     beta'
        constant_expressions = '${R}    ${T}    ${chi12}    ${chi13}    ${chi23}    ${N1}   ${N2}   ${N3}    ${s}    ${beta}'
        expression = 's*(R*T*(c1*log(c1)/N1 + c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2)) + beta*(1/c1 + 1/c2 + 1/(1-c1-c2)))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c2'
        sum_materials = 'f_mix'
        derivative_order = 2
    []
[]

[Postprocessors]
    # Calculate total free energy at each timestep
    [total_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_density
        execute_on = 'initial timestep_end'
    []
    [nodes] # Number of nodes in mesh
        type = NumNodes
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-8
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e4 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    scaling_group_variables = 'c1 c2; w1 w2'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/3phase_Mdecay_3_0.1
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3phase_Mdecay_3_0.1
    []
    # print_linear_residuals = true
[]

# [Debug]
#   show_var_residual_norms = true
# []