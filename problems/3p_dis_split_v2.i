nx = 100     # number of elements in x
ny = 120     # number of elements in y
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
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

R = 1  # Universal gas constant
T = 1 # Temperature in Kelvin
beta = 1e-3*R*T
delta = 0

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
        combinatorial_geometry = 'y > 100/${ny}'
        input = 2d
    []
[]

[Variables]
    # polymer volume fraction
    [c1]
        order = FIRST
        family = LAGRANGE
        scaling = 1e-10
    []
    # Chemical potential (nJ/mol)
    [w1]
        order = FIRST
        family = LAGRANGE
        scaling = 1e10
    []
    # polymer volume fraction
    [c2]
        order = FIRST
        family = LAGRANGE
        scaling = 1e-10
    []
    # Chemical potential (nJ/mol)
    [w2]
        order = FIRST
        family = LAGRANGE
        scaling = 1e10
    []
[]

[ICs]
    # [c1]
    #     type = SolutionIC
    #     from_variable = 'c1_rescale'
    #     solution_uo = 2phase
    #     variable = c1
    #     block = 0
    # []
    # [c2]
    #     type = SolutionIC
    #     from_variable = 'c2_rescale'
    #     solution_uo = 2phase
    #     variable = c2
    #     block = 0
    # []
    [c1]
        type = SolutionIC
        from_variable = 'c'
        solution_uo = 2phase
        variable = c1
        block = 0
    []
    [c2]
        type = CoupledValueFunctionIC
        function = c_2phase
        variable = c2
        v = c1
        block = 0
    []
    [top_c1]
        type = ConstantIC
        value = ${delta}
        variable = c1
        block = 1
    []
    [top_c2]
        type = ConstantIC
        value = ${delta}
        variable = c2
        block = 1
    []
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
    #     v = 'c1 c2'
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
#   [2phase]
#     type = SolutionUserObject
#     mesh = 'output/2phase_copy_0.4.e'
#     system_variables = 'c1_rescale c2_rescale'
#     timestep = LATEST
#   []
  [2phase]
    type = SolutionUserObject
    mesh = 'output/2phase.e'
    system_variables = 'c'
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
    [f_density_0]
        order = CONSTANT
        family = MONOMIAL
        block = 0
    []
    [f_density_1]
        order = CONSTANT
        family = MONOMIAL
        block = 1
    []
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    [f_int_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [f_int_density_0]
        order = CONSTANT
        family = MONOMIAL
        block = 0
    []
    [f_int_density_1]
        order = CONSTANT
        family = MONOMIAL
        block = 1
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
    [coupled_parsed1_0]
        type = SplitCHParsed
        variable = c1
        coupled_variables = 'c2'
        f_name = f_mix_0
        kappa_name = kappa
        w = w1
        block = 0
    []
    [coupled_parsed1_1]
        type = SplitCHParsed
        variable = c1
        coupled_variables = 'c2'
        f_name = f_mix_1
        kappa_name = kappa
        w = w1
        block = 1
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
    [coupled_parsed2_0]
        type = SplitCHParsed
        variable = c2
        coupled_variables = 'c1'
        f_name = f_mix_0
        kappa_name = kappa
        w = w2
        block = 0
    []
    [coupled_parsed2_1]
        type = SplitCHParsed
        variable = c2
        coupled_variables = 'c1'
        f_name = f_mix_1
        kappa_name = kappa
        w = w2
        block = 1
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density_0]
        type = TotalFreeEnergy
        variable = f_density_0
        f_name = 'f_tot_0'
        kappa_names = 'kappa kappa'
        interfacial_vars = 'c1 c2'
        block = 0
    []
    [f_density_1]
        type = TotalFreeEnergy
        variable = f_density_1
        f_name = 'f_tot_1'
        kappa_names = 'kappa kappa'
        interfacial_vars = 'c1 c2'
        block = 1
    []
    # calculate c3
    [c3]
        type = ParsedAux
        variable = c3
        coupled_variables = 'c1 c2'
        expression = '1 - c1 - c2'
    []
    # calculate interfacial energy density
    [f_int_density_0]
        type = ParsedAux
        variable = f_int_density_0
        coupled_variables = 'f_density_0'
        material_properties = 'f_tot_0'
        expression = 'f_density_0 - f_tot_0'
        block = 0
    []
    [f_int_density_1]
        type = ParsedAux
        variable = f_int_density_1
        coupled_variables = 'f_density_1'
        material_properties = 'f_tot_1'
        expression = 'f_density_1 - f_tot_1'
        block = 1
    []
[]

# [BCs]
#     [top1]
#         type = DirichletBC
#         variable = c1
#         boundary = 2
#         # block = 1
#         value = ${delta}
#     []
#     [top2]
#         type = DirichletBC
#         variable = c2
#         boundary = 2
#         # block = 1
#         value = ${delta}
#     []
# []

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
        coupled_variables = 'c1 c2'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        expression = '(M*(1-c1)^2)/s'
        # derivative_order = 2
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        expression = '(M*(1-c2)^2)/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy_0]
        type = DerivativeParsedMaterial
        property_name = f_mix_0
        coupled_variables = 'c1 c2'
        constant_names = 'R      T       chi12      chi13       chi23     N1        N2      N3       s     beta'
        constant_expressions = '${R}    ${T}    ${chi12}    ${chi13}    ${chi23}    ${N1}   ${N2}   ${N3}    ${s}    ${beta}'
        expression = 'if(1-c1-c2>0, s*(R*T*(c1*log(c1)/N1 +
                     c2*log(c2)/N2 + (1 - c1 - c2)*log(1 - c1 - c2)/N3 +
                     chi12*c1*c2 + chi13*c1*(1 - c1 - c2) + chi23*c2*(1 - c1 - c2) +
                    beta*(1/c1 + 1/c2 + 1/(1 - c1 - c2)))), s*(R*T*(c1*log(c1)/N1 + c2*log(c2)/N2 + chi12*c1*c2)))'
        derivative_order = 2
        block = 0
    []
    [mixing_energy_1]
        type = DerivativeParsedMaterial
        property_name = f_mix_1
        coupled_variables = 'c1 c2'
        constant_names = 'R      T       chi12      chi13       chi23     N1        N2      N3       s     beta'
        constant_expressions = '${R}    ${T}    ${chi12}    ${chi13}    ${chi23}    ${N1}   ${N2}   ${N3}    ${s}    ${beta}'
        expression = 'if(c1>0, if(c2>0, s*(R*T*(c1*log(c1)/N1 +
                     c2*log(c2)/N2 + (1 - c1 - c2)*log(1 - c1 - c2)/N3 +
                     chi12*c1*c2 + chi13*c1*(1 - c1 - c2) + chi23*c2*(1 - c1 - c2) +
                    beta*(1/c1 + 1/c2 + 1/(1 - c1 - c2)))), 0), 0)'
        derivative_order = 2
        block = 1
    []
    # Total free energy
    # Sum of all the parts
    [free_energy_0]
        type = DerivativeSumMaterial
        property_name = f_tot_0
        coupled_variables = 'c1 c2'
        sum_materials = 'f_mix_0'
        derivative_order = 2
        block = 0
    []
    [free_energy_1]
        type = DerivativeSumMaterial
        property_name = f_tot_1
        coupled_variables = 'c1 c2'
        sum_materials = 'f_mix_1'
        derivative_order = 2
        block = 1
    []
[]

[Preconditioning]
    [coupled]
      type = SMP
      full = true
    []
[]

[Postprocessors]
    # Calculate total free energy at each timestep
    [total_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_density
        execute_on = 'initial timestep_end'
    []
    [interfacial_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_int_density
        execute_on = 'initial timestep_end'
    []
    [./elapsed]
        type = PerfGraphData
        section_name = "Root"
        data_type = total
    [../]
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 30
    nl_max_its = 100
    nl_abs_tol = 1e-6

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-10
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e0 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'c1 c2; w1 w2'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/3p_dis_split_v2_t1
        time_step_interval = 1
        execute_on = 'TIMESTEP_BEGIN INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3p_dis_split_v2_t1
    []
[]

[Debug]
  show_var_residual_norms = true
[]