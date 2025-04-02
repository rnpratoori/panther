n = 100     # number of elements per side
d = 1       # ND size of the side
# D = 1e+0    # actual size of the side in 0.1 mum
# l = 1e-3    # Kuhn statistical length in 0.1 mum
a = 0.3     # type A monomer density
chi = 2.0   # Flory-Huggins parameter
N = 5       # Degree of polymerisation
M = 1       # Initial mobility, depends on swell ratio
s = 1e+4    # Scaling factor
k = 1e-0    # gradient energy coefficient

R = 10  # Universal gas constant
T = 300 # Temperature in Kelvin

[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}
    ymax = ${d}
    # uniform_refine = 2
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
    [pvfIC_1]
        type = RandomIC
        variable = c1
        seed = 123
        min = '${fparse a-0.04}'
        max = '${fparse a+0.04}'
    []
    [pvfIC_2]
        type = RandomIC
        variable = c2
        seed = 12
        min = '${fparse a-0.04}'
        max = '${fparse a+0.04}'
    []
[]

[AuxVariables]
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[Kernels]
    [w1_dot]
        type = CoupledTimeDerivative
        variable = w1
        v = c1
    []
    # adding nonlocal term to the energy
    [coupled_res1]
        type = SplitCHWRes
        variable = w1
        mob_name = M
    []
    [coupled_parsed1]
        type = SplitCHParsed
        variable = c1
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    [w2_dot]
        type = CoupledTimeDerivative
        variable = w2
        v = c2
    []
    # adding nonlocal term to the energy
    [coupled_res2]
        type = SplitCHWRes
        variable = w2
        mob_name = M
    []
    [coupled_parsed2]
        type = SplitCHParsed
        variable = c2
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
[]

# [BCs]
#     [Periodic]
#         [all]
#             auto_direction = 'x y'
#         []
#     []
# []

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'M   kappa'
        prop_values = '${fparse M/s} ${fparse k*s}'
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c2'
        constant_names = 'R      T       chi     N       s'
        constant_expressions = '${R}    ${T}    ${chi}  ${N}    ${s}'
        expression = 's*(R*T*(c1*log(c1)/N + c2*log(c2)/N + (1-c1-c2)*log(1-c1-c2) + chi*c1*c2 + chi*c1*(1-c1-c2) + chi*c2*(1-c1-c2)))'
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

    petsc_options_iname = '-pc_type -sub_ksp_type
                            -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm          preonly
                               lu           2'

    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'

    l_tol = 1e-6
    l_abs_tol = 1e-9
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-9

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e-2 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'u w'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output_cl/3phase
        time_step_interval = 1
        execute_on = 'TIMESTEP_END FINAL'
    []
    [csv]
        type = CSV
        file_base = output_cl/3phase
    []
[]