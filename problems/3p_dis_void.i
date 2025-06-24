nx = 100     # number of elements in x
ny = 102     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.00       # ND size of the side in y
# a = 0.5     # type A monomer density
M1 = 1       # Initial mobility, depends on swell ratio
M3 = 1e-0    # Initial mobility, depends on swell ratio
M4 = 0e-3    # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn3 = 1e-2  # Cahn number
k4 = ${fparse Cn3^2}    # gradient energy coefficient
k_void = 1e2
epsilon = 1e-2

# chi12 = 1.0   # Flory-Huggins parameter
# chi13 = 10.0   # Flory-Huggins parameter
# chi23 = 10.0   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 100     # Penalty term for void
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
beta = 1e-3
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
        combinatorial_geometry = 'y > ${nx}/${ny}'
        input = 2d
        # uniform_refine = 2
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
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w3]
        order = FIRST
        family = LAGRANGE
    []
    [c4]
        order = FIRST
        family = SCALAR
    []
[]

[ICs]
    [c1]
        type = SolutionIC
        from_variable = 'c1'
        solution_uo = 2phase_void
        variable = c1
        block = 0
    []
    [c3]
        type = ConstantIC
        value = ${delta}
        variable = c3
        block = 0
    []
    [c4]
        type = SolutionIC
        from_variable = 'c3'
        solution_uo = 2phase_void
        variable = c4
        block = 0
    []
    [top_c1]
        type = ConstantIC
        value = ${delta}
        variable = c1
        block = 1
    []
    [top_c3]
        type = ConstantIC
        value = ${fparse 1-delta}
        variable = c3
        block = 1
    []
[]

[UserObjects]
    [2phase_void]
      type = SolutionUserObject
      mesh = 'output/2p_void_spline.e'
      system_variables = 'c1 c3'
      timestep = LATEST
    []
[]

[AuxVariables]
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [f_int_density]
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
        mob_name = M1
    []
    [coupled_parsed1]
        type = SplitCHParsed
        variable = c1
        coupled_variables = 'c3 c4'
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    [w3_dot]
        type = CoupledTimeDerivative
        variable = w3
        v = c3
    []
    # adding nonlocal term to the energy
    [coupled_res3]
        type = SplitCHWRes
        variable = w3
        mob_name = M3
    []
    [coupled_parsed3]
        type = SplitCHParsed
        variable = c3
        coupled_variables = 'c1 c4'
        f_name = f_mix
        kappa_name = kappa
        w = w3
    []
[]

[ScalarKernels]
    [dc4]
        type = ODETimeDerivative
        variable = c4
    []
    [void_decay]
        type = ParsedODEKernel
        variable = c4
        constant_names = 'k_void epsilon'
        constant_expressions = '${k_void} ${epsilon}'
        expression = '-k_void * c4 * tanh(1000 * (c3 - 0.5))'
        coupled_variables = 'c3'
    []
    [solvent_gain]
        type = ParsedODEKernel
        variable = c3
        constant_names = 'k_void epsilon'
        constant_expressions = '${k_void} ${epsilon}'
        expression = 'k_void * c4 * tanh(1000 * (c3 - 0.5))'
        coupled_variables = 'c4'
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa    kappa    kappa4'
        interfacial_vars = 'c1  c3  c4'
    []
    # calculate interfacial energy density
    [f_int_density]
        type = ParsedAux
        variable = f_int_density
        coupled_variables = 'f_density'
        material_properties = 'f_tot'
        expression = 'f_density - f_tot'
    []
[]

[BCs]
    [top1]
        type = DirichletBC
        variable = c1
        boundary = 2
        # block = 1
        value = ${delta}
    []
    [top3]
        type = DirichletBC
        variable = c3
        boundary = 2
        # block = 1
        value = ${fparse 1-delta}
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa     kappa4'
        prop_values = '${fparse k*s}    ${fparse k4*s}'
    []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c4'
        constant_names = 'M1     s'
        constant_expressions = '${M1} ${s}'
        expression = '(M1*16*(c1^2*(1-c1)^2)*(1-c4))/s'
        # derivative_order = 2
    []
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c3 c4'
        constant_names = 'M3     s'
        constant_expressions = '${M3} ${s}'
        expression = '(M3*16*(c3^2*(1-c3)^2))/s'
        # derivative_order = 2
    []
    [mobility4]
        type = DerivativeParsedMaterial
        property_name = M4
        coupled_variables = 'c4'
        constant_names = 'M4     s'
        constant_expressions = '${M4} ${s}'
        expression = '(M4*16*c4^2*(1-c4)^2)/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeSpline2Material
        triangle_file = 'triangles.csv'      
        coefficient_file = 'coefficients.csv'      
        property_name = 'f_mix'           
        coupled_variables = 'c1 c3'          
        derivative_order = 2                 
    []
    # beta penalty term
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = 'f_beta'
        coupled_variables = 'c1 c3'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'if (c1>0, if(c3>0, if(1-c1-c3>0, beta*(1/c1 + 1/c3 + 1/(1 - c1 - c3)), beta*(1/c1 + 1/c3)),
                    if(1-c1-c3>0, beta*(1/c1 + 1/(1 - c1 - c3)), beta*(1/c1))),
                    if(c3>0, if(1-c1-c3>0, beta*(1/c3 + 1/(1 - c1 - c3)), beta*(1/c3)),
                    if(1-c1-c3>0, beta*(1/(1-c1-c3)), 0)))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c3'
        sum_materials = 'f_mix f_beta'
        derivative_order = 2
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
    [step_size]
        type = TimestepSize
    []

[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    # petsc_options = '-pc_svd_monitor -ksp_view'
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    # petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    # petsc_options = '-pc_svd_monitor -ksp_view'
    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'svd'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
    nl_max_its = 100
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-6
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e0 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    # off_diagonal_
    scaling_group_variables = 'c1 c3 c4; w1 w3 w4'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/3p_dis_void_spline
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3p_dis_void_spline
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []