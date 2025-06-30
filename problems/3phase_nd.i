# nx = 400     # number of elements per side
# ny = 300      # number of elements per side
# dx = 4       # ND size of the side
# dy = 3        # ND size of the side
# c = 0.05     # solvent density
M = 1.0       # Initial mobility, depends on swell ratio
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
s = 1e-0    # Scaling factor
    
# 1 - drug
# 2 - polymer
# 3 - solvent   
chi12 = 0.42   # Flory-Huggins parameter
chi13 = 0.42   # Flory-Huggins parameter
chi23 = 1.21   # Flory-Huggins parameter
N1 = 10        # Degree of polymerisation
N2 = 100        # Degree of polymerisation
N3 = 1        # Degree of polymerisation
beta = 1e-3
delta = 0.025


[Mesh]
  [main]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 400
    ny = 300
    xmax = 400
    ymax = 300
  []
  [top_block]
    type = ParsedSubdomainMeshGenerator
    input = main
    block_id = 1
    combinatorial_geometry = 'y > 100'
  []
[]

[Variables]
    # polymer volume fraction
    [c1]
        order = FIRST
        family = LAGRANGE
        # scaling = 1e-30
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
        # scaling = 1e-30
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
        from_variable = 'c'
        solution_uo = 2phase
        variable = c1
        block = 0
    []
    [c2]
        type = SolutionIC
        from_variable = 'c2'
        solution_uo = 2phase
        variable = c2
        block = 0
    []
    [top_c1]
        type = RandomIC
        variable = c1
        seed = 123
        min = '${fparse delta*0.5}'
        max = '${fparse delta*1.5}'
        block = 1
    []
    [top_c2]
        type = RandomIC
        variable = c2
        seed = 12
        min = '${fparse delta*0.5}'
        max = '${fparse delta*1.5}'
        block = 1
    []
[]

[UserObjects]
    [2phase]
      type = SolutionUserObject
      mesh = 'output/2phase_nd_2.e'
      system_variables = 'c c2'
      timestep = LATEST
    []
    [./normal_noise]
        type = ConservedNormalNoise
    [../]
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
    # adding nonlocal term to the energy
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
    [./conserved_langevin1]
        type = ConservedLangevinNoise
        amplitude = 0.02
        variable = c1
        noise = normal_noise
    []
    [./conserved_langevin2]
        type = ConservedLangevinNoise
        amplitude = 0.02
        variable = c2
        noise = normal_noise
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
    # calculate interfacial energy density
    [f_int_density]
        type = ParsedAux
        variable = f_int_density
        coupled_variables = 'f_density'
        material_properties = 'f_tot'
        expression = 'f_density - f_tot'
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa'
        prop_values = '${fparse k*s}'
    []
    # [mat]
    #     type = GenericFunctionMaterial
    #     prop_names = 'kappa M1 M2'
    #     prop_values = '${fparse k*s}   ${fparse M/s}   ${fparse M/s}'
    # []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c2'
        constant_names = 'M        s'
        constant_expressions = '${M}   ${s}'
        # expression = 'if(1-c1-c2>0, M*c1^2*(1-c1-c2)^2/s, 0)'
        expression = '(M*10^(5*(1-c1-c2)-1)/1.58)/s'
        # expression = 'if (c1>0, if(c1<1, (M)/s, 0), 0)'
        # derivative_order = 2
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2'
        constant_names = 'M        s'
        constant_expressions = '${M}   ${s}'
        expression = '(M*10^(15*(1-c1-c2)-3)/(0.79))/s'
        # expression = 'if(1-c1-c2>0, M*c2^2*(1-c1-c2)^2/s, 0)'
        # expression = '(M)/s'
        # expression = 'if (c1>0, if(c1<1, (M)/s, 0), 0)'
        # derivative_order = 2
    []
    # mixing energy based on    
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c2'
        constant_names =        'chi12      chi13       chi23     N1        N2      N3       s     beta'
        constant_expressions = '${chi12}    ${chi13}    ${chi23}    ${N1}   ${N2}   ${N3}    ${s}    ${beta}'
        expression = 'if(1-c1-c2>0, s*((c1*log(c1)/N1 + c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2) + beta*(1/c1^2 + 1/c2^2 + 1/(1-c1-c2)^2))), s*((c1*log(c1)/N1 + c2*log(c2)/N2 + chi12*c1*c2 + beta*(1/c1^2 + 1/c2^2))))'
        # expression = 'if(c2>0, if(c1>0, if(1-c1-c2>0, s*((c1*log(c1)/N1 + c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2) + beta*(1/c1^2 + 1/c2^2 + 1/(1-c1-c2)^2))), s*((c1*log(c1)/N1 + c2*log(c2)/N2 + chi12*c1*c2 + beta*(1/c1^2 + 1/c2^2)))), if(1-c2>0, s*(((1-c2)*log(1-c2)/N3 + c2*log(c2)/N2 + chi23*c2*(1-c2) + beta*(1/(1-c2)^2 + 1/c2^2))), s*((c2*log(c2)/N2 + beta*(1/c2^2))))), if(c1>0, if(1-c1>0, s*((c1*log(c1)/N1 + (1-c1)*log(1-c1)/N3 + chi13*c1*(1-c1) + beta*(1/c1^2 + 1/(1-c1)^2))), s*(((1-c1)*log(1-c1)/N3 + beta*(1/(1-c1)^2)))), 0))'
        derivative_order = 2
    []
    # beta penalty term
    # [beta_penalty]
    #     type = DerivativeParsedMaterial
    #     property_name = 'f_beta'
    #     coupled_variables = 'c1 c2'
    #     constant_names = 'beta'
    #     constant_expressions = '${beta}'
    #     expression = 'if (c1>0, if(c2>0, if(1-c1-c2>0, beta*(1/c1 + 1/c2 + 1/(1 - c1 - c2)), beta*(1/c1 + 1/c2)),
    #                 if(1-c1-c2>0, beta*(1/c1 + 1/(1 - c1 - c2)), beta*(1/c1))),
    #                 if(c2>0, if(1-c1-c2>0, beta*(1/c2 + 1/(1 - c1 - c2)), beta*(1/c2)),
    #                 if(1-c1-c2>0, beta*(1/(1-c1-c2)), 0)))'
    #     derivative_order = 2
    # []
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

    # petsc_options = '-pc_svd_monitor -ksp_view'
    # petsc_options = '-ksp_converged_reason -snes_converged_reason'
    petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31                  preonly      ilu          1'

    petsc_options_iname = '-pc_type -ksp_type'
    petsc_options_value = 'gamg      gmres'

    # line_search = 'basic'

    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'lu'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 100
    nl_max_its = 50
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-8
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e2 # seconds

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
        file_base = output/3phase_wp
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3phase_wp
    []
[]

[Debug]
    show_var_residual_norms = true
[]