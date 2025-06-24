nx = 100     # number of elements in x
ny = 102     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.00       # ND size of the side in y
M1 = 1       # Initial mobility, depends on swell ratio
M3 = 1e-0    # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn3 = 1e-2  # Cahn number
k4 = ${fparse Cn3^2}    # gradient energy coefficient

beta = 1e-3
delta = 0

# NEW PARAMETERS for void filling
filling_rate = 1e3     # Rate of void to solvent conversion
activation_threshold = 0.05  # Minimum c3 needed to trigger conversion

[Mesh]
    [2d]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
    []
    [c3_domain]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        combinatorial_geometry = 'y > ${nx}/${ny}'
        input = 2d
    []
[]

[Variables]
    [c1]
        order = FIRST
        family = LAGRANGE
    []
    [w1]
        order = FIRST
        family = LAGRANGE
    []
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    [w3]
        order = FIRST
        family = LAGRANGE
    []
    [c4]
        order = FIRST
        family = LAGRANGE
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
    # Initialize void in middle region
    [c4_void_region]
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
    [top_c4]
        type = ConstantIC
        value = 0.0
        variable = c4
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

[Kernels]
    # Standard Cahn-Hilliard kernels for c1
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
        coupled_variables = 'c3 c4'
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    
    # Standard Cahn-Hilliard kernels for c3
    [w3_dot]
        type = CoupledTimeDerivative
        variable = w3
        v = c3
    []
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
        kappa_name = kappa4
        w = w3
    []
    
    # NEW: c3 production from void filling reaction
    [c3_void_reaction]
        type = CoupledForce
        variable = w3
        v = void_to_solvent_rate
    []
    
    # NEW: Time evolution for c4 (void)
    [c4_dot]
        type = TimeDerivative
        variable = c4
    []
    
    # NEW: Void consumption reaction
    [c4_void_consumption]
        type = CoupledForce
        variable = c4
        v = void_to_solvent_rate
        coef = -1.0  # Negative coefficient to consume void
    []
    
    # Optional: Small diffusion for c4 to smooth interfaces
    [c4_diffusion]
        type = MatDiffusion
        variable = c4
        diffusivity = D4
    []
[]

# Add auxiliary variables
[AuxVariables]
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    [void_to_solvent_rate]
        order = FIRST
        family = LAGRANGE
    []
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [f_int_density]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    # Calculate c2 = 1 - c1 - c3 - c4
    [c2_calculation]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c1 c3 c4'
        expression = 'max(0.0, 1.0 - c1 - c3 - c4)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    # NEW: Calculate void filling reaction rate
    [void_reaction_rate]
        type = ParsedAux
        variable = void_to_solvent_rate
        coupled_variables = 'c3 c4'
        constant_names = 'k_fill threshold'
        constant_expressions = '${filling_rate} ${activation_threshold}'
        # Reaction occurs when c3 is present and c4 > 0
        # Rate proportional to both concentrations
        expression = 'if(c3 > threshold & c4 > threshold, k_fill * c3 * c4, 0.0)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa4 kappa4'
        interfacial_vars = 'c1 c3 c4'
    []
    
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
        prop_names = 'kappa kappa4 D4'
        prop_values = '${fparse k*s} ${fparse k4*s} 1e-6'  # Small diffusivity for c4
    []
    
    # Modified mobility for c1 - cannot diffuse into void regions
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c4'
        constant_names = 'M1 s'
        constant_expressions = '${M1} ${s}'
        # Reduce mobility in void regions
        expression = '(M1*16*(c1^2*(1-c1)^2)*(1-c4))/s'
    []
    
    # Enhanced mobility for c3 when interacting with voids
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c3 c4'
        constant_names = 'M3 s invasion_boost'
        constant_expressions = '${M3} ${s} 10.0'
        # Enhanced mobility when solvent contacts void
        expression = 'if(c4 > 0.01 & c3 > 0.01, 
                         (M3*invasion_boost*(16*(c3^2*(1-c3)^2)) + M3*c4)/s,
                         (M3*(16*(c3^2*(1-c3)^2)) + M3*c4)/s)'
    []
    
    # Modified mixing energy including void interactions
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = 'f_mix'
        coupled_variables = 'c1 c3 c4'
        constant_names = 'chi13 void_energy'
        constant_expressions = '2.0 -1.0'  # Negative void_energy favors c3-c4 contact
        # Mixing energy + void interaction term
        expression = 'chi13*c1*c3*(1-c4) + void_energy*c3*c4'
        derivative_order = 2
    []
    
    # Beta penalty term (logarithmic barriers)
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = 'f_beta'
        coupled_variables = 'c1 c3 c4'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'if(c1>1e-6, if(c3>1e-6, if(c4>1e-6, if(1-c1-c3-c4>1e-6, 
                         beta*(log(c1) + log(c3) + log(c4) + log(1-c1-c3-c4)), 
                         beta*(log(c1) + log(c3) + log(c4))), 
                         if(1-c1-c3-c4>1e-6, beta*(log(c1) + log(c3) + log(1-c1-c3-c4)), 
                         beta*(log(c1) + log(c3)))), 
                      if(c4>1e-6, if(1-c1-c3-c4>1e-6, beta*(log(c1) + log(c4) + log(1-c1-c3-c4)), 
                         beta*(log(c1) + log(c4))), 
                         if(1-c1-c3-c4>1e-6, beta*(log(c1) + log(1-c1-c3-c4)), beta*log(c1)))), 0)'
        derivative_order = 2
    []
    
    # Total free energy
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c3 c4'
        sum_materials = 'f_mix f_beta'
        derivative_order = 2
    []
[]

[BCs]
    [top1]
        type = DirichletBC
        variable = c1
        boundary = 2
        value = ${delta}
    []
    [top3]
        type = DirichletBC
        variable = c3
        boundary = 2
        value = ${fparse 1-delta}
    []
    [top4]
        type = DirichletBC
        variable = c4
        boundary = 2
        value = 0.0
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
    nl_max_its = 100
    nl_abs_tol = 1e-10

    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1.0e-10
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e0

    automatic_scaling = true
    scaling_group_variables = 'c1 c3 c4; w1 w3'
[]

[Postprocessors]
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
    [total_void]
        type = ElementIntegralVariablePostprocessor
        variable = c4
        execute_on = 'initial timestep_end'
    []
    [total_solvent]
        type = ElementIntegralVariablePostprocessor
        variable = c3
        execute_on = 'initial timestep_end'
    []
    [total_c1]
        type = ElementIntegralVariablePostprocessor
        variable = c1
        execute_on = 'initial timestep_end'
    []
    # NEW: Monitor reaction rate
    [max_reaction_rate]
        type = ElementExtremeValue
        variable = void_to_solvent_rate
        value_type = max
        execute_on = 'timestep_end'
    []
    [avg_reaction_rate]
        type = ElementAverageValue
        variable = void_to_solvent_rate
        execute_on = 'timestep_end'
    []
    [step_size]
        type = TimestepSize
    []
[]

[Outputs]
    exodus = true
    csv = true
    console = true
    [pgraph]
        type = PerfGraphOutput
        execute_on = 'final'  # print performance info
        level = 1
    []
[]