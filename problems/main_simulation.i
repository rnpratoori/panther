nx = 100     # number of elements in x
ny = 102     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.00       # ND size of the side in y
# a = 0.5     # type A monomer density
M1 = 1       # Initial mobility, depends on swell ratio
# M2 = 1e-0    # Initial mobility, depends on swell ratio
M3 = 1e-0    # Initial mobility, depends on swell ratio
# M4 = 0e-3    # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn3 = 1e-2  # Cahn number
k4 = ${fparse Cn3^2}    # gradient energy coefficient

# chi12 = 1.0   # Flory-Huggins parameter
# chi13 = 10.0   # Flory-Huggins parameter
# chi23 = 10.0   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 100     # Penalty term for void
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
beta = 1e-3
# alpha = 1e10
delta = 0

# Approach 3: Post-processor Method for Void-to-Solvent Conversion
# This approach uses a two-step process:
# 1. Normal Cahn-Hilliard evolution
# 2. Post-processing conversion step

# MAIN INPUT FILE
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
    [c4]  # Void as full variable
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
    [c4_void_region]
        type = FunctionIC
        variable = c4
        function = 'if(x > 0.3 & x < 0.7 & y > 0.3 & y < 0.7, 1.0, 0.0)'
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
    [interface_detector]
        type = LayeredSideDiffusiveFluxAverage
        variable = c3
        direction = x
        num_layers = ${nx}
        boundary = 'left right top bottom'
        diffusivity = M3           # <-- Add this line, or use your actual property name
    []
[]

[AuxVariables]
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    [interface_flag]  # Flag to mark interface regions
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

[Kernels]
    # Standard Cahn-Hilliard for c1
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
    
    # Standard Cahn-Hilliard for c3
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
        kappa_name = kappa
        w = w3
    []
    
    # Simple evolution for c4 (void) - will be modified by post-processor
    [c4_dot]
        type = TimeDerivative
        variable = c4
    []
    [c4_diffusion]
        type = Diffusion
        variable = c4
    []
[]

[AuxKernels]
    # Calculate c2 (implicit component)
    [c2_calculation]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c1 c3 c4'
        expression = 'max(0.0, 1.0 - c1 - c3 - c4)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    # NEW: Mark interface regions where conversion should occur
    [interface_detection]
        type = ParsedAux
        variable = interface_flag
        coupled_variables = 'c3 c4'
        constant_names = 'threshold'
        constant_expressions = '0.05'
        # Flag regions where c3 and c4 coexist above threshold
        expression = 'if(c3 > threshold & c4 > threshold, 1.0, 0.0)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    # Energy density calculations
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa kappa'
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

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa'
        prop_values = '${fparse k*s}'
    []
    
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c4'
        constant_names = 'M1 s'
        constant_expressions = '${M1} ${s}'
        expression = '(M1*16*(c1^2*(1-c1)^2)*(1-c4))/s'
    []
    
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c3 c4'
        constant_names = 'M3 s'
        constant_expressions = '${M3} ${s}'
        expression = '(M3*(16*(c3^2*(1-c3)^2)) + 100*c4)/s'
    []
    
    # Simple mixing energy (replace with your spline if needed)
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = 'f_mix'
        coupled_variables = 'c1 c3 c4'
        constant_names = 'chi13'
        constant_expressions = '2.0'
        expression = 'chi13*c1*c3*(1-c4)'
        derivative_order = 2
    []
    
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = 'f_beta'
        coupled_variables = 'c1 c3 c4'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'if(c1>0, if(c3>0, if(c4>0, if(1-c1-c3-c4>0, 
                         beta*(1/c1 + 1/c3 + 1/c4 + 1/(1-c1-c3-c4)), 
                         beta*(1/c1 + 1/c3 + 1/c4)), 
                         if(1-c1-c3-c4>0, beta*(1/c1 + 1/c3 + 1/(1-c1-c3-c4)), 
                         beta*(1/c1 + 1/c3))), 
                      if(c4>0, if(1-c1-c3-c4>0, beta*(1/c1 + 1/c4 + 1/(1-c1-c3-c4)), 
                         beta*(1/c1 + 1/c4)), 
                         if(1-c1-c3-c4>0, beta*(1/c1 + 1/(1-c1-c3-c4)), beta*(1/c1)))), 0)'
        derivative_order = 2
    []
    
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c3 c4'
        sum_materials = 'f_mix f_beta'
        derivative_order = 2
    []
[]

# NEW: Post-processing conversion using MultiApp
[MultiApps]
    [conversion_app]
        type = TransientMultiApp
        input_files = 'conversion_postprocess.i'
        execute_on = 'timestep_end'
        max_procs_per_app = 1
    []
[]

# NEW: Transfer data to/from conversion app
[Transfers]
    [send_to_conversion]
        type = MultiAppCopyTransfer
        to_multi_app = conversion_app
        source_variable = 'c1 c3 c4 interface_flag'
        variable = 'c1 c3 c4 interface_flag'
    []
    
    [receive_from_conversion]
        type = MultiAppCopyTransfer
        from_multi_app = conversion_app
        source_variable = 'c3_new c4_new'
        variable = 'c3 c4'
    []
[]

[Postprocessors]
    [total_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_density
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
    [interface_area]
        type = ElementIntegralVariablePostprocessor
        variable = interface_flag
        execute_on = 'initial timestep_end'
    []
    [step_size]
        type = TimestepSize
    []
[]

[Preconditioning]
    [coupled]
      type = SMP
      full = true
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

[Outputs]
    [ex]
        type = Exodus
        file_base = output/3p_dis_void_postprocess
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3p_dis_void_postprocess
    []
[]