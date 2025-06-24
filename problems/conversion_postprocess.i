# conversion_postprocess.i
# This file handles the void-to-solvent conversion logic
# It receives c1, c3, c4, interface_flag from main app
# Returns modified c3_new, c4_new

[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 100
    ny = 102
    xmax = 1.0
    ymax = 1.0
[]

[Variables]
    # Receive these from main app
    [c1]
        order = FIRST
        family = LAGRANGE
    []
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    [c4]
        order = FIRST
        family = LAGRANGE
    []
    [interface_flag]
        order = FIRST
        family = LAGRANGE
    []
    
    # New variables to send back
    [c3_new]
        order = FIRST
        family = LAGRANGE
    []
    [c4_new]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    # Initialize with received values
    [c1_ic]
        type = ConstantIC
        variable = c1
        value = 0.0
    []
    [c3_ic]
        type = ConstantIC
        variable = c3
        value = 0.0
    []
    [c4_ic]
        type = ConstantIC
        variable = c4
        value = 0.0
    []
    [interface_flag_ic]
        type = ConstantIC
        variable = interface_flag
        value = 0.0
    []
    # [c3_new_ic]
    #     type = ConstantIC
    #     variable = c3_new
    #     value = 0.0
    # []
    [c4_new_ic]
        type = ConstantIC
        variable = c4_new
        value = 0.0
    []
[]

[AuxVariables]
    [conversion_rate]
        order = FIRST
        family = LAGRANGE
    []
    [total_converted]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    # Dummy kernels to make variables active
    [c1_dummy]
        type = Diffusion
        variable = c1
    []
    [c3_dummy]
        type = Diffusion
        variable = c3
    []
    [c4_dummy]
        type = Diffusion
        variable = c4
    []
    [interface_flag_dummy]
        type = Diffusion
        variable = interface_flag
    []
    [c3_new_dummy]
        type = Diffusion
        variable = c3_new
    []
    # [c4_new_dummy]
    #     type = Diffusion
    #     variable = c4_new
    # []
[]

[AuxKernels]
    # Calculate conversion rate based on interface detection
    [conversion_rate_calc]
        type = ParsedAux
        variable = conversion_rate
        coupled_variables = 'c3 c4 interface_flag'
        constant_names = 'max_conversion_rate'
        constant_expressions = '0.1'  # Convert up to 10% per timestep
        # Conversion rate proportional to interface strength
        expression = 'if(interface_flag > 0.5, max_conversion_rate * min(c4, c3), 0)'
        execute_on = 'timestep_end'
    []
    
    # Calculate total amount converted
    [total_converted_calc]
        type = ParsedAux
        variable = total_converted
        coupled_variables = 'conversion_rate'
        expression = 'conversion_rate'
        execute_on = 'timestep_end'
    []
    
    # # Calculate new c4 (reduced by conversion)
    # [c4_new_calc]
    #     type = ParsedAux
    #     variable = c4_new
    #     coupled_variables = 'c4 conversion_rate'
    #     expression = 'max(0, c4 - conversion_rate)'
    #     execute_on = 'timestep_end'
    # []
    
    # # Calculate new c3 (increased by conversion)
    # [c3_new_calc]
    #     type = ParsedAux
    #     variable = c3_new
    #     coupled_variables = 'c3 conversion_rate'
    #     expression = 'min(1, c3 + conversion_rate)'
    #     execute_on = 'timestep_end'
    # []
[]

[BCs]
    # No flux boundary conditions
    [c1_bc]
        type = NeumannBC
        variable = c1
        boundary = 'left right top bottom'
        value = 0
    []
    [c3_bc]
        type = NeumannBC
        variable = c3
        boundary = 'left right top bottom'
        value = 0
    []
    [c4_bc]
        type = NeumannBC
        variable = c4
        boundary = 'left right top bottom'
        value = 0
    []
    [interface_flag_bc]
        type = NeumannBC
        variable = interface_flag
        boundary = 'left right top bottom'
        value = 0
    []
    [c3_new_bc]
        type = NeumannBC
        variable = c3_new
        boundary = 'left right top bottom'
        value = 0
    []
    [c4_new_bc]
        type = NeumannBC
        variable = c4_new
        boundary = 'left right top bottom'
        value = 0
    []
[]

[Postprocessors]
    [total_void_before]
        type = ElementIntegralVariablePostprocessor
        variable = c4
        execute_on = 'initial timestep_end'
    []
    [total_void_after]
        type = ElementIntegralVariablePostprocessor
        variable = c4_new
        execute_on = 'initial timestep_end'
    []
    [total_solvent_before]
        type = ElementIntegralVariablePostprocessor
        variable = c3
        execute_on = 'initial timestep_end'
    []
    [total_solvent_after]
        type = ElementIntegralVariablePostprocessor
        variable = c3_new
        execute_on = 'initial timestep_end'
    []
    [total_conversion_amount]
        type = ElementIntegralVariablePostprocessor
        variable = total_converted
        execute_on = 'initial timestep_end'
    []
[]

[Executioner]
    type = Steady  # Single solve per timestep
    solve_type = 'PJFNK'
    
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
    
    nl_max_its = 10
    nl_abs_tol = 1e-12
[]

[Outputs]
    [csv]
        type = CSV
        file_base = output/conversion_postprocess
    []
[]