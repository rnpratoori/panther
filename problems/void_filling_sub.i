# Sub App - void_filling_sub.i
# Handles rapid void filling kinetics when c3 contacts c4

# Parameters for void filling
filling_rate = 1e4      # Fast reaction rate
activation_threshold = 0.05
completion_threshold = 1e-4  # Stop when c4 < this value

[Mesh]
    # Same mesh as main app
    [2d]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 100
        ny = 102
        xmax = 1.00
        ymax = 1.00
    []
    [c3_domain]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        combinatorial_geometry = 'y > 100/102'
        input = 2d
    []
[]

[Variables]
    # Only the variables that participate in void filling
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    [c4]
        order = FIRST
        family = LAGRANGE
    []
    
    # c1 is passive in this subapp - just receives values
    [c1]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    # Time evolution for c3 (gains from void filling)
    [c3_dot]
        type = TimeDerivative
        variable = c3
    []
    [c3_production]
        type = CoupledForce
        variable = c3
        v = void_filling_rate
    []
    # Small diffusion to smooth the conversion
    [c3_diffusion]
        type = MatDiffusion
        variable = c3
        diffusivity = D_fill
    []
    
    # Time evolution for c4 (consumed by void filling)
    [c4_dot]
        type = TimeDerivative
        variable = c4
    []
    [c4_consumption]
        type = CoupledForce
        variable = c4
        v = void_filling_rate
        coef = -1.0  # Consume c4
    []
    # Small diffusion for c4
    [c4_diffusion]
        type = MatDiffusion
        variable = c4
        diffusivity = D_fill
    []
    
    # c1 doesn't evolve in subapp - just hold its value
    # No kernels for c1
[]

[AuxVariables]
    [void_filling_rate]
        order = FIRST
        family = LAGRANGE
    []
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    [reaction_complete]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    # Calculate void filling reaction rate
    [filling_rate_calc]
        type = ParsedAux
        variable = void_filling_rate
        coupled_variables = 'c3 c4'
        constant_names = 'k_fill threshold'
        constant_expressions = '${filling_rate} ${activation_threshold}'
        # Fast reaction when both phases are present
        expression = 'if(c3 > threshold & c4 > threshold, k_fill * c3 * c4, 0.0)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    # Calculate c2 for completeness
    [c2_calculation]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c1 c3 c4'
        expression = 'max(0.0, 1.0 - c1 - c3 - c4)'
        execute_on = 'timestep_begin timestep_end'
    []
    
    # Check if reaction is complete
    [completion_check]
        type = ParsedAux
        variable = reaction_complete
        coupled_variables = 'c4'
        constant_names = 'threshold'
        constant_expressions = '${completion_threshold}'
        expression = 'if(c4 < threshold, 1.0, 0.0)'
        execute_on = 'timestep_end'
    []
[]

[Materials]
    [diffusion_props]
        type = GenericConstantMaterial
        prop_names = 'D_fill'
        prop_values = '1e-4'  # Small diffusion for smoothing
    []
[]

[BCs]
    # Same boundary conditions as main app
    [top_c3]
        type = DirichletBC
        variable = c3
        boundary = 2
        value = 1.0  # Solvent reservoir
    []
    [top_c4]
        type = DirichletBC
        variable = c4
        boundary = 2
        value = 0.0  # No void at top
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    # Tighter tolerances for accurate void filling
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    petsc_options_iname = '-pc_type -ksp_gmres_restart'
    petsc_options_value = 'lu       31'

    l_tol = 1e-12
    l_abs_tol = 1e-12
    l_max_its = 100
    nl_max_its = 50
    nl_abs_tol = 1e-12

    # Fast, small timesteps for reaction
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1e-8  # Start small
        cutback_factor = 0.5
        growth_factor = 1.2
        optimal_iterations = 8
    []

    # Run until void filling is complete or max time
    end_time = 1e-2  # Short time for fast reaction
    
    # Stop early if reaction is complete everywhere
    steady_state_detection = true
    steady_state_tolerance = 1e-8
[]

[Postprocessors]
    [total_void_remaining]
        type = ElementIntegralVariablePostprocessor
        variable = c4
        execute_on = 'timestep_end'
    []
    [total_solvent]
        type = ElementIntegralVariablePostprocessor
        variable = c3
        execute_on = 'timestep_end'
    []
    [max_reaction_rate]
        type = ElementExtremeValue
        variable = void_filling_rate
        value_type = max
        execute_on = 'timestep_end'
    []
    [avg_completion]
        type = ElementAverageValue
        variable = reaction_complete
        execute_on = 'timestep_end'
    []
    
    # Termination criteria - stop when void is mostly gone
    [void_completion_fraction]
        type = ParsedPostprocessor
        pp_names = 'total_void_remaining'
        # Assuming initial void was ~0.16 (0.4x0.4 region in 1x1 domain)
        expression = 'max(0, 1.0 - total_void_remaining/0.16)'
    []
    
    [sub_step_size]
        type = TimestepSize
    []
[]

# Termination when void filling is essentially complete
[UserObjects]
    [terminator]
        type = Terminator
        expression = 'void_completion_fraction > 0.99'
        execute_on = 'timestep_end'
        message = 'Void filling complete - terminating subapp'
    []
[]

[Outputs]
    exodus = true
    csv = true
    console = true
    # Reduced output frequency for subapp
    interval = 10
[]