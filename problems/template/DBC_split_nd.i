n = 100     # number of elements per side
d = 1       # ND size of the side
D = 1e+0    # actual size of the side in 0.1 mum
l = 1e-3    # Kuhn statistical length in 0.1 mum
a = a_val   # type A monomer density
chi = x_val # Flory-Huggins parameter
N = N_val   # Degree of polymerisation
M_in = 1    # Initial mobility, depends on swell ratio
s = 1e+4    # Scaling factor
seed = r_val# Seed value for IC

[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}  # nd
    ymax = ${d}  # nd
    uniform_refine = 2
[]
  
[Variables]
    # difference in the volume fractions of the 2 phases
    [./u]
        order = FIRST
        family = LAGRANGE
        [./InitialCondition]
            type = RandomIC
            seed = ${seed}
            min = ${fparse 2*(a-0.5)-0.1}
            max = ${fparse 2*(a-0.5)+0.1}
        [../]
    [../]
    # Chemical potential (nJ/mol)
    [./w]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[Functions]
    # A ParsedFunction to define time dependent Mobility
    [./mobility_func]
        type = ParsedFunction
        expression = '${M_in}'
    [../]
[]
  
[AuxVariables]
    # polymer volume fraction
    [./pvf]
    [../]
    # segment number fraction
    [./m]
        order = CONSTANT
        family = MONOMIAL
    [../]
    # used to describe the exponential func to be used in ParsedMaterial
    [./mobility_temp]
    [../]
    # Local free energy density (nJ/mol)
    [./f_density]
        order = CONSTANT
        family = MONOMIAL
    [../]
[]
  
[Kernels]
    [./w_dot]
        type = CoupledTimeDerivative
        variable = w
        v = u
    [../]
    # adding nonlocal term to the energy
    [./coupled_res]
        type = SplitCHPhaseSep
        variable = w
        mob_name = M
        c = u
        m = m
        sigma_name = sigma
        kappa_name = kappa
    [../]
    [./coupled_parsed]
        type = SplitCHParsed
        variable = u
        f_name = f_loc
        kappa_name = kappa
        w = w
    [../]
[]
  
[AuxKernels]
    # calculate polymer volume fraction from difference in volume fractions
    [./pvf]
        type = ParsedAux
        variable = pvf
        coupled_variables = 'u'
        expression = '(u+1)/2'
    [../]
    # assign segment number fraction
    [./m]
        type = ParsedAux
        variable = m
        expression = '2*${a} - 1'
    [../]
    # calculate M
    [./mobility]
        type = FunctionAux
        variable = mobility_temp
        function = 'mobility_func'
        execute_on = timestep_begin
    [../]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [./f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_loc'
        kappa_names = 'kappa'
        interfacial_vars = u
    [../]
[]
  
[BCs]
    [./Periodic]
        [./all]
        auto_direction = 'x y'
        [../]
    [../]
[]
  
[Materials]
    # All properties are non-dimensional
    # In the ND formulation, M is 1.0
    [./mat]
        type = GenericFunctionMaterial
        prop_names  = 'kappa    sigma'
        prop_values = '${fparse ((l^2)/(3*a*(1-a)*chi*D^2))*s}
                        ${fparse ((36*D^2)/((l*a*(1-a)*N)^2))*s}'
    [../]
    [./mobility]
        type = ParsedMaterial
        property_name  = M
        coupled_variables = mobility_temp
        constant_names = 'M0'
        constant_expressions = '1e-0'
        expression = 'if((M0 * mobility_temp / ${s})<=0, 0, (M0 * mobility_temp / ${s}))'
    [../]
    # free energy density function (nJ/mol/nm^2)
    # same as in CHMath
    [./local_energy]
        type = DerivativeParsedMaterial
        property_name = f_loc
        coupled_variables = u
        constant_names = 'W1'
        constant_expressions = '1/4'
        expression = 'W1*${s}*(u^2 - 1)^2'
        derivative_order = 2
    [../]
[]
  
[Postprocessors]
    # Calculate total free energy at each timestep
    [./total_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_density
        execute_on = 'initial timestep_end'
    [../]
[]
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2
  
    # Preconditioning using the additive Schwartz method and LU decomposition
    petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm          preonly         lu      2'
  
    l_tol = 1e-3
    l_abs_tol = 1e-9
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-8
    
    [./TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    [../]
  
    end_time = 1.0 # seconds

    [./Adaptivity]
      coarsen_fraction = 0.1
      refine_fraction = 0.7
      max_h_level = 2
    [../]
[]
  
[Outputs]
    [ex]
        type = Exodus
        file_base = nd_a${a}_x${chi}_N${N}
        time_step_interval = 10
        execute_on = 'TIMESTEP_END FINAL'
    []
    [csv]
        type = CSV
        file_base = nd_a${a}_x${chi}_N${N}_e
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []
  
