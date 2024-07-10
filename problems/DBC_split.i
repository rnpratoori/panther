n = 400
d = 8

[Mesh]
    # generate a 2D, 1 mum x 1 mum mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    # nx = 25
    # ny = 25
    xmax = ${d}  # 0.1 mum
    ymax = ${d}  # 0.1 mum
    # uniform_refine = 2
[]
  
[Variables]
    # difference in the volume fractions of the 2 phases
    [./u]
        order = FIRST
        family = LAGRANGE
        [./InitialCondition]
            type = RandomIC
            seed = 123
            min = -0.1
            max =  0.1
        [../]
    [../]
    # Chemical potential (J/mol)
    [./w]
        order = FIRST
        family = LAGRANGE
    [../]
[]
  
[AuxVariables]
    # polymer volume fraction
    [./pvf]
        # order = FIRST
        # family = LAGRANGE
    [../]
    # average u
    [./bar_u]
        order = CONSTANT
        family = MONOMIAL
    [../]
    # Local free energy density (J/mol)
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
        bar_c = bar_u
        alpha_name = alpha
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
    # calculate bar_u
    [./bar_u]
        type = ProjectionAux
        v = u
        variable = bar_u
    [../]
    # # calculate energy density from local and gradient energies (J/mol/mum^2)
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
    # Units of M are m^2 mol / (J s)
    # Units of kappa are J m^2 / mol
    # consider adding a scaling factor in the future
    [./mat]
        type = GenericConstantMaterial
        prop_names  = 'M   kappa    alpha'
        prop_values = '1e-02    5e-03   1e-01'
                        # M*mum_m^2/eV_J
                        # kappa*eV_J*mu_m^2
    [../]
    # free energy density function (J/mol/mum^2)
    # same as in CHMath
    [./local_energy]
        type = DerivativeParsedMaterial
        property_name = f_loc
        coupled_variables = u
        constant_names = 'W1'
        constant_expressions = '1/4'
        expression = 'W1*(u^2 - 1)^2'
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
    petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type'
    petsc_options_value = 'asm      preonly       lu          '
  
    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'
  
    l_tol = 1e-3
    l_max_its = 100
    nl_max_its = 20
    nl_abs_tol = 1e-8
    
    [./TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 7
    [../]
  
    end_time = 1000.0 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'u w'
  
    # [./Adaptivity]
    #   coarsen_fraction = 0.1
    #   refine_fraction = 0.7
    #   max_h_level = 2
    # [../]
[]
  
[Outputs]
    [ex]
        type = Exodus
        file_base = ${n}_${d}
    []
    [csv]
        type = CSV
        file_base = ${n}_${d}_e
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []
  
