n = 200     # number of elements per side
d = 1000    # size of the sample in nm
eps = 10    # interface width in nm


[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    # nx = 25
    # ny = 25
    xmax = ${d}  # nm
    ymax = ${d}  # nm
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
    # Chemical potential (nJ/mol)
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
        bar_c = bar_u
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
    # Units of M are nm^3 / (kg s)
    # Units of kappa are nm^2
    # Units of sigma are / s
    # scaling factor d is multiplied to f_loc, kappa and sigma
    # and divided to 
    # d = 1e-1
    [./mat]
        type = GenericFunctionMaterial
        prop_names  = 'M   kappa    sigma'
        prop_values = '1e-0/1e-1 ${fparse (eps^2)*1e-1}  7.4e-2*1e-1'

    [../]
    # free energy density function (nJ/mol/nm^2)
    # same as in CHMath
    [./local_energy]
        type = DerivativeParsedMaterial
        property_name = f_loc
        coupled_variables = u
        constant_names = 'W1    d'
        constant_expressions = '1/4 1e-1'
        expression = 'W1*d*(u^2 - 1)^2'
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
    l_abs_tol = 1e-9
    l_max_its = 100
    nl_max_its = 30
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
        file_base = ad_nm_${n}_${d}
    []
    [csv]
        type = CSV
        file_base = ad_nm_${n}_${d}_e
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []
  
