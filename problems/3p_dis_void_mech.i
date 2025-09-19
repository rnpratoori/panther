timestep = 010

nx = 202     # number of elements in x
ny = 100     # number of elements in y
# dx = 1.02       # ND size of the side in x
dy = 1.00       # ND size of the side in y

[GlobalParams]
    displacements = 'disp_x disp_y'
    large_kinematics = true
[]

[Mesh]
    # add_subdomain_ids = '1'
    # uniform_refine = 2
    [2d]
        # generate a 2D mesh
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmin = -0.01
        xmax = 2.01
        ymax = ${dy}
    []
    # Subdomain for polymer block
    [active_domain]
        type = SubdomainBoundingBoxGenerator
        input = 2d
        block_id = 1
        bottom_left = '0 0 0'
        top_right = '2.00 1.00 0'
    []
[]

[Variables]
    [disp_x]
    []
    [disp_y]
    []
[]

[AuxVariables]
    [c1]
        order = FIRST
        family = LAGRANGE
        [InitialCondition]
            type = SolutionIC
            from_variable = 'c1'
            solution_uo = 3phase
            variable = c1
            block = 1
            from_subdomains = 0
        []
    []
    [c2]
        order = FIRST
        family = LAGRANGE
        [InitialCondition]
            type = SolutionIC
            from_variable = 'c2'
            solution_uo = 3phase
            variable = c2
            block = 1
            from_subdomains = 0
        []
    []
    [stress_xx]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [stress_xx]
      type = RankTwoAux
      variable = stress_xx
      rank_two_tensor = cauchy_stress
      index_i = 0
      index_j = 0
      execute_on = 'timestep_end'
    []
[]

[ICs]
    [0c1]
        type = ConstantIC
        value = 1
        variable = c1
        block = 0
    []
    [0c2]
        type = ConstantIC
        value = 0
        variable = c2
        block = 0
    []
[]

[UserObjects]
    [3phase]
      type = SolutionUserObject
      mesh = 'output/mechic/3pv_0.3_0.05_0.2_mech_ic.e-s${timestep}'
      system_variables = 'c1 c2'
      # timestep = ${timestep}
    []
[]

[Kernels]
    [sdx]
      type = TotalLagrangianStressDivergence
      variable = disp_x
      component = 0
    []
    [sdy]
      type = TotalLagrangianStressDivergence
      variable = disp_y
      component = 1
    []
[]

[Functions]
    [pullx]
      type = ParsedFunction
      expression = '1e0*t'
    []
  #   [pully]
  #     type = ParsedFunction
  #     expression = '-2000 * t'
  #   []
    # [mu_function]
    #   type = ParsedFunction
    #   expression = '1e6*exp(-10*(1-c1-c2))'
    # []
[]

[BCs]
    [./left_x_fixed]
      type = DirichletBC
      preset = true
      variable = disp_x
      boundary = 'left'
      value = 0
    [../]
    [./bottom_y_fixed]
      type = DirichletBC
      preset = true
      variable = disp_y
      boundary = 'bottom'
      value = 0
    [../]
    [./right_x_strain]
      type = FunctionDirichletBC
      variable = disp_x
      boundary = 'right'
      function = pullx
    [../]
[]
  
  
[Materials]
    [compute_stress]
      type = ComputeNeoHookeanStress
      lambda = lambda
      mu = mu
    []
    [mu]
      type = ParsedMaterial
        property_name = mu
        coupled_variables = 'c1 c2'
        expression = '(1e5*c1 + 1e6*c2)*exp(-10*(1-c1-c2-0.5))'
        output_properties = 'mu'
        outputs = 'ex'
    []
    [lambda]
      type = ParsedMaterial
      property_name = lambda
      coupled_variables = 'c1 c2'
      expression = '5e10'
      output_properties = 'lambda'
      outputs = 'ex'
    []
    [compute_strain]
      type = ComputeLagrangianStrain
    []
[]

[Postprocessors]
    [reaction_force_x]
        type        = SideIntegralVariablePostprocessor
        variable    = stress_xx
        boundary    = right
        execute_on  = 'timestep_end'
    []
  
    [avg_disp_right]
      type = SideAverageValue
      variable = disp_x
      boundary = right
      execute_on = 'timestep_end'
    []
  []
  
  
[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]  
  
[Executioner]
    type = Transient
  
    solve_type = 'newton'
    # line_search = none
  
    # petsc_options = ''
    petsc_options_iname = '-pc_type -ksp_type'
    petsc_options_value = 'lu gmres'

    # petsc_options = '-pc_svd_monitor -ksp_view'
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    # petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'
    
  
    reuse_preconditioner = false
    reuse_preconditioner_max_linear_its = 20
  
    nl_abs_tol = 1e-6
    nl_rel_tol = 1e-6
    nl_max_its = 20
    l_tol = 1e-6
    l_max_its = 50

    [TimeStepper]
        type = ConstantDT
        dt = 1e-2
      []
  
    # [TimeStepper]
    #   # Turn on time stepping
    #   type = IterationAdaptiveDT
    #   dt = 1.0e-4
    #   cutback_factor = 0.8
    #   growth_factor = 1.5
    #   optimal_iterations = 10
    # []
    end_time = 1.0
  
    [Predictor]
          type = SimplePredictor
          scale = 1
    []

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]
  
  
  [Outputs]
    [ex]
        type = Exodus
        file_base = 'output/mech_void/mech_${timestep}'
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = 'output/mech_void/mech_${timestep}'
        execute_on = 'timestep_end'
        show = 'reaction_force_x avg_disp_right'
      []
  []

