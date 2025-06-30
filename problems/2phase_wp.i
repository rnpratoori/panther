nx = 400     # number of elements per side
ny = 100      # number of elements per side
dx = 400       # ND size of the side
dy = 100       # ND size of the side
a = 0.3     # type A monomer density
# M = 1e0     # Initial mobility, depends on swell ratio
# Cn = 5e-2   # Cahn number
# k = ${fparse Cn^2}    # gradient energy coefficient

# 1 - drug
# 2 - polymer
chi12_d = 0.42   # Flory-Huggins parameter
chi12_o = 2.42   # Flory-Huggins parameter
N1 = 10        # Degree of polymerisation
N2 = 100        # Degree of polymerisation
R = 8.314
T = 298
V = 40e-6
beta = ${fparse 1e-3*R*T/V}
L1 = 150e6
L2 = 150e6
Tm1 = 400
Tm2 = 320
lambda = 1.0e-9
cs0 = 0.2


[Mesh]
    [2p]
        # generate a 2D mesh
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
        # uniform_refine = 2
    []
[]

[Variables]
    # polymer volume fraction
    [c]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w]
        order = FIRST
        family = LAGRANGE
    []
    [eta]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    [pvfIC]
        type = RandomIC
        variable = c
        seed = 123
        min = '${fparse a-0.2}'
        max = '${fparse a+0.2}'
    []
    [etaIC]
        type = MultiBoundingBoxIC
        variable = eta
        corners = '327 6 0  309 42 0    354 18 0    16 50 0
                65 93 0     239 0 0     283 53 0    25 66 0
                5 31 0      325 38 0    387 5 0     318 28 0
                265 39 0    203 28 0    94 60 0     13 62 0'
        opposite_corners = '332 11 0    314 47 0    359 23 0    21 55 0
                70 98 0     244 5 0     288 58 0    30 71 0
                10 36 0     330 43 0    392 10 0    323 33 0
                270 44 0    208 33 0    99 65 0     18 67 0'
        inside = 0
        outside = 1
    []
[]

# [Functions]
#   [cs_time]
#     type = ParsedFunction
#     expression = '${cs0}*(1-t/200)'
#   []
# []

[AuxVariables]  
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    [f_int_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [cs_time]
    []
[]

[Kernels]
    [w_dot]
        type = CoupledTimeDerivative
        variable = w
        v = c
    []
    # adding nonlocal term to the energy
    [coupled_res]
        type = SplitCHWRes
        variable = w
        mob_name = M
    []
    [coupled_parsed]
        type = SplitCHParsed
        variable = c
        f_name = f_tot
        coupled_variables = 'eta'
        kappa_name = kappa_c
        w = w
    []
    [eta_dot]
        type = TimeDerivative
        variable = eta
    []
    [AC_res]
        type = AllenCahn
        variable = eta
        coupled_variables = 'c'
        f_name = f_tot
        mob_name = L
    []
    [AC_interface]
        type = ACInterface
        variable = eta
        mob_name = L
        kappa_name = kappa_eta
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa_c kappa_eta'
        interfacial_vars = 'c eta'
    []
    # calculate c2 from c
    [c2]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c'
        expression = '1 - c'
    []
    # calculate interfacial energy density
    [f_int_density]
        type = ParsedAux
        variable = f_int_density
        coupled_variables = 'f_density'
        material_properties = 'f_tot'
        expression = 'f_density - f_tot'
    []
    [cs_time]
        type = ParsedAux
        variable = cs_time
        expression = '${cs0}*(1-t/200)'
        use_xyzt = true
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa_c kappa_eta'
        prop_values = '${fparse 7.5e-3/lambda^2} ${fparse 2e-10/lambda^2}'
    []
    [cs]
        type = DerivativeParsedMaterial
        property_name = cs  
        coupled_variables = 'cs_time'
        expression = 'cs_time'
        outputs = 'ex'
    []  
    [mobility_c]
        type = DerivativeParsedMaterial
        property_name = M
        coupled_variables = 'c'
        material_property_names = 'cs h'
        expression = 'h*10^(5*cs-1)/1.58e8'
    []
    [mobility_eta]
        type = DerivativeParsedMaterial
        property_name = L
        coupled_variables = 'c'
        material_property_names = 'cs   kappa_eta'
        constant_names = 'lambda'
        constant_expressions = '${lambda}'
        expression = '(0.3*10^(5*cs-1) + (1-0.3)*10^(15*cs-3) + cs)*lambda^2/kappa_eta'
    []
    [f_constants]
        type = GenericFunctionMaterial
        prop_names = 'f_o_pure_1 f_o_pure_2 f_d_pure_1 f_d_pure_2'
        prop_values = '${fparse L1*(T-Tm1)/Tm1}  ${fparse L2*(T-Tm2)/Tm2}   0   0'
    []
    [f_ordered]
        type = DerivativeParsedMaterial
        property_name = f_o_pure
        coupled_variables = 'c'
        material_property_names = 'f_o_pure_1 f_o_pure_2'
        expression = 'c*f_o_pure_1 + (1-c)*f_o_pure_2'
        derivative_order = 2
    []
    [f_o_mixing]
        type = DerivativeParsedMaterial
        property_name = f_o_mix
        coupled_variables = 'c'
        constant_names = 'chi12 N1 N2 R T V'
        constant_expressions = '${chi12_o} ${N1} ${N2} ${R} ${T} ${V}'
        expression = '(R*T/V)*(c*log(c)/N1 + (1-c)*log(1-c)/N2 + chi12*c*(1-c))'
        derivative_order = 2
    []
    [f_disordered]
        type = DerivativeParsedMaterial
        property_name = f_d_pure
        coupled_variables = 'c'
        material_property_names = 'f_d_pure_1 f_d_pure_2'
        expression = 'c*f_d_pure_1 + (1-c)*f_d_pure_2'
        derivative_order = 2
    []
    [f_d_mixing]
        type = DerivativeParsedMaterial
        property_name = f_d_mix
        coupled_variables = 'c'
        constant_names = 'chi12 N1 N2 R T V'
        constant_expressions = '${chi12_d} ${N1} ${N2} ${R} ${T} ${V}'
        expression = '(R*T/V)*(c*log(c)/N1 + (1-c)*log(1-c)/N2 + chi12*c*(1-c))'
        derivative_order = 2
    []
    [f_o]
        type = DerivativeSumMaterial
        property_name = f_o
        coupled_variables = 'c'
        sum_materials = 'f_o_pure f_o_mix'
        derivative_order = 2
    []
    [f_d]
        type = DerivativeSumMaterial
        property_name = f_d
        coupled_variables = 'c'
        sum_materials = 'f_d_pure f_d_mix'
        derivative_order = 2
    []
    [switching_function]
        type = SwitchingFunctionMaterial
        eta = eta
        h_order = HIGH
    []
    [barrier_function]
        type = BarrierFunctionMaterial
        eta = eta
        g_order = SIMPLE
    []
    [free_energy]
        type = DerivativeTwoPhaseMaterial
        property_name = f
        eta = eta
        fa_name = f_o
        fb_name = f_d
        coupled_variables = 'c'
        derivative_order = 2
        W = 30
    []
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = f_beta
        coupled_variables = 'c'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'beta*(1/c + 1/(1-c))'
        derivative_order = 2
    []
    [total_free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c eta'
        sum_materials = 'f f_beta'
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
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    # petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31                  preonly      ilu          1'


    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'

    line_search = 'basic'

    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'

    l_tol = 1e-10
    l_abs_tol = 1e-12
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-12

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e2 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    scaling_group_variables = 'c w'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/2phase_wp
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/2phase_wp
    []
[]