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
# s - solvent
chi12_d = 0.42   # Flory-Huggins parameter
chi12_o = 2.42   # Flory-Huggins parameter
chi1s_d = 0.08   # Flory-Huggins parameter
chi1s_o = 1.61   # Flory-Huggins parameter
chi2s_d = 0.40   # Flory-Huggins parameter
chi2s_o = 2.42   # Flory-Huggins parameter
N1 = 10        # Degree of polymerisation
N2 = 100        # Degree of polymerisation
Ns = 1        # Degree of polymerisation
R = 8.314e14
T = 298
V = 4e22
beta = ${fparse 1e-3*R*T/V}
L1 = 1.5e-5
L2 = 1.5e-5
Ls = 1.5e-5
Tm1 = 400
Tm2 = 320
Tms = 180
lambda = 1.0
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
    [eta]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    [pvfIC1]
        type = RandomIC
        variable = c1
        seed = 123
        min = '${fparse a-0.2}'
        max = '${fparse a+0.2}'
    []
    [pvfIC2]
        type = RandomIC
        variable = c2
        seed = 123
        min = '${fparse (1-a)-0.2}'
        max = '${fparse (1-a)+0.2}'
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

[AuxVariables]  
    [f_density]
        order = CONSTANT
        family = MONOMIAL
    []
    # [c2]
    #     order = FIRST
    #     family = LAGRANGE
    # []
    [f_int_density]
        order = CONSTANT
        family = MONOMIAL
    []
    [cs_time]
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
        coupled_variables = 'c2 eta'
        f_name = f_tot
        kappa_name = kappa_c
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
        coupled_variables = 'c1 eta'
        f_name = f_tot
        kappa_name = kappa_c
        w = w2
    []
    [eta_dot]
        type = TimeDerivative
        variable = eta
    []
    [AC_res]
        type = AllenCahn
        variable = eta
        coupled_variables = 'c1 c2'
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
        kappa_names = 'kappa_c kappa_c kappa_eta'
        interfacial_vars = 'c1 c2 eta'
    []
    # # calculate c2 from c
    # [c2]
    #     type = ParsedAux
    #     variable = c2
    #     coupled_variables = 'c'
    #     expression = '1 - c'
    # []
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
        prop_names = 'kappa_c   kappa_eta'
        prop_values = '${fparse 2.48e1}    ${fparse 2.48e-5}'
    []
    [d2fdc1]
        type = DerivativeParsedMaterial
        property_name = d2fdc1
        # coupled_variables = 'c1 c2'
        constant_names = 'chi1s_d   chi1s_o'
        constant_expressions = '${chi1s_d}  ${chi1s_o}'
        material_property_names = 'h    cs'
        expression = '6.19e-6*(h*(-2*chi1s_d) + (1-h)*(-2*chi1s_o))'
    []
    [d2fdc2]
        type = DerivativeParsedMaterial
        property_name = d2fdc2
        # coupled_variables = 'c1 c2'
        constant_names = 'chi2s_d   chi2s_o'
        constant_expressions = '${chi2s_d}  ${chi2s_o}'
        material_property_names = 'h    cs'
        expression = '6.19e-6*(h*(-2*chi2s_d) + (1-h)*(-2*chi2s_o))'
    []
    [diffusivity_c1]
        type = DerivativeParsedMaterial
        property_name = D1
        # coupled_variables = 'c1 c2'
        material_property_names = 'cs'
        expression = 'if(cs<0.2, 10^(5*cs-1), 1)'
    []
    [diffusivity_c2]
        type = DerivativeParsedMaterial
        property_name = D2
        # coupled_variables = 'c1 c2'
        material_property_names = 'cs'
        expression = 'if(cs<0.2, 10^(15*cs-3), 1)'
    []
    [cs]
        type = DerivativeParsedMaterial
        property_name = cs  
        coupled_variables = 'cs_time'
        expression = 'cs_time'
        outputs = 'ex'
    []  
    [mobility_c1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c2'
        material_property_names = 'cs   h   d2fdc1  D1'
        expression = 'h*D1/d2fdc1'
    []
    [mobility_c2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2'
        material_property_names = 'cs   h   d2fdc2  D2'
        expression = 'h*D2/d2fdc2'
    []
    [mobility_eta]
        type = DerivativeParsedMaterial
        property_name = L
        coupled_variables = 'c1 c2'
        material_property_names = 'cs   kappa_eta'
        constant_names = 'lambda    a'
        constant_expressions = '${lambda}   ${a}'
        expression = '(a*(1-cs)*10^(5*cs-1) + (1-a)*(1-cs)*10^(15*cs-3) + cs)/1'
    []
    [f_constants]
        type = GenericFunctionMaterial
        prop_names = 'f_o_pure_1    f_o_pure_2  f_o_pure_s  f_d_pure_1  f_d_pure_2  f_d_pure_s'
        prop_values = '${fparse L1*(T-Tm1)/Tm1} ${fparse L2*(T-Tm2)/Tm2}    ${fparse Ls*(T-Tms)/Tms}    0   0   0'
    []
    [f_ordered]
        type = DerivativeParsedMaterial
        property_name = f_o_pure
        coupled_variables = 'c1 c2'
        material_property_names = 'f_o_pure_1   f_o_pure_2  f_o_pure_s  cs'
        expression = 'c1*(1-cs)*f_o_pure_1 + c2*(1-cs)*f_o_pure_2 + cs*f_o_pure_s'
        derivative_order = 2
    []
    [f_o_mixing]
        type = DerivativeParsedMaterial
        property_name = f_o_mix
        coupled_variables = 'c1 c2'
        material_property_names = 'cs'
        constant_names = 'chi12 chi1s   chi2s   N1  N2  Ns  R   T   V'
        constant_expressions = '${chi12_o}  ${chi1s_o}  ${chi2s_o}  ${N1}   ${N2}   ${Ns}   ${R}    ${T}    ${V}'
        expression = '(R*T/V)*(c1*(1-cs)*log(c1*(1-cs))/N1 + c2*(1-cs)*log(c2*(1-cs))/N2 + cs*log(cs)/Ns + chi12*c1*c2*(1-cs)^2 + chi1s*c1*cs*(1-cs) + chi2s*c2*cs*(1-cs))'
        derivative_order = 2
    []
    [f_disordered]
        type = DerivativeParsedMaterial
        property_name = f_d_pure
        coupled_variables = 'c1 c2'
        material_property_names = 'f_d_pure_1   f_d_pure_2  f_d_pure_s  cs'
        expression = 'c1*(1-cs)*f_d_pure_1 + c2*(1-cs)*f_d_pure_2 + cs*f_d_pure_s'
        derivative_order = 2
    []
    [f_d_mixing]
        type = DerivativeParsedMaterial
        property_name = f_d_mix
        coupled_variables = 'c1 c2'
        material_property_names = 'cs'
        constant_names = 'chi12 chi1s   chi2s   N1  N2  Ns  R   T   V'
        constant_expressions = '${chi12_d}  ${chi1s_d}  ${chi2s_d}  ${N1}   ${N2}   ${Ns}   ${R}    ${T}    ${V}'
        expression = '(R*T/V)*(c1*(1-cs)*log(c1*(1-cs))/N1 + c2*(1-cs)*log(c2*(1-cs))/N2 + cs*log(cs)/Ns + chi12*c1*c2*(1-cs)^2 + chi1s*c1*cs*(1-cs) + chi2s*c2*cs*(1-cs))'
        derivative_order = 2
    []
    [f_o]
        type = DerivativeSumMaterial
        property_name = f_o
        coupled_variables = 'c1 c2'
        sum_materials = 'f_o_pure f_o_mix'
        derivative_order = 2
    []
    [f_d]
        type = DerivativeSumMaterial
        property_name = f_d
        coupled_variables = 'c1 c2'
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
        coupled_variables = 'c1 c2'
        derivative_order = 2
        W = 3.2e-6
    []
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = f_beta
        coupled_variables = 'c1 c2'
        material_property_names = 'cs'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'beta*(1/(c1*cs) + 1/(c2*cs) + 1/cs)'
        derivative_order = 2
    []
    [total_free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c2  eta'
        sum_materials = 'f  f_beta'
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


    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'

    line_search = 'basic'

    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
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

    # # Automatic scaling for u and w
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
        file_base = output/3phase_wp_p1
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3phase_wp_p1
    []
[]

[Debug]
    show_var_residual_norms = true
    show_execution_order = 'XFEM_MARK FORWARD ADJOINT HOMOGENEOUS_FORWARD ADJOINT_TIMESTEP_BEGIN ADJOINT_TIMESTEP_END NONE INITIAL LINEAR NONLINEAR_CONVERGENCE NONLINEAR POSTCHECK TIMESTEP_END TIMESTEP_BEGIN MULTIAPP_FIXED_POINT_END MULTIAPP_FIXED_POINT_BEGIN FINAL FAILED CUSTOM ALWAYS TRANSFER'
[]