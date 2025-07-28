n = 100     # number of elements per side
d = 1       # ND size of the side
a = 0.5     # type A monomer density
M1 = 1       # Initial mobility, depends on swell ratio
M3 = 0e-3    # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn3 = 1e-2  # Cahn number
k3 = ${fparse Cn3^2}    # gradient energy coefficient

# chi12 = 1.0   # Flory-Huggins parameter
# chi13 = 10.0   # Flory-Huggins parameter
# chi23 = 10.0   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 100     # Penalty term for void
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
beta = 1e-3
delta = 0

[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}
    ymax = ${d}
    add_subdomain_ids = '0 1'
[]

[MeshModifiers]
    [void]
      type = CoupledVarThresholdElementSubdomainModifier
      coupled_var = c3
      criterion_type = ABOVE
      subdomain_id = 1
      threshold = 1e-6
    []
[]

[Variables]
    # polymer volume fraction
    [c1]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w1]
        order = FIRST
        family = LAGRANGE
    []
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w3]
        order = FIRST
        family = LAGRANGE   
    []
[]

[ICs]
    [pvfIC_1]
        type = RandomConstraintIC
        variable = c1
        min = '${fparse a-0.05}'
        max = '${fparse a+0.05}'
        seed = 123
        # distribution = Normal_a
        coupled = c3
    []
    [voidIC]
        type = LatticeSmoothCircleIC
        variable = c3
        invalue = ${fparse 1.0-delta}
        outvalue = ${delta}
        circles_per_side = '2 2'
        pos_variation = 0.1
        radius = 0.1
        int_width = 0.01
        radius_variation_type = uniform
        avoid_bounds = true
    []
[]

[Distributions]
    [Normal_a]
        type = Normal
        mean = ${a}
        standard_deviation = 0.02
    []
[]

[AuxVariables]
    [c2]
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
        coupled_variables = 'c3'
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    [w3_dot]
        type = CoupledTimeDerivative
        variable = w3
        v = c3
    []
    # adding nonlocal term to the energy
    [coupled_res3]
        type = SplitCHWRes
        variable = w3
        mob_name = M3
    []
    [coupled_parsed3]
        type = SplitCHParsed
        variable = c3
        coupled_variables = 'c1'
        f_name = f_mix
        kappa_name = kappa3
        w = w3
    []
[]

[AuxKernels]
    [c2]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c1 c3'
        expression = '1 - c1 - c3'
    []
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa    kappa3'
        interfacial_vars = 'c1  c3'
    []
    # calculate interfacial energy density
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
        prop_names = 'kappa     kappa3'
        prop_values = '${fparse k*s}    ${fparse k3*s}'
    []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c3'
        constant_names = 'M1     s'
        constant_expressions = '${M1} ${s}'
        expression = '(M1*16*(c1^2*(1-c1)^2)*(1-c3))/s'
        # derivative_order = 2
    []
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c3'
        constant_names = 'M3     s'
        constant_expressions = '${M3} ${s}'
        expression = '(M3*16*c3^2*(1-c3)^2)/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1'
        # constant_names =        'p      q       r       s
        #                         t      u        v       S
        #                         c0      beta    z'
        # constant_expressions = '${p}    ${q}    ${r}    ${s}
        #                         ${t}    ${u}    ${v}    ${S}
        #                         ${c0}   ${beta} ${z}'
        expression = 'if(c1<1.0000000000000000e-08, -2.5057439080344793e+13*(c1-0.0000000000000000e+00)^3 + 7.8674294099482000e+06*(c1-0.0000000000000000e+00)^2 + -2.9561685501914474e+00*(c1-0.0000000000000000e+00) + 0.0000000000000000e+00,
            if(c1<9.9999999999999995e-08, -2.5057439080335312e+13*(c1-1.0000000000000000e-08)^3 + 7.1157062375379289e+06*(c1-1.0000000000000000e-08)^2 + -2.8063371937165869e+00*(c1-1.0000000000000000e-08) + -2.8800000000000000e-08,
            if(c1<9.9999999999999995e-07, -1.1069255004891592e+11*(c1-9.9999999999999995e-08)^3 + 3.5019768584737869e+05*(c1-9.9999999999999995e-08)^2 + -2.1344058406119077e+00*(c1-9.9999999999999995e-08) + -2.4200000000000002e-07,
            if(c1<1.0000000000000001e-05, -1.7368391407276378e+09*(c1-9.9999999999999995e-07)^3 + 5.1327800715307094e+04*(c1-9.9999999999999995e-07)^2 + -1.7730329027054919e+00*(c1-9.9999999999999995e-07) + -1.9599999999999999e-06,
            if(c1<1.0000000000000000e-04, -1.4666624999236980e+07*(c1-1.0000000000000001e-05)^3 + 4.4331439156609822e+03*(c1-1.0000000000000001e-05)^2 + -1.2711844010267801e+00*(c1-1.0000000000000001e-05) + -1.5025900000000000e-05,
            if(c1<1.0000000000000000e-03, -1.5752393557192417e+05*(c1-1.0000000000000000e-04)^3 + 4.7315516586699800e+02*(c1-1.0000000000000000e-04)^2 + -8.2961748368926191e-01*(c1-1.0000000000000000e-04) + -1.0421600000000001e-04,
            if(c1<1.7977659999999999e-03, 4.4670927656117537e+03*(c1-1.0000000000000000e-03)^3 + 4.7840539822803173e+01*(c1-1.0000000000000000e-03)^2 + -3.6072134856844118e-01*(c1-1.0000000000000000e-03) + -5.8245100000000000e-04,
            if(c1<3.5945321000000001e-03, -6.5593815167842477e+03*(c1-1.7977659999999999e-03)^3 + 5.8531624004556605e+01*(c1-1.7977659999999999e-03)^2 + -2.7586125292054398e-01*(c1-1.7977659999999999e-03) + -8.3750700000000003e-04,
            if(c1<5.3912980999999997e-03, -9.5089857212632137e+02*(c1-3.5945321000000001e-03)^3 + 2.3174600965582986e+01*(c1-3.5945321000000001e-03)^2 + -1.2905427773522352e-01*(c1-3.5945321000000001e-03) + -1.1822520000000000e-03,
            if(c1<7.1880642000000002e-03, -9.9597377835417899e+02*(c1-5.3912980999999997e-03)^3 + 1.8048974294047621e+01*(c1-5.3912980999999997e-03)^2 + -5.4985159310278087e-02*(c1-5.3912980999999997e-03) + -1.3448320000000000e-03,
            if(c1<8.9848302000000001e-03, -4.6066282646900726e+02*(c1-7.1880642000000002e-03)^3 + 1.2680378529740516e+01*(c1-7.1880642000000002e-03)^2 + 2.2830011844371736e-04*(c1-7.1880642000000002e-03) + -1.3911360000000001e-03,
            if(c1<1.0781595999999999e-02, -3.3947886385268129e+02*(c1-8.9848302000000001e-03)^3 + 1.0197268617550295e+01*(c1-8.9848302000000001e-03)^2 + 4.1334078672692807e-02*(c1-8.9848302000000001e-03) + -1.3524609999999999e-03,
            if(c1<1.2578361999999999e-02, -2.8120077574184234e+02*(c1-1.0781595999999999e-02)^3 + 8.3673765803702587e+00*(c1-1.0781595999999999e-02)^2 + 7.4690398253450627e-02*(c1-1.0781595999999999e-02) + -1.2472420000000000e-03,
            if(c1<1.4375127999999999e-02, -2.8659960397195075e+01*(c1-1.2578361999999999e-02)^3 + 6.8516206012905370e+00*(c1-1.2578361999999999e-02)^2 + 1.0203537494355461e-01*(c1-1.2578361999999999e-02) + -1.0876589999999999e-03,
            if(c1<1.6171893999999999e-02, -7.1972940960162930e+02*(c1-1.4375127999999999e-02)^3 + 6.6971348740813728e+00*(c1-1.4375127999999999e-02)^2 + 1.2637931812401684e-01*(c1-1.4375127999999999e-02) + -8.8237199999999997e-04,
            if(c1<1.1293752000000000e-01, -1.1001590632023886e+01*(c1-1.6171893999999999e-02)^3 + 2.8175788769645367e+00*(c1-1.6171893999999999e-02)^2 + 1.4347503229162858e-01*(c1-1.6171893999999999e-02) + -6.3785199999999995e-04,
            if(c1<2.0970314000000001e-01, 1.4703666032807175e-01*(c1-1.1293752000000000e-01)^3 + -3.7614853654604224e-01*(c1-1.1293752000000000e-01)^2 + 3.7972156751761715e-01*(c1-1.1293752000000000e-01) + 2.9659995000000001e-02,
            if(c1<3.0646876000000001e-01, -7.6783767459934260e-01*(c1-2.0970314000000001e-01)^3 + -3.3346425574791871e-01*(c1-2.0970314000000001e-01)^2 + 3.1105544571136101e-01*(c1-2.0970314000000001e-01) + 6.3015114999999997e-02,
            if(c1<4.0323438000000000e-01, -8.6313817159341000e-02*(c1-3.0646876000000001e-01)^3 + -5.5636512167380936e-01*(c1-3.0646876000000001e-01)^2 + 2.2495055431093347e-01*(c1-3.0646876000000001e-01) + 8.9296449999999999e-02,
            if(c1<5.0000000000000000e-01, -8.2896614040814401e-02*(c1-4.0323438000000000e-01)^3 + -5.8142175176978028e-01*(c1-4.0323438000000000e-01)^2 + 1.1485190207430300e-01*(c1-4.0323438000000000e-01) + 1.0577615100000000e-01,
            if(c1<5.9676562000000000e-01, 8.2896408401759783e-02*(c1-5.0000000000000000e-01)^3 + -6.0548637853046061e-01*(c1-5.0000000000000000e-01)^2 + 9.6275940691947781e-10*(c1-5.0000000000000000e-01) + 1.1137056400000001e-01,
            if(c1<6.9353123999999999e-01, 8.6314845354622069e-02*(c1-5.9676562000000000e-01)^3 + -5.8142181146615268e-01*(c1-5.9676562000000000e-01)^2 + -1.1485190592534063e-01*(c1-5.9676562000000000e-01) + 1.0577615100000000e-01,
            if(c1<7.9029685999999999e-01, 7.6783376745728615e-01*(c1-6.9353123999999999e-01)^3 + -5.5636488288832131e-01*(c1-6.9353123999999999e-01)^2 + -2.2495054083230170e-01*(c1-6.9353123999999999e-01) + 8.9296449999999999e-02,
            if(c1<8.8706247999999999e-01, -1.4702205995513731e-01*(c1-7.9029685999999999e-01)^3 + -3.3346515119349757e-01*(c1-7.9029685999999999e-01)^2 + -3.1105549577485037e-01*(c1-7.9029685999999999e-01) + 6.3015114999999997e-02,
            if(c1<9.8382811000000003e-01, 1.1001536771074347e+01*(c1-8.8706247999999999e-01)^3 + -3.7614519354920373e-01*(c1-8.8706247999999999e-01)^2 + -3.7972138074229178e-01*(c1-8.8706247999999999e-01) + 2.9659995000000001e-02,
            if(c1<9.8562486999999999e-01, 7.1980887143528207e+02*(c1-9.8382811000000003e-01)^3 + 2.8175667163143059e+00*(c1-9.8382811000000003e-01)^2 + -1.4347568899636559e-01*(c1-9.8382811000000003e-01) + -6.3785199999999995e-04,
            if(c1<9.8742163999999999e-01, 2.8509378376754860e+01*(c1-9.8562486999999999e-01)^3 + 6.6975380798343593e+00*(c1-9.8562486999999999e-01)^2 + -1.2637932930283793e-01*(c1-9.8562486999999999e-01) + -8.8237199999999997e-04,
            if(c1<9.8921840000000005e-01, 2.8134008569972605e+02*(c1-9.8742163999999999e-01)^3 + 6.8512124671923988e+00*(c1-9.8742163999999999e-01)^2 + -1.0203534078245668e-01*(c1-9.8742163999999999e-01) + -1.0876589999999999e-03,
            if(c1<9.9101516999999995e-01, 3.3937772036966066e+02*(c1-9.8921840000000005e-01)^3 + 8.3677143043379090e+00*(c1-9.8921840000000005e-01)^2 + -7.4690581916440804e-02*(c1-9.8921840000000005e-01) + -1.2472420000000000e-03,
            if(c1<9.9281193999999995e-01, 4.6073393659106887e+02*(c1-9.9101516999999995e-01)^3 + 1.0197065424223586e+01*(c1-9.9101516999999995e-01)^2 + -4.1333942643555374e-02*(c1-9.9101516999999995e-01) + -1.3524609999999999e-03,
            if(c1<9.9460870000000001e-01, 9.9586756994319614e+02*(c1-9.9281193999999995e-01)^3 + 1.2680564169969795e+01*(c1-9.9281193999999995e-01)^2 + -2.2810411759646189e-04*(c1-9.9281193999999995e-01) + -1.3911360000000001e-03,
            if(c1<9.9640547000000002e-01, 9.5115450783761219e+02*(c1-9.9460870000000001e-01)^3 + 1.8048569214883361e+01*(c1-9.9460870000000001e-01)^2 + 5.4984773582974322e-02*(c1-9.9460870000000001e-01) + -1.3448320000000000e-03,
            if(c1<9.9820222999999997e-01, 6.5588680457606251e+03*(c1-9.9640547000000002e-01)^3 + 2.3175586870025548e+01*(c1-9.9640547000000002e-01)^2 + 1.2905510051165620e-01*(c1-9.9640547000000002e-01) + -1.1822520000000000e-03,
            if(c1<9.9900000000000000e-01, -4.4654743087282859e+03*(c1-9.9820222999999997e-01)^3 + 5.8529722119727467e+01*(c1-9.9820222999999997e-01)^2 + 2.7585993149208049e-01*(c1-9.9820222999999997e-01) + -8.3750700000000003e-04,
            if(c1<9.9990000000000001e-01, 1.5752289473127548e+05*(c1-9.9900000000000000e-01)^3 + 4.7842457801904345e+01*(c1-9.9900000000000000e-01)^2 + 3.6072046546816416e-01*(c1-9.9900000000000000e-01) + -5.8245100000000000e-04,
            if(c1<9.9999000000000005e-01, 1.4666629972857999e+07*(c1-9.9990000000000001e-01)^3 + 4.7315427357635099e+02*(c1-9.9990000000000001e-01)^2 + 8.2961752370860264e-01*(c1-9.9990000000000001e-01) + -1.0421600000000001e-04,
            if(c1<9.9999899999999997e-01, 1.7368390875549436e+09*(c1-9.9999000000000005e-01)^3 + 4.4331443662495221e+03*(c1-9.9999000000000005e-01)^2 + 1.2711844012931011e+00*(c1-9.9999000000000005e-01) + -1.5025900000000000e-05,
            if(c1<9.9999990000000005e-01, 1.1069255652491737e+11*(c1-9.9999899999999997e-01)^3 + 5.1327799729845981e+04*(c1-9.9999899999999997e-01)^2 + 1.7730328981538199e+00*(c1-9.9999899999999997e-01) + -1.9599999999999999e-06,
            if(c1<9.9999998999999995e-01, 2.5057438123848500e+13*(c1-9.9999990000000005e-01)^3 + 3.5019770237413980e+05*(c1-9.9999990000000005e-01)^2 + 2.1344058500800980e+00*(c1-9.9999990000000005e-01) + -2.4200000000000002e-07,
            if(c1<1.0000000000000000e+00, 2.5057438123840047e+13*(c1-9.9999998999999995e-01)^3 + 7.1157059880792024e+06*(c1-9.9999998999999995e-01)^2 + 2.8063371814527889e+00*(c1-9.9999998999999995e-01) + -2.8800000000000000e-08,
            0.0))))))))))))))))))))))))))))))))))))))))'
        derivative_order = 2
    []
    # beta penalty term
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = 'f_beta'
        coupled_variables = 'c1'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'if (c1>0, beta*(1/c1 + 1/(1 - c1)), 0)'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1'
        sum_materials = 'f_mix f_beta'
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
    [step_size]
        type = TimestepSize
    []

[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    # petsc_options = '-pc_svd_monitor -ksp_view'
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    # petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    # petsc_options = '-pc_svd_monitor -ksp_view'
    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'svd'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
    nl_max_its = 100
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e0 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    # off_diagonal_
    scaling_group_variables = 'c1 c3; w1 w3'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/2p_void_spline
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/2p_void_spline
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []