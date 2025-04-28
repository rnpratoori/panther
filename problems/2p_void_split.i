n = 100     # number of elements per side
d = 1       # ND size of the side
a = 0.5     # type A monomer density
chi12 = 2.0   # Flory-Huggins parameter
chi13 = 10.0   # penalty term for void
chi23 = 10.0   # penalty term for void
N1 = 5       # Degree of polymerisation
N2 = 5       # Degree of polymerisation
N3 = 100     # Penalty term for void
M1 = 1       # Mobility for c1
M3 = 1       # Mobility for c3
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn3 = 1e-2  # Cahn number for c3
k3 = ${fparse Cn3^2}    # gradient energy coefficient

R = 1  # Universal gas constant
T = 1 # Temperature in Kelvin
beta = 1e-3*R*T
delta = 1e-4

[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}
    ymax = ${d}
    # uniform_refine = 2
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
        min = '${fparse a-0.04}'
        max = '${fparse a+0.04}'
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
        pos_variation = 0.2
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
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa3'
        interfacial_vars = 'c1 c3'
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
        prop_names = 'kappa      kappa3'
        prop_values = '${fparse k*s}   ${fparse k3*s}'
    []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1'
        constant_names = 'M1     s'
        constant_expressions = '${M1} ${s}'
        expression = '(M1*(1-c1)^2)/s'
        # derivative_order = 2
    []
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c3'
        constant_names = 'M3     s'
        constant_expressions = '${M3} ${s}'
        expression = 'if (c3>0.9, 0, (M3*(1-c3)^2)/s)'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c3'
        constant_names = 'R      T       chi12      chi13       chi23     N1        N2      N3       s      beta'
        constant_expressions = '${R}    ${T}    ${chi12}    ${chi13}    ${chi23}  ${N1} ${N2}    ${N3}    ${s}    ${beta}'
        expression = 'if(c3>1e-4, s*(R*T*(c1*log(c1)/N1 + c3*log(c3)/N3 + (1-c1-c3)*log(1-c1-c3)/N2 + chi12*c1*(1-c1-c3) + chi13*c1*c3 + chi23*c3*(1-c1-c3)) + beta*(1/c1 + 1/c3 + 1/(1-c1-c3))), s*(R*T*(c1*log(c1)/N1 + (1-c1-c3)*log(1-c1-c3)/N2 + chi12*c1*(1-c1-c3))))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c3'
        sum_materials = 'f_mix'
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

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'lu'
    
    # line_search = 'basic'

    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 50
    nl_max_its = 30
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-10
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e-2 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'c1 w1; c2 w2'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/2p_void_split_t1
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/2p_void_split_t1
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []