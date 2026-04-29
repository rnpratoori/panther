nx = 100     # number of elements in x
ny = 200     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 2.00       # ND size of the side in y
M = 1e-0       # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

chi12 = 1.0   # Flory-Huggins parameter
chi13 = 0.1   # Flory-Huggins parameter
chi23 = 0.1   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 1       # Degree of polymerisation
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
A00 = -0.485203
A10 = -0.384112
A20 = 1.4
A30 = 0.133333
A40 = 1.73333
A50 = -1.76
A60 = 7.89333
A01 = -0.384112
A11 = 2.0
A21 = 2.0
A31 = 2.66667
A41 = 4.0
A51 = 6.4
A61 = 10.6667
A02 = 1.4
A12 = 2.0
A22 = 4.0
A32 = 8.0
A42 = 16.0
A52 = 32.0
A62 = 64.0
A03 = 0.133333
A13 = 2.66667
A23 = 8.0
A33 = 21.3333
A43 = 53.3333
A53 = 128.0
A63 = 298.667
A04 = 1.73333
A14 = 4.0
A24 = 16.0
A34 = 53.3333
A44 = 160.0
A54 = 448.0
A64 = 1194.67
A05 = -1.76
A15 = 6.4
A25 = 32.0
A35 = 128.0
A45 = 448.0
A55 = 1433.6
A65 = 4300.8
A06 = 7.89333
A16 = 10.6667
A26 = 64.0
A36 = 298.667
A46 = 1194.67
A56 = 4300.8
A66 = 14336.0
c1_0 = 0.25
c2_0 = 0.25
beta = 1.0e-3       # Stability parameter
delta = 0.025

[Mesh]
    [2d]
        # generate a 2D mesh
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
        # uniform_refine = 2
    []
    # Subdomain for ramp
    [c3_domain]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        combinatorial_geometry = 'y > 1'
        input = 2d
        # uniform_refine = 2
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
    # polymer volume fraction
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w2]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    [c1]
        type = SolutionIC
        from_variable = 'c'
        solution_uo = 2phase
        variable = c1
        block = 0
    []
    [c2]
        type = SolutionIC
        from_variable = 'c2'
        solution_uo = 2phase
        variable = c2
        block = 0
    []
    [top_c1]
        type = ConstantIC
        value = ${delta}
        variable = c1
        block = 1
    []
    [top_c2]
        type = ConstantIC
        value = ${delta}
        variable = c2
        block = 1
    []
[]

[UserObjects]
  [2phase]
    type = SolutionUserObject
    mesh = 'ic_2p/2phase_0.3.e'
    system_variables = 'c c2'
    timestep = LATEST
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
    [coupled_res1]
        type = SplitCHWRes
        variable = w1
        mob_name = M1
    []
    [coupled_parsed1]
        type = SplitCHParsed
        variable = c1
        coupled_variables = 'c2'
        f_name = f_mix
        kappa_name = kappa
        w = w1
    []
    [w2_dot]
        type = CoupledTimeDerivative
        variable = w2
        v = c2
    []
    [coupled_res2]
        type = SplitCHWRes
        variable = w2
        mob_name = M2
    []
    [coupled_parsed2]
        type = SplitCHParsed
        variable = c2
        coupled_variables = 'c1'
        f_name = f_mix
        kappa_name = kappa
        w = w2
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa'
        interfacial_vars = 'c1 c2'
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

[BCs]
    [top1]
        type = DirichletBC
        variable = c1
        boundary = 2
        # block = 1
        value = ${delta}
    []
    [top2]
        type = DirichletBC
        variable = c2
        boundary = 2
        # block = 1
        value = ${delta}
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa'
        prop_values = '${fparse k*s}'
    []
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c2'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        # expression = 'M/s'
        expression = '(M*exp((15*(1-c1-c2)-3)))/s'
        # expression = '(M*16*c1^2*(1-c1)^2)/s'
        # expression = 'if (c1>0, if(c1<1, (M)/s, 0), 0)'
        # derivative_order = 2
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        # expression = 'M/s'
        expression = '(M*exp((15*(1-c1-c2)-3)))/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = 'f_mix'           
        coupled_variables = 'c1 c2'
        constant_names =      'A00    A10    A20    A30    A40    A50    A60
                                A01    A11    A21    A31    A41    A51    A61
                                A02    A12    A22    A32    A42    A52    A62
                                A03    A13    A23    A33    A43    A53    A63
                                A04    A14    A24    A34    A44    A54    A64
                                A05    A15    A25    A35    A45    A55    A65
                                A06    A16    A26    A36    A46    A56    A66
                                s      c1_0   c2_0  chi12   chi13   chi23'
        constant_expressions = '${A00} ${A10} ${A20} ${A30} ${A40} ${A50} ${A60}
                                ${A01} ${A11} ${A21} ${A31} ${A41} ${A51} ${A61}
                                ${A02} ${A12} ${A22} ${A32} ${A42} ${A52} ${A62}
                                ${A03} ${A13} ${A23} ${A33} ${A43} ${A53} ${A63}
                                ${A04} ${A14} ${A24} ${A34} ${A44} ${A54} ${A64}
                                ${A05} ${A15} ${A25} ${A35} ${A45} ${A55} ${A65}
                                ${A06} ${A16} ${A26} ${A36} ${A46} ${A56} ${A66}
                                ${s}   ${c1_0} ${c2_0} ${chi12} ${chi13} ${chi23}'
        expression = 's*(A00 +
                    A10*(c1-c1_0) + A01*(c2-c2_0) +
                    A20*(c1-c1_0)^2 + A11*(c1-c1_0)*(c2-c2_0) + A02*(c2-c2_0)^2 +
                    A30*(c1-c1_0)^3 + A21*(c1-c1_0)^2*(c2-c2_0) + A12*(c1-c1_0)*(c2-c2_0)^2 + A03*(c2-c2_0)^3 +
                    A40*(c1-c1_0)^4 + A31*(c1-c1_0)^3*(c2-c2_0) + A22*(c1-c1_0)^2*(c2-c2_0)^2 + A13*(c1-c1_0)*(c2-c2_0)^3 + A04*(c2-c2_0)^4 +
                    A50*(c1-c1_0)^5 + A41*(c1-c1_0)^4*(c2-c2_0) + A32*(c1-c1_0)^3*(c2-c2_0)^2 + A23*(c1-c1_0)^2*(c2-c2_0)^3 + A14*(c1-c1_0)*(c2-c2_0)^4 + A05*(c2-c2_0)^5 +
                    A60*(c1-c1_0)^6 + A51*(c1-c1_0)^5*(c2-c2_0) + A42*(c1-c1_0)^4*(c2-c2_0)^2 + A33*(c1-c1_0)^3*(c2-c2_0)^3 + A24*(c1-c1_0)^2*(c2-c2_0)^4 + A15*(c1-c1_0)*(c2-c2_0)^5 + A06*(c2-c2_0)^6 +
                    chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2))'
        #      
        derivative_order = 2                 
    []
    # beta penalty term
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = f_beta
        coupled_variables = 'c1 c2'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'beta*(1/c1 + 1/c2 + 1/(1 - c1 - c2))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c2'
        sum_materials = 'f_mix  f_beta'
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

    petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew '

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    # petsc_options = '-pc_svd_monitor -ksp_view'
    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'svd'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
    nl_max_its = 30
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-7
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e-4 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    # off_diagonal_
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
        file_base = output/3phase_5
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3phase_5
    []
[]

# [Debug]
#   show_var_residual_norms = true
# []