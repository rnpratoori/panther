nx = 100     # number of elements in x
ny = 102     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.00       # ND size of the side in y
# a = 0.3     # type A monomer density
# b = 0.3     # type B monomer density
M = 1e-0       # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

# chi12 = 2.0   # Flory-Huggins parameter
# chi13 = 0.1   # Flory-Huggins parameter
# chi23 = 0.1   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 1       # Degree of polymerisation
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
A00 = 3.35203e-1    # changes
A10 = 1.15888e-1
A20 = 1.3           # changes
A30 = 1.33333e-1
A40 = 1.73333
A50 = -1.76
A60 = 7.89333
A01 = 1.15888e-1
A11 = 3.8           # changes
A21 = 2
A31 = 2.66667
A41 = 4
A51 = 6.4
A61 = 1.06667e1
A02 = 1.3           # changes
A12 = 2
A22 = 4
A32 = 8
A42 = 1.6e1
A52 = 3.2e1
A62 = 6.4e1
A03 = 1.33333e-1
A13 = 2.66667
A23 = 8
A33 = 2.13333e1
A43 = 5.33333e1
A53 = 1.28e2
A63 = 2.98667e2
A04 = 1.73333
A14 = 4
A24 = 1.6e1
A34 = 5.33333e1
A44 = 1.60e2
A54 = 4.48e2
A64 = 1.19467e3
A05 = -1.76
A15 = 6.4
A25 = 3.2e1
A35 = 1.28e2
A45 = 4.48e2
A55 = 1.4336e3
A65 = 4.3008e3
A06 = 7.89333
A16 = 1.06667e1
A26 = 6.4e1
A36 = 2.98667e2
A46 = 1.19467e3
A56 = 4.3008e3
A66 = 1.4336e4
c1_0 = 0.25
c3_0 = 0.25
beta = 1.0e-3       # Stability parameter
delta = 0

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
        combinatorial_geometry = 'y > 100/${ny}'
        input = 2d
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
    [c1]
        type = SolutionIC
        from_variable = 'c'
        solution_uo = 2phase
        variable = c1
        block = 0
    []
    [c3]
        type = ConstantIC
        value = ${delta}
        variable = c3
        block = 0
    []
    [top_c1]
        type = ConstantIC
        value = ${delta}
        variable = c1
        block = 1
    []
    [top_c3]
        type = ConstantIC
        value = ${fparse 1-delta}
        variable = c3
        block = 1
    []
[]

[UserObjects]
  [2phase]
    type = SolutionUserObject
    mesh = 'output/2phase_taylor.e'
    system_variables = 'c'
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
        kappa_name = kappa
        w = w3
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa kappa'
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

[BCs]
    [top1]
        type = DirichletBC
        variable = c1
        boundary = 2
        # block = 1
        value = ${delta}
    []
    [top3]
        type = DirichletBC
        variable = c3
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
        coupled_variables = 'c1 c3'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        # expression = '(M*4*c1*(1-c1))/s'
        expression = '(M)/s'
        # derivative_order = 2
    []
    [mobility3]
        type = DerivativeParsedMaterial
        property_name = M3
        coupled_variables = 'c1 c3'
        constant_names = 'M     s'
        constant_expressions = '${M} ${s}'
        expression = '(M)/s'
        # expression = '(M*4*c3*(1-c3))/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c3'
        constant_names =       'A00    A10    A20    A30    A40    A50    A60
                               A01    A11    A21    A31    A41    A51    A61
                               A02    A12    A22    A32    A42    A52    A62
                               A03    A13    A23    A33    A43    A53    A63
                               A04    A14    A24    A34    A44    A54    A64
                               A05    A15    A25    A35    A45    A55    A65
                               A06    A16    A26    A36    A46    A56    A66    
                               s      c1_0   c3_0   beta'
        constant_expressions = '${A00} ${A10} ${A20} ${A30} ${A40} ${A50} ${A60}
                               ${A01} ${A11} ${A21} ${A31} ${A41} ${A51} ${A61}
                               ${A02} ${A12} ${A22} ${A32} ${A42} ${A52} ${A62}
                               ${A03} ${A13} ${A23} ${A33} ${A43} ${A53} ${A63}
                               ${A04} ${A14} ${A24} ${A34} ${A44} ${A54} ${A64}
                               ${A05} ${A15} ${A25} ${A35} ${A45} ${A55} ${A65}
                               ${A06} ${A16} ${A26} ${A36} ${A46} ${A56} ${A66}
                               ${s}   ${c1_0} ${c3_0} ${beta}'
        expression = 's*(A00 + 
                    A10*(c1-c1_0) + A01*(c3-c3_0) + 
                    A20*(c1-c1_0)^2 + A11*(c1-c1_0)*(c3-c3_0) + A02*(c3-c3_0)^2 +
                    A30*(c1-c1_0)^3 + A21*(c1-c1_0)^2*(c3-c3_0) + A12*(c1-c1_0)*(c3-c3_0)^2 + A03*(c3-c3_0)^3 + 
                    A40*(c1-c1_0)^4 + A31*(c1-c1_0)^3*(c3-c3_0) + A22*(c1-c1_0)^2*(c3-c3_0)^2 + A13*(c1-c1_0)*(c3-c3_0)^3 + A04*(c3-c3_0)^4 + 
                    A50*(c1-c1_0)^5 + A41*(c1-c1_0)^4*(c3-c3_0) + A32*(c1-c1_0)^3*(c3-c3_0)^2 + A23*(c1-c1_0)^2*(c3-c3_0)^3 + A14*(c1-c1_0)*(c3-c3_0)^4 + A05*(c3-c3_0)^5 + 
                    A60*(c1-c1_0)^6 + A51*(c1-c1_0)^5*(c3-c3_0) + A42*(c1-c1_0)^4*(c3-c3_0)^2 + A33*(c1-c1_0)^3*(c3-c3_0)^3 + A24*(c1-c1_0)^2*(c3-c3_0)^4 + A15*(c1-c1_0)*(c3-c3_0)^5 + A06*(c3-c3_0)^6
                    )'
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

    line_search = 'basic'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-8
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    # dt = 1.0e-8

    end_time = 1e0 # seconds

    # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'c1 c2; w1 w2'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 2
    # []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = output/3p_dis_taylor
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3p_dis_taylor
    []
[]

# [Debug]
#   show_var_residual_norms = true
# []