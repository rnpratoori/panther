# Three-phase void transport/dissolution case with adaptive refinement.
# c1/c2 evolve by Cahn-Hilliard while eta tracks the void phase.
nx = 50     # number of elements in x
ny = 75     # number of elements in y
dx = 2.00       # ND size of the side in x
dy = 3.00       # ND size of the side in y
M = 1e0       # Initial mobility, depends on swell ratio
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient
Cn_eta = 1e-2  # Cahn number for AC variable
k_eta = ${fparse Cn_eta^2}    # gradient energy coefficient for AC

# Flory-Huggins constants used to build the Taylor approximation below.
chi12 = 1.0   # Flory-Huggins parameter
chi13 = 0.004   # Flory-Huggins parameter
chi23 = 0.1   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 1       # Degree of polymerisation
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
# Coefficients for the two-variable Taylor-expanded ternary free energy.
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
# delta_c = 0.025
# delta_eta = 0

W = 100
A = 5

[Mesh]
    uniform_refine = 3
    add_subdomain_ids = '1'
    [2d]
        # generate a 2D mesh
        # type = DistributedRectilinearMeshGenerator
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
    []
    [c3_subdomain]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        combinatorial_geometry = 'y > 1'
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
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    # Chemical potential (nJ/mol)
    [w2]
        order = FIRST
        family = LAGRANGE
    []
    # void variable
    [eta]
        order = FIRST
        family = LAGRANGE
    []
    [weta]
        order = FIRST
        family = LAGRANGE
    []
    [ct]
        order = FIRST
        family = LAGRANGE
    []
[]

[ICs]
    [c1]
        type = SolutionIC
        from_variable = 'c1'
        solution_uo = 2phase
        variable = c1
        from_subdomains = '0 1 2'
    []
    [c2]
        type = SolutionIC
        from_variable = 'c2'
        solution_uo = 2phase
        variable = c2
        from_subdomains = '0 1 2'
    []
    [eta]
        type = SolutionIC
        from_variable = 'eta'
        solution_uo = 2phase
        variable = eta
        from_subdomains = '0 1 2'
    []
    [ct]
        type = ConstantIC
        value = 0
        variable = ct
    []
[]

[UserObjects]
    [2phase]
        type = SolutionUserObject
        mesh = /work/mech-ai/rnp/MOOSE/projects/panther/problems/ic/square/3pv_0.3_ic_0.10_0.2.e
        system_variables = 'c1 c2 eta'
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
    [c3]
        order = FIRST
        family = LAGRANGE
    []
    [voids]
        order = CONSTANT
        family = MONOMIAL
    []
    [ct]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    # Cahn-Hilliard equation for polymer 1
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
        coupled_variables = 'c2 eta'
        f_name = f_tot
        kappa_name = kappa
        w = w1
    []
    # Cahn-Hilliard equation for polymer 2
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
        coupled_variables = 'c1 eta'
        f_name = f_tot
        kappa_name = kappa
        w = w2
    []
    # Cahn-Hilliard equation for void filling
    [weta_dot]
        type = CoupledTimeDerivative
        variable = weta
        v = eta
    []
    [coupled_reseta]
        type = SplitCHWRes
        variable = weta
        mob_name = Meta
    []
    [coupled_parsedeta]
        type = SplitCHParsed
        variable = eta
        coupled_variables = 'c3'
        f_name = f_sol
        kappa_name = kappa
        w = weta
    []
    # # Allen-Cahn for neotissue formation at c1
    # [c1t_bulk]
    #     type = AllenCahn
    #     f_name = f_t
    #     variable = c1
    # []
    # [c1t_interface]
    #     type = ACInterface
    #     variable = c1
    #     kappa_name = kappa_c1t
    # []
    # [c1t_time]
    #     type = TimeDerivative
    #     variable = c1
    # []
    # # Allen-Cahn for neotissue formation at c2
    # [c2t_bulk]
    #     type = AllenCahn
    #     f_name = f_t
    #     variable = c2
    # []
    # [c2t_interface]
    #     type = ACInterface
    #     variable = c2
    #     kappa_name = kappa_c2t
    # []
    # [c2t_time]
    #     type = TimeDerivative
    #     variable = c2
    # []
    # Allen-Cahn for neotissue formation
    [ct_bulk]
        type = AllenCahn
        f_name = f_t
        variable = ct
    []
    [ct_interface]
        type = ACInterface
        variable = ct
        kappa_name = kappa_ct
    []
    [ct_time]
        type = TimeDerivative
        variable = ct
    []
    # Sink term
    [c1_sink]
        type = ADReaction
        variable = c1
        reaction_rate = kr
    []
    [c2_sink]
        type = ADReaction
        variable = c2
        reaction_rate = kr
    []
    # source term
    [ct_source]
        type = ADReaction
        variable = ct
        reaction_rate = -kr
        coupled_variables = 'c1 c2'
    []
[]

[AuxKernels]
    # # calculate energy density from local and gradient energies (J/mol/mum^2)
    # [f_density]
    #     type = TotalFreeEnergy
    #     variable = f_density
    #     f_name = 'f_tot'
    #     kappa_names = 'kappa kappa'
    #     interfacial_vars = 'c1 c2'
    # []
    # # calculate interfacial energy density
    # [f_int_density]
    #     type = ParsedAux
    #     variable = f_int_density
    #     coupled_variables = 'f_density'
    #     material_properties = 'f_tot'
    #     expression = 'f_density - f_tot'
    # []
    # calculate solvent phase volume fraction
    [c3_polymer]
        type = ParsedAux
        variable = c3
        coupled_variables = 'c1 c2 eta'
        expression = 'if(eta<0, 1 - c1 - c2, 1 - c1 - c2 - (1+eta)/2)'
        execute_on = 'INITIAL TIMESTEP_END'
    []
    # count number of active voids
    [voids]
        type = FeatureFloodCountAux
        variable =  voids
        flood_counter = 'voids'
        execute_on = 'INITIAL TIMESTEP_END'
    []
[]

[Materials]
    # gradient energy coefficients
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa kappa_eta kappa_ct'
        prop_values = '${fparse k} ${fparse k_eta} ${fparse k}'
    []
    # mobility for void filling
    [mobility_eta]
        type = DerivativeParsedMaterial
        property_name = Meta
        coupled_variables = 'c3  eta'
        constant_names = 'M'
        constant_expressions = '${M}'
        expression = 'if(c3>0.05, M*1e5*(c3-0.05), 0)'
    []
    # mobility for polymers
    [mobility1]
        type = DerivativeParsedMaterial
        property_name = M1
        coupled_variables = 'c1 c2 c3 eta'
        constant_names = 'M'
        constant_expressions = '${M}'
        expression = 'if(eta<0, min((M*exp((15*c3-3)))*(1-eta)/2,1), 0)'
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2 c3 eta'
        constant_names = 'M'
        constant_expressions = '${M}'
        # expression = 'if(eta<0, if(c3<0.2, (M*exp((11*c3-2.2)))*(1-eta)/2, (M*exp((15*c3-3)))*(1-eta)/2), 0)'
        expression = 'if(eta<0, min((M*exp((11*c3-2.2)))*(1-eta)/2,1), 0)'
    []
    # mobility for tissue formation
    [mobility_ct]
        type = GenericConstantMaterial
        prop_names = L
        prop_values = 1
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = 'f_mix'           
        coupled_variables = 'c1 c2 c3'
        constant_names =      'A00    A10    A20    A30    A40    A50    A60
                                A01    A11    A21    A31    A41    A51    A61
                                A02    A12    A22    A32    A42    A52    A62
                                A03    A13    A23    A33    A43    A53    A63
                                A04    A14    A24    A34    A44    A54    A64
                                A05    A15    A25    A35    A45    A55    A65
                                A06    A16    A26    A36    A46    A56    A66
                                c1_0   c2_0  chi12   chi13   chi23'
        constant_expressions = '${A00} ${A10} ${A20} ${A30} ${A40} ${A50} ${A60}
                                ${A01} ${A11} ${A21} ${A31} ${A41} ${A51} ${A61}
                                ${A02} ${A12} ${A22} ${A32} ${A42} ${A52} ${A62}
                                ${A03} ${A13} ${A23} ${A33} ${A43} ${A53} ${A63}
                                ${A04} ${A14} ${A24} ${A34} ${A44} ${A54} ${A64}
                                ${A05} ${A15} ${A25} ${A35} ${A45} ${A55} ${A65}
                                ${A06} ${A16} ${A26} ${A36} ${A46} ${A56} ${A66}
                                ${c1_0} ${c2_0} ${chi12} ${chi13} ${chi23}'
        expression = 'A00 +
                    A10*(c1-c1_0) + A01*(c2-c2_0) +
                    A20*(c1-c1_0)^2 + A11*(c1-c1_0)*(c2-c2_0) + A02*(c2-c2_0)^2 +
                    A30*(c1-c1_0)^3 + A21*(c1-c1_0)^2*(c2-c2_0) + A12*(c1-c1_0)*(c2-c2_0)^2 + A03*(c2-c2_0)^3 +
                    A40*(c1-c1_0)^4 + A31*(c1-c1_0)^3*(c2-c2_0) + A22*(c1-c1_0)^2*(c2-c2_0)^2 + A13*(c1-c1_0)*(c2-c2_0)^3 + A04*(c2-c2_0)^4 +
                    A50*(c1-c1_0)^5 + A41*(c1-c1_0)^4*(c2-c2_0) + A32*(c1-c1_0)^3*(c2-c2_0)^2 + A23*(c1-c1_0)^2*(c2-c2_0)^3 + A14*(c1-c1_0)*(c2-c2_0)^4 + A05*(c2-c2_0)^5 +
                    A60*(c1-c1_0)^6 + A51*(c1-c1_0)^5*(c2-c2_0) + A42*(c1-c1_0)^4*(c2-c2_0)^2 + A33*(c1-c1_0)^3*(c2-c2_0)^3 + A24*(c1-c1_0)^2*(c2-c2_0)^4 + A15*(c1-c1_0)*(c2-c2_0)^5 + A06*(c2-c2_0)^6 +
                    chi12*c1*c2 + chi13*c1*c3 + chi23*c2*c3'
        derivative_order = 2
    []
    # beta penalty term
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = f_beta
        coupled_variables = 'c1 c2 c3 eta'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'if(eta<0, if(c3>1e-3, beta*(1/c1 + 1/c2 + 1/c3), beta*(1/c1 + 1/c2)), 0)'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c2 c3 eta'
        sum_materials = 'f_mix'
        derivative_order = 2
    []
    # CH free energy for void filling
    [ac_energy]
        type = DerivativeParsedMaterial
        property_name = f_sol
        coupled_variables = 'c3 eta'
        constant_names = 'W A'
        constant_expressions = '${W} ${A}'
        expression = '(1-c3)*(1+eta)^4'
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
    # Calculate interfacial energy at each timestep
    [interfacial_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_int_density
        execute_on = 'initial timestep_end'
    []
    # Calculate total elapsed time
    [elapsed]
        type = PerfGraphData
        section_name = "Root"
        data_type = total
    []
    [step_size]             # Size of the time step
        type = TimestepSize
    []
    [nodes]                 # Number of nodes in mesh
        type = NumNodes
    []
    # Calculate number of voids
    [voids]
        type = FeatureFloodCount
        variable = eta
        threshold = 0
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2

    petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew -ksp_monitor_cancel'

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    l_tol = 1e-6
    l_max_its = 2000
    l_abs_tol = 1e-10
    nl_max_its = 30
    nl_abs_tol = 1e-10

    dtmax = 2e-9
    
    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-10
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 3e-4 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    scaling_group_variables = 'c1 c2; w1 w2'

    # [Adaptivity]
    #     coarsen_fraction = 0.1
    #     refine_fraction = 0.7
    #     max_h_level = 3
    # []
[]

[Adaptivity]
    marker = marker_c3
    max_h_level = 2
    stop_time = 1e-7
    [Indicators]
        [indicator_c3]
            type = GradientJumpIndicator
            variable = c3
            block = 1
        []
        [indicator_eta]
            type = GradientJumpIndicator
            variable = eta
        []
    []
    [Markers]
        [marker_c3]
            type = ErrorFractionMarker
            indicator = indicator_c3
            block = 1
            coarsen = 0.2
            refine = 0.7
        []
        [marker_eta]
            type = ErrorFractionMarker
            indicator = indicator_eta
            coarsen = 0.2
            refine = 0.7
        []
        [combined_marker]
            type = ComboMarker
            markers = 'marker_c3 marker_eta'
        []
    []
[]

[Times]
    [out_times]
        type = CSVFileTimes
        files = ../../out_times.csv
    []
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = 3pv_0.3_0.10_0.2_out
        time_step_interval = 100
        execute_on = 'INITIAL FINAL TIMESTEP_END'
    []
    [csv]
        type = CSV
        file_base = 3pv_0.3_0.10_0.2_out
    []
    [ex_mech]
        type = Exodus
        file_base = /work/mech-ai-scratch/rnp/output/polydeg_1/square/ic/3pv_0.3_0.10_0.2_mechic
        execute_on = 'TIMESTEP_END'
        sync_only = true
        sync_times_object = out_times
    []
[]
