nx = 100     # number of elements in x
ny = 102     # number of elements in y
dx = 1.00       # ND size of the side in x
dy = 1.02       # ND size of the side in y
chi12 = 2.0   # Flory-Huggins parameter
chi13 = 0.1   # Flory-Huggins parameter
chi23 = 0.1   # Flory-Huggins parameter
N1 = 5       # Degree of polymerisation
N2 = 5       # Degree of polymerisation
N3 = 1       # Degree of polymerisation
M = 1e-3       # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

R = 1  # Universal gas constant
T = 1 # Temperature in Kelvin
beta = ${fparse 0e-3*R*T}
delta = 1e-2

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
        type = CoupledValueFunctionIC
        function = c_2phase
        variable = c2
        v = c1
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

[Functions]
    [c_2phase]
      type = ParsedFunction
      expression = '1 - x - ${delta}'
    []
[]

[UserObjects]
  [2phase]
    type = SolutionUserObject
    mesh = 'output/2phase_fh.e'
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
        constant_names = 'M     s       delta'
        constant_expressions = '${M} ${s}   ${delta}'
        expression = '(M*exp(1*(1-c1-c2-delta)-0.2))/s'
        # derivative_order = 2
    []
    [mobility2]
        type = DerivativeParsedMaterial
        property_name = M2
        coupled_variables = 'c1 c2'
        constant_names = 'M     s   delta'
        constant_expressions = '${M} ${s}   ${delta}'
        expression = '(M*exp(1*(1-c1-c2-delta)-0.2))/s'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c1 c2'
        constant_names = 'R      T       chi12      chi13       chi23     N1        N2      N3       s     beta'
        constant_expressions = '${R}    ${T}    ${chi12}    ${chi13}    ${chi23}    ${N1}   ${N2}   ${N3}    ${s}    ${beta}'
        expression = 's*(R*T*(c1*log(c1)/N1 + c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2)))'
        # expression = 'if (c1>0, if(c2>0, if(1-c1-c2>0, s*(R*T*(c1*log(c1)/N1 + c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi12*c1*c2 + chi13*c1*(1-c1-c2) + chi23*c2*(1-c1-c2))), s*(R*T*(c1*log(c1)/N1 + c2*log(c2)/N2 + chi12*c1*c2))),
            # if(1-c1-c2>0, s*(R*T*(c1*log(c1)/N1 + (1-c1-c2)*log(1-c1-c2)/N3 + chi13*c1*(1-c1-c2))), s*R*T*(c1*log(c1)/N1))),
            # if(c2>0, if(1-c1-c2>0, s*(R*T*(c2*log(c2)/N2 + (1-c1-c2)*log(1-c1-c2)/N3 + chi23*c2*(1-c1-c2))), s*R*T*(c2*log(c2)/N2)),
            # if(1-c1-c2>0, s*R*T*((1-c1-c2)*log(1-c1-c2)/N3), 0)))'
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
        # expression = 'if (c1>0, if(c2>0, if(1-c1-c2>0, beta*(1/c1 + 1/c2 + 1/(1 - c1 - c2)), beta*(1/c1 + 1/c2)),
                    # if(1-c1-c2>0, beta*(1/c1 + 1/(1 - c1 - c2)), beta*(1/c1))),
                    # if(c2>0, if(1-c1-c2>0, beta*(1/c2 + 1/(1 - c1 - c2)), beta*(1/c2)),
                    # if(1-c1-c2>0, beta*(1/(1-c1-c2)), 0)))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c1 c2'
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
        file_base = output/3p_dis_fh
        time_step_interval = 1
        execute_on = 'TIMESTEP_BEGIN INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/3p_dis_fh
    []
[]

# [Debug]
#   show_var_residual_norms = true
# []