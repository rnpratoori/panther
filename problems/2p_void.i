n = 100     # number of elements per side
d = 1       # ND size of the side
a = 0.5     # type A monomer density
M = 1       # Initial mobility, depends on swell ratio
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

# Flory-Huggins approximation
chi12 = 1.0   # Flory-Huggins parameter
# chi13 = 10.0   # Flory-Huggins parameter
# chi23 = 10.0   # Flory-Huggins parameter
# N1 = 5       # Degree of polymerisation
# N2 = 5       # Degree of polymerisation
# N3 = 100     # Penalty term for void
# R = 1  # Universal gas constant
# T = 1 # Temperature in Kelvin
p = -1.38629e-1      # 0th coefficient of taylor function
q = 0               # 1st coefficient of taylor function
r = 0.4            # 2nd coefficient of taylor function
s = 0               # 3rd coefficient of taylor function
t = 2.66667e-1      # 4th coefficient of taylor function
u = 0               # 5th coefficient of taylor function
v = 4.26667e-1      # 6th coefficient of taylor function
c0 = 0.5
beta = 1e-3
delta = 0

[GlobalParams]
  block = 0
[]

[Mesh]
    # [2p]
        # generate a 2D mesh
        type = GeneratedMesh
        dim = 2
        nx = ${n}
        ny = ${n}
        xmax = ${d}
        ymax = ${d}
        add_subdomain_ids = '0  1'
    # []
[]

[MeshModifiers]
  [void]
    type = CoupledVarThresholdElementSubdomainModifier
    coupled_var = eta
    criterion_type = ABOVE
    subdomain_id = 1
    threshold = 1e-6
    execute_on = 'INITIAL TIMESTEP_BEGIN'
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
[]

[ICs]
    [pvfIC_1]
        type = RandomIC
        variable = c
        min = '${fparse a-0.04}'
        max = '${fparse a+0.04}'
        seed = 123
        # # distribution = Normal_a
        # coupled = eta
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
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    # AC variable
    [eta]
        order = FIRST
        family = LAGRANGE
        [InitialCondition]
          type = LatticeSmoothCircleIC
            variable = eta
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
        kappa_name = kappa
        w = w
    []
[]

[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa'
        interfacial_vars = 'c'
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
        prop_names = 'kappa'
        prop_values = '${fparse k}'
    []
    [mobility]
        type = DerivativeParsedMaterial
        property_name = M
        coupled_variables = 'c  eta'
        constant_names = 'M'
        constant_expressions = '${M}'
        expression = '(M*(1-c)^2)*(1-eta)'
        # derivative_order = 2
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c'
        constant_names =        'p      q       r       s
                                t      u        v
                                c0     chi'
        constant_expressions = '${p}    ${q}    ${r}    ${s}
                                ${t}    ${u}    ${v}
                                ${c0}   ${chi12}'
        expression = 'p + q*(c-c0) + r*(c-c0)^2 + s*(c-c0)^3 + t*(c-c0)^4 + u*(c-c0)^5 + v*(c-c0)^6 + chi*c*(1-c)'
        derivative_order = 2
    []
    # beta penalty term
    [beta_penalty]
        type = DerivativeParsedMaterial
        property_name = f_beta
        coupled_variables = 'c'
        constant_names = 'beta'
        constant_expressions = '${beta}'
        expression = 'beta*(1/c + 1/(1-c))'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c'
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

    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly      ilu          1'

    line_search = 'basic'

    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'lu'

    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    l_max_its = 200
    nl_max_its = 100
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-8
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e0 # seconds

    # # Automatic scaling for u and w
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
        file_base = output/2p_void
        time_step_interval = 1
        execute_on = 'TIMESTEP_END INITIAL FINAL'
    []
    [csv]
        type = CSV
        file_base = output/2p_void
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []