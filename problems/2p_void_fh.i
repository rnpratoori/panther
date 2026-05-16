rc = 0.10
dc = 0.4

nx = 200     # number of elements per side
ny = 100     # number of elements per side
dx = 2       # ND size of the side
dy = 1       # ND size of the side
a = 0.5     # type A monomer density
M = 1       # Initial mobility, depends on swell ratio
Cn = 5e-2  # Cahn number
k = ${fparse Cn^2}    # gradient energy coefficient

# Flory-Huggins approximation
chi = 1.0   # Flory-Huggins parameter
N1 = 5       # Degree of polymerisation
N2 = 5       # Degree of polymerisation

[GlobalParams]
  block = 0
[]

[Mesh]
    add_subdomain_ids = '1'
    [2p]
        # generate a 2D mesh
        type = DistributedRectilinearMeshGenerator
        dim = 2
        nx = ${nx}
        ny = ${ny}
        xmax = ${dx}
        ymax = ${dy}
    []
[]

[MeshModifiers]
    [void]
        type = CoupledVarThresholdElementSubdomainModifier
        coupled_var = eta
        criterion_type = ABOVE
        subdomain_id = 1
        complement_subdomain_id = 0
        threshold = 0
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
[]

[Variables]
    # polymer volume fraction
    [c]
        order = FIRST
        family = LAGRANGE
        block = '0  1'
    []
    # Chemical potential (nJ/mol)
    [w]
        order = FIRST
        family = LAGRANGE
    []
    # void variable
    [eta]
        order = FIRST
        family = LAGRANGE
        block = '0  1'
    []
[]

[ICs]
    [c]
        type = RandomIC
        variable = c
        seed = 123
        distribution = Normal_a
    []
    [eta]
        type = SpecifiedSmoothCircleIC
        radii =         '0.10 0.10 0.10
                        0.10 0.10 0.10'
        x_positions =   '0.20 0.80 1.40
                        0.20 0.80 1.40'
        y_positions =   '0.80 0.80 0.80
                        0.20 0.20 0.20'
        z_positions =   '0 0 0
                        0 0 0'
        variable = eta
        invalue = ${fparse 1.0}
        outvalue = ${fparse -1.0}
        int_width = 0
    []
[]

[Distributions]
    [Normal_a]
        type = Normal
        mean = ${a}
        standard_deviation = 0.04
    []
[]

[AuxVariables]
    [c2]
        order = FIRST
        family = LAGRANGE
    []
    # Variables to be read in dissolution simulation
    [c1_total]
        order = FIRST
        family = LAGRANGE
        block = '0  1'
    []
    [c2_total]
        order = FIRST
        family = LAGRANGE
        block = '0  1'
    []
    [cv_total]
        order = FIRST
        family = LAGRANGE
        block = '0  1'
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
    [null]
        type = NullKernel
        variable = eta
        block = '0  1'
    []
    [null_c]
        type = NullKernel
        variable = c
        block = 1
    []
[]

[AuxKernels]
    # calculate c2 from c
    [c2]
        type = ParsedAux
        variable = c2
        coupled_variables = 'c'
        expression = '1 - c'
        block = 0
        execute_on = 'INITIAL'
    []
    [c1_total]
        type = ParsedAux
        variable = c1_total
        coupled_variables = 'c  eta'
        expression = 'if(eta<0., c, 0)'
        block = '0  1'
    []
    [c2_total]
        type = ParsedAux
        variable = c2_total
        coupled_variables = 'c  eta'
        expression = 'if(eta<0., 1 - c, 0)'
        block = '0  1'
    []
    [cv_total]
        type = ParsedAux
        variable = cv_total
        coupled_variables = 'eta'
        expression = 'if(eta>=0., 1, 0)'
        block = '0  1'
    []
[]

[Materials]
    [mat]
        type = GenericFunctionMaterial
        prop_names = 'kappa'
        prop_values = '${fparse k}'
        block = '0  1'
    []
    [mobility]
        type = DerivativeParsedMaterial
        property_name = M
        coupled_variables = 'c  eta'
        constant_names = 'M'
        constant_expressions = '${M}'
        expression = '(M*16*c^2*(1-c)^2)*(1-eta)/2'
    []
    # mixing energy based on
    # Flory-Huggins theory
    [mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c'
        constant_names =        'chi    N1     N2'
        constant_expressions = '${chi}    ${N1}    ${N2}'
        expression = 'c*log(c)/N1 + (1-c)*log(1-c)/N2 + chi*c*(1-c)'
        derivative_order = 2
    []
    # Total free energy
    # Sum of all the parts
    [free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c'
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

    petsc_options = '-ksp_converged_reason -snes_converged_reason -snes_ksp_ew -ksp_monitor_cancel'

    # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31                  preonly      ilu          1'

    petsc_options_iname = '-pc_type -ksp_type -pc_factor_mat_solver_type'
    petsc_options_value = 'lu       preonly   mumps'

    line_search = 'basic'

    l_tol = 1e-10
    l_abs_tol = 1e-10
    nl_max_its = 30
    nl_abs_tol = 1e-10

    [TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    []

    end_time = 1e0 # seconds

    # Automatic scaling for u and w
    automatic_scaling = true
    scaling_group_variables = 'c w'
[]

[Outputs]
    [ex]
        type = Exodus
        file_base = ic_2pv/2pv_${a}_ic_${rc}_${dc}_sharp
        time_step_interval = 10
        execute_on = 'INITIAL TIMESTEP_END FINAL'
    []
[]