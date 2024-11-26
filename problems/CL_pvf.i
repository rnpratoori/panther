# Unit conversions
# ev_J = 6.24e18    # Coversion of energy
ev_J = 1
d_f = 1e6           # factor to convert to m from the chosen units
s = 1e-0            # Scaling factor

# Simulation parameters
n = 050     # number of elements per side
d = ${fparse d_f*1e-3}          # size of the side

# System parameters
a = 0.81    # type A monomer density
chi = 2.4   # Flory-Huggins parameter
N1 = 23     # Segment number for polymer
N2 = 5      # Segment number for additive
R = 8.314   # Universal gas constant
T = 453     # Temperature in Kelvin
M = ${fparse d_f^5*6.52e-17}    # Initial mobility, depends on swell ratio
k = ${fparse 4.5e-5/d_f}        # gradient energy coefficient
eps = ${fparse d_f*1e-4}        # interface width
E = 8e6                         # Elastic modulus

[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}
    ymax = ${d}
    uniform_refine = 2
[]
  
[Variables]
    # polymer volume fraction
    [./c]
        order = FIRST
        family = LAGRANGE
    [../]
    # Chemical potential (nJ/mol)
    [./w]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[ICs]
    [pvfIC]
        type = RandomIC
        variable = c
        seed = 123
        min = ${fparse a-0.05}
        max = ${fparse a+0.05}
    []
[]

[AuxVariables]
    [./f_density]
        order = CONSTANT
        family = MONOMIAL
    [../]
[]
  
[Kernels]
    [./w_dot]
        type = CoupledTimeDerivative
        variable = w
        v = c
    [../]
    # adding nonlocal term to the energy
    [./coupled_res]
        type = SplitCHWRes
        variable = w
        mob_name = M
    [../]
    [./coupled_parsed]
        type = SplitCHParsed
        variable = c
        f_name = f_mix
        kappa_name = kappa
        w = w
    [../]
[]
  
[AuxKernels]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [./f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa'
        interfacial_vars = c
    [../]
[]
  
[BCs]
    [./Periodic]
        [./all]
            auto_direction = 'x y'
        [../]
    [../]
[]
  
[Materials]
    [./mat]
        type = GenericFunctionMaterial
        prop_names  = 'M   kappa'
        prop_values = '${fparse M/ev_J/s} ${fparse k*ev_J*s}'
    [../]
    # free energy density function (nJ/mol/nm^2)
    # local energy as a double well potential
    [./local_energy]
        type = DerivativeParsedMaterial
        property_name = f_loc
        coupled_variables = 'c'
        constant_names =        'W1    eps'
        constant_expressions =  '1/4    ${eps}'
        expression = '(W1/eps^2)*${ev_J}*${s}*c^2*(c - 1)^2'
        derivative_order = 2
    [../]
    # mixing energy based on 
    # Flory-Huggins theory
    [./mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'c'
        constant_names =        'R      T       N1      N2      s       sw
                                chi     ev_J    d_f'
        constant_expressions = '${R}    ${T}   ${N1}    ${N2}   ${s}    1.0
                                ${chi}  ${ev_J} ${d_f}'
        expression = 's*ev_J*(R*T/d_f^3)*((sw*c*log(sw*c))/N1
                    +((1-sw*c)*log(1-sw*c))/N2+(chi*sw*c*(1-sw*c)))'
        derivative_order = 2
    [../]
    # elastic energy
    [./elastic_energy]
        type = DerivativeParsedMaterial
        property_name = f_el
        coupled_variables = 'c'
        constant_names =        'E      s'
        constant_expressions =  '${E}   ${s}'
        expression = 's*(E/3)*(1/(c^(2/3))-1)'
        derivative_order = 2
    [../]
    # Total free energy
    # Sum of all the parts
    [./free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'c'
        sum_materials = 'f_mix'
        derivative_order = 2
    [../]
[]
  
[Postprocessors]
    # Calculate total free energy at each timestep
    [./total_energy]
        type = ElementIntegralVariablePostprocessor
        variable = f_density
        execute_on = 'initial timestep_end'
    [../]
    [./nodes]                 # Number of nodes in mesh
        type = NumNodes
    [../]
[]
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2
  
    petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type
                            -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly
                               lu           2'
  
    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'
  
    l_tol = 1e-6
    l_abs_tol = 1e-9
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-9
    
    [./TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    [../]
  
    end_time = 31536000 # seconds

    # # Automatic scaling for c and w
    # automatic_scaling = true
    # scaling_group_variables = 'c w'
  
    [./Adaptivity]
      coarsen_fraction = 0.1
      refine_fraction = 0.7
      max_h_level = 2
    [../]
[]
  
[Outputs]
    [ex_keya_2]
        type = Exodus
        time_step_interval = 1
    []
    [csv_keya_2]
        type = CSV
    []
[]