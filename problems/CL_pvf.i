# Unit conversions
# ev_J = 6.24e18    # Coversion of energy
ev_J = 1
d_f = 1e0           # factor to convert to m from the chosen units
s = 1e-5            # Scaling factor

# Simulation parameters
n = 100     # number of elements per side
d = ${fparse d_f*5e-3}          # size of the side

# System parameters
a = 0.88    # type A monomer density
chi = 8e-3  # Flory-Huggins parameter
N1 = 1000   # Segment number for elastomer matrix (very large number)
N2 = 87.7   # Segment number for silicon fluid (DPDM-005-088)
R = 8.314   # Universal gas constant
T = 300     # Temperature in Kelvin
E = 0.5e6   # Elastic modulus (V31-151)
# nc = 3e-4   # crosslink density (V31-151)
M = ${fparse d_f^5*2.74e-15}    # Initial mobility, depends on swell ratio
k = ${fparse 1e-8/d_f}          # gradient energy coefficient
eps = ${fparse d_f*1e-4}        # interface width
vs = ${fparse (d_f*1e-2)^3*81.2}    # volume of repetition unit
V2 = ${fparse vs*N2}            # volume of solvent

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
        f_name = 'f_mix'
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
    # Flory-Huggins parameter
    # dependent on polymer volume fraction
    # [./chi]
    #     type = DerivativeParsedMaterial
    #     property_name = chi
    #     coupled_variables = 'c'
    #     constant_names =        'N2      vs     nc      s'
    #     constant_expressions = '${N2}   ${vs}   ${nc}   ${s}'
    #     expression = '(-1/(N2*c^2))*(log(1-c)+c+nc*vs*((1/c)-(c/2)))'
    #     derivative_order = 2
    #     # outputs = ex
    # [../]
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
        # material_property_names = 'chi'
        # constant_names =        'R      T       N1      N2      s       sw
        #                         ev_J    d_f     chi'
        constant_names =        'R      T       N1      V2      s       vs
                                ev_J    d_f     chi'
        constant_expressions = '${R}    ${T}   ${N1}    ${V2}   ${s}    ${vs}
                                ${ev_J} ${d_f}  ${chi}'
        expression = 's*ev_J*(R*T/(c*d_f^3))*((c*0*log(c))/N1
                    +((1-c)*log(1-c))/V2+(chi*c*(1-c)/vs))'
        # expression = 's*ev_J*(R*T/d_f^3)*((sw*c*log(sw*c))/N1
        #             +((1-sw*c)*log(1-sw*c))/N2+(chi*sw*c*(1-sw*c)))'
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
        sum_materials = 'f_loc f_mix f_el'
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
        dt = 1.0e0
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    [../]
  
    end_time = 31536000 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'u w'
  
    [./Adaptivity]
      coarsen_fraction = 0.1
      refine_fraction = 0.7
      max_h_level = 2
    [../]
[]
  
[Outputs]
    [ex_noPBC_new]
        type = Exodus
        time_step_interval = 1
    []
    [csv_noPBC_new]
        type = CSV
    []
[]
