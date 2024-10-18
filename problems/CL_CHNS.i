n = 400     # number of elements per side
d = 1       # ND size of the side
D = 1e+0    # actual size of the side in 1 mum
l = 1e-3    # Kuhn statistical length in 1 mum
a = 0.5     # type A monomer density
# chi = 0.077 # Flory-Huggins parameter
N = 300     # Degree of polymerisation
M_in = 4    # Initial mobility, depends on swell ratio
s = 1e+0    # Scaling factor

N2 = 30.1
N1 = 3000
E = 0.4e6   #
nc = 3e-4   # V33
vr = 78.5   # Considering V33 and DM-100-030 are similar
vs = 78.5   # DM-100-030
R = 8.314   # Universal gas constant
T = 300     # Temperature in Kelvin

V1 = ${fparse N1*vr}  # Volume
V2 = ${fparse N2*vs}


[Mesh]
    # generate a 2D mesh
    type = GeneratedMesh
    dim = 2
    nx = ${n}
    ny = ${n}
    xmax = ${d}  # nd
    ymax = ${d}  # nd
    # uniform_refine = 2
[]
  
[Variables]
    # difference in the volume fractions of the 2 phases
    [./u]
        order = FIRST
        family = LAGRANGE
        [./InitialCondition]
            type = RandomIC
            seed = 123
            min = ${fparse 2*(a-0.5)-0.1}
            max = ${fparse 2*(a-0.5)+0.1}
        [../]
    [../]
    # Chemical potential (nJ/mol)
    [./w]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    # # polymer volume fraction
    # [./pvf]
    # [../]
    # Local free energy density (nJ/mol)
    [./f_density]
        order = CONSTANT
        family = MONOMIAL
    [../]
[]
  
[Kernels]
    [./w_dot]
        type = CoupledTimeDerivative
        variable = w
        v = u
    [../]
    # adding nonlocal term to the energy
    [./coupled_res]
        type = SplitCHWRes
        variable = w
        mob_name = M
    [../]
    [./coupled_parsed]
        type = SplitCHParsed
        variable = u
        f_name = f_tot
        kappa_name = kappa
        w = w
    [../]
[]
  
[AuxKernels]
    # # calculate polymer volume fraction from difference in volume fractions
    # [./pvf]
    #     type = ParsedAux
    #     variable = phi
    #     coupled_variables = 'u'
    #     expression = '(u+1)/2'
    # [../]
    # calculate energy density from local and gradient energies (J/mol/mum^2)
    [./f_density]
        type = TotalFreeEnergy
        variable = f_density
        f_name = 'f_tot'
        kappa_names = 'kappa'
        interfacial_vars = u
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
    # All properties are non-dimensional
    # In the ND formulation, M is 1.0
    [./mat]
        type = GenericFunctionMaterial
        prop_names  = 'M   kappa'
        prop_values = '${fparse 1.0/s} ${fparse 0.1/s}'
    [../]
    # polymer volume fraction
    # defined for convenience
    [./pvf]
        type = DerivativeParsedMaterial
        property_name = phi
        coupled_variables = 'u'
        expression = '(u+1)/2'
        derivative_order = 2
    [../]
    # Flory-Huggins parameter
    # dependent on polymer volume fraction
    [./chi]
        type = DerivativeParsedMaterial
        property_name = chi
        coupled_variables = 'u'
        constant_names =        'N2     nc      vs      s'
        constant_expressions = '${N2}   ${nc}   ${vs}   ${s}'
        expression = '(-1/(N2*phi^2))*(log(1-phi)+phi+nc*vs*((1/phi)-(phi/2)))'
        derivative_order = 2
    [../]
    # free energy density function (nJ/mol/nm^2)
    # local energy as a double well potential
    [./local_energy]
        type = DerivativeParsedMaterial
        property_name = f_loc
        coupled_variables = 'u'
        constant_names = 'W1'
        constant_expressions = '1/4'
        expression = 'W1*${s}*(u^2 - 1)^2'
        derivative_order = 2
    [../]
    # mixing energy based on 
    # Flory-Huggins theory
    [./mixing_energy]
        type = DerivativeParsedMaterial
        property_name = f_mix
        coupled_variables = 'u'
        constant_names =        'R      T       V1      V2      vr      s'
        constant_expressions = '${R}    ${T}   ${V1}    ${V2}   ${vr}   ${s}'
        expression = 's*(R*T/phi)*((phi*log(phi))/V1+((1-phi)*log(1-phi))/V2+
                        (chi*phi*(1-phi))/vr)'
        derivative_order = 2
    [../]
    # elastic energy
    [../elastic_energy]
        type = DerivativeParsedMaterial
        property_name = f_el
        coupled_variables = 'u'
        constant_names =        'E      s'
        constant_expressions =  '${E}   ${s}'
        expression = 's*(E/2)*(1/(phi^(2/3))-1)'
        derivative_order = 2
    [../]
    # total free energy
    # sum of local, mixing and elastic energy
    [../free_energy]
        type = DerivativeSumMaterial
        property_name = f_tot
        coupled_variables = 'u'
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
[]
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    scheme = bdf2
  
    # Preconditioning using the additive Schwartz method and LU decomposition
    petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm                  preonly       lu           2'
  
    # # Alternative preconditioning options using Hypre (algebraic multi-grid)
    # petsc_options_iname = '-pc_type -pc_hypre_type'
    # petsc_options_value = 'hypre    boomeramg'
  
    l_tol = 1e-3
    l_abs_tol = 1e-9
    l_max_its = 30
    nl_max_its = 30
    nl_abs_tol = 1e-8
    
    [./TimeStepper]
        # Turn on time stepping
        type = IterationAdaptiveDT
        dt = 1.0e-4
        cutback_factor = 0.8
        growth_factor = 1.5
        optimal_iterations = 10
    [../]
  
    end_time = 1.0 # seconds

    # # Automatic scaling for u and w
    # automatic_scaling = true
    # scaling_group_variables = 'u w'
  
    # [./Adaptivity]
    #   coarsen_fraction = 0.1
    #   refine_fraction = 0.7
    #   max_h_level = 2
    # [../]
[]
  
[Outputs]
    # file_base = /Users/rnp/panther/test/
    [ex]
        type = Exodus
        # file_base = /Users/rnp/panther/M/100/${M_in}/nd_${N}_${a}
        time_step_interval = 1
    []
    [csv]
        type = CSV
        # file_base = /Users/rnp/panther/M/100/${M_in}/nd_${N}_${a}_e
    []
[]

# [Debug]
#     show_var_residual_norms = true
# []