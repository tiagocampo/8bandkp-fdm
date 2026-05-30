module test_parameters
  use funit
  use definitions
  use parameters
  implicit none

contains

  !@test
  subroutine test_gaas_parameters()
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

    ! Check key GaAs parameters
#line 18 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.067_dp, params(1)%meff, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 19 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(1.519_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 20 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.341_dp, params(1)%deltaSO, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 20) )
  if (anyExceptions()) return
#line 21 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 21 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(6.98_dp, params(1)%gamma1, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 21) )
  if (anyExceptions()) return
#line 22 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 22 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.06_dp, params(1)%gamma2, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 22) )
  if (anyExceptions()) return
#line 23 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 23 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.93_dp, params(1)%gamma3, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 23) )
  if (anyExceptions()) return
#line 24 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 24 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(28.8_dp, params(1)%EP, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 24) )
  if (anyExceptions()) return
#line 25 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_gaas_parameters

  !@test
  subroutine test_inas_parameters()
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

#line 35 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.417_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 36 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.39_dp, params(1)%deltaSO, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 37 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(20.0_dp, params(1)%gamma1, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 38 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(8.5_dp, params(1)%gamma2, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 38) )
  if (anyExceptions()) return
#line 39 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_inas_parameters

  !@test
  subroutine test_ec_ev_consistency()
    ! EC should equal EV + Eg for all materials
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)
#line 49 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"

    material(1) = "InAs"
    call paramDatabase(material, 1, params)
#line 53 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"

    material(1) = "AlSb"
    call paramDatabase(material, 1, params)
#line 57 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_ec_ev_consistency

  !@test
  subroutine test_p_from_ep()
    ! P = sqrt(EP * const) for all materials
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)
#line 68 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(sqrt(params(1)%EP * const), params(1)%P, tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_p_from_ep

  !@test
  subroutine test_multi_material()
    type(paramStruct) :: params(3)
    character(len=255) :: material(3)

    material(1) = "AlSb"
    material(2) = "GaSb"
    material(3) = "InAs"
    call paramDatabase(material, 3, params)

    ! AlSb Eg = 2.386
#line 82 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.386_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! GaSb Eg = 0.812
#line 84 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.812_dp, params(2)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InAs Eg = 0.417
#line 86 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.417_dp, params(3)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_multi_material

  !@test
  subroutine test_dielectric_constants()
    ! Verify eps0 values against Vurgaftman 2001 Table III
    type(paramStruct) :: params(4)
    character(len=255) :: material(4)

    material(1) = "GaAs"
    material(2) = "InAs"
    material(3) = "AlSb"
    material(4) = "InSb"
    call paramDatabase(material, 4, params)

    ! GaAs eps0 = 12.90 (Vurgaftman 2001)
#line 102 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(12.90_dp, params(1)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InAs eps0 = 15.15 (Vurgaftman 2001)
#line 104 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(15.15_dp, params(2)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! AlSb eps0 = 12.04 (Vurgaftman 2001)
#line 106 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(12.04_dp, params(3)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InSb eps0 = 16.80 (Vurgaftman 2001)
#line 108 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(16.80_dp, params(4)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_dielectric_constants

  !@test
  subroutine test_strain_elastic_constants()
    ! Verify elastic constants from Vurgaftman 2001 Table XII
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

#line 120 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(1221.0_dp, params(1)%C11, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 121 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(566.0_dp, params(1)%C12, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 122 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(599.0_dp, params(1)%C44, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_elastic_constants

  !@test
  subroutine test_strain_lattice_constants()
    ! Verify lattice constants from Vurgaftman 2001 Table I
    type(paramStruct) :: params(2)
    character(len=255) :: material(2)

    material(1) = "GaAs"
    material(2) = "InSb"
    call paramDatabase(material, 2, params)

#line 135 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(5.65325_dp, params(1)%a0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 136 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(6.4794_dp, params(2)%a0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_lattice_constants

  !@test
  subroutine test_strain_deformation_potentials()
    ! Verify deformation potentials from Vurgaftman 2001 Table XIII
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

#line 148 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(-5.08_dp, params(1)%ac, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 149 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(1.00_dp, params(1)%av, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 150 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(-1.8_dp, params(1)%b_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 151 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(-3.6_dp, params(1)%d_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_deformation_potentials

  !@test
  subroutine test_strain_w_variant_same_as_base()
    ! W-variant materials should have same elastic/strain params as base
    type(paramStruct) :: params(2)
    character(len=255) :: material(2)

    material(1) = "GaAs"
    material(2) = "GaAsW"
    call paramDatabase(material, 2, params)

#line 164 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%C11, params(2)%C11, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 165 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%C12, params(2)%C12, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 166 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%C44, params(2)%C44, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 167 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%a0, params(2)%a0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 168 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%ac, params(2)%ac, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 169 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%av, params(2)%av, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 170 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%b_dp, params(2)%b_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 171 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%d_dp, params(2)%d_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_w_variant_same_as_base

  !@test
  subroutine test_strain_defaults_zero()
    ! Vacuum should have zero strain parameters (default)
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "Vacuum"
    call paramDatabase(material, 1, params)

#line 183 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%C11, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 184 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%C12, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 185 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%C44, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 186 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%a0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 186) )
  if (anyExceptions()) return
#line 187 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 187 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%ac, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 187) )
  if (anyExceptions()) return
#line 188 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 188 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%av, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 189 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%b_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 190 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.0_dp, params(1)%d_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 190) )
  if (anyExceptions()) return
#line 191 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_defaults_zero

  !@test
  subroutine test_strain_insb_corrected()
    ! Verify InSb values match Vurgaftman 2001 Table XII/XIII
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "InSb"
    call paramDatabase(material, 1, params)

#line 202 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(373.5_dp, params(1)%C12, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 203 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(311.1_dp, params(1)%C44, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 203) )
  if (anyExceptions()) return
#line 204 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 204 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(-2.0_dp, params(1)%b_dp, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 204) )
  if (anyExceptions()) return
#line 205 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_strain_insb_corrected

  !@test
  subroutine test_alloy_strain_interpolation()
    ! Verify Vegard interpolation for Al20Ga80As
    ! C11 = 1221*0.8 + 1250*0.2 = 1226.8
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "Al20Ga80As"
    call paramDatabase(material, 1, params)

#line 217 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(1226.8_dp, params(1)%C11, tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 218 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(5.65482_dp, params(1)%a0, tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_alloy_strain_interpolation

  ! ==================================================================
  ! validate_semantic tests — passing cases only.
  ! error stop cannot be caught by pFUnit 4.x, so failure cases
  ! are tested manually.
  ! ==================================================================

  !@test
  subroutine test_validate_semantic_bandStructure_ok()
    ! bandStructure has no additional semantic constraints — always passes.
    type(simulation_config) :: cfg
    ! cfg is default-initialized; no fields need setting.
    call validate_semantic(cfg, 'bandStructure')
    ! If we reach here, the test passes (no error stop).
  end subroutine test_validate_semantic_bandStructure_ok

  !@test
  subroutine test_validate_semantic_gfactor_nsteps0_ok()
    ! gfactor with nsteps=0 should pass regardless of mode.
    type(simulation_config) :: cfg
    cfg%wave_vector%nsteps = 0
    cfg%wave_vector%mode = 'sweep'
    call validate_semantic(cfg, 'gfactor')
  end subroutine test_validate_semantic_gfactor_nsteps0_ok

  !@test
  subroutine test_validate_semantic_gfactor_mode_k0_ok()
    ! gfactor with mode='k0' should pass even if nsteps /= 0.
    type(simulation_config) :: cfg
    cfg%wave_vector%nsteps = 100
    cfg%wave_vector%mode = 'k0'
    call validate_semantic(cfg, 'gfactor')
  end subroutine test_validate_semantic_gfactor_mode_k0_ok

  !@test
  subroutine test_validate_semantic_optics_enabled_ok()
    ! opticalProperties with optics%enabled=.true. passes.
    type(simulation_config) :: cfg
    cfg%optics%enabled = .true.
    call validate_semantic(cfg, 'opticalProperties')
  end subroutine test_validate_semantic_optics_enabled_ok

  !@test
  subroutine test_validate_semantic_topology_enabled_ok()
    ! topologicalAnalysis with topo%enabled=.true. and mode set passes.
    type(simulation_config) :: cfg
    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qhe'
    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_enabled_ok

  ! ------------------------------------------------------------------
  ! S1: wire gfactor with valid bandIdx in range [1, num_cb-1]
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_gfactor_bandidx_wire_ok()
    type(simulation_config) :: cfg

    cfg%confinement = 'wire'
    cfg%wave_vector%nsteps = 0
    cfg%wave_vector%mode = 'k0'
    cfg%band_idx = 1
    cfg%bands%num_cb = 4

    call validate_semantic(cfg, 'gfactor')
  end subroutine test_validate_semantic_gfactor_bandidx_wire_ok

  ! ------------------------------------------------------------------
  ! S2: QSHE Z2 with QW confinement accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_qshe_qw_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qshe'
    cfg%confinement = 'qw'

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_qshe_qw_ok

  ! ------------------------------------------------------------------
  ! S3/S4: BdG mode with enabled bdg section and wire confinement
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_bdg_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'bdg'
    cfg%bdg%enabled = .true.
    cfg%confinement = 'wire'

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_bdg_ok

  ! ------------------------------------------------------------------
  ! S5/S6/S7: spectral function with valid parameters
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_spectral_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qhe'
    cfg%confinement = 'bulk'
    cfg%topo%compute_spectral = .true.
    cfg%topo%spectral_eta = 0.001_dp
    cfg%topo%spectral_nk = 100
    cfg%topo%spectral_nE = 200

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_spectral_ok

  ! ------------------------------------------------------------------
  ! S8: sweep_model='bhz_analytic' accepted (confinement-agnostic)
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_sweep_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qhe'
    cfg%confinement = 'qw'
    cfg%topo%sweep_model = 'bhz_analytic'

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_sweep_ok

  ! ------------------------------------------------------------------
  ! S9: conductance with recognized method accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_conductance_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qhe'
    cfg%confinement = 'bulk'
    cfg%topo%compute_conductance = .true.
    cfg%topo%conductance_method = 'kubo_chern'

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_conductance_ok

  ! ------------------------------------------------------------------
  ! S10: topology mode='qhe' recognized (redundant with existing test
  ! but explicitly exercises S10 acceptance path)
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_semantic_topology_mode_ok()
    type(simulation_config) :: cfg

    cfg%topo%enabled = .true.
    cfg%topo%mode = 'qhe'
    cfg%confinement = 'bulk'

    call validate_semantic(cfg, 'topologicalAnalysis')
  end subroutine test_validate_semantic_topology_mode_ok

end module test_parameters

module Wraptest_parameters
   use FUnit
   use test_parameters
   implicit none
   private

contains


end module Wraptest_parameters

function test_parameters_suite() result(suite)
   use FUnit
   use test_parameters
   use Wraptest_parameters
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_parameters_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_gaas_parameters', &
      test_gaas_parameters))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_inas_parameters', &
      test_inas_parameters))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_ec_ev_consistency', &
      test_ec_ev_consistency))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_p_from_ep', &
      test_p_from_ep))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_multi_material', &
      test_multi_material))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_dielectric_constants', &
      test_dielectric_constants))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_elastic_constants', &
      test_strain_elastic_constants))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_lattice_constants', &
      test_strain_lattice_constants))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_deformation_potentials', &
      test_strain_deformation_potentials))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_w_variant_same_as_base', &
      test_strain_w_variant_same_as_base))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_defaults_zero', &
      test_strain_defaults_zero))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_strain_insb_corrected', &
      test_strain_insb_corrected))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_alloy_strain_interpolation', &
      test_alloy_strain_interpolation))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_bandStructure_ok', &
      test_validate_semantic_bandStructure_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_gfactor_nsteps0_ok', &
      test_validate_semantic_gfactor_nsteps0_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_gfactor_mode_k0_ok', &
      test_validate_semantic_gfactor_mode_k0_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_optics_enabled_ok', &
      test_validate_semantic_optics_enabled_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_enabled_ok', &
      test_validate_semantic_topology_enabled_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_gfactor_bandidx_wire_ok', &
      test_validate_semantic_gfactor_bandidx_wire_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_qshe_qw_ok', &
      test_validate_semantic_topology_qshe_qw_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_bdg_ok', &
      test_validate_semantic_topology_bdg_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_spectral_ok', &
      test_validate_semantic_topology_spectral_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_sweep_ok', &
      test_validate_semantic_topology_sweep_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_conductance_ok', &
      test_validate_semantic_topology_conductance_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_semantic_topology_mode_ok', &
      test_validate_semantic_topology_mode_ok))
   call suite%addTest(t)


end function test_parameters_suite

