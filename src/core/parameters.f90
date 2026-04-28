module parameters

  use definitions

  implicit none

  contains

    subroutine paramDatabase(material,nlayers,params)

      integer, intent(in) :: nlayers
      character(len = 255), intent(in) :: material(nlayers)
      type(paramStruct), intent(inout) :: params(nlayers)

      real(kind=dp) :: gamma0, EP, A, P, gamma1, gamma2, gamma3
      integer :: i

      do i = 1, nlayers, 1
        select case(material(i))

          case ("GaAs")

            params(i)%meff    = 0.067_dp
            params(i)%EP      = 28.8_dp
            params(i)%Eg      = 1.519_dp
            params(i)%deltaSO = 0.341_dp
            params(i)%gamma1  = 6.98_dp
            params(i)%gamma2  = 2.06_dp
            params(i)%gamma3  = 2.93_dp
            params(i)%EV      = -0.8
            params(i)%EC      = 0.719
            params(i)%eps0    = 12.90_dp
            ! Strain parameters (Vurgaftman 2001)
            params(i)%C11   = 1221.0_dp
            params(i)%C12   = 566.0_dp
            params(i)%C44   = 599.0_dp
            params(i)%a0    = 5.65325_dp
            params(i)%ac    = -7.17_dp
            params(i)%av    = 1.16_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.8_dp

          case ("GaAsW")

            params(i)%meff    = 0.0665_dp
            params(i)%EP      = 28.89_dp
            params(i)%Eg      = 1.519_dp
            params(i)%deltaSO = 0.341_dp
            params(i)%gamma1  = 6.85_dp
            params(i)%gamma2  = 2.10_dp
            params(i)%gamma3  = 2.90_dp
            params(i)%EV      = -0.8_dp
            params(i)%EC      = 0.719_dp
            params(i)%eps0    = 12.90_dp
            ! Strain: same as GaAs
            params(i)%C11   = 1221.0_dp
            params(i)%C12   = 566.0_dp
            params(i)%C44   = 599.0_dp
            params(i)%a0    = 5.65325_dp
            params(i)%ac    = -7.17_dp
            params(i)%av    = 1.16_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.8_dp

          case ("Al20Ga80As")

            params(i)%meff    = 0.0836_dp
            params(i)%EP      = 27.26_dp
            params(i)%Eg      = 1.835_dp
            params(i)%deltaSO = 0.3288_dp
            params(i)%gamma1  = 6.336_dp
            params(i)%gamma2  = 1.812_dp
            params(i)%gamma3  = 2.628_dp
            params(i)%EV      = -0.906
            params(i)%EC      = 0.929
            params(i)%eps0    = 12.332_dp
            ! Strain: Vegard interpolation GaAs/AlAs, x=0.20
            params(i)%C11   = 1221.0_dp*0.8_dp + 1250.0_dp*0.2_dp
            params(i)%C12   = 566.0_dp*0.8_dp + 534.0_dp*0.2_dp
            params(i)%C44   = 599.0_dp*0.8_dp + 542.0_dp*0.2_dp
            params(i)%a0    = 5.65325_dp*0.8_dp + 5.6611_dp*0.2_dp
            params(i)%ac    = -7.17_dp*0.8_dp + (-5.64_dp)*0.2_dp
            params(i)%av    = 1.16_dp*0.8_dp + 2.47_dp*0.2_dp
            params(i)%b_dp  = -2.0_dp*0.8_dp + (-2.3_dp)*0.2_dp
            params(i)%d_dp  = -4.8_dp*0.8_dp + (-3.4_dp)*0.2_dp

          case ("Al15Ga85As")

            params(i)%meff    = 0.0795_dp
            params(i)%EP      = 27.645_dp
            params(i)%Eg      = 1.756_dp
            params(i)%deltaSO = 0.3319_dp
            params(i)%gamma1  = 6.497_dp
            params(i)%gamma2  = 1.874_dp
            params(i)%gamma3  = 2.7035_dp
            params(i)%EV      = -0.8795
            params(i)%EC      = 0.8765
            params(i)%eps0    = 12.474_dp
            ! Strain: Vegard interpolation GaAs/AlAs, x=0.15
            params(i)%C11   = 1221.0_dp*0.85_dp + 1250.0_dp*0.15_dp
            params(i)%C12   = 566.0_dp*0.85_dp + 534.0_dp*0.15_dp
            params(i)%C44   = 599.0_dp*0.85_dp + 542.0_dp*0.15_dp
            params(i)%a0    = 5.65325_dp*0.85_dp + 5.6611_dp*0.15_dp
            params(i)%ac    = -7.17_dp*0.85_dp + (-5.64_dp)*0.15_dp
            params(i)%av    = 1.16_dp*0.85_dp + 2.47_dp*0.15_dp
            params(i)%b_dp  = -2.0_dp*0.85_dp + (-2.3_dp)*0.15_dp
            params(i)%d_dp  = -4.8_dp*0.85_dp + (-3.4_dp)*0.15_dp

          case ("Al30Ga70As")

            params(i)%meff    = 0.093_dp
            params(i)%EP      = 26.32_dp
            params(i)%Eg      = 1.977_dp
            params(i)%deltaSO = 0.353_dp
            params(i)%gamma1  = 6.107_dp
            params(i)%gamma2  = 1.773_dp
            params(i)%gamma3  = 2.543_dp
            params(i)%EV      = -0.959_dp
            params(i)%EC      = 1.018_dp
            params(i)%eps0    = 12.048_dp
            ! Strain: Vegard interpolation GaAs/AlAs, x=0.30
            params(i)%C11   = 1221.0_dp*0.70_dp + 1250.0_dp*0.30_dp
            params(i)%C12   = 566.0_dp*0.70_dp + 534.0_dp*0.30_dp
            params(i)%C44   = 599.0_dp*0.70_dp + 542.0_dp*0.30_dp
            params(i)%a0    = 5.65325_dp*0.70_dp + 5.6611_dp*0.30_dp
            params(i)%ac    = -7.17_dp*0.70_dp + (-5.64_dp)*0.30_dp
            params(i)%av    = 1.16_dp*0.70_dp + 2.47_dp*0.30_dp
            params(i)%b_dp  = -2.0_dp*0.70_dp + (-2.3_dp)*0.30_dp
            params(i)%d_dp  = -4.8_dp*0.70_dp + (-3.4_dp)*0.30_dp

          case ("Ga47In53AsW")

            params(i)%meff    = 0.038_dp
            params(i)%EP      = 25.26_dp
            params(i)%Eg      = 0.8166_dp
            params(i)%deltaSO = 0.362_dp
            params(i)%gamma1  = 11.97_dp
            params(i)%gamma2  = 4.36_dp
            params(i)%gamma3  = 5.15_dp
            params(i)%EV      = -0.689_dp
            params(i)%EC      = 0.128_dp
            params(i)%eps0    = 14.0925_dp
            ! Strain: Vegard interpolation GaAsW/InAsW, x_Ga=0.47, x_In=0.53
            params(i)%C11   = 1221.0_dp*0.47_dp + 832.9_dp*0.53_dp
            params(i)%C12   = 566.0_dp*0.47_dp + 452.6_dp*0.53_dp
            params(i)%C44   = 599.0_dp*0.47_dp + 395.9_dp*0.53_dp
            params(i)%a0    = 5.65325_dp*0.47_dp + 6.0583_dp*0.53_dp
            params(i)%ac    = -7.17_dp*0.47_dp + (-5.08_dp)*0.53_dp
            params(i)%av    = 1.16_dp*0.47_dp + 1.00_dp*0.53_dp
            params(i)%b_dp  = -2.0_dp*0.47_dp + (-1.8_dp)*0.53_dp
            params(i)%d_dp  = -4.8_dp*0.47_dp + (-3.6_dp)*0.53_dp

          case ("Al47In53AsW")

            params(i)%meff    = 0.0779_dp
            params(i)%EP      = 21.7_dp
            params(i)%Eg      = 1.693_dp
            params(i)%deltaSO = 0.342_dp
            params(i)%gamma1  = 6.17_dp
            params(i)%gamma2  = 1.62_dp
            params(i)%gamma3  = 2.31_dp
            params(i)%EV      = -0.938_dp
            params(i)%EC      = 0.755_dp
            params(i)%eps0    = 12.7577_dp
            ! Strain: Vegard interpolation AlAsW/InAsW, x_Al=0.47, x_In=0.53
            params(i)%C11   = 1250.0_dp*0.47_dp + 832.9_dp*0.53_dp
            params(i)%C12   = 534.0_dp*0.47_dp + 452.6_dp*0.53_dp
            params(i)%C44   = 542.0_dp*0.47_dp + 395.9_dp*0.53_dp
            params(i)%a0    = 5.6611_dp*0.47_dp + 6.0583_dp*0.53_dp
            params(i)%ac    = -5.64_dp*0.47_dp + (-5.08_dp)*0.53_dp
            params(i)%av    = 2.47_dp*0.47_dp + 1.00_dp*0.53_dp
            params(i)%b_dp  = -2.3_dp*0.47_dp + (-1.8_dp)*0.53_dp
            params(i)%d_dp  = -3.4_dp*0.47_dp + (-3.6_dp)*0.53_dp

          case ("InAs")

            params(i)%meff    = 0.026_dp
            params(i)%EP      = 21.5_dp
            params(i)%Eg      = 0.417_dp
            params(i)%deltaSO = 0.39_dp
            params(i)%gamma1  = 20.0_dp
            params(i)%gamma2  = 8.5_dp
            params(i)%gamma3  = 9.2_dp
            params(i)%EV      = -0.59
            params(i)%EC      = -0.173
            params(i)%eps0    = 15.15_dp
            ! Strain parameters (Vurgaftman 2001)
            params(i)%C11   = 832.9_dp
            params(i)%C12   = 452.6_dp
            params(i)%C44   = 395.9_dp
            params(i)%a0    = 6.0583_dp
            params(i)%ac    = -5.08_dp
            params(i)%av    = 1.00_dp
            params(i)%b_dp  = -1.8_dp
            params(i)%d_dp  = -3.6_dp

          case ("InAsW")

            params(i)%meff    = 0.0229_dp
            params(i)%EP      = 22.2_dp
            params(i)%Eg      = 0.418_dp
            params(i)%deltaSO = 0.38_dp
            params(i)%gamma1  = 20.4_dp
            params(i)%gamma2  = 8.3_dp
            params(i)%gamma3  = 9.1_dp
            params(i)%EV      = -0.59_dp
            params(i)%EC      = -0.172_dp
            params(i)%eps0    = 15.15_dp
            ! Strain: same as InAs
            params(i)%C11   = 832.9_dp
            params(i)%C12   = 452.6_dp
            params(i)%C44   = 395.9_dp
            params(i)%a0    = 6.0583_dp
            params(i)%ac    = -5.08_dp
            params(i)%av    = 1.00_dp
            params(i)%b_dp  = -1.8_dp
            params(i)%d_dp  = -3.6_dp

          case ("AlAs")

            params(i)%meff    = 0.15_dp
            params(i)%EP      = 21.1_dp
            params(i)%Eg      = 3.099_dp
            params(i)%deltaSO = 0.28_dp
            params(i)%gamma1  = 3.76_dp
            params(i)%gamma2  = 0.82_dp
            params(i)%gamma3  = 1.42_dp
            params(i)%EV      = -1.33
            params(i)%EC      = 1.769
            params(i)%eps0    = 10.06_dp
            ! Strain parameters (Vurgaftman 2001)
            params(i)%C11   = 1250.0_dp
            params(i)%C12   = 534.0_dp
            params(i)%C44   = 542.0_dp
            params(i)%a0    = 5.6611_dp
            params(i)%ac    = -5.64_dp
            params(i)%av    = 2.47_dp
            params(i)%b_dp  = -2.3_dp
            params(i)%d_dp  = -3.4_dp

          case ("AlAsW")

            params(i)%meff    = 0.15_dp
            params(i)%EP      = 21.12_dp
            params(i)%Eg      = 3.13_dp
            params(i)%deltaSO = 0.3_dp
            params(i)%gamma1  = 3.25_dp
            params(i)%gamma2  = 0.65_dp
            params(i)%gamma3  = 1.21_dp
            params(i)%EV      = -1.33_dp
            params(i)%EC      = 1.800_dp
            params(i)%eps0    = 10.06_dp
            ! Strain: same as AlAs
            params(i)%C11   = 1250.0_dp
            params(i)%C12   = 534.0_dp
            params(i)%C44   = 542.0_dp
            params(i)%a0    = 5.6611_dp
            params(i)%ac    = -5.64_dp
            params(i)%av    = 2.47_dp
            params(i)%b_dp  = -2.3_dp
            params(i)%d_dp  = -3.4_dp

          case ("GaP")

            params(i)%meff    = 0.13_dp
            params(i)%EP      = 31.4_dp
            params(i)%Eg      = 2.35_dp
            params(i)%deltaSO = 0.08_dp
            params(i)%gamma1  = 4.05_dp
            params(i)%gamma2  = 0.49_dp
            params(i)%gamma3  = 2.93_dp
            params(i)%EV      = -1.27_dp
            params(i)%EC      = 1.08_dp
            params(i)%eps0    = 11.11_dp
            ! Strain parameters (Vurgaftman 2001 Table XII/XIII)
            ! av uses positive sign convention: P_eps = -av*Tr(eps)
            params(i)%C11   = 1405.0_dp
            params(i)%C12   = 620.3_dp
            params(i)%C44   = 703.3_dp
            params(i)%a0    = 5.4505_dp
            params(i)%ac    = -8.2_dp
            params(i)%av    = 1.7_dp
            params(i)%b_dp  = -1.6_dp
            params(i)%d_dp  = -4.6_dp

          case ("AlP")

            params(i)%meff    = 0.22_dp
            params(i)%EP      = 17.7_dp
            params(i)%Eg      = 2.52_dp
            params(i)%deltaSO = 0.07_dp
            params(i)%gamma1  = 3.35_dp
            params(i)%gamma2  = 0.71_dp
            params(i)%gamma3  = 1.23_dp
            params(i)%EV      = -1.74_dp
            params(i)%EC      = 0.78_dp
            params(i)%eps0    = 9.02_dp
            ! Strain parameters (Vurgaftman 2001 Table XII/XIII)
            ! av uses positive sign convention: P_eps = -av*Tr(eps)
            ! Note: many AlP values are estimated (Vurgaftman 2001)
            params(i)%C11   = 1330.0_dp
            params(i)%C12   = 630.3_dp
            params(i)%C44   = 615.2_dp
            params(i)%a0    = 5.4672_dp
            params(i)%ac    = -5.7_dp
            params(i)%av    = 3.0_dp
            params(i)%b_dp  = -1.5_dp
            params(i)%d_dp  = -4.6_dp

          case ("InP")

            params(i)%meff    = 0.0795_dp
            params(i)%EP      = 20.7_dp
            params(i)%Eg      = 1.4236_dp
            params(i)%deltaSO = 0.108_dp
            params(i)%gamma1  = 5.08_dp
            params(i)%gamma2  = 1.60_dp
            params(i)%gamma3  = 2.10_dp
            params(i)%EV      = -0.70_dp
            params(i)%EC      = 0.7236_dp
            params(i)%eps0    = 12.61_dp
            ! Strain parameters (Vurgaftman 2001 Table XII/XIII)
            ! av uses positive sign convention: P_eps = -av*Tr(eps)
            params(i)%C11   = 1011.0_dp
            params(i)%C12   = 561.0_dp
            params(i)%C44   = 456.0_dp
            params(i)%a0    = 5.8687_dp
            params(i)%ac    = -6.0_dp
            params(i)%av    = 0.6_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -5.0_dp

          case ("InPW")

            params(i)%meff    = 0.0803_dp
            params(i)%EP      = 20.56_dp
            params(i)%Eg      = 1.423_dp
            params(i)%deltaSO = 0.11_dp
            params(i)%gamma1  = 4.95_dp
            params(i)%gamma2  = 1.65_dp
            params(i)%gamma3  = 2.35_dp
            params(i)%EV      = -0.70_dp
            params(i)%EC      = 0.723_dp
            params(i)%eps0    = 12.61_dp
            ! Strain: same as InP
            params(i)%C11   = 1011.0_dp
            params(i)%C12   = 561.0_dp
            params(i)%C44   = 456.0_dp
            params(i)%a0    = 5.8687_dp
            params(i)%ac    = -6.35_dp
            params(i)%av    = 1.70_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -5.6_dp

          case ("GaSb")

            params(i)%meff    = 0.039_dp
            params(i)%EP      = 27.0_dp
            params(i)%Eg      = 0.812_dp
            params(i)%deltaSO = 0.76_dp
            params(i)%gamma1  = 13.4_dp
            params(i)%gamma2  = 4.7_dp
            params(i)%gamma3  = 6.0_dp
            params(i)%eps0    = 15.70_dp
            ! Strain parameters (Vurgaftman 2001)
            params(i)%C11   = 884.2_dp
            params(i)%C12   = 402.4_dp
            params(i)%C44   = 432.3_dp
            params(i)%a0    = 6.0959_dp
            params(i)%ac    = -7.5_dp
            params(i)%av    = 0.8_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.7_dp

          case ("GaSbW")

            params(i)%meff    = 0.041_dp
            params(i)%EP      = 22.37_dp
            params(i)%Eg      = 0.812_dp
            params(i)%deltaSO = 0.76_dp
            params(i)%gamma1  = 13.4_dp
            params(i)%gamma2  = 4.7_dp
            params(i)%gamma3  = 5.7_dp
            params(i)%EV      = -0.03_dp
            params(i)%EC      = 0.782_dp
            params(i)%eps0    = 15.70_dp
            ! Strain: same as GaSb
            params(i)%C11   = 884.2_dp
            params(i)%C12   = 402.4_dp
            params(i)%C44   = 432.3_dp
            params(i)%a0    = 6.0959_dp
            params(i)%ac    = -7.5_dp
            params(i)%av    = 0.8_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.7_dp

          case ("AlSb")

            params(i)%meff    = 0.14_dp
            params(i)%EP      = 18.7_dp
            params(i)%Eg      = 2.386_dp
            params(i)%deltaSO = 0.676_dp
            params(i)%gamma1  = 5.18_dp
            params(i)%gamma2  = 1.19_dp
            params(i)%gamma3  = 1.97_dp
            params(i)%EV      = -0.41_dp
            params(i)%EC      = 1.976_dp
            params(i)%eps0    = 12.04_dp
            ! Strain parameters (Vurgaftman 2001)
            params(i)%C11   = 876.5_dp
            params(i)%C12   = 434.1_dp
            params(i)%C44   = 407.6_dp
            params(i)%a0    = 6.1355_dp
            params(i)%ac    = -4.5_dp
            params(i)%av    = 1.4_dp
            params(i)%b_dp  = -1.35_dp
            params(i)%d_dp  = -4.3_dp

          case ("AlSbW")

            params(i)%meff    = 0.12_dp
            params(i)%EP      = 18.8_dp
            params(i)%Eg      = 2.384_dp
            params(i)%deltaSO = 0.673_dp
            params(i)%gamma1  = 4.15_dp
            params(i)%gamma2  = 1.01_dp
            params(i)%gamma3  = 1.71_dp
            params(i)%EV      = -0.41_dp
            params(i)%EC      = 1.974_dp
            params(i)%eps0    = 12.04_dp
            ! Strain: same as AlSb
            params(i)%C11   = 876.5_dp
            params(i)%C12   = 434.1_dp
            params(i)%C44   = 407.6_dp
            params(i)%a0    = 6.1355_dp
            params(i)%ac    = -4.5_dp
            params(i)%av    = 1.4_dp
            params(i)%b_dp  = -1.35_dp
            params(i)%d_dp  = -4.3_dp

          case ("InSb")

            params(i)%meff    = 0.0135_dp
            params(i)%EP      = 23.3_dp
            params(i)%Eg      = 0.235_dp
            params(i)%deltaSO = 0.81_dp
            params(i)%gamma1  = 34.8_dp
            params(i)%gamma2  = 15.5_dp
            params(i)%gamma3  = 16.5_dp
            params(i)%EV      = 0
            params(i)%EC      = 0.235
            params(i)%eps0    = 16.80_dp
            ! Strain parameters (Vurgaftman 2001 Table XII/XIII)
            ! av uses positive sign convention: P_eps = -av*Tr(eps)
            params(i)%C11   = 684.7_dp
            params(i)%C12   = 373.5_dp
            params(i)%C44   = 311.1_dp
            params(i)%a0    = 6.4794_dp
            params(i)%ac    = -6.94_dp
            params(i)%av    = 0.36_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.7_dp

          case ("InSbW")

            params(i)%meff    = 0.0139_dp
            params(i)%EP      = 24.4_dp
            params(i)%Eg      = 0.237_dp
            params(i)%deltaSO = 0.81_dp
            params(i)%gamma1  = 37.1_dp
            params(i)%gamma2  = 16.5_dp
            params(i)%gamma3  = 17.7_dp
            params(i)%EV      = 0.0_dp
            params(i)%EC      = 0.237_dp
            params(i)%eps0    = 16.80_dp
            ! Strain: same as InSb (Vurgaftman 2001)
            params(i)%C11   = 684.7_dp
            params(i)%C12   = 373.5_dp
            params(i)%C44   = 311.1_dp
            params(i)%a0    = 6.4794_dp
            params(i)%ac    = -6.94_dp
            params(i)%av    = 0.36_dp
            params(i)%b_dp  = -2.0_dp
            params(i)%d_dp  = -4.7_dp

          case ("Vacuum")

            params(i)%meff    = 1.0_dp
            params(i)%EP      = 0.0_dp
            params(i)%Eg      = 3.0_dp
            params(i)%deltaSO = 0.7021_dp
            params(i)%gamma1  = -1.0_dp
            params(i)%gamma2  = 0.0_dp
            params(i)%gamma3  = 0.0_dp
            params(i)%EV      = -0.1326_dp
            params(i)%EC      = 2.8674_dp
            params(i)%eps0    = 1.0_dp

          case ("Al")

            params(i)%meff    = 1.0_dp
            params(i)%EP      = 21.999_dp
            params(i)%Eg      = 0.0_dp
            params(i)%EV      = -11.6_dp
            params(i)%EC      = -11.6_dp
            params(i)%deltaSO = 0.0_dp
            params(i)%gamma1  = -1.0_dp
            params(i)%gamma2  = 0.0_dp
            params(i)%gamma3  = 0.0_dp
            params(i)%eps0    = 1.0_dp

          case ("InAs10Sb90")

            params(i)%meff    = 0.0116_dp
            params(i)%EP      = 23.12_dp
            params(i)%Eg      = 0.1749_dp
            params(i)%EV      = -0.504_dp
            params(i)%EC      = -0.1835_dp
            params(i)%deltaSO = 0.66_dp
            params(i)%gamma1  = 33.32_dp
            params(i)%gamma2  = 14.8_dp
            params(i)%gamma3  = 15.77_dp
            params(i)%eps0    = 16.635_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.10
            params(i)%C11   = 832.9_dp*0.10_dp + 684.7_dp*0.90_dp
            params(i)%C12   = 452.6_dp*0.10_dp + 373.5_dp*0.90_dp
            params(i)%C44   = 395.9_dp*0.10_dp + 311.1_dp*0.90_dp
            params(i)%a0    = 6.0583_dp*0.10_dp + 6.4794_dp*0.90_dp
            params(i)%ac    = -5.08_dp*0.10_dp + (-6.94_dp)*0.90_dp
            params(i)%av    = 1.00_dp*0.10_dp + 0.36_dp*0.90_dp
            params(i)%b_dp  = -1.8_dp*0.10_dp + (-2.0_dp)*0.90_dp
            params(i)%d_dp  = -3.6_dp*0.10_dp + (-4.7_dp)*0.90_dp

          case ("InAs20Sb80")

            params(i)%meff    = 0.0104_dp
            params(i)%EP      = 22.94_dp
            params(i)%Eg      = 0.1322_dp
            params(i)%EV      = -0.424_dp
            params(i)%EC      = -0.1826_dp
            params(i)%deltaSO = 0.534_dp
            params(i)%gamma1  = 31.84_dp
            params(i)%gamma2  = 14.1_dp
            params(i)%gamma3  = 15.04_dp
            params(i)%eps0    = 16.47_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.20
            params(i)%C11   = 832.9_dp*0.20_dp + 684.7_dp*0.80_dp
            params(i)%C12   = 452.6_dp*0.20_dp + 373.5_dp*0.80_dp
            params(i)%C44   = 395.9_dp*0.20_dp + 311.1_dp*0.80_dp
            params(i)%a0    = 6.0583_dp*0.20_dp + 6.4794_dp*0.80_dp
            params(i)%ac    = -5.08_dp*0.20_dp + (-6.94_dp)*0.80_dp
            params(i)%av    = 1.00_dp*0.20_dp + 0.36_dp*0.80_dp
            params(i)%b_dp  = -1.8_dp*0.20_dp + (-2.0_dp)*0.80_dp
            params(i)%d_dp  = -3.6_dp*0.20_dp + (-4.7_dp)*0.80_dp

          case ("InAs30Sb70")

            params(i)%gamma1  = 30.36_dp
            params(i)%gamma2  = 13.4_dp
            params(i)%gamma3  = 14.31_dp
            params(i)%Eg      = 0.1069_dp
            params(i)%meff    = 0.0099_dp
            params(i)%deltaSO = 0.432_dp
            params(i)%EP      = 22.76_dp
            params(i)%EV      = -0.35_dp
            params(i)%EC      = -0.1703_dp
            params(i)%eps0    = 16.305_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.30
            params(i)%C11   = 832.9_dp*0.30_dp + 684.7_dp*0.70_dp
            params(i)%C12   = 452.6_dp*0.30_dp + 373.5_dp*0.70_dp
            params(i)%C44   = 395.9_dp*0.30_dp + 311.1_dp*0.70_dp
            params(i)%a0    = 6.0583_dp*0.30_dp + 6.4794_dp*0.70_dp
            params(i)%ac    = -5.08_dp*0.30_dp + (-6.94_dp)*0.70_dp
            params(i)%av    = 1.00_dp*0.30_dp + 0.36_dp*0.70_dp
            params(i)%b_dp  = -1.8_dp*0.30_dp + (-2.0_dp)*0.70_dp
            params(i)%d_dp  = -3.6_dp*0.30_dp + (-4.7_dp)*0.70_dp

          case ("InAs40Sb60")

            params(i)%gamma1  = 28.88_dp
            params(i)%gamma2  = 12.7_dp
            params(i)%gamma3  = 13.58_dp
            params(i)%Eg      = 0.099_dp
            params(i)%meff    = 0.0101_dp
            params(i)%deltaSO = 0.354_dp
            params(i)%EP      = 22.58_dp
            params(i)%EV      = -0.282_dp
            params(i)%EC      = -0.1466_dp
            params(i)%eps0    = 16.14_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.40
            params(i)%C11   = 832.9_dp*0.40_dp + 684.7_dp*0.60_dp
            params(i)%C12   = 452.6_dp*0.40_dp + 373.5_dp*0.60_dp
            params(i)%C44   = 395.9_dp*0.40_dp + 311.1_dp*0.60_dp
            params(i)%a0    = 6.0583_dp*0.40_dp + 6.4794_dp*0.60_dp
            params(i)%ac    = -5.08_dp*0.40_dp + (-6.94_dp)*0.60_dp
            params(i)%av    = 1.00_dp*0.40_dp + 0.36_dp*0.60_dp
            params(i)%b_dp  = -1.8_dp*0.40_dp + (-2.0_dp)*0.60_dp
            params(i)%d_dp  = -3.6_dp*0.40_dp + (-4.7_dp)*0.60_dp

          case ("InAs50Sb50")

            params(i)%meff    = 0.011_dp
            params(i)%EP      = 22.4_dp
            params(i)%Eg      = 0.1085_dp
            params(i)%EV      = -0.22_dp
            params(i)%EC      = -0.1115_dp
            params(i)%deltaSO = 0.3_dp
            params(i)%gamma1  = 27.4_dp
            params(i)%gamma2  = 12.0_dp
            params(i)%gamma3  = 12.85_dp
            params(i)%eps0    = 15.975_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.50
            params(i)%C11   = 832.9_dp*0.50_dp + 684.7_dp*0.50_dp
            params(i)%C12   = 452.6_dp*0.50_dp + 373.5_dp*0.50_dp
            params(i)%C44   = 395.9_dp*0.50_dp + 311.1_dp*0.50_dp
            params(i)%a0    = 6.0583_dp*0.50_dp + 6.4794_dp*0.50_dp
            params(i)%ac    = -5.08_dp*0.50_dp + (-6.94_dp)*0.50_dp
            params(i)%av    = 1.00_dp*0.50_dp + 0.36_dp*0.50_dp
            params(i)%b_dp  = -1.8_dp*0.50_dp + (-2.0_dp)*0.50_dp
            params(i)%d_dp  = -3.6_dp*0.50_dp + (-4.7_dp)*0.50_dp

          case ("InAs60Sb40")

            params(i)%gamma1  = 25.92_dp
            params(i)%gamma2  = 11.30_dp
            params(i)%gamma3  = 12.12_dp
            params(i)%Eg      = 0.1354_dp
            params(i)%meff    = 0.0126_dp
            params(i)%deltaSO = 0.27_dp
            params(i)%EP      = 22.22_dp
            params(i)%EV      = -0.164_dp
            params(i)%EC      = -0.065_dp
            params(i)%eps0    = 15.81_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.60
            params(i)%C11   = 832.9_dp*0.60_dp + 684.7_dp*0.40_dp
            params(i)%C12   = 452.6_dp*0.60_dp + 373.5_dp*0.40_dp
            params(i)%C44   = 395.9_dp*0.60_dp + 311.1_dp*0.40_dp
            params(i)%a0    = 6.0583_dp*0.60_dp + 6.4794_dp*0.40_dp
            params(i)%ac    = -5.08_dp*0.60_dp + (-6.94_dp)*0.40_dp
            params(i)%av    = 1.00_dp*0.60_dp + 0.36_dp*0.40_dp
            params(i)%b_dp  = -1.8_dp*0.60_dp + (-2.0_dp)*0.40_dp
            params(i)%d_dp  = -3.6_dp*0.60_dp + (-4.7_dp)*0.40_dp

          case ("InAs70Sb30")

            params(i)%gamma1  = 24.44_dp
            params(i)%gamma2  = 10.60_dp
            params(i)%gamma3  = 11.39_dp
            params(i)%Eg      = 0.1797_dp
            params(i)%meff    = 0.0149_dp
            params(i)%deltaSO = 0.264_dp
            params(i)%EP      = 22.04_dp
            params(i)%EV      = -0.114_dp
            params(i)%EC      = -0.0071_dp
            params(i)%eps0    = 15.645_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.70
            params(i)%C11   = 832.9_dp*0.70_dp + 684.7_dp*0.30_dp
            params(i)%C12   = 452.6_dp*0.70_dp + 373.5_dp*0.30_dp
            params(i)%C44   = 395.9_dp*0.70_dp + 311.1_dp*0.30_dp
            params(i)%a0    = 6.0583_dp*0.70_dp + 6.4794_dp*0.30_dp
            params(i)%ac    = -5.08_dp*0.70_dp + (-6.94_dp)*0.30_dp
            params(i)%av    = 1.00_dp*0.70_dp + 0.36_dp*0.30_dp
            params(i)%b_dp  = -1.8_dp*0.70_dp + (-2.0_dp)*0.30_dp
            params(i)%d_dp  = -3.6_dp*0.70_dp + (-4.7_dp)*0.30_dp

          case ("InAs80Sb20")

            params(i)%gamma1  = 22.96_dp
            params(i)%gamma2  = 9.90_dp
            params(i)%gamma3  = 10.66_dp
            params(i)%Eg      = 0.2414_dp
            params(i)%meff    = 0.0179_dp
            params(i)%deltaSO = 0.282_dp
            params(i)%EP      = 21.86_dp
            params(i)%EV      = -0.07_dp
            params(i)%EC      = 0.0622_dp
            params(i)%eps0    = 15.48_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.80
            params(i)%C11   = 832.9_dp*0.80_dp + 684.7_dp*0.20_dp
            params(i)%C12   = 452.6_dp*0.80_dp + 373.5_dp*0.20_dp
            params(i)%C44   = 395.9_dp*0.80_dp + 311.1_dp*0.20_dp
            params(i)%a0    = 6.0583_dp*0.80_dp + 6.4794_dp*0.20_dp
            params(i)%ac    = -5.08_dp*0.80_dp + (-6.94_dp)*0.20_dp
            params(i)%av    = 1.00_dp*0.80_dp + 0.36_dp*0.20_dp
            params(i)%b_dp  = -1.8_dp*0.80_dp + (-2.0_dp)*0.20_dp
            params(i)%d_dp  = -3.6_dp*0.80_dp + (-4.7_dp)*0.20_dp

          case ("InAs90Sb10")

            params(i)%gamma1  = 21.48_dp
            params(i)%gamma2  = 9.20_dp
            params(i)%gamma3  = 9.93_dp
            params(i)%Eg      = 0.3205_dp
            params(i)%meff    = 0.0216_dp
            params(i)%deltaSO = 0.324_dp
            params(i)%EP      = 21.68_dp
            params(i)%EV      = -0.032_dp
            params(i)%EC      = 0.1429_dp
            params(i)%eps0    = 15.315_dp
            ! Strain: Vegard interpolation InAs/InSb, x_As=0.90
            params(i)%C11   = 832.9_dp*0.90_dp + 684.7_dp*0.10_dp
            params(i)%C12   = 452.6_dp*0.90_dp + 373.5_dp*0.10_dp
            params(i)%C44   = 395.9_dp*0.90_dp + 311.1_dp*0.10_dp
            params(i)%a0    = 6.0583_dp*0.90_dp + 6.4794_dp*0.10_dp
            params(i)%ac    = -5.08_dp*0.90_dp + (-6.94_dp)*0.10_dp
            params(i)%av    = 1.00_dp*0.90_dp + 0.36_dp*0.10_dp
            params(i)%b_dp  = -1.8_dp*0.90_dp + (-2.0_dp)*0.10_dp
            params(i)%d_dp  = -3.6_dp*0.90_dp + (-4.7_dp)*0.10_dp

          case ("Al63In37Sb")

            params(i)%meff    = 0.0137_dp
            params(i)%EP      = 21.598_dp
            params(i)%Eg      = 0.9306_dp
            params(i)%EV      = -0.1326_dp
            params(i)%EC      = 0.7981_dp
            params(i)%deltaSO = 0.7021_dp
            params(i)%gamma1  = 21.7616_dp
            params(i)%gamma2  = 9.0398_dp
            params(i)%gamma3  = 10.0529_dp
            params(i)%eps0    = 13.8012_dp
            ! Strain: Vegard interpolation AlSb/InSb, x_Al=0.63
            params(i)%C11   = 876.5_dp*0.63_dp + 684.7_dp*0.37_dp
            params(i)%C12   = 434.1_dp*0.63_dp + 373.5_dp*0.37_dp
            params(i)%C44   = 407.6_dp*0.63_dp + 311.1_dp*0.37_dp
            params(i)%a0    = 6.1355_dp*0.63_dp + 6.4794_dp*0.37_dp
            params(i)%ac    = -4.5_dp*0.63_dp + (-6.94_dp)*0.37_dp
            params(i)%av    = 1.4_dp*0.63_dp + 0.36_dp*0.37_dp
            params(i)%b_dp  = -1.35_dp*0.63_dp + (-2.0_dp)*0.37_dp
            params(i)%d_dp  = -4.3_dp*0.63_dp + (-4.7_dp)*0.37_dp

          case default
            stop "material not defined"

        end select

        ! Compute derived parameters from raw material data
        params(i)%P = sqrt(params(i)%EP*const)
        params(i)%A = 1.0_dp/params(i)%meff

        if (renormalization .or. params(i)%A < 0) then

          gamma0 = params(i)%A
          EP = (gamma0-1.0_dp)*params(i)%Eg*(params(i)%Eg+params(i)%deltaSO)/(params(i)%Eg+(2.0_dp/3.0_dp)*params(i)%deltaSO)
          P = sqrt(EP*const)
          A = - ((params(i)%Eg+(2.0_dp/3.0_dp)*params(i)%deltaSO)/(params(i)%Eg+params(i)%deltaSO))*(EP/params(i)%Eg) + gamma0

          gamma1 = const*(params(i)%gamma1 - EP/(3.0_dp*params(i)%Eg))
          gamma2 = const*(params(i)%gamma2 - EP/(6.0_dp*params(i)%Eg))
          gamma3 = const*(params(i)%gamma3 - EP/(6.0_dp*params(i)%Eg))

          print *, "Material: ", trim(material(i))
          print *, "Old parameters / Renormalized parameters"
          print *, "EP :", params(i)%EP, EP
          print *, "P  :", params(i)%P, P
          print *, "A  :", params(i)%A, A
          print *, "gamma1 :", params(i)%gamma1, gamma1
          print *, "gamma2 :", params(i)%gamma2, gamma2
          print *, "gamma3 :", params(i)%gamma3, gamma3
          print *, " "

          params(i)%EP = EP
          params(i)%P = P
          params(i)%A = const*A
          params(i)%gamma1 = gamma1
          params(i)%gamma2 = gamma2
          params(i)%gamma3 = gamma3

        else

          print *, "Material: ", trim(material(i))
          print *, "Parameters"
          print *, "EP :", params(i)%EP
          print *, "P  :", params(i)%P
          print *, "A  :", params(i)%A
          print *, "gamma1 :", params(i)%gamma1
          print *, "gamma2 :", params(i)%gamma2
          print *, "gamma3 :", params(i)%gamma3
          print *, " "

        end if
      end do



    end subroutine paramDatabase

end module parameters
