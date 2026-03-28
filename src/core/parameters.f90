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
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.519_dp
            params(i)%deltaSO = 0.341_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.98_dp
            params(i)%gamma2  = 2.06_dp
            params(i)%gamma3  = 2.93_dp
            params(i)%EV      = -0.8
            params(i)%EC      = 0.719

          case ("GaAsW")

            params(i)%meff    = 0.0665_dp
            params(i)%EP      = 28.89_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.519_dp
            params(i)%deltaSO = 0.341_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.85_dp
            params(i)%gamma2  = 2.10_dp
            params(i)%gamma3  = 2.90_dp

          case ("Al20Ga80As")

            params(i)%meff    = 0.0836_dp
            params(i)%EP      = 27.26_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.835_dp
            params(i)%deltaSO = 0.3288_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.336_dp
            params(i)%gamma2  = 1.812_dp
            params(i)%gamma3  = 2.628_dp
            params(i)%EV      = -0.906
            params(i)%EC      = 0.929


          case ("Al15Ga85As")

            params(i)%meff    = 0.0795_dp
            params(i)%EP      = 27.645_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.756_dp
            params(i)%deltaSO = 0.3319_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.497_dp
            params(i)%gamma2  = 1.874_dp
            params(i)%gamma3  = 2.7035_dp
            params(i)%EV      = -0.8795
            params(i)%EC      = 0.8765

          case ("Al30Ga70As")

            params(i)%meff    = 0.093_dp
            params(i)%EP      = 26.32_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.977_dp
            params(i)%deltaSO = 0.353_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.107_dp
            params(i)%gamma2  = 1.773_dp
            params(i)%gamma3  = 2.543_dp

          case ("Ga47In53AsW")

            params(i)%meff    = 0.038_dp
            params(i)%EP      = 25.26_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.8166_dp
            params(i)%deltaSO = 0.362_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 11.97_dp
            params(i)%gamma2  = 4.36_dp
            params(i)%gamma3  = 5.15_dp

          case ("Al47In53AsW")

            params(i)%meff    = 0.0779_dp
            params(i)%EP      = 21.7_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.693_dp
            params(i)%deltaSO = 0.342_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 6.17_dp
            params(i)%gamma2  = 1.62_dp
            params(i)%gamma3  = 2.31_dp

          case ("InAs")

            params(i)%meff    = 0.026_dp
            params(i)%EP      = 21.5_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.417_dp
            params(i)%deltaSO = 0.39_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 20.0_dp
            params(i)%gamma2  = 8.5_dp
            params(i)%gamma3  = 9.2_dp
            params(i)%EV      = -0.59
            params(i)%EC      = -0.173

          case ("InAsW")

            params(i)%meff    = 0.0229_dp
            params(i)%EP      = 22.2_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.418_dp
            params(i)%deltaSO = 0.38_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 20.4_dp
            params(i)%gamma2  = 8.3_dp
            params(i)%gamma3  = 9.1_dp

          case ("AlAs")

            params(i)%meff    = 0.15_dp
            params(i)%EP      = 21.1_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 3.099_dp
            params(i)%deltaSO = 0.28_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 3.76_dp
            params(i)%gamma2  = 0.82_dp
            params(i)%gamma3  = 1.42_dp
            params(i)%EV      = -1.33
            params(i)%EC      = 1.769

          case ("AlAsW")

            params(i)%meff    = 0.15_dp
            params(i)%EP      = 21.12_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 3.13_dp
            params(i)%deltaSO = 0.3_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 3.25_dp
            params(i)%gamma2  = 0.65_dp
            params(i)%gamma3  = 1.21_dp

          case ("GaP")

            params(i)%meff    = 0.13_dp
            params(i)%EP      = 31.4_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 2.35_dp
            params(i)%deltaSO = 0.08_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 4.05_dp
            params(i)%gamma2  = 0.49_dp
            params(i)%gamma3  = 2.93_dp

          case ("AlP")

            params(i)%meff    = 0.22_dp
            params(i)%EP      = 17.7_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 2.52_dp
            params(i)%deltaSO = 0.07_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 3.35_dp
            params(i)%gamma2  = 0.71_dp
            params(i)%gamma3  = 1.23_dp

          case ("InP")

            params(i)%meff    = 0.0795_dp
            params(i)%EP      = 20.7_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.4236_dp
            params(i)%deltaSO = 0.108_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 5.08_dp
            params(i)%gamma2  = 1.60_dp
            params(i)%gamma3  = 2.10_dp

          case ("InPW")

            params(i)%meff    = 0.0803_dp
            params(i)%EP      = 20.56_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 1.423_dp
            params(i)%deltaSO = 0.11_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 4.95_dp
            params(i)%gamma2  = 1.65_dp
            params(i)%gamma3  = 2.35_dp

          case ("GaSb")

            params(i)%meff    = 0.039_dp
            params(i)%EP      = 27.0_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.812_dp
            params(i)%deltaSO = 0.76_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 13.4_dp
            params(i)%gamma2  = 4.7_dp
            params(i)%gamma3  = 6.0_dp

          case ("AlSb")

            params(i)%meff    = 0.14_dp
            params(i)%EP      = 18.7_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 2.386_dp
            params(i)%deltaSO = 0.676_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 5.18_dp
            params(i)%gamma2  = 1.19_dp
            params(i)%gamma3  = 1.97_dp
            params(i)%EV      = -0.41_dp
            params(i)%EC      = 1.976_dp

          case ("AlSbW")

            params(i)%meff    = 0.12_dp
            params(i)%EP      = 18.8_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 2.384_dp
            params(i)%deltaSO = 0.673_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 4.15_dp
            params(i)%gamma2  = 1.01_dp
            params(i)%gamma3  = 1.71_dp

          case ("InSb")

            params(i)%meff    = 0.0135_dp
            params(i)%EP      = 23.3_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.235_dp
            params(i)%deltaSO = 0.81_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 34.8_dp
            params(i)%gamma2  = 15.5_dp
            params(i)%gamma3  = 16.5_dp
            params(i)%EV      = 0
            params(i)%EC      = 0.235

          case ("InSbW")

            params(i)%meff    = 0.0139_dp
            params(i)%EP      = 24.4_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.237_dp
            params(i)%deltaSO = 0.81_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 37.1_dp
            params(i)%gamma2  = 16.5_dp
            params(i)%gamma3  = 17.7_dp

          case ("Vacuum")

            params(i)%meff    = 1!0.14_dp
            params(i)%EP      = 0!18.7_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 3.0_dp
            params(i)%deltaSO = 0.7021_dp
            gamma0         = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = -1!5.18_dp
            params(i)%gamma2  = 0!1.19_dp
            params(i)%gamma3  = 0!1.97_dp
            params(i)%EV      = -0.1326_dp
            params(i)%EC      = 2.8674_dp

          case ("Al")

            params(i)%meff    = 1
            params(i)%EP      = 21.999_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0
            params(i)%EV      = -11.6
            params(i)%EC      = -11.6
            params(i)%deltaSO = 0
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = -1
            params(i)%gamma2  = 0
            params(i)%gamma3  = 0

          case ("InAs10Sb90")

            params(i)%meff    = 0.0116_dp
            params(i)%EP      = 23.12_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.1749_dp
            params(i)%EV      = -0.504_dp
            params(i)%EC      = -0.1835_dp
            params(i)%deltaSO = 0.66_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 33.32_dp
            params(i)%gamma2  = 14.8_dp
            params(i)%gamma3  = 15.77_dp

          case ("InAs20Sb80")

            params(i)%meff    = 0.0104_dp
            params(i)%EP      = 22.94_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.1322_dp
            params(i)%EV      = -0.424_dp
            params(i)%EC      = -0.1826_dp
            params(i)%deltaSO = 0.534_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 31.84_dp
            params(i)%gamma2  = 14.1_dp
            params(i)%gamma3  = 15.04_dp

          case ("InAs30Sb70")

            params(i)%gamma1  = 30.36_dp
            params(i)%gamma2  = 13.4_dp
            params(i)%gamma3  = 14.31_dp
            params(i)%Eg      = 0.1069_dp
            params(i)%meff    = 0.0099_dp
            params(i)%deltaSO = 0.432_dp
            params(i)%EP      = 22.76_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.35
            params(i)%EC      = -0.1703
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("InAs40Sb60")

            params(i)%gamma1  = 28.88_dp
            params(i)%gamma2  = 12.7_dp
            params(i)%gamma3  = 13.58_dp
            params(i)%Eg      = 0.099_dp
            params(i)%meff    = 0.0101_dp
            params(i)%deltaSO = 0.354_dp
            params(i)%EP      = 22.58_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.282
            params(i)%EC      = -0.1466
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("InAs50Sb50")

            params(i)%meff    = 0.011_dp
            params(i)%EP      = 22.4_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.1085_dp
            params(i)%EV      = -0.22
            params(i)%EC      = -0.1115
            params(i)%deltaSO = 0.3_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 27.4_dp
            params(i)%gamma2  = 12.0_dp
            params(i)%gamma3  = 12.85_dp

          case ("InAs60Sb40")

            params(i)%gamma1  = 25.92_dp
            params(i)%gamma2  = 11.30_dp
            params(i)%gamma3  = 12.12_dp
            params(i)%Eg      = 0.1354_dp
            params(i)%meff    = 0.0126_dp
            params(i)%deltaSO = 0.27_dp
            params(i)%EP      = 22.22_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.164
            params(i)%EC      = -0.065
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("InAs70Sb30")

            params(i)%gamma1  = 24.44_dp
            params(i)%gamma2  = 10.60_dp
            params(i)%gamma3  = 11.39_dp
            params(i)%Eg      = 0.1797_dp
            params(i)%meff    = 0.0149_dp
            params(i)%deltaSO = 0.264_dp
            params(i)%EP      = 22.04_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.114
            params(i)%EC      = -0.0071
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("InAs80Sb20")

            params(i)%gamma1  = 22.96_dp
            params(i)%gamma2  = 9.90_dp
            params(i)%gamma3  = 10.66_dp
            params(i)%Eg      = 0.2414_dp
            params(i)%meff    = 0.0179_dp
            params(i)%deltaSO = 0.282_dp
            params(i)%EP      = 21.86_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.07
            params(i)%EC      = 0.0622
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("InAs90Sb10")

            params(i)%gamma1  = 21.48_dp
            params(i)%gamma2  = 9.20_dp
            params(i)%gamma3  = 9.93_dp
            params(i)%Eg      = 0.3205_dp
            params(i)%meff    = 0.0216_dp
            params(i)%deltaSO = 0.324_dp
            params(i)%EP      = 21.68_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%EV      = -0.032
            params(i)%EC      = 0.1429
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0

          case ("Al63In37Sb")

            params(i)%meff    = 0.0137_dp
            params(i)%EP      = 21.598_dp
            params(i)%P       = dsqrt(params(i)%EP*const)
            params(i)%Eg      = 0.9306_dp
            params(i)%EV      = -0.1326
            params(i)%EC      = 0.7981
            params(i)%deltaSO = 0.7021_dp
            gamma0            = 1.0_dp/params(i)%meff
            params(i)%A       = gamma0
            params(i)%gamma1  = 21.7616_dp
            params(i)%gamma2  = 9.0398_dp
            params(i)%gamma3  = 10.0529_dp

          case default
            stop "material not defined"

        end select


        if (renormalization .eqv. .True. .or. params(i)%A < 0) then

          EP = (gamma0-1.0_dp)*params(i)%Eg*(params(i)%Eg+params(i)%deltaSO)/(params(i)%Eg+(2.0_dp/3.0_dp)*params(i)%deltaSO)
          P = dsqrt(EP*const)
          A = - ((params(i)%Eg+(2.0_dp/3.0_dp)*params(i)%deltaSO)/(params(i)%Eg+params(i)%deltaSO))*(EP/params(i)%Eg) + gamma0

          gamma1 = const*(params(i)%gamma1 - EP/(3.*params(i)%Eg))
          gamma2 = const*(params(i)%gamma2 - EP/(6.*params(i)%Eg))
          gamma3 = const*(params(i)%gamma3 - EP/(6.*params(i)%Eg))

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
