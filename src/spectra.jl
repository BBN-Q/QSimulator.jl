export #Types
       SpectrumApproximation,
       #Methods
       approx_spectrum,
       #Values
       mathieu_spectrum_EJoverEC_100_400

type SpectrumApproximation
    label::String
    max_dim::Int
    lower_bound::Float64
    upper_bound::Float64
    energy::Vector{Function}
end

function approx_spectrum( dim::Int, E_J::Float64, E_C::Float64, approx::SpectrumApproximation )
    x = E_J/E_C

    if x < approx.lower_bound || x > approx.upper_bound
        error("Mathieu spectrum approximation only implemented for $(approx.lower_bound) < Ej/Ec < $(approx.upper_bound)")
    end

    if dim > approx.max_dim
        error("Spectrum approximation $(approx.label) only implemented up to $(approx.max_dim) levels.")
    end

    return Float64[ E_C*energy[i](x) for i in 1:dim ]
end

#Computing the Mathieu characteristic function numerically requires
# care, so a simpler approach is to generate Chebyshev polynomial
# approximations to the transmon spectrum with Mathematica. The
# broader the EJ/EC ratio these approximation need to cover, the
# higher the order of the polynomial, so it is simple to choose the
# smallest range that is covered by the asymmetry of the transmon (if
# it is runable) or by the expected parameter range of the transmon
# (if it is fixed frequency). For efficiency's sake, it is also a good
# idea to write the polynomials in Horner form.
transmon_EoverEC_0_100_400(x) =
    3.27079 + 
     x*(-0.825172 + 
        x*(-0.00126499 + 
           x*(9.51911E-6 + 
              x*(-5.46144E-8 + 
                 x*(2.28516E-10 + 
                    x*(-6.82832E-13 + 
                       x*(1.41587E-15 + 
                        x*(-1.93291E-18 + (1.56096E-21 - 
                        5.64582E-25*x)*x))))))))

transmon_EoverEC_1_100_400(x) =
    9.22171 + 
     x*(-0.473818 + 
        x*(-0.00381804 + 
           x*(0.0000287684 + 
              x*(-1.65179E-7 + 
                 x*(6.91475E-10 + 
                    x*(-2.06691E-12 + 
                       x*(4.28689E-15 + 
                        x*(-5.85347E-18 + (4.7278E-21 - 
                        1.7102E-24*x)*x))))))))

transmon_EoverEC_2_100_400(x) =
    13.8518 + 
     x*(-0.115993 + 
        x*(-0.0064614 + 
           x*(0.0000488557 + 
              x*(-2.81089E-7 + 
                 x*(1.17829E-9 + 
                    x*(-3.52542E-12 + 
                       x*(7.31704E-15 + 
                        x*(-9.99632E-18 + (8.07738E-21 - 
                        2.92285E-24*x)*x))))))))


transmon_EoverEC_3_100_400(x) =
    16.7595 + 
     x*(0.258476 + 
        x*(-0.00935046 + 
           x*(0.0000713021 + 
              x*(-4.12412E-7 + 
                 x*(1.7351E-9 + 
                    x*(-5.20532E-12 + 
                       x*(1.08259E-14 + 
                        x*(-1.48139E-17 + (1.19858E-20 - 
                        4.34182E-24*x)*x))))))))


transmon_EoverEC_4_100_400(x) =
    17.4525 + 
     x*(0.659739 + 
        x*(-0.0126029 + 
           x*(0.0000969444 + 
              x*(-5.62932E-7 + 
                 x*(2.37253E-9 + 
                    x*(-7.12232E-12 + 
                       x*(1.48141E-14 + 
                        x*(-2.02669E-17 + (1.63916E-20 - 
                        5.93503E-24*x)*x))))))))

mathieu_spectrum_EJoverEC_100_400 = SpectrumApproximation("Chebyshev polynomial for EJ/EC in [100,400]",
                                                          5, 100., 400., 
                                                          [transmon_EoverEC_0_100_400, 
                                                           transmon_EoverEC_1_100_400,
                                                           transmon_EoverEC_2_100_400,
                                                           transmon_EoverEC_3_100_400,
                                                           transmon_EoverEC_4_100_400] )


