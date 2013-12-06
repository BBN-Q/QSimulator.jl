export FFTransmon,
       Resonator,
       TunableTransmon,
       Qubit

#Resonator 
type Resonator <: QSystem
    label::String
    freq::Float64
    dim::Int
end 
hamiltonian(r::Resonator) = r.freq*number(r)

#Helper functions to approximate exact Mathieu spectrum for 
# 400 > Ej/Ec > 100, n_g ~ .499, and at most 4 levels. Going to
# these lengths for the "self" Hamiltonian without making similar
# efforts in the dipole interation is likely to be questionable.
mathieu_spectrum0(x) =  4.3076 + 
                        x*(-0.867497 + 
                          x*(-0.000513746 + 
                            x*(1.84278E-6   + 
                              x*(-4.27568E-9 + 
                                (5.57896E-12 - 3.09745E-15 x)*x))));
mathieu_spectrum1(x) =  2.3556 + 
                        x*(-0.601775 + 
                          x*(-0.00154639  + 
                            x*(5.55068E-6   + 
                              x*(-1.28839E-8 + 
                                (1.68153E-11 - 9.33749E-15 x)*x))));
mathieu_spectrum2(x) = 17.471  + 
                       x*(-0.274884 + 
                         x*(-0.003428    + 
                           x*(0.0000156315 + 
                             x*(-4.92999E-8 + 
                               x*(9.86038E-11 + 
                                 (-1.12111E-13 +5.5065E-17 x)*x)))));
mathieu_spectrum3(x) = 20.0562 + 
                       x*(0.104463  + 
                         x*(-0.00621289  + 
                           x*(0.000034579  + 
                             x*(-1.38748E-7 + 
                               x*(3.77111E-10 + 
                                 x*(-6.5644E-13  + 
                                   (6.58996E-16 - 2.89781E-19 x)*x))))));
mathieu_spectrum4(x) = 17.6273 + 
                       x*(0.651625  + 
                         x*(-0.0124386   + 
                           x*(0.0000950317 + 
                             x*(-5.48747E-7 + 
                               x*(2.30243E-9  + 
                                 x*(-6.88843E-12 + 
                                   x*(1.42932E-14 + 
                                     x*(-1.95252E-17 + 
                                       (1.57814E-20 - 5.71476E-24 x)*x))))))));

function approx_mathieu_spectrum( dim::Int, E_J::Float64, E_C::Float64 )
    x = E_J/E_C

    if x < 100. || x > 400.
        error("Mathieu spectrum approximation only implemented for 400. > Ej/Ec > 100.")
    end

    if dim == 1
        error("A one level transmon ... really?")
    else if dim == 2
        [mathieu_spectrum0(x), mathieu_spectrum1(x)]
    else if dim == 3
        [mathieu_spectrum0(x), mathieu_spectrum1(x), mathieu_spectrum2(x)]
    else if dim == 4
        [mathieu_spectrum0(x), mathieu_spectrum1(x), mathieu_spectrum2(x), mathieu_spectrum3(x)]
    else if dim == 5
        [mathieu_spectrum0(x), mathieu_spectrum1(x), mathieu_spectrum2(x), mathieu_spectrum3(x), mathieu_spectrum4(x)]
    else
        error("Mathieu spectrum approximation only implemented up to 5 levels.")
    end

end

#Fixed frequency transmon 
type FFTransmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64
    dim::Int
    hamiltonian::Matrix{Float64}
    #Inner constructor needed to make sure the hamiltonian is properly innitalized
    function FFTransmon(label::String, E_C::Float64, E_J::Float64, dim::Int)
        new(label, E_C, E_J, dim, diagm(E_C*approx_mathieu_spectrum(dim,E_J,E_C)))
    end
end 

function hamiltonian(tt::FFTransmon, t::Float64=0.0)
    return tt.hamiltonian
end

#Tunable transmon 
type TunableTransmon <: QSystem
    label::String
    E_C::Float64
    E_J::Float64 #sum of junction E_J's
    d::Float64 #asymmetry parameter
    dim::Int
    fluxBias::Float64 # flux bias in units of Phi_0
    flux::Float64 #total flux in units of Phi_0
end

TunableTransmon(label::String, 
                E_C::Float64, 
                E_J::Float64, 
                d::Float64, 
                dim::Int, 
                fluxBias::Float64) = TunableTransmon(label, E_C, E_J, d, dim, fluxBias, fluxBias)

#Helper function to calculate effective EJ for a transmon
scale_EJ(flux::Float64, d::Float64) = cos(pi*flux)*sqrt(1 + d^2*(tan(pi*flux)^2))

function hamiltonian(tt::TunableTransmon, t::Float64=0.0)
    myE_J = tt.E_J*scale_EJ(tt.flux, tt.d)
    return diagm(tt.E_C*approx_mathieu_spectrum(tt.dim,myE_J,tt.E_C))
end

#Basic two-level qubit
type Qubit <: QSystem
    label::String
    freq::Float64
end
dim(q::Qubit) = 2
hamiltonian(q::Qubit) = q.freq*number(q)
