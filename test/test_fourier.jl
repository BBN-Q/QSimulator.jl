using QuadGK: quadgk
using PyPlot
const plt = PyPlot
plt.ioff()

EC, EJ1, EJ2 = 0.172, 12.71, 3.69
freq = 0.1
amp = 0.6
harmonics = collect(-50:50)
parks = [0.0, 0.25]

@testset "FourierSeries" begin
	times = collect(0.0:.01:50)
	for park in parks
		ϕ(t) = park + amp * sin(2π * freq * t)
		ω(t) = perturbative_transmon_freq(EC, EJ1, EJ2, ϕ(t))
		fs = FourierSeries(ω, freq, harmonics)
		ans = eval_series(fs, times)
		@test isapprox(ω.(times), ans, rtol=1e-9)
		if "plot" in ARGS
			plt.plot(times, ω.(times))
			plt.plot(times, ans, "--")
			plt.show()
		end
	end
end

@testset "rotating_frame_series" begin
	times = collect(0.0:.05:10)
	for park in parks
		ϕ(t) = park + amp * sin(2π * freq * t)
		ω(t) = perturbative_transmon_freq(EC, EJ1, EJ2, ϕ(t))
		fs = FourierSeries(ω, freq, harmonics)
		Φ(t) = exp(1im * quadgk(ω, 0.0, t)[1])
		cfs = rotating_frame_series(fs, harmonics)
		# note that we must multiply by the term corresponding to fs.terms[0]
		# as this was ignored in rotating_frame_series
		ans = exp.(1im * fs.terms[0] * times) .* eval_series(cfs, times)
		@test isapprox(Φ.(times), ans, rtol=1e-6)
		if "plot" in ARGS
			plt.plot(times, real.(Φ.(times)))
			plt.plot(times, real.(ans), "--")
			plt.show()
		end
	end
end
