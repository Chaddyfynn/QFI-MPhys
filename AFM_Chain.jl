using Sunny, GLMakie
include("QFI.jl")

units = Units(:meV, :angstrom)
a = 5 # (Å)
# Make the axes we don't care about large
latvecs = lattice_vectors(a, 2a, 5a, 90, 90, 90) 
positions = [[0 0 0]]
cryst = Crystal(latvecs, positions; )
view_crystal(cryst)

dims = (2,1,1)
sys = System(cryst, [1 => Moment(s=1, g=1)], :dipole;dims)
J = 1.0 # (>0 => AFM)

# Weak 3D Coupling
J2 = 1
J3 = 1
set_exchange!(sys, J2, Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, J3, Bond(1, 1, [0, 0, 1]))

# Anisotropy
# D = 1
# set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])

measure = ssf_trace(sys; )
swt = SpinWaveTheory(sys; measure) # constructs an object to do LSWT on 

kernel = lorentzian(fwhm=0.1) 
qs = [[0, 0, 0], [1/2, 0, 0], [1 0 0]]
path = q_space_path(cryst, qs, 500)
energies = range(0, 3, 300)

res = intensities(swt, path; energies, kernel)
plot_intensities(res;)

res_bands = intensities_bands(swt, path)
plot_intensities(res_bands; units)

# Intensities is like an experimental measurement
# Intensities Bands is like a hand calculation (Matrix)

nQFI_numerical = nqfi(res, cryst;)
nQFI_bands = nqfi(res_bands, cryst;)

max_nqfi(nQFI_numerical, path)
max_nqfi(nQFI_bands, path)

fig, ax = plot_nqfi(path, nQFI_numerical;)
fig, ax = plot_nqfi(path, nQFI_bands;)

# THIS ONE DOES DSSF
begin
    fig = Figure()
    ax1 = Axis(fig[1,1],xlabel="q", ylabel="Intensity (DSSF)")
    ax2 = Axis(fig[1,2], xlabel="q", ylabel="Log Intensity (DSSF)")
    xs = [q[1] for q in path.qs]
    ys = (res.data[1, :] +res.data[2, :] )
    logys = log10.(res.data[1, :] +res.data[2, :])
    lines!(ax1,xs, ys; )
    lines!(ax2,xs, logys; )
    fig    
end

res_int = [sum(res.data[:,i])/size(res.data,2) for i ∈ range(1,size(res.data,2))]

# Numerical Integration of res
begin
    fig = Figure()
    ax1 = Axis(fig[1,1],xlabel="q", ylabel="res Integral")
    xs = [q[1] for q in path.qs]
    ys = res_int
    lines!(ax1,xs, ys; )
    fig    
end

# Using intensities bands
begin
    fig = Figure()
    kT = 0
    # Extract the first component in the q path (x val)
    xs = [q[1] for q in path.qs]
    # Is this numerical integration of the res_bands ?
    # ys1 = res_bands.data[1,:] + res_bands.data[2,:]
    ys1 = vec(sum(res_bands.data; dims=1))
    ys2 = vec(sum(res_bands.disp; dims=1)) # ħω/2β ?
    # dot syntax is for vectorised operations
    ħω_part = tanh.(ys2 ./ 2kT).*(1 .- exp.(-ys2 ./ kT))
    S2 = 1 ./(4ys1)
    QFI = ys1 .* ħω_part
    ax = Axis(fig[1, 1],
    xlabel = "Q_x",
    ylabel = "nQFI",
)
    lines!(ax, xs, QFI)
    fig
end

# Using intensities
begin
    kT = 0
    # Extract the first component in the q path (x val)
    xs = [q[1] for q in path.qs]
    ys1 = [sum(res.data[:,i])/size(res.data,2) for i ∈ range(1,size(res.data,2))]
    ys2 = res_bands.disp[1,:] + res_bands.disp[2,:]
    # dot syntax is for vectorised operations
    ħω_part = tanh.(ys2 ./ 2kT).*(1 .- exp.(-ys2 ./ kT))
    S2 = 1 ./(4ys1)
    QFI = ys1 .* ħω_part
    fig = lines(xs, QFI; axis=(xlabel="q", ylabel="QFI"))
end