using Sunny, GLMakie # import Sunny and a plotting package GLMakie. 
include("QFI.jl")

# Set up a cubic crystal with 2 of the axes much longer than the x-axis (which will be our chain direction)
# the length of these other axes are not important, we just want them > a
units = Units(:meV, :angstrom)
a = 5 # (Å)
latvecs = lattice_vectors(a, 2a, 5a, 90, 90, 90) 
positions = [[0 0 0]] #put atoms of the vertices of this lattice
cryst = Crystal(latvecs, positions; )
view_crystal(cryst) #we can plot this and view the atoms/bonds

# Now we define a spin system. This means telling the code what spins each atom has
# We need to tell Sunny what crystal to use (the one we created above), the magnetic moment of the atom (here we choose site one to have an 
# S=1) and the mode of calculation, here we choose dipole. 
# Dipole mode is standard LSWT i.e. we consider only a single flavor of Holstein Primakoff boson. This is akin to treating each site as a 2-level system.
# Note that g is the g-factor which is 2 for a spin only magnetic ion.
# We also specify the kwarg "dims" which is the size of the supercell. Since we have an AFM, the supercell should be (2,1,1)
dims = (2,1,1)
sys = System(cryst, [1 => Moment(s=1, g=1)], :dipole;dims)
# Now we assign bonds to this system. The labeling mirrors that of the view_crystal function
J = - 1.0 # (<0 => FM)
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
# we randomize the spins and then perform an energy minimization of the 
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])


# S^αβ(q,w) is a 3x3 matrix with the elements corresponding to α,β ∈ [x,y,z]. We need to specify how we wish to contract this matrix to plot.
# In conventional neutron scattering they measure the component perpendicular to  the momentum transfer. To specify this we use ssf_perp(sys). 
# The trace is specificed by ssf_trace
measure = ssf_trace(sys; )
# measure = ssf_custom((q, ssf) -> real(ssf[3,3]), sys; apply_g=false) # this is the Lehmann representation of the dynamical susceptibility.
swt = SpinWaveTheory(sys; measure) # constructs an object to do LSWT on 

kernel = lorentzian(fwhm=0.01) # chooses a kernel. In the Lehmann Representation of the dynamical susceptibilty we saw a delta function. This here approximates this with a finite width Lorentzian
# Choose some reciprocal space path. Since it is 1-d only 1 direction makes sense.
qs = [[0, 0, 0], [1/2, 0, 0], [1 0 0]]
path = q_space_path(cryst, qs, 500)
energies = range(0, 3, 1000) # choose energy range
# calculate ∑ⱼSʲʲ(q,ω) and plot.
res = intensities(swt, path; energies, kernel)
plot_intensities(res;)

res_bands = intensities_bands(swt, path)
plot_intensities(res_bands; units)

nQFI_numerical = nqfi(res, cryst; )
nQFI_bands = nqfi(res_bands, cryst; )

max_nqfi(nQFI_numerical, path)
max_nqfi(nQFI_bands, path)

plot_nqfi(path, nQFI_numerical;)
fig, ax = plot_nqfi(path, nQFI_bands;)
ylims!(ax, 0, 1.1)
fig

# Intensities is like an experimental measurement
# Intensities Bands is like a hand calculation (Matrix)

begin
    fig = Figure()
    ax1 = Axis(fig[1,1],xlabel="q", ylabel="Intensity")
    ax2 = Axis(fig[1,2], xlabel="q", ylabel="Log Intensity")
    xs = [q[1] for q in path.qs]
    ys = (res.data[1, :] +res.data[2, :] )
    logys = log10.(res.data[1, :] +res.data[2, :])
    lines!(ax1,xs, ys; )
    lines!(ax2,xs, logys; )
    fig    
end

res_int = [sum(res.data[:,i])/size(res.data,2) for i ∈ range(1,size(res.data,2))]

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
    # Extract the first component in the q path (x val)
    xs = [q[1] for q in path.qs]
    # Is this numerical integration of `the res_bands ?
    ys1 = res_bands.data[1,:] + res_bands.data[2,:]
    ys2 = res_bands.disp[1,:] + res_bands.disp[2,:] # ħω/2β ?
    # dot syntax is for vectorised operations
    ħω_part = tanh.(ys2).*(1 .- exp.(-ys2))
    S2 = 1 ./(4ys1)
    QFI = ys1 .* ħω_part
    # QFI =  ħω_part
    fig = lines(xs, QFI; axis=(xlabel="q", ylabel="QFI"))
end

# Using intensities
begin
    # Extract the first component in the q path (x val)
    xs = [q[1] for q in path.qs]
    ys1 = [sum(res.data[:,i])/size(res.data,2) for i ∈ range(1,size(res.data,2))]
    ys2 = res_bands.disp[1,:] + res_bands.disp[2,:]
    # dot syntax is for vectorised operations
    ħω_part = tanh.(ys2).*(1 .- exp.(-ys2))
    S2 = 1 ./(4ys1)
    QFI = ys1 .* ħω_part
    #QFI =  ħω_part
    fig = lines(xs, QFI; axis=(xlabel="q", ylabel="QFI"))
end