using Sunny, GLMakie
include("QFI.jl")

# Define lattice constants for YbMgGaO4
a = 3.403  
c = 25.145 

latvecs = lattice_vectors(a, a, c, 90, 90, 120)  

positions = [[0,0,0]] 

cryst = Crystal(latvecs, positions, 1)
view_crystal(cryst; ndims=2)

dims=(10,10,1)
spininfos = [1 => Moment(; s=1, g=1)] #g=2
sys = System(cryst, spininfos, :dipole; dims)

J1pm   = -0.236 # (meV)
J1pmpm = -0.161
J1zpm  = -0.261
J2pm   = 0.026
J3pm   = 0.166
J′0pm  = 0.037
J′1pm  = 0.013
J′2apm = 0.068

J1zz   = -0.236
J2zz   = 0.113
J3zz   = 0.211
J′0zz  = -0.036
J′1zz  = 0.051
J′2azz = 0.073

J1xx = J1pm + J1pmpm
J1yy = J1pm - J1pmpm
J1yz = J1zpm

set_exchange!(sys, [J1xx   0.0    0.0;
                    0.0    J1yy   J1yz;
                    0.0    J1yz   J1zz], Bond(1,1,[1,0,0]))
set_exchange!(sys, [J2pm   0.0    0.0;
                    0.0    J2pm   0.0;
                    0.0    0.0    J2zz], Bond(1,1,[1,2,0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)


measure = ssf_trace(sys; apply_g = false)
#measure = ssf_custom((q, ssf) -> real(ssf[1,1]), sys)
swt = SpinWaveTheory(sys; measure) 

kernel = lorentzian(fwhm=0.1)

print_irreducible_bz_paths(cryst)

Z = [0, 0, 1/2]
V₂ = [1/2, -1/2, 0]
U₂ = [-1/2, 0, 1/2]
Γ = [0, 0, 0]
Y = [0, 1/2, 0]
X = [1/2, 0, 0]
R₂ = [-1/2, -1/2, 1/2]
T₂ = [0, -1/2, 1/2]

qs_1 = [Γ, X, Γ]
qs_2 = [Y, Γ, Z, Y]
qs_3 = [R₂, Γ, T₂, R₂]
qs_4 = [U₂, Γ, V₂, U₂]

qs = [[0, 0, 0], [1/2, 0, 0], [1 0 0]]
path = q_space_path(cryst, qs, 500)
energies = range(0, 1.0, 4000) 
res = intensities(swt, path; energies, kernel)
plot_intensities(res;)

# res_bands = intensities_bands(swt, path)
# plot_intensities(res_bands; )


nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
fig = plot_nqfi(path, nQFI_Crystal; scale=1.0)
save("High Symmetry R2 Γ T2 Path YbMgGaO4.png", fig)

max_nqfi(nQFI_Crystal, path)

sys_inhom = to_inhomogeneous(repeat_periodically(sys, (2, 2, 1)))
for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
    noise = randn()/3
    set_exchange_at!(sys_inhom, [J1xx   0.0    0.0;
                                 0.0    J1yy   0.0;
                                0.0    0.0   J1zz] + [1 0 0; 0 1 0; 0 0 1] .* noise, site1, site2; offset)
end

randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom, maxiters=5_000)
plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)

measure = ssf_trace(sys_inhom; apply_g = false)
#measure = ssf_custom((q, ssf) -> real(ssf[1,1]), sys)
swt_dis = SpinWaveTheoryKPM(sys_inhom; measure, tol=0.01) 
swt_dis1 = SpinWaveTheory(sys_inhom; measure) 

kernel = lorentzian(fwhm=0.1)
qs = [[0, 0, 0], [1/2, 0, 0], [1 0 0]]
path = q_space_path(cryst, qs, 500)
energies = range(0, 3.5, 4000) 
res_dis1 = intensities(swt_dis1, path; energies, kernel)
plot_intensities(res_dis1;)

res_dis = intensities_bands(swt_dis1, path)
plot_intensities(res_dis; )

nQFI_Crystal = nqfi(res_dis; prefactor=1.0, cryst)

plot_nqfi(path, nQFI_Crystal; scale=1.0)

maximum_QFI(nQFI_Crystal)