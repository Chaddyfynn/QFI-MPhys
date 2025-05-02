using Sunny, GLMakie
include("QFI.jl")

function YbMgGaO4(noise_amp)
    # Define lattice constants for YbMgGaO4
    a = 3.403  
    c = 25.145 
    # Configure Crystal
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)  
    positions = [[0,0,0]] 
    cryst = Crystal(latvecs, positions, 1)
    # view_crystal(cryst; ndims=2)
    # Setup Spins
    dims=(10,10,1)
    spininfos = [1 => Moment(; s=1, g=1)] #g=2
    sys = System(cryst, spininfos, :dipole; dims)

    # Hamiltonian
    J1pm = 0.1094
    J1pmpm = 0.013
    J1zz = 0.1264
    J2pm = 0.024
    J2zz = 0.0278
    J1xx = J1pm + J1pmpm
    J1yy = J1pm - J1pmpm

    set_exchange!(sys, [J1xx   0.0    0.0;
                        0.0    J1yy   0.0;
                        0.0    0.0   J1zz], Bond(1,1,[1,0,0]))
    set_exchange!(sys, [J2pm   0.0    0.0;
                        0.0    J2pm   0.0;
                        0.0    0.0    J2zz], Bond(1,1,[1,2,0]))

    # Noise
    sys_inhom = to_inhomogeneous(repeat_periodically(sys, (1, 1, 1)))
    for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
        noise = randn()/3
        set_exchange_at!(sys_inhom, [J1xx   0.0    0.0;
                                    0.0    J1yy   0.0;
                                    0.0    0.0   J1zz] + [1 0 0; 0 1 0; 0 0 1] .* noise .* noise_amp, site1, site2; offset)
    end

    # Ground Sate
    randomize_spins!(sys_inhom)
    minimize_energy!(sys_inhom, maxiters=5_000)
    # plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)
    return sys_inhom, cryst
end

# Iterates over val_range for cryst_func, plots max nQFI vs val_range
function examine_effect(cryst_func::Function, val_range, qs; y_text::Float64=1.0)
    max_nqfis = []
    locations = []
    for val ∈ val_range
        print("Calculaing cryst_funct at $val\n")
        sys, cryst = cryst_func(val)
        measure = ssf_trace(sys; apply_g = false)
        swt = SpinWaveTheory(sys; measure) 
        kernel = lorentzian(fwhm=0.1)
        path = q_space_path(cryst, qs, 500)
        energies = range(0, 1, 500) 
        res = intensities(swt, path; energies, kernel)
        nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
        max_val = max_nqfi(nQFI_Crystal, path)
        push!(max_nqfis, max_val[1])
        push!(locations, max_val[2])
    end
    fig = Figure()
    ax = Axis(
    fig[1, 1], 
    xlabel="Parameter",
    ylabel="Max nQFI"
    )
    lines!(ax, val_range, max_nqfis)
    text!(val_range[1], max_nqfis[1] + findmax(max_nqfis)[1]/33; text=string(locations[1]))
    for i in range(2, size(val_range,1))
        if locations[i] != locations[i-1]
            text!(val_range[i], max_nqfis[i] + findmax(max_nqfis)/33 * y_text; text=string.(locations[i]))
        end
    end
    ylims!(ax, 0, findmax(max_nqfis)[1]*1.2)
    fig
    return fig, ax
end

sys_inhom, cryst = YbMgGaO4(0)

measure = ssf_trace(sys_inhom; apply_g = false)
#measure = ssf_custom((q, ssf) -> real(ssf[1,1]), sys)
swt_dis_kpm = SpinWaveTheoryKPM(sys_inhom; measure, tol=0.01) 
swt_dis = SpinWaveTheory(sys_inhom; measure) 

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
qs_2 = [Y, Γ, Z]
qs_3 = [R₂, Γ, T₂]
qs_4 = [U₂, Γ, V₂]

qs = [[0, 0, 0], [1/2, 0, 0], [1 0 0]]
path = q_space_path(cryst, qs_3, 500)
energies = range(0, 1, 500) 
res_dis = intensities(swt_dis, path; energies, kernel)
plot_intensities(res_dis;)
# res_dis_bands = intensities_bands(swt_dis, path;)
# plot_intensities(res_dis_bands;)

# res_di_kpms = intensities(swt_dis_kpm, path; energies, kernel)
# plot_intensities(res_dis_kpm; )
# res_di_kpms_bands = intensities_bands(swt_dis_kpm, path)
# plot_intensities(res_dis_kpm_bands; )

nQFI_Crystal = nqfi(res_dis, cryst; prefactor=1.0)

plot_nqfi(path, nQFI_Crystal; scale=1.0)

max_nqfi(nQFI_Crystal, path)

fig, ax = examine_effect(YbMgGaO4, range(1e-3,0.5,5), qs_3;)
fig


sys, cryst = YbMgGaO4(1)
measure = ssf_trace(sys_inhom; apply_g = false)
swt = SpinWaveTheory(sys_inhom; measure) 
kernel = lorentzian(fwhm=0.1)
path = q_space_path(cryst, qs_3, 500)
energies = range(0, 1, 500) 
res = intensities(swt, path; energies, kernel)
plot_intensities(res;)
nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
plot_nqfi(path, nQFI_Crystal)
max_val = max_nqfi(nQFI_Crystal, path)