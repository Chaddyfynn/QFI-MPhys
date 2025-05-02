using Sunny, GLMakie, Random, LinearAlgebra
include("QFI.jl")

function NiPS3_crystal(;dims=(2,2,1))
    a = 5.811
    b = 10.064
    c = 6.594
    S=1
    γ = 1 - 1/(2S)
    latvecs = lattice_vectors(a, b, c, 90, 106.997, 90)
    positions = [[0.,    0.33335,    0. ]]
    cryst = Crystal(latvecs, positions,12)
    spininfos = [1 => Moment(; s=1, g=2)] 
    sys = System(cryst, spininfos, :dipole; dims)
    J1a = -2.7
    J1b = −2.0
    J2 = 0.2
    J3 = 13.9
    J4 = −0.38
    Ax = −0.010/γ
    Az = 0.21/γ
    set_exchange!(sys,J1a ,Bond(1,2,[0,0,0]))
    set_exchange!(sys,J1b ,Bond(2,3,[0,0,0]))
    set_exchange!(sys,J2,Bond(1,3,[0,0,0]))
    set_exchange!(sys,J2,Bond(1,1,[1,0,0]))
    set_exchange!(sys,J3,Bond(2,3,[1,0,0]))
    set_exchange!(sys,J3,Bond(1,4,[0,0,0]))
    set_exchange!(sys,J4 ,Bond(1,1,[0,0,1]))
    set_exchange!(sys,J4 ,Bond(1,2,[0,0,-1])) #not sure if they included this 
    set_onsite_coupling!(sys, S -> Ax*S[1]^2 + Az*S[3]^2, 1)
    randomize_spins!(sys)
    minimize_energy!(sys)
    return sys, cryst
end

sys, cryst = NiPS3_crystal(dims=(1,1,1)) # If you change values play around with dims to check the GS still fits inside the crystallographic unit cell
print_wrapped_intensities(sys)
plot_spins(sys;color = [s[3] for s ∈ sys.dipoles ])

# measure = ssf_trace(sys; ) # swap to ssf_trace for our calculations
measure = ssf_custom((q, ssf) -> real(ssf[3,3]), sys;apply_g=false) # individual components
swt = SpinWaveTheory(sys; measure)
kernel = lorentzian(fwhm=1.25)
qs = [[1, 3, 0], [1/2, 3, 0], [1/2, 5/2, 0],[1,5/2,0], [1, 3, 0]]

path = q_space_path(cryst, qs,800)
res_bands = intensities_bands(swt, path; )
# plot_intensities(res_bands; )

kernel = lorentzian(fwhm=1)
energies = range(0, 60, 300);  # 0 < ω < 10 (meV)

# res = intensities(swt, path; energies, kernel )
# plot_intensities(res;colormap = :viridis)

# res_bands = intensities_bands(swt, path)
# plot_intensities(res_bands;)

# q_eval = [1/2,5/2,0]
# q_index = findfirst(item -> item == q_eval, path.qs)

# nQFI_NiPS3 = nqfi(res, cryst; prefactor=0.5)
nQFI_NiPS3_bands = nqfi(res_bands, cryst; prefactor=0.5)

# labels = ["Γ", "Y", "C", "Z", "Γ"]
# plot_nqfi(path, nQFI_NiPS3; xlabels=labels)
# plot_nqfi(path, nQFI_NiPS3_bands; xlabels=labels)

max_nqfi(nQFI_NiPS3_bands, path)