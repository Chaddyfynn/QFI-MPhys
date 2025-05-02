# # SW15 - Ba₃NbFe₃Si₂O₁₄
#
# This is a Sunny port of [SpinW Tutorial
# 15](https://spinw.org/tutorials/15tutorial), originally authored by Sandor
# Toth. It calculates the linear spin wave theory spectrum of Ba₃NbFe₃Si₂O₁₄.
# The ground state is an incommensurate spiral, which can be directly studied
# using the functions [`minimize_spiral_energy!`](@ref) and
# [`SpinWaveTheorySpiral`](@ref).

# Load packages 

using Sunny, GLMakie
include("QFI.jl")

# Specify the Ba₃NbFe₃Si₂O₁₄ [`Crystal`](@ref) cell following [Marty et al.,
# Phys. Rev. Lett. **101**, 247201
# (2008)](http://dx.doi.org/10.1103/PhysRevLett.101.247201).

units = Units(:meV, :angstrom)
a = b = 8.539 # (Å)
c = 5.2414
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
types = ["Fe", "Nb", "Ba", "Si", "O", "O", "O"]
positions = [[0.24964,0,0.5], [0,0,0], [0.56598,0,0], [2/3,1/3,0.5220],
             [2/3,1/3,0.2162], [0.5259,0.7024,0.3536], [0.7840,0.9002,0.7760]]
langasite = Crystal(latvecs, positions, 150; types)
cryst = subcrystal(langasite, "Fe")
view_crystal(cryst)

# Create a [`System`](@ref) and set exchange interactions as parametrized in
# [Loire et al., Phys. Rev. Lett. **106**, 207201
# (2011)](http://dx.doi.org/10.1103/PhysRevLett.106.207201).

sys = System(cryst, [1 => Moment(s=5/2, g=2)], :dipole)
J₁ = 0.85
J₂ = 0.24
J₃ = 0.053
J₄ = 0.017
J₅ = 0.24
set_exchange!(sys, J₁, Bond(3, 2, [1,1,0]))
set_exchange!(sys, J₄, Bond(1, 1, [0,0,1]))
set_exchange!(sys, J₂, Bond(1, 3, [0,0,0]))

# The final two exchanges are set according to the desired chirality ``ϵ_T`` of
# the magnetic structure.

ϵT = -1
if ϵT == -1
    set_exchange!(sys, J₃, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₅, Bond(3, 2, [1,1,1]))
elseif ϵT == 1
    set_exchange!(sys, J₅, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₃, Bond(3, 2, [1,1,1]))
else
    error("Chirality must be ±1")
end

# This compound is known to have a spiral order with approximate propagation
# wavevector ``𝐤 ≈ [0, 0, 1/7]``. Search for this magnetic order with
# [`minimize_spiral_energy!`](@ref). Due to reflection symmetry, one of two
# possible propagation wavevectors may appear, ``𝐤 = ± [0, 0, 0.1426…]``.
# Note that ``k_z = 0.1426…`` is very close to ``1/7 = 0.1428…``.

axis = [0, 0, 1]
randomize_spins!(sys)
k = minimize_spiral_energy!(sys, axis)

# We can visualize the full magnetic cell using [`repeat_periodically_as_spiral`](@ref),
# which includes 7 rotated copies of the chemical cell.

sys_enlarged = repeat_periodically_as_spiral(sys, (1, 1, 7); k, axis)
plot_spins(sys_enlarged; color=[S[1] for S in sys_enlarged.dipoles])

# One could perform a spin wave calculation using either
# [`SpinWaveTheory`](@ref) on `sys_enlarged`, or [`SpinWaveTheorySpiral`](@ref)
# on the original `sys`. The latter has some restrictions on the interactions,
# but allows for our slightly incommensurate wavevector ``𝐤``.

measure = ssf_perp(sys)
swt = SpinWaveTheorySpiral(sys; measure, k, axis)

# Calculate broadened intensities for a path ``[0, 1, L]`` through reciprocal
# space

print_irreducible_bz_paths(cryst)

Sunny.irreducible_bz_paths(cryst)[2]

qs = [[0, 1, -1], [0, 1, -1+1], [0, 1, -1+2], [0, 1, -1+3]]
path = q_space_path(cryst, qs, 400)
energies = range(0, 6, 400)
res = intensities(swt, path; energies, kernel=gaussian(fwhm=0.25))
plot_intensities(res; units, saturation=0.7, colormap=:jet, title="Scattering intensities")

# Use [`ssf_custom_bm`](@ref) to calculate the imaginary part of
# ``\mathcal{S}^{2, 3}(𝐪, ω) - \mathcal{S}^{3, 2}(𝐪, ω)``. In polarized
# neutron scattering, it is conventional to express the 3×3 structure factor
# matrix ``\mathcal{S}^{α, β}(𝐪, ω)`` in the Blume-Maleev polarization axis
# system. Specify the scattering plane ``[0, K, L]`` via the spanning vectors
# ``𝐮 = [0, 1, 0]`` and ``𝐯 = [0, 0, 1]``.
measure = ssf_custom_bm(sys; u=[0, 1, 0], v=[0, 0, 1]) do q, ssf
    imag(ssf[2,3] - ssf[3,2])
end
swt = SpinWaveTheorySpiral(sys; measure, k, axis)
res = intensities(swt, path; energies, kernel=gaussian(fwhm=0.25))
plot_intensities(res; units, saturation=0.8, allpositive=false,
                 title="Im[S²³(q, ω) - S³²(q, ω)]")

nqfi_crystal = nqfi(res, cryst; prefactor=0.5)
plot_nqfi(path, nqfi_crystal)