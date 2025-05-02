using Pkg
Pkg.activate("Inhom")
using Sunny, GLMakie, DelimitedFiles

PATH = "./Figs/Production/YMGO/"


include("QFI.jl")
include("Crystals.jl")
include("dom_wrapped_intensity.jl")

# Simple Crystal and System
sys, cryst = YbMgGaO4(0.5)
#print_irreducible_bz_paths(cryst)
print_wrapped_intensities1(sys)
# 0 Noise: [1/2,0,0] (collinear in y=x), [0, 1/2, 0] (collinear in y - vertical), [1/2, 1/2, 0] (collinear in y=-x)
# plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

# High Symmetry Paths/Points

M = [1/2, 0, 0]
A = [0, 0, 1/2]
H = [1/3, 1/3, 1/2]
K = [1/3, 1/3, 0]
Γ = [0, 0, 0]
L = [1/2, 0, 1/2]

pt_1 = [0.2777, 0.2777, 0]

qs_1 = [Γ, M, K, Γ, A, L, H, A]
qs_1_half_1 = [Γ, M, K]
qs_1_half_2 = [Γ, L, H, A] # This one seems to contain enough information ordered
qs_2 = [L, M]
qs_3 = [H, K]

width_val = 1/30 # Must be small otherwise increase path_length arg in examine_effect
width = [width_val,width_val,width_val]
qs_fast = [L-width, L, L+width, K-width, K, K+width, pt_1-width, pt_1, pt_1+width]

#---------------------------
#       SINGLE QFI
#--------------------------- 

sys, cryst = YbMgGaO4(0.3, parameter="J1zz")#args=("anneal")
print_wrapped_intensities1(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

measure = ssf_trace(sys; apply_g = false)
swt = SpinWaveTheory(sys; measure) # SpinWaveTheory
# swt_kpm = SpinWaveTheoryKPM(sys; measure, tol=0.01)
kernel = lorentzian(fwhm=0.1)
path = q_space_path(cryst, qs_1_half_2, 500)
energies=range(0, 1, 300)
res = intensities(swt, path; energies, kernel)
# res_bands = intensities_bands(swt, path)

plot_intensities(res)

nQFI_Crystal = nqfi(res, cryst; prefactor=1/3)
fig, ax = plot_nqfi(path, nQFI_Crystal)

#---------------------------
#     PROBE 1 PARAMETER
#---------------------------

# Change variables to fig, ax, fig2, ax2 if method is "default"

locations, max_nqfis, gnd_states, nqfi_std_devs, dom_vecs = examine_effect(YbMgGaO4, range(0,0.5,50), qs_fast; stochastic_sample=24, model="swt",
 intens_type="intensities", parameter="J1xx", energies=[0,1,100], path_length=50, method="data", verbosity="quiet",
 average=false, args=["domain"])

# Figure plotting and modifications ----------------

# ax.title = L"\text{Max nQFI vs Lower Energy Cutoff}"
# ax.xlabel = L"\text{Lower energy cutoff, meV}"
# xlims!(ax, -0.005, 0.5)
# ylims!(ax, -1, 5)
# vlines!(ax, 0.237, color = (:blue, 0.5), linestyle = :dot)

# fig
# fig2

# Save stuff --------------

# GLMakie.save("fig.png", fig) # If method is "default"
dom_save = []
weight_save = []
for elem in dom_vecs
  push!(dom_save, elem[1])
end
for elem in dom_vecs
  push!(weight_save, elem[2])
end
NAME = "YMGO_Jxx_domains_high_res"
writedlm(PATH*NAME*".txt", max_nqfis) # If method is "data
writedlm(PATH*NAME*"_vecs.txt", dom_save) # If method is "data
writedlm(PATH*NAME*"_weights.txt", weight_save) # If method is "data
# Sims to do: Averages, Jxx, Jyy, Jzz. Find seeds for each domain. Then using fixed seed do each noise again Jxx, Jyy, Jzz.


#---------------------------
#     PROBE 2 PARAMETERS
#--------------------------- 
data = []
for noise in range(0,0.2,5)
  push!(data, examine_effect(XXZ_triangular, range(0.5,1.5,100), qs_fast; stochastic_sample=12, model="swt",
  intens_type="intensities", parameter="D", energies=[0,1,100], graph_detail="default", path_length=30,
  args=(noise, "J1zz"), method=data))
end

begin
  fig3 = Figure()
  ax3 = Axis(
          fig3[1,1],
          title=L"\text{Max QFI vs } D ",
          xlabel=L"D=J_{zz}",
          ylabel="Max nQFI"
      )
  scatter_fig = []
  for points in data
    push!(scatter_fig,scatter!(ax3, points[1].content[1].scene.plots[1][1][]))
    errorbars!(ax3, points[1].content[1].scene.plots[2][1][], whiskerwidth = 10)
  end
  Legend(fig3[1, 2],
  [i for i in scatter_fig],
  ["0.00","0.05","0.10","0.15", "0.20"])
  fig3
end



# Assertion: Points at (j,k,C) for const C have the same intensities 
#---------------------------
#    TESTING INTENSITIES
#--------------------------- 
qs_test = []
for i in [0,1/2,1/3]
  for j in range(-1,1,5)
    push!(qs_test, [i,i,j])
  end
end

sys, cryst = YbMgGaO4(0, parameter="J1zz")

measure = ssf_trace(sys; apply_g = false)
swt = SpinWaveTheory(sys; measure) # SpinWaveTheory
swt_kpm = SpinWaveTheoryKPM(sys; measure, tol=0.01)
kernel = lorentzian(fwhm=0.1)
path = q_space_path(cryst, qs_test, 1000)
energies=range(0, 0.6, 300)
res = intensities(swt_kpm, path; energies, kernel)

plot_intensities(res)
# Assertion is true

#---------------------------
#    SIMULATED ANNEALING
#--------------------------- 

a = 3.403  
c = 25.145 
latvecs = lattice_vectors(a, a, c, 90, 90, 120)  
positions = [[0,0,0]] 
cryst = Crystal(latvecs, positions) # Used to be 
# view_crystal(cryst; ndims=2)
dims=(16,16,1) # Default dims
spininfos = [1 => Moment(; s=1, g=1)] #g=2
sys = System(cryst, spininfos, :dipole; dims)

J1pm = 0.1094
J1pmpm = 0.013
J1zz = 0.1264
J2pm = 0.024
J2zz = 0.0278
J1xx = J1pm + J1pmpm
J1yy = J1pm - J1pmpm
Js = [J1zz, J1xx, J1yy]
J_max = findmax(Js)[1]
J1 = [J1xx 0 0; 0 J1yy 0; 0 0 J1zz]

set_exchange!(sys, J1, Bond(1,1,[1,0,0]))
set_exchange!(sys, [J2pm   0.0    0.0;
                    0.0    J2pm   0.0;
                    0.0    0.0    J2zz], Bond(1,1,[1,2,0]))

# sys = to_inhomogeneous(repeat_periodically(sys, (1, 1, 1)))

# for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
#     # The random noise is scaled to the maximum J so that val=1.0 is noise of order J_max (depending on randn distribution)
#     noise = randn()*J_max
#     J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
#     R = J_rot*inv(J1)
#     J_noise = R*(J1 + [1 0 0; 0 1 0; 0 0 1] .* noise .* val)
#     set_exchange_at!(sys, J_noise, site1, site2; offset)
# end

randomize_spins!(sys)

kT_upper = 10J_max  # where J is the largest term in the Hamiltonian
kT_target = 0.1 # some low but finite temperature
kTs = [kT_target + (kT_upper-kT_target) * 0.9^k  for k in 0:100] # define a temperature schedule from high to low T logarithmically
langevin = Langevin(; damping=0.2, kT=kT_upper)
suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.2473; # put in the output of the step above
Es = anneal!(sys,langevin,kTs,500) #play around with the last number 
# plot(Es) # show the energies in this simulated anneal


minimize_energy!(sys, maxiters=10_000)
plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

