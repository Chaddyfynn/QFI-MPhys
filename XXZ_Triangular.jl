using Pkg
Pkg.activate("Inhom")
using Sunny, GLMakie, DelimitedFiles

PATH = "./Figs/Production/XXZ/"

# Include Files
include("QFi.jl")
include("Crystals.jl")

# Simple Crystal and System
sys, cryst = XXZ_triangular(0)
print_irreducible_bz_paths(cryst)

# High Symmetry Points/Paths
M = [1/2, 0, 0]
A = [0, 0, 1/2]
H = [1/3, 1/3, 1/2]
K = [1/3, 1/3, 0]
Γ = [0, 0, 0]
L = [1/2, 0, 1/2]

qs_1 = [Γ, M, K, Γ, A, L, H, A]
qs_2 = [L, M]
qs_3 = [H, K]

qs_neat = [Γ, M, K]
# Small path through max QFI position [-1/3,-1/3,1/3]
width_val = 1/50 # Must be small otherwise increase path_length arg in examine_effect
width = [width_val,width_val,width_val]
qs_fast = [K-width, K, K+width]

#---------------------------
#       SINGLE QFI
#--------------------------- 
sys, cryst = XXZ_triangular(0, parameter="J1xx")
# plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

measure = ssf_trace(sys; apply_g = false)
swt = SpinWaveTheory(sys; measure) # SpinWaveTheory
# kpm = SpinWaveTheoryKPM(sys; measure, tol=0.01)
kernel = lorentzian(fwhm=0.1)
path = q_space_path(cryst, qs_fast, 30)
energies=range(0, 20, 100)
res = intensities(swt, path; energies, kernel)
res_bands = intensities_bands(swt, path)

plot_intensities(res)
# plot_intensities(res_bands)

nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
plot_nqfi(path, nQFI_Crystal)
max_nqfi(nQFI_Crystal, path) # about 8.5 at high res and low res (good)
  # 663 when ultra high res at 0 disorder, 
  # 2.247
path.qs[1056]
# max happens at [-1/3,-1/3,1/3] *approximately

#---------------------------
#     PROBE 1 PARAMETER
#---------------------------
# See QFI.jl (or Hover over examine_effect to see Docstring) for information on how to use examine_effect
locations, max_nqfis, gnd_states, nqfi_std_devs = examine_effect(XXZ_triangular, range(0,1,50), qs_fast; stochastic_sample=24, model="swt",
 intens_type="intensities", parameter="J1zz", energies=[0,10,200], path_length=25, method="data", verbosity="quiet",
  args=["domain"], average=false)

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
NAME = "XXZ_Jzz_domains_high_res"
writedlm(PATH*NAME*".txt", max_nqfis) # If method is "data
writedlm(PATH*NAME*"_vecs.txt", dom_save) # If method is "data
writedlm(PATH*NAME*"_weights.txt", weight_save) # If method is "data

#---------------------------
#     PROBE 2 PARAMETERS
#--------------------------- 
data = []
for noise in range(0,0.2,5)
  push!(data, examine_effect(XXZ_triangular, range(0.5,1.5,100), qs_fast; stochastic_sample=12, model="swt",
  intens_type="intensities", parameter="D", energies=[0,4,100], graph_detail="default", path_length=30,
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

#---------------------------
#       0 NOISE SCALE
#---------------------------
clean_qfis = []
noise_0_qfis = []
Threads.@threads for i in range(0,0,36)
  println("Calculating $i")
  sys, cryst = XXZ_triangular(0, parameter="energy-cutoff") # This disables the sys_inhom section
  measure = ssf_trace(sys; apply_g = false)
  swt = SpinWaveTheory(sys; measure) # SpinWaveTheory
  kernel = lorentzian(fwhm=0.1)
  path = q_space_path(cryst, qs_fast, 25)
  energies=range(0, 20, 200)
  res = intensities(swt, path; energies, kernel)
  nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
  push!(clean_qfis,max_nqfi(nQFI_Crystal, path)[1]) # about 8.5 at high res and low res (good)
end
Threads.@threads for i in range(0,0,36)
  println("Calculating $i")
  sys, cryst = XXZ_triangular(0, parameter="J1zz") # This uses sys_inhom but sets noise to 0
  measure = ssf_trace(sys; apply_g = false)
  swt = SpinWaveTheory(sys; measure) # SpinWaveTheory
  kernel = lorentzian(fwhm=0.1)
  path = q_space_path(cryst, qs_fast, 25)
  energies=range(0, 20, 200)
  res = intensities(swt, path; energies, kernel)
  nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
  push!(noise_0_qfis,max_nqfi(nQFI_Crystal, path)[1]) # about 8.5 at high res and low res (good)
end
begin
  fig = Figure()
  ax = Axis(
    fig[1, 1], 
        title="Clean vs 0 Noise (Random Seed)",
        xlabel="Test",
        xticks = (0:1, ["Clean", "Noise=0"]),
        ylabel="Max nQFI",
  )
  scatter!([0,1],[mean(clean_qfis),mean(noise_0_qfis)])
  errorbars!([0,1], [mean(clean_qfis),mean(noise_0_qfis)], [std(clean_qfis),std(noise_0_qfis)], [std(clean_qfis),std(noise_0_qfis)], whiskerwidth = 10)
  ylims!(ax, 0, findmax([mean(clean_qfis), mean(noise_0_qfis)])[1]*1.4)

  fig
end