using Sunny, GLMakie
units = Units(:meV, :angstrom)
latvecs = lattice_vectors(3, 8, 8, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(1, 1, 1))
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
axis = [0, 1, 0]
k = [1/2,0,0]
measure = ssf_custom((q, ssf) -> real(ssf[1,1]), sys)
swt_rot = SpinWaveTheorySpiral(sys; measure, k, axis)
res = intensities(swt_rot, path; energies, kernel=gaussian(fwhm=0.25))
plot_intensities(res;)

res = intensities_bands(swt_rot, path)
plot_intensities(res; units)