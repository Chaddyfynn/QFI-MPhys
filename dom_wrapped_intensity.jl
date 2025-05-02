using Sunny, GLMakie,Printf, FFTW, LinearAlgebra

import Pkg
Pkg.add("FFTW")

function print_wrapped_intensities1(sys::System{N}; nmax=10) where N
    #sys.crystal == orig_crystal(sys) || error("Cannot perform this analysis on reshaped system.")

    s = reinterpret(reshape, Float64, sys.dipoles)
    V = prod(sys.dims) # number of spins in sublattice
    cell_dims = (2,3,4)
    sk = FFTW.fft(s, cell_dims) / √V
    Sk = real.(conj.(sk) .* sk)
    dims = (1,5) # sum over spin index and atom (sublattice) index
    # In Julia 1.9 this becomes: sum(eachslice(dat; dims)))
    Sk = dropdims(sum(Sk; dims); dims)
    @assert sum(Sk) ≈ norm(sys.dipoles)^2

    weight = Sk .* (100 / norm(sys.dipoles)^2)
    p = sortperm(-weight[:])

    #println("Dominant wavevectors for spin sublattices:\n")
    for (i, m) in enumerate(CartesianIndices(sys.dims)[p])
        k = (Tuple(m) .- 1) ./ sys.dims
        k = [ki > 1/2 ? ki-1 : ki for ki in k]
        kstr = Sunny.fractional_vec3_to_string(k)
        
        if weight[m] < 0.01
            break
        end

        wstr = @sprintf "%6.2f" weight[m]
        nspaces = 20
        spacing = " " ^ max(0, nspaces - length(kstr))
        #print("    $kstr $spacing $wstr%")
        #println(i == 1 ? " weight" : "")
        if i > nmax
            spacing = " " ^ (nspaces-1)
            println("    ... $spacing ...")
            break
        end
        return k, parse(Float64, wstr)
    end
end