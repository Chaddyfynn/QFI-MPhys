using Pkg
Pkg.activate("Inhom")
using Sunny, GLMakie, Statistics

include("dom_wrapped_intensity.jl")

##----------
## nQFI Type 
#-----------
struct nQFI 
    data::Vector{Float64} # nQFI values at each q (s) along QPath (SPath)
end
#-----------------
## nQFI functions
#-----------------

# TODO: Update ħω functions
ħω(res::Sunny.Intensities) = [sum(res.data[:,i])/size(res.data,1) for i in range(1,size(res.data,2))]
ħω(res::Sunny.BandIntensities) = vec(sum(res.disp; dims=1))

## Structure Factors
struct_fac(res::Sunny.Intensities) = [sum(res.data[:,i])/size(res.data,1) for i in range(1,size(res.data,2))] .* (findmax(res.energies)[1] - findmin(res.energies)[1])
struct_fac(res::Sunny.BandIntensities) = vec(sum(res.data, dims=1))


"""Calculates nQFI from the given 'res', at temp 'T', with prefactor 'prefactor'"""
function nqfi(res::Sunny.BandIntensities, cryst::Sunny.Crystal; T::Float64=0.0, prefactor::Float64=1.0, normatm::Bool=true, measure::String="trace")
    qfi = prefactor .* struct_fac(res) .* tanh.(ħω(res)/(2*1.38e-23*T)) .* (1 .- exp.(-ħω(res) ./ (1.38e-23*T)))
    # if measure == "trace"
    #     qfi ./= 3
    # elseif measure == "component"
    #     qfi .*= 1
    # else
    #     print("Unexpected measure. Use either 'trace' or 'component'")
    # end
    normatm ? (return nQFI(qfi ./ Sunny.natoms(cryst))) : (return nQFI(qfi))
end

"""
Calculates the nQFI for an Intensities or BandIntensities Type and given Crystal
Multiplies nQFI by prefactor. Default prefactor=1.0
Normalised by number of atoms in unit cell. To turn off set normatm=false
Returns nQFI Type variable
"""
function nqfi(res::Sunny.Intensities, cryst::Sunny.Crystal; T::Float64=0.0, prefactor::Float64=1.0, normatm::Bool=true, measure::String="trace")
    qfi = prefactor .* struct_fac(res) .* tanh.(ħω(res)/(2*1.38e-23*T)) .* (1 .- exp.(-ħω(res) ./ (1.38e-23*T)))
    # if measure == "trace"
    #     qfi ./= 3
    # elseif measure == "component"
    #     qfi .*= 1
    # else
    #     print("Unexpected measure. Use either 'trace' or 'component'")
    # end
    normatm ? (return nQFI(qfi ./ Sunny.natoms(cryst))) : (return nQFI(qfi))
end

"""Finds maximum value of nQFI. Returns tuple of (max nQFI, Q Point, index)"""
function max_nqfi(nqfi::nQFI, path::Sunny.QPath)
    max_data = findmax(nqfi.data)
    return (max_data[1], path.qs[max_data[2]], max_data[2])
end


##---------
## Plotting
##---------

"""Converts QPath to 1D 'S' Path parameterised by the length along the QPath
Preserves the distance between Q points"""
function spath(path::Sunny.QPath) 
    s = []
    magnitude(Q) = sqrt(Q[1]^2+Q[2]^2+Q[3]^2)
    push!(s, magnitude(path.qs[1])) # Set initial s to initial |Q|
    for i in range(1,size(path.qs,1)-1)
        δ = abs(magnitude(path.qs[i]) - magnitude(path.qs[i+1])) # Distance between Q points
        push!(s, s[i]+δ)
    end

    s_ticks = [s[i] for i ∈ path.xticks[1]] # Get the values of s at the Q tick indices
    return s, s_ticks
end

"""Plots nQFI vs Q Path"""
function plot_nqfi(qpath::Sunny.QPath, nqfi::nQFI; scale::Float64=1.0, xlabels=qpath.xticks[2])
    xs, ticks = spath(qpath) # Parameterise QPath as 1D SPath
    fig = Figure()
    ax = Axis(
    fig[1, 1], 
    xticks = (ticks, xlabels),
    xticklabelrotation = (xlabels == qpath.xticks[2]) ? 45 : 0, # Rotate the labels if they're Q Points not symbols
    xlabel="Q",
    ylabel="nQFI"
    )
    lines!(ax, xs, scale .* nqfi.data)
    fig
    return fig, ax # fig is returned so you can mess with the attributes if you need
end


#--------------------------
# Automation and Cool Stuff
#--------------------------

# ----------------------------------------------------------------------------------------------------
"""
    Examine Effect Function
    Used to probe the effect of different parameters on nQFI
    For a given cryst_func (a function of the parameter to be tested that returns sys and cryst)
    finds the max nQFI given by cryst_func at each val in iterable val_range on the given qpath qs

    Parameters
    cryst_func: Function of two parameters (parameter value, parameter name). returns sys, cryst
    val_range: The values to be tested. Iterable. e.g. range(start,end,num_points)
    qs: The points in the qpath. e.g. [R₂, Γ, T₂] or [[-1/2, -1/2, 1/2], [0, 0, 0], [0, -1/2, 1/2]]
    
    Optional Parameters
    plot_title: Self explanatory. Can be String or latex string
        default - "Max nQFI vs Parameter"
    x_label: As above
        default - "Parameter"
    intens_type: Method of calculating res (intensities). Choose from "intensities" or "intensities_bands"
        default - "intensities"
    model: Choose from "swt" or "kpm". "kpm" tends to be a bit faster.
        default - "swt"
    y_text: The additional y displacement of the Q labels (Float64)
        default - 0.0
    energies: Energies at which to be evaluated. Tuple or Array. (min_energy, max_energy, num_points::Int64). 
        default - [0,1,100]
        caveat - First entry is irrelevant if parameter="energy-cutoff"
    round_to: Number of decimal places at which to round the Q point labels. Int64.
        default - 2
    parameter: The parameter to be probed (String). This depends almost entirely on which parameters are defined in cryst_func
        default - "default" (see cryst_func)
        optional(s) - "energy-cutoff" (The minimum energy in energies)
    stochastic_sample - Number of calculations at each val (to be averaged over) (Int64).
        default - The number of available threads (1 if unconfigured)
    verbosity - The frequency of print and debug statements. Choose from "loud", "quiet", or "none"
    graph_detail - The scale of additional detail on the graphs (e.g. Q Points, red lines, error bars)
        default - "default", single plot with error bars
        optional(s) - "minimal", "default", "all"
    path_length - The length of the QPath. Int64
    method - The return variable, Str
        default - "fig", returns figures of the computations
        optional(s) - "data"
"""
function examine_effect(cryst_func::Function, val_range, qs; plot_title="Max nQFI vs Parameter", x_label="Parameter",
     intens_type::String="intensities", model::String="swt", y_text::Float64=0.0, energies=[0,1,100], round_to::Int64=2,
     parameter::String="default", stochastic_sample::Int64=Threads.nthreads(), verbosity::String="quiet", graph_detail::String="default",
     path_length::Int64=250, method="fig", average=true, args=0)
    if !(verbosity in ["loud", "quiet", "none"])
        print("How can verbosity be $verbosity? Setting it to quiet by default.")
        verbosity = "quiet"
    end

    # If the parameter is anything but the cutoff parameter, then fix energies as normal
    if !(parameter=="energy-cutoff")
        energies = range(energies[1],energies[2],Int(energies[3]))
    end

    # This gets the res, cryst, and path but is its own function so that we can error catch if there is instability
    function get_res_cryst_path(val, iters, ens) # ens is the same as energies. Different name because of nested function
        iters+=1
        if iters > 10
            print("Exceeded maximum iterations. Ending.\n")
            return -1
        end
        try
            sys, cryst = cryst_func(val; parameter=parameter, args=0) # REMOVE ARGS=0 IF ARGS IN CRYST
            measure = ssf_trace(sys; apply_g = false)
            if model == "swt"
                mod = SpinWaveTheory(sys; measure) # SpinWaveTheory (mod=swt)
            elseif model == "kpm"
                mod = SpinWaveTheoryKPM(sys; measure, tol=0.01) 
            else
                print("The model you requested '$model' is not available.")
                return -1
            end
            kernel = lorentzian(fwhm=0.1)
            path = q_space_path(cryst, qs, path_length) # Change number of points in Q path here.
            if intens_type == "intensities"
                res = intensities(mod, path; energies=ens, kernel) # Intensities (try with band_intensities)
            elseif intens_type == "intensities_bands"
                res = intensities_bands(mod, path;)
            else
                print("Intensity type '$intens_type' not recognised!\n")
                return -1
            end
            return res, cryst, path, sys, iters
        catch y
            println(y)
            if verbosity in ["loud", "quiet"]
                print("I ran into an error... But I'll try again. You sit there and look pretty.\n")
            end
            return get_res_cryst_path(val, iters, ens)
        end
    end

    # Empty arrays to store final data
    max_nqfis = []
    locations = [] # Q values at the max nQFI
    gnd_states = []
    nqfi_std_devs = []
    dom_vecs = []

    if verbosity in ["loud", "quiet"]
        nthreads = Threads.nthreads()
        if nthreads > 1
            println("Using $nthreads threads.")
        else
            println("Using 1 thread. You can increase this in your environment settings for better performance.")
        end
    end
    # Test system at each val (will parallel process over JULIA_NUM_THREADS)
    for val ∈ val_range
        if verbosity in ["loud", "quiet"]
            print("Calculaing cryst_funct at $val\n")
        end
        
        # Prepare memory
        inter_maxes=Array{Float64}(undef, stochastic_sample)
        inter_locations=Array{Vector{Float64}}(undef, stochastic_sample)
        inter_gnd_states=Array{Float64}(undef, stochastic_sample)
        inter_dom_vec=Array{Tuple{Vector{Float64}, Float64}}(undef, stochastic_sample)

        # Take n samples where n is stochastic_sample
        Threads.@threads for i =1:stochastic_sample
            if (stochastic_sample != 1) && (verbosity == "loud")
                # If only doing 1 point per val then no need to care about printing
                    print("Calculaing stochastic point $i\n")
            end

            # Try getting the res, cryst, path (will iterate up to 10 times if unstable)
            if !(parameter=="energy-cutoff")
                # Default
                res, cryst, path, sys, iters = get_res_cryst_path(val, 0, energies)
            elseif parameter=="energy-cutoff"
                # Set lower energy cutoff to val in energies
                ens = range(val, energies[2], energies[3])
                res, cryst, path, sys, iters = get_res_cryst_path(val, 0, ens)
            else
                # Failure
                print("'$parameter' is not a valid parameter.")
                return -1
            end
            # Calculate the max nQFI, Q point, and append to intermediate arrays
            nQFI_Crystal = nqfi(res, cryst; prefactor=1.0)
            max_val = max_nqfi(nQFI_Crystal, path)

            inter_maxes[Int(i)] = max_val[1]
            inter_locations[Int(i)] = max_val[2]
            inter_gnd_states[Int(i)] = energy(sys)


            if "domain" in args
                inter_dom_vec[Int(i)] = print_wrapped_intensities1(sys; nmax=10)
            end

        end

        # Once all samples made at val, check that Q points match, average the max nQFI and append'
        for elem in inter_locations
            if (elem != inter_locations[1]) && (verbosity in ["loud"])
                print("Max nQFI is not stable at "*string(elem)*" in the path!\n")
                break
            end
        end

        if stochastic_sample != 1
            push!(nqfi_std_devs, std(inter_maxes))
        end

        if average
            push!(max_nqfis, mean(inter_maxes))
            push!(locations, inter_locations[1])
            push!(gnd_states, mean(inter_gnd_states))
        else
            for i in range(1, stochastic_sample)
                push!(max_nqfis, inter_maxes[i])
                push!(locations, inter_locations[i])
                push!(gnd_states, inter_gnd_states[i])
                push!(dom_vecs, inter_dom_vec[i])
            end
        end
    end

    if method == "data"
        if !("domain" in args)
            return locations, max_nqfis, gnd_states, nqfi_std_devs
        else
            return locations, max_nqfis, gnd_states, nqfi_std_devs, dom_vecs
        end
    else
        # TODO: Combine plots side by side (option)
        # nQFI Plot
        fig = Figure()
        ax = Axis(
        fig[1, 1], 
            title=plot_title,
            xlabel=x_label,
            ylabel="Max nQFI"
        )
        scatter!(ax, val_range, max_nqfis)
        if graph_detail == "all"
            text!(val_range[1], max_nqfis[1] + findmax(max_nqfis)[1]/33; text=string(locations[1]), rotation=3/2 * pi)
            for i in range(2, size(val_range,1))
                if locations[i] != locations[i-1]
                    text!(val_range[i], max_nqfis[i] + findmax(max_nqfis)[1]/33 + y_text; text=string(round.(locations[i], digits = round_to)), rotation=3/2 * pi)
                    vlines!(ax, val_range[i], color = (:red, 0.5), linestyle = :dash)
                end
            end
        end
        if (stochastic_sample != 1) && (graph_detail != "minimal")
            errorbars!(val_range, max_nqfis, nqfi_std_devs, nqfi_std_devs, whiskerwidth = 10)
        end
        ylims!(ax, 0, findmax(max_nqfis)[1]*1.2)

        # Ground States Plot
        fig2 = Figure()
        ax2 = Axis(
            fig2[1,1],
            title="Ground State Energies",
            xlabel=x_label,
            ylabel="Energy, meV"
        )
        ground_states = lines!(ax2, val_range, gnd_states)
        Legend(fig2[1, 2],
        [ground_states],
        ["Ground States"])

        return fig, ax, fig2, ax2
    end
end

function plot_examine_effect(locations, max_nqfis, gnd_states, nqfi_std_devs)
    fig = Figure()
    ax = Axis(
    fig[1, 1], 
        title="Max nQFI vs Variable",
        xlabel="Variable",
        ylabel="Max nQFI"
    )
    scatter!(ax, val_range, max_nqfis)
    if graph_detail == "all"
        text!(val_range[1], max_nqfis[1] + findmax(max_nqfis)[1]/33; text=string(locations[1]), rotation=3/2 * pi)
        for i in range(2, size(val_range,1))
            if locations[i] != locations[i-1]
                text!(val_range[i], max_nqfis[i] + findmax(max_nqfis)[1]/33 + y_text; text=string(round.(locations[i], digits = round_to)), rotation=3/2 * pi)
                vlines!(ax, val_range[i], color = (:red, 0.5), linestyle = :dash)
            end
        end
    end
    if (stochastic_sample != 1) && (graph_detail != "minimal")
        errorbars!(val_range, max_nqfis, nqfi_std_devs, nqfi_std_devs, whiskerwidth = 10)
    end
    ylims!(ax, 0, findmax(max_nqfis)[1]*1.2)

    # Ground States Plot
    fig2 = Figure()
    ax2 = Axis(
        fig2[1,1],
        title="Ground State Energies",
        xlabel="Variable",
        ylabel="Energy, meV"
    )
    ground_states = lines!(ax2, val_range, gnd_states)
    Legend(fig2[1, 2],
    [ground_states],
    ["Ground States"])

    return fig, ax, fig2, ax2
end

