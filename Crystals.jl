using Pkg
Pkg.activate("Inhom")
using Sunny, GLMakie, LinearAlgebra

function YbMgGaO4(val; parameter::String="default", args=0)
    # Hamiltonian
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

    # Unpack args (useful for testing two or more parameters)
    forced_val = 0
        forced_parameter = "no_parameter" 

    if args == "anneal"
        kT_upper = 10J_max  # where J is the largest term in the Hamiltonian
        kT_target = 0.1 # some low but finite temperature
        kTs = [kT_target + (kT_upper-kT_target) * 0.9^k  for k in 0:100] # define a temperature schedule from high to low T logarithmically
        
    elseif args[1] in ["J1zz", "J1yy", "J1xx", "J1yz", "xy_block"]
        forced_val = args[1] # Define args[1] as a forced value of noise
        forced_parameter = args[2] # Define args[2] as the corresponding parameter
    end

    # Choose parameters from "default" (J_ii), "J1zz", "J1yy", "J1xx", "J1yz", "dims"
    # Define lattice constants for YbMgGaO4
    a = 3.403  
    c = 25.145 
    # Configure Crystal
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)  
    positions = [[0,0,0]] 
    cryst = Crystal(latvecs, positions) # Used to be 
    # view_crystal(cryst; ndims=2)
    # Setup Spins
    if parameter == "dims"
        dims=(Int(3*val),Int(3*val),1) # Multiples of three for this symmetry
    else
        dims=(16,16,1) # Default dims
    end
    spininfos = [1 => Moment(; s=1, g=1)] #g=2
    sys = System(cryst, spininfos, :dipole; dims) # REMOVE SEED FOR RANDOM GENERATION ------------------------------------------------------------------------------------

    set_exchange!(sys, J1, Bond(1,1,[1,0,0]))
    set_exchange!(sys, [J2pm   0.0    0.0;
                        0.0    J2pm   0.0;
                        0.0    0.0    J2zz], Bond(1,1,[1,2,0]))

    # Noise
    if !(parameter in ["dims", "energy-cutoff"])
        sys = to_inhomogeneous(repeat_periodically(sys, (1, 1, 1)))
    end

   if (parameter == "default") || (forced_parameter == "default")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            # The random noise is scaled to the maximum J so that val=1.0 is noise of order J_max (depending on randn distribution)
            noise = randn()*J_max
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [1 0 0; 0 1 0; 0 0 1] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1zz") || (forced_parameter == "J1zz")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*J1zz
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 0 0; 0 0 1] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1yy") || (forced_parameter == "J1yy")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*J1yy
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 1 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1xx") || (forced_parameter == "J1xx")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*J1xx
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = *(J1 + [1 0 0; 0 0 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1yz") || (forced_parameter == "J1yz")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 0 1; 0 1 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "xy_block") || (forced_parameter == "xy_block")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [1 1 0; 1 1 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif parameter == "dims"
        0 
    elseif parameter == "energy-cutoff"
        0
    else
        print("The parameter you wish to probe ('$parameter') was not expected.\n")
    end
    # Ground Sate
    randomize_spins!(sys)
    if args == "anneal"
        langevin = Langevin(; damping=0.2, kT=kT_upper)
        # suggest_timestep(sys, langevin; tol=1e-2)
        langevin.dt = 0.2473; # put in the output of the step above
        Es = anneal!(sys,langevin,kTs,500) #play around with the last number 
        plot(Es) # show the energies in this simulated anneal
    end
    minimize_energy!(sys, maxiters=10_000)
    # plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)
    return sys, cryst
end

function XXZ_triangular(val; parameter::String="default", args=0)
    # Unpack args (useful for testing two or more parameters)
    if args == 0
        forced_val = 0
        forced_parameter = "no_parameter" 
    elseif args[2] in ["J1zz", "J1yy", "J1xx", "J1yz", "xy_block"]
        forced_val = args[1] # Define args[1] as a forced value of noise
        forced_parameter = args[2] # Define args[2] as the corresponding parameter
    elseif "anneal" in args
        kT_upper = 10J_max  # where J is the largest term in the Hamiltonian
        kT_target = 0.1 # some low but finite temperature
        kTs = [kT_target + (kT_upper-kT_target) * 0.9^k  for k in 0:100] # define a temperature schedule from high to low T logarithmically
    end
    a = 3.403  
    c = 25.145 
    
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)  
    
    positions = [[0,0,0]] 
    
    cryst = Crystal(latvecs, positions)
    if parameter == "dims"
        dims=(Int(val),Int(val),1) # Multiples of three for this symmetry
    else
        dims=(15,15,1) # Default dims
    end
    spininfos = [1 => Moment(; s=1, g=1)] #g=2
    sys = System(cryst, spininfos, :dipole; dims) # REMOVE SEED FOR RANDOM GENERATION ------------------------------------------------------------------------------------
    
    J = 1.0
    if parameter == "D"
        D = val
    else
        D = 1.5
        # 2.1775 is the value I expect to get max QFI at noise of 0.5
    end

    #------- THIS LINE MARKS THE LAST PLACE WHERE val AND forced_val ARE DISTINCT --------------
    # TODO: Make this a bit more functional rather than procedural
    if args != 0
        val = forced_val
    end

    J1 = [1.0 0 0; 0 1.0 0; 0 0 D]
    set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
    
    if !(parameter in ["dims", "energy-cutoff"]) || (forced_parameter in ["default", "J1zz", "J1yy", "J1xx", "J1yz"])
        sys = to_inhomogeneous(repeat_periodically(sys, (1, 1, 1)))
    end
    if (parameter == "default") || (forced_parameter == "default")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            # The random noise is scaled to the maximum J so that val=1.0 is noise of order J_max (depending on randn distribution)
            noise = randn()*maximum([1,D])
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [1 0 0; 0 1 0; 0 0 1] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1zz") || (forced_parameter == "J1zz")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*D
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 0 0; 0 0 1] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1yy") || (forced_parameter == "J1yy")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 1 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1xx") || (forced_parameter == "J1xx")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = *(J1 + [1 0 0; 0 0 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "J1yz") || (forced_parameter == "J1yz")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [0 0 0; 0 0 1; 0 1 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif  (parameter == "xy_block") || (forced_parameter == "xy_block")
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys, Bond(1,1,[1,0,0]))
            noise = randn()*1
            J_rot = Sunny.get_exchange_at(sys, site1, site2; offset)
            R = J_rot*inv(J1)
            J_noise = R*(J1 + [1 1 0; 1 1 0; 0 0 0] .* noise .* val)
            set_exchange_at!(sys, J_noise, site1, site2; offset)
        end
    elseif parameter in ["dims", "energy-cutoff", "D"]
        0 
    else
        print("The parameter you wish to probe ('$parameter') was not expected.\n")
    end
    
    
    randomize_spins!(sys)
    if "anneal" in args
        langevin = Langevin(; damping=0.2, kT=kT_upper)
        suggest_timestep(sys, langevin; tol=1e-2)
        # langevin.dt = ; # put in the output of the step above
        Es = anneal!(sys,sampler,kTs,500) #play around with the last number 
        plot(Es) # show the energies in this simulated anneal
        # then do minimize_energy!()
        # compare the energy you get with the energies you've been getting by just doing minimize_energy! straight after randomize_spins! 
    end
    minimize_energy!(sys, maxiters=25_000)
    return sys, cryst
end

function anneal!(sys, sampler, kTs,nsweeps)
    Es = zeros(length(kTs))        # Buffer for saving energy as we proceed
    for (i, kT) in enumerate(kTs)
        sampler.kT = kT
        for j ∈ 1:nsweeps
            step!(sys, sampler)
        end                
        Es[i] = energy(sys)   # Query the energy
    end
    return Es    # Return the energy values collected during annealing
end