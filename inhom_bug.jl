using Pkg
Pkg.activate("Inhom")
using Sunny, GLMakie,LinearAlgebra 

function set_up_sys(dims;Jxx=1,Jyy=1,Jzz=1)
    a = 3  
    c = 90
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)  
    positions = [[0,0,0]] 
    cryst = Crystal(latvecs, positions)
    spininfos = [1 => Moment(; s=1, g=1)] 
    sys = System(cryst, spininfos, :dipole; dims,seed=2 )
    set_exchange!(sys, [Jxx   0.0    0.0;
                        0.0    Jyy   0.0;
                        0.0    0.0   Jzz], Bond(1,1,[1,0,0]))
    randomize_spins!(sys)
    minimize_energy!(sys;maxiters=2_000)
    return sys
end

dims = (8,8,1)
Jxx = 1
Jyys =  0:0.1:1
Jzzs = 0:0.1:1

energies_clean = zeros(Float64,length(Jyys),length(Jzzs))
energies_dis = zeros(Float64,length(Jyys),length(Jzzs))
energies_diff = zeros(Float64,length(Jyys),length(Jzzs))

#Jyys=[0.5]
#Jzzs=[0.25]

for jjz in 1:length(Jzzs)
    for jjy in 1:length(Jyys)
        Jxx = 1
        Jyy = Jyys[jjy]
        Jzz = Jzzs[jjz]
        sys = set_up_sys(dims;Jxx=Jxx,Jyy,Jzz)
        sys_inhom = to_inhomogeneous(sys)
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
            Jmat = Sunny.get_exchange_at(sys_inhom, site1, site2; offset)
            set_exchange_at!(sys_inhom, Jmat, site1, site2; offset)
        end
        randomize_spins!(sys_inhom)
        minimize_energy!(sys_inhom)
        randomize_spins!(sys)
        minimize_energy!(sys)
        E_dis=energy_per_site(sys_inhom)
        E_clean= energy_per_site(sys)
        energies_clean[jjy,jjz] = E_clean
        energies_dis[jjy,jjz] = E_dis
        energies_diff[jjy,jjz] = E_clean-E_dis
    end
end

begin
    fig = Figure()
    crange = (-2.25,-0.75)
    crange2 = (0,0.001)
    ax1 = Axis(fig[1,1];title="Clean Energies",xlabel="Jyy",ylabel="Jzz")
    ax2 = Axis(fig[1,2];title="Disordered Energies",xlabel="Jyy",ylabel="Jzz")
    ax3 = Axis(fig[1,4];title="|Difference|",xlabel="Jyy",ylabel="Jzz")
    heatmap!(ax1,Jyys,Jzzs,energies_clean;colorrange=crange)
    heatmap!(ax2,Jyys,Jzzs,energies_dis;colorrange=crange)
    Colorbar(fig[1,3];colorrange=crange)
    heatmap!(ax3,Jyys,Jzzs,norm.(energies_diff);colorrange=crange2)
    Colorbar(fig[1,5];colorrange=crange2)
    fig
end
