using Sunny, GLMakie, Statistics

include("QFI.jl")

a = 3.403  
c = 25.145 

latvecs = lattice_vectors(a, a, c, 90, 90, 120)  

positions = [[0,0,0]] 

cryst = Crystal(latvecs, positions)
view_crystal(cryst; ndims=2)


dims=(14,14,1)
spininfos = [1 => Moment(; s=1, g=1)] #g=2
sys = System(cryst, spininfos, :dipole; dims, seed = 2)

"""J = 1.0
D = 0.5
J1 = [1.0 0 0; 0 1.0 0; 0 0 D]
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))""";

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

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)


randomize_spins!(sys)
minimize_energy!(sys,maxiters = 50000)
plot_spins(sys; ndims=2)

function vacancies_system(sys::Sunny.System, array, add_vac)
    sys_inhom = to_inhomogeneous(sys) #repeat_periodically(sys, (1, 1, 1)))
    for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
        J_rot = Sunny.get_exchange_at(sys_inhom, site1, site2; offset)
        set_exchange_at!(sys_inhom, J_rot, site1, site2; offset)
    end


    if add_vac == true
        for site in eachsite(sys_inhom)
            i = 1 
            i += rand(array)
            if iseven(i)
                set_vacancy_at!(sys_inhom, site)
            end
        end       
    end
    return sys_inhom
end


function vacancy_func(sys_inhom::Sunny.System, path)    
    randomize_spins!(sys_inhom)
    minimize_energy!(sys_inhom, maxiters=10_000)

    measure = ssf_trace(sys_inhom; apply_g = false)
    swt_dis = SpinWaveTheory(sys_inhom; measure)
    res_dis = intensities_bands(swt_dis, path)

    nQFI = nqfi(res_dis; prefactor=1.0, cryst)
    maxQFI = maximum_QFI(nQFI)
    return  maxQFI

end

M = [1/2, 0, 0]
A = [0, 0, 1/2]
H = [1/3, 1/3, 1/2]
K = [1/3, 1/3, 0]
Γ = [0, 0, 0]
L = [1/2, 0, 1/2]

qs = [[0.49,0,0],M, K,[0.34,0.34,0]] 

path = q_space_path(cryst, qs, 30)


scale = 0.0
vacancy_scale = []
array1 = [2,2,2,2,2,2,2,2,2,2,1]
errors = []
qfi_list = []
full_list = []
vac_num = []

add_vac = false
at_half = false

for each in range(0,15)
    some = []
    if each == 1
        add_vac = true
    end     
    sys_inhom = vacancies_system(sys, array1, add_vac)
    for things in range(1, 50)
        a = true
        while a == true
            try
                d = vacancy_func(sys_inhom, path)
                append!(some, d[1])

                a = false
            catch d
               print("Instability error")
            end
        end
    print(1)
    end
    push!(full_list, some)
    error = std(some)
    s = sum(some) ./size(some)
    append!(qfi_list, s[1])
    append!(errors, error)
    if each == 0
        append!(vacancy_scale, 0)
        print(0)
    elseif each < 10
        append!(vacancy_scale, (1 ./(10-(each - 1))))
        print((1 ./(10-(each - 1))))
    else
        append!(vacancy_scale, (1 - 1 ./(size(array1)[1])))
        print((1 - 1 ./(size(array1)[1])))
    end
    
    if size(array1)[1] == 2
        at_half = true
    end
    if at_half == true
        append!(array1,1)
    end
    if at_half == false
        deleteat!(array1,1)
    end
               
end

qfi_list

fig = Figure(size=(900,600))
    ax = Axis(
    fig[1, 1], 
    title = "QFI vs Vacancy Fraction (YbMgGaO4)",
    xlabel="Vacancy Fraction",
    ylabel="QFI"
    )
    scatter!(ax, vacancy_scale, qfi_list, color = :red, label = "Vacancy system")
    errorbars!(ax, vacancy_scale, qfi_list, errors, color = :black)
    lines!(ax, vacancy_scale, qfi_list)
    axislegend()
fig