using DelimitedFiles # for writing the csv file
using DataStructures # for DefaultDict

L = 2
N = L ^ 2

spins = Iterators.product(fill((-1,1),N)...)
spins = map(x -> reshape(collect(x), (L,L)), spins)[:]

state_summary = []
degeneracies  = DefaultDict(0)

Sigmas = [0, 0, 0, 0]

for spin in spins
    positive_spins = count(i -> i == 1, spin)
    E_h = spin .* circshift(spin, (1,0))
    E_v = spin .* circshift(spin, (0,1))
    E_s = -sum(E_h + E_v)
    M_s = sum(spin)

    Sigmas[1] += E_s
    Sigmas[2] += E_s ^ 2
    Sigmas[3] += abs(M_s)
    Sigmas[4] += M_s ^ 2

    degeneracies[E_s] += 1

    push!(state_summary, [positive_spins, E_s, M_s])
end

for state in state_summary
    push!(state, degeneracies[state[2]])
end

pushfirst!(state_summary, ["1+ spins", "E(s)", "M(s)", "Degeneracies"])
writedlm("output/state_summary.csv",  state_summary, ',')


println("⟨ϵ⟩:   ", Sigmas[1] / N)
println("⟨ϵ²⟩:  ", Sigmas[2] / (N ^ 2))
println("⟨|m|⟩: ", Sigmas[3] / N)
println("⟨m²⟩:  ", Sigmas[4] / (N ^ 2))
