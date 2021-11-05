using DelimitedFiles # for writing the csv file
using DataStructures # for DefaultDict

N = 2

spins = Iterators.product(fill((-1,1),N^2)...)
spins = map(x -> reshape(collect(x), (2,2)), spins)[:]

state_summary = []
degeneracies  = DefaultDict(0)

for spin in spins
    positive_spins = count(i -> i == 1, spin)
    E_h = spin .* circshift(spin, (1,0))
    E_v = spin .* circshift(spin, (0,1))
    E_s = -sum(E_h + E_v)
    M_s = sum(spin)

    degeneracies[E_s] += 1

    push!(state_summary, [positive_spins, E_s, M_s])
end

for state in state_summary
    push!(state, degeneracies[state[2]])
end

pushfirst!(state_summary, ["1+ spins", "E(s)", "M(s)", "Degeneracies"])
writedlm("output/state_summary.csv",  state_summary, ',')
