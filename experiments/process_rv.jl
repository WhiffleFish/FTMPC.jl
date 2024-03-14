using DataFrames
using JLD2

hists_unit = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "hists")
hists_consensus = load(joinpath(@__DIR__, "results/rendezvous_threaded_consensus.jld2"), "hists")
simtime = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "simtime")
nvals = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "nvals")
failtimes = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "failtimes")
ndelays = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "ndelays")

df_unit = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_consensus = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])

histcnt = 1
dlytimes = 0:ndelays
for n in nvals
    for ft in failtimes
        for dt in dlytimes
            #println("n: $n, failtime: $ft, delaytime: $dt, histcount: $histcnt")
            push!(df_unit, (n, ft, dt, size(hists_unit[histcnt].x, 2) == simtime))
            push!(df_consensus, (n, ft, dt, size(hists_consensus[histcnt].x, 2) == simtime))
            histcnt += 1
        end
    end
end