using DataFrames
using JLD2

hists = load(joinpath(@__DIR__, "results/rendezvous_threaded.jld2"), "hists")
simtime = load(joinpath(@__DIR__, "results/rendezvous_threaded.jld2"), "simtime")
nvals = load(joinpath(@__DIR__, "results/rendezvous_threaded.jld2"), "nvals")
failtimes = load(joinpath(@__DIR__, "results/rendezvous_threaded.jld2"), "failtimes")
ndelays = load(joinpath(@__DIR__, "results/rendezvous_threaded.jld2"), "ndelays")

# Create the DataFrame
df = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
# Populate the DataFrame
histcnt = 1
for n in nvals
    for ft in failtimes
        dlytimes = ft:ft + ndelays
        for dt in dlytimes
            println("n: $n, failtime: $ft, delaytime: $dt, histcount: $histcnt")
            push!(df, (n, ft, dt, size(hists[histcnt].x, 2) == simtime))
            histcnt += 1
        end
    end
end