using DataFrames
using JLD2

hists_unit = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "hists")
hists_consensus = load(joinpath(@__DIR__, "results/hex_threaded_consensus.jld2"), "hists")
simtime = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "simtime")
pos2d = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "pos2d")
failtimes = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "failtimes")
ndelays = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "ndelays")

df_unit = DataFrame(pos = Vector{Float64}[], failtime = Int[], delaytime = Int[], result = Bool[])
df_consensus = DataFrame(pos = Vector{Float64}[], failtime = Int[], delaytime = Int[], result = Bool[])

histcnt = 1
dlytimes = 0:ndelays
for pos in pos2d
    for ft in failtimes
        for dt in dlytimes
            #println("pos: $pos, failtime: $ft, delaytime: $dt, histcount: $histcnt")
            push!(df_unit, (pos, ft, dt, size(hists_unit[histcnt].x, 2) == simtime))
            push!(df_consensus, (pos, ft, dt, size(hists_consensus[histcnt].x, 2) == simtime))
            histcnt += 1
        end
    end
end