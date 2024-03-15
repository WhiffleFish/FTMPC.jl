using DataFrames
using JLD2

hists_unit = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "hists")
hists_max = load(joinpath(@__DIR__, "results/hex_threaded_max.jld2"), "hists")
hists_nonrobust = load(joinpath(@__DIR__, "results/hex_threaded_nonrobust.jld2"), "hists")
hists_consensus = load(joinpath(@__DIR__, "results/hex_threaded_consensus.jld2"), "hists")

simtime = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "simtime")
pos2d = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "pos2d")
failtimes = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "failtimes")
ndelays = load(joinpath(@__DIR__, "results/hex_threaded_unit.jld2"), "ndelays")

df_unit = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_max = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_nonrobust = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_consensus = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])

histcnt = 1
dlytimes = 0:ndelays
for pos in pos2d
    for ft in failtimes
        for dt in dlytimes
            #println("n: $n, failtime: $ft, delaytime: $dt, histcount: $histcnt")
            if isassigned(hists_unit, histcnt)
                push!(df_unit, (pos, ft, dt, size(hists_unit[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            if isassigned(hists_max, histcnt)
                push!(df_max, (pos, ft, dt, size(hists_max[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            if isassigned(hists_nonrobust, histcnt)
                push!(df_nonrobust, (pos, ft, dt, size(hists_nonrobust[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            if isassigned(hists_consensus, histcnt)
                push!(df_consensus, (pos, ft, dt, size(hists_consensus[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            histcnt += 1
        end
    end
end

println("Unit planner success (%): ", sum(df_unit[:,end])/length(df_unit[:,end]))
println("Max planner success (%): ", sum(df_max[:,end])/length(df_max[:,end]))
println("Nonrobust planner success (%): ", sum(df_nonrobust[:,end])/length(df_nonrobust[:,end]))
println("Consensus planner success (%): ", sum(df_consensus[:,end])/length(df_consensus[:,end]))
