using DataFrames
using JLD2
using Statistics

function process_costs(hist, x_ref)
    costs = zeros(size(hist.x)[2])
    for i in 1:size(hist.x)[2]
        costs[i] = norm(hist.x[:,i]-x_ref, 2)
    end
    return costs
end

function process_costs(hists::AbstractArray, x_ref)
    println("Processing costs")
    costs = zeros(length(hists))
    for i in 1:length(hists)
        costs[i] = last(process_costs(hists[i], x_ref))
    end
    return costs
end

hists_unit = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "hists")
hists_max = load(joinpath(@__DIR__, "results/rendezvous_threaded_max.jld2"), "hists")
hists_nonrobust = load(joinpath(@__DIR__, "results/rendezvous_threaded_nonrobust.jld2"), "hists")
hists_consensus = load(joinpath(@__DIR__, "results/rendezvous_threaded_consensus.jld2"), "hists")

simtime = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "simtime")
nvals = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "nvals")
failtimes = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "failtimes")
ndelays = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "ndelays")
x_ref = load(joinpath(@__DIR__, "results/rendezvous_threaded_unit.jld2"), "x_ref")

df_unit = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_max = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_nonrobust = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])
df_consensus = DataFrame(nval = Float64[], failtime = Int[], delaytime = Int[], result = Bool[])

histcnt = 1
dlytimes = 0:ndelays
for n in nvals
    for ft in failtimes
        for dt in dlytimes
            #println("n: $n, failtime: $ft, delaytime: $dt, histcount: $histcnt")
            if isassigned(hists_unit, histcnt)
                push!(df_unit, (n, ft, dt, size(hists_unit[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            if isassigned(hists_max, histcnt)
                push!(df_max, (n, ft, dt, size(hists_max[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            if isassigned(hists_nonrobust, histcnt)
                push!(df_nonrobust, (n, ft, dt, size(hists_nonrobust[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
                deleteat!(hists_nonrobust, histcnt) 
            end
            if isassigned(hists_consensus, histcnt)
                push!(df_consensus, (n, ft, dt, size(hists_consensus[histcnt].x, 2) == simtime))
            else
                @warn "histcnt: $histcnt is not assigned"
            end
            histcnt += 1
        end
    end
end
println("TOTAL TRIALS: ", length(df_unit[:,end]))
println("-------------------------------------------------------")
println("Unit planner successes: ", sum(df_unit[:,end])," | success (%): ", sum(df_unit[:,end])/length(df_unit[:,end]), " | final cost avg: ", mean(process_costs(hists_unit, x_ref)))
println("Max planner successes: ", sum(df_max[:,end])," | success (%): ", sum(df_max[:,end])/length(df_max[:,end]), " | final cost avg: ", mean(process_costs(hists_max, x_ref)))
println("Nonrobust planner successes: ", sum(df_nonrobust[:,end])," | success (%): ", sum(df_nonrobust[:,end])/length(df_nonrobust[:,end]), " | final cost avg: ", mean(process_costs(hists_nonrobust, x_ref)))
println("Consensus planner successes: ", sum(df_consensus[:,end])," | success (%): ", sum(df_consensus[:,end])/length(df_consensus[:,end]), " | final cost avg: ", mean(process_costs(hists_consensus, x_ref)))