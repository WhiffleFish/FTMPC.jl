@recipe function plot(sim::SimHist)
    layout := (1, 2)
    @series begin
        subplot := 1
        labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
        sim.t[1:size(sim.x,2)], trans_states(flip_z(sim.x))'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t[1:size(sim.u,2)], sim.u'
    end
end

@recipe function plot(sim::ModeChangeSimHist, rv=true)
    layout := (1, 2)
    @series begin
        subplot := 1
        labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
        sim.t[1:size(sim.x,2)], sim.x[1:6, :]'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t[1:size(sim.u,2)], sim.u'
    end
end

@recipe function plot(sim::ModeChangeSimHist)
    layout := (1, 2)
    @series begin
        subplot := 1
        labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
        sim.t[1:size(sim.x,2)], trans_states(flip_z(sim.x))'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t[1:size(sim.u,2)], sim.u'
    end
end


rectangle(w, h, x, y) = Plots.Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

@recipe function plot(hist::ModeChangeSimHist, Δt, side, ground)
    xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    tsim = length(x)
    tvec = 0.0:Δt:tsim*Δt-Δt
    recend = tvec[end] + 1
    barrcolor = "lightgray"

    xlims --> (0.0,tsim*Δt-Δt)
    layout := (3, 1)
    size --> (600,600)
    linewidth --> 2
    #linecolor --> "gray"
    legend --> false
    xguidefontsize --> 15
    yguidefontsize --> 15
    fontfamily --> "Computer Modern"


    # Subplot 1
    @series begin
        subplot := 1
        ylims --> [-side,side] + [-1,1]
        tvec, x
    end
    @series begin
        subplot := 1
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-5,0,-side)
    end
    @series begin
        subplot := 1
        opacity --> .5
        color --> barrcolor
        ylabel --> "x(m)"
        rectangle(recend,5,0,side)
    end

    # Subplot 2
    @series begin
        subplot := 2
        ylims --> [-side,side] + [-1,1]
        tvec, y
    end
    @series begin
        subplot := 2
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-5,0,-side)
    end
    @series begin
        subplot := 2
        opacity --> .5
        color --> barrcolor
        ylabel --> "y(m)"
        rectangle(recend,5,0,side)
    end

    # Subplot 3
    @series begin
        subplot := 3
        ylims --> [-ground,ground] + [-1,1]
        xlabel --> "time(s)"
        tvec, z
    end
    @series begin
        subplot := 3
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-5,0,-ground)
    end
    @series begin
        subplot := 3
        opacity --> .5
        color --> barrcolor
        xlabel --> "time(s)"
        ylabel --> "z(m)"
        rectangle(recend,5,0,ground)
    end

end

function crashfilter!(states, ground, side)
    
    crashind = findfirst(x -> any(abs.(x).>[side,side,ground]), [states[i,:] for i in eachindex(states[:,1])])
    if isnothing(crashind)
        return states
    else
        states = states[1:crashind-1,:]
    end
    return states
end

@recipe function plot(hists::Vector{ModeChangeSimHist}, Δt, side, ground, xref, failtime, delaytime)
    x = Vector{Float64}[]
    y = Vector{Float64}[]
    z = Vector{Float64}[]
    tvec = Vector{Float64}[]

    for hist ∈ hists
        xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
        #xyz = crashfilter!(xyz, ground, side)
        xh=xyz[:,1];yh=xyz[:,2];zh=xyz[:,3]
        push!(x,xh) ; push!(y,yh) ; push!(z,zh)
        tsim = length(xh)
        push!(tvec, Δt:Δt:tsim*Δt)
    end

    tsimmax = length(x[end])
    recend = tvec[end][end] + 1
    barrcolor = "lightgray"

    #labels = ["\n  Non-Robust\n" "\n  Unitary-Consensus\n" "\n  Max-Consensus\n" "\n  Reference\n"]
    #labels = ["Non-Robust" "Unitary-Consensus             " "Feasibility-Guided MPC" "" "" "" "" ""]#
    labels = ["" "" "" "Rotor Fail" "IMM Delay" "Reference" "" ""]
    layout := (3,1) #@layout [grid(3, 1) a{0.25w}]#(3,1)
    xlims --> (Δt/2,tsimmax*Δt+0.01)
    size --> (700,600)
    linewidth --> 2
    xguidefontsize --> 15
    yguidefontsize --> 15
    #legend --> false#:outertop
    #legend_columns --> -1 =#
    fontfamily --> "Computer Modern"
    label --> labels

    #= @series begin
        subplot := 4
        seriestype := scatter
        lims --> (0,0.5)
        legendfontsize --> 7
        legend --> :left
        fg_color_legend --> nothing
        label --> labels
        mc --> colors
        frame --> :none
        #rectangle(recend,-2,0,-side)
        (-n:-1)',(-n:-1)'#[1,2,3,4]#,[1,2,3,4] 
    end =#

    # Subplot 1
    for i in eachindex(hists)
        @series begin
            subplot := 1
            ylims --> [-side,side] + [-1,1]
            linewidth --> 4
            tvec[i], x[i]
        end
    end
    @series begin # fail
        subplot := 1
        seriestype:= :vline
        linestyle --> :dashdot
        color --> :black 
        [tvec[1][failtime]]
    end
    @series begin # delay
        subplot := 1
        seriestype:= :vline
        linestyle --> :dashdotdot
        color --> :blue
        [tvec[1][failtime + delaytime]]
    end
    @series begin
        subplot := 1
        seriestype:= :hline
        linestyle --> :dash
        color --> :firebrick
        [xref[1]]
    end
    @series begin
        subplot := 1
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,-side)
    end
    @series begin
        subplot := 1
        opacity --> .5
        color --> barrcolor
        ylabel --> "x(m)"
        ylims --> (-side-1,side+1)
        #legend --> false
        #legend --> :outertop
        legendfontsize --> 10
        #fg_legend --> :transparent
        #legend_columns --> -1
        rectangle(recend,2,0,side)
    end

    # Subplot 2
    for i in eachindex(hists)
        @series begin
            subplot := 2
            ylims --> [-side,side] + [-1,1]
            tvec[i], y[i]
        end
    end
    @series begin # fail
        subplot := 2
        seriestype:= :vline
        linestyle --> :dashdot
        color --> :black 
        [tvec[2][failtime]]
    end
    @series begin # delay
        subplot := 2
        seriestype:= :vline
        linestyle --> :dashdotdot
        color --> :blue
        [tvec[2][failtime + delaytime]]
    end
    @series begin
        subplot := 2
        seriestype:= :hline
        linestyle --> :dash
        color --> :firebrick
        [xref[2]]
    end
    @series begin
        subplot := 2
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,-side)
    end
    @series begin
        subplot := 2
        opacity --> .5
        color --> barrcolor
        ylabel --> "y(m)"
        ylims --> (-side-1,side+1)
        legend --> false
        rectangle(recend,2,0,side)
    end

    # Subplot 3
    for i in eachindex(hists)
        @series begin
            subplot := 3
            ylims --> [-ground,ground] + [-1,1]
            xlabel --> "Time(s)"
            tvec[i], z[i]
        end
    end
    @series begin # fail
        subplot := 3
        seriestype:= :vline
        linestyle --> :dashdot
        color --> :black 
        [tvec[1][failtime]]
    end
    @series begin # delay
        subplot := 3
        seriestype:= :vline
        linestyle --> :dash
        color --> :blue
        [tvec[1][failtime + delaytime]]
    end
    @series begin
        subplot := 3
        seriestype:= :hline
        linestyle --> :dash
        color --> :firebrick
        [-xref[3]]
    end
    @series begin
        subplot := 3
        opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,-ground)
    end
    @series begin
        subplot := 3
        opacity --> .5
        color --> barrcolor
        xlabel --> "Time(s)"
        ylabel --> "z(m)"
        ylims --> (-ground-2,ground+2)
        legend --> false
        rectangle(recend,2,0,ground)
    end
end