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

@recipe function plot(sim::ModeChangeSimHist, rv::Bool)
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
        xyz = crashfilter!(xyz, ground, side)
        xh=xyz[:,1];yh=xyz[:,2];zh=xyz[:,3]
        push!(x,xh) ; push!(y,yh) ; push!(z,zh)
        tsim = length(xh)
        push!(tvec, 0.0:Δt:tsim*Δt-Δt)
    end

    tmax_ind = findmax([length(x[i]) for i in eachindex(hists)])[2]
    tmax = tvec[tmax_ind]
    tsimmax = length(x[tmax_ind])
    recend = tmax[end] + 1
    barrcolor = colorant"#FFE6E6"

    #labels = ["\n  Non-Robust\n" "\n  Unitary-Consensus\n" "\n  Max-Consensus\n" "\n  Reference\n"]
    #labels = ["Non-Robust" "Unitary-Consensus             " "Feasibility-Guided MPC" "" "" "" "" ""]#
    labels = ["" "" "" "Rotor Fail" "IMM Delay" "Reference" "" ""]
    layout := (3,1) #@layout [grid(3, 1) a{0.25w}]#(3,1)
    xlims --> (0.0,tsimmax*Δt-Δt)
    size --> (700,600)
    #linewidth -->2
    xguidefontsize --> 15
    yguidefontsize --> 15
    legend --> false#:outertop
    #legend_columns --> -1 =#
    fontfamily --> "Computer Modern"
    label --> labels

    hist_colors = [:cyan3, :maroon3, :blue, :gold2]
    num_plots = 3
    
    # Subplot 1
    for i in eachindex(hists)
        feasvec = [info.feas for info ∈ hists[i].info]
        infeas_ind = findfirst(x -> x==false, feasvec)
        if !isnothing(infeas_ind)
            @series begin
                subplot := 1
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                opacity --> .9
                color --> hist_colors[i]
                tvec[i][1:infeas_ind], x[i][1:infeas_ind]
            end
            @series begin
                subplot := 1
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                opacity --> .2
                color --> hist_colors[i]
                tvec[i][infeas_ind:end], x[i][infeas_ind:end]
            end
        else
            @series begin
                subplot := 1
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                color --> hist_colors[i]
                tvec[i], x[i]
            end
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :vline
            linestyle --> :dashdot
            color --> :black 
            [tmax[failtime]]
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :vline
            linestyle --> :dashdotdot
            color --> :blue
            [tmax[failtime + delaytime]]
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :hline
            linestyle --> :dash
            color --> :firebrick
            [xref[i]]
        end
    end
    @series begin
        subplot := 1
        color --> barrcolor
        rectangle(recend,-2,0,-side)
    end
    @series begin
        subplot := 1
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
        feasvec = [info.feas for info ∈ hists[i].info]
        infeas_ind = findfirst(x -> x==false, feasvec)
        if !isnothing(infeas_ind)
            @series begin
                subplot := 2
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                opacity --> .9
                color --> hist_colors[i]
                tvec[i][1:infeas_ind], y[i][1:infeas_ind]
            end
            @series begin
                subplot := 2
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                opacity --> .3
                color --> hist_colors[i]
                tvec[i][infeas_ind:end], y[i][infeas_ind:end]
            end
        else
            @series begin
                subplot := 2
                ylims --> [-side,side] + [-1,1]
                linewidth --> 4
                color --> hist_colors[i]
                tvec[i], y[i]
            end
        end
    end
    @series begin
        subplot := 2
        color --> barrcolor
        rectangle(recend,-2,0,-side)
    end
    @series begin
        subplot := 2
        color --> barrcolor
        ylabel --> "y(m)"
        ylims --> (-side-1,side+1)
        legend --> false
        rectangle(recend,2,0,side)
    end

    # Subplot 3
    for i in eachindex(hists)
        feasvec = [info.feas for info ∈ hists[i].info]
        infeas_ind = findfirst(x -> x==false, feasvec)
        if !isnothing(infeas_ind)
            @series begin
                subplot := 3
                ylims --> [-ground,ground] + [-1,1]
                linewidth --> 4
                opacity --> .9
                color --> hist_colors[i]
                tvec[i][1:infeas_ind], z[i][1:infeas_ind]
            end
            @series begin
                subplot := 3
                ylims --> [-ground,ground] + [-1,1]
                linewidth --> 4
                opacity --> .2
                color --> hist_colors[i]
                tvec[i][infeas_ind:end], z[i][infeas_ind:end]
            end
        else
            @series begin
                subplot := 3
                ylims --> [-ground,ground] + [-1,1]
                linewidth --> 4
                xlabel --> "Time(s)"
                color --> hist_colors[i]
                tvec[i], z[i]
            end
        end
    end
    @series begin
        subplot := 3
        color --> barrcolor
        rectangle(recend,-2,0,-ground)
    end
    @series begin
        subplot := 3
        color --> barrcolor
        xlabel --> "Time(s)"
        ylabel --> "z(m)"
        ylims --> (-ground-2,ground+2)
        legend --> false
        rectangle(recend,2,0,ground)
    end
end

@recipe function plot(hists::Vector{ModeChangeSimHist}, Δt, barriers, xref, failtime, delaytime)
    x = Vector{Float64}[]
    y = Vector{Float64}[]
    z = Vector{Float64}[]
    tvec = Vector{Float64}[]

    side, ylower, yupper, ground = barriers[1], barriers[2], barriers[3], barriers[4]

    for hist ∈ hists
        xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
        #xyz = crashfilter!(xyz, ground, side)
        xh=xyz[:,1];yh=xyz[:,2];zh=xyz[:,3]
        push!(x,xh) ; push!(y,yh) ; push!(z,zh)
        tsim = length(xh)
        push!(tvec, 0.0:Δt:tsim*Δt-Δt)
    end
    tmax_ind = findmax([length(x[i]) for i in eachindex(hists)])[2]
    tmax = tvec[tmax_ind]
    tsimmax = length(x[tmax_ind])
    recend = tmax[end] + 1
    barrcolor = colorant"#FFE6E6"#RGB(255.0, 230.0, 230.0)#"lightsalmon1"#"lightgray"

    #labels = ["\n  Non-Robust\n" "\n  Unitary-Consensus\n" "\n  Max-Consensus\n" "\n  Reference\n"]
    #labels = ["Non-Robust" "Unitary-Consensus             " "Feasibility-Guided MPC" "" "" "" "" ""]#
    labels = ["" "" "" "Rotor Fail" "IMM Delay" "Reference" "" ""]
    layout := (3,1) #@layout [grid(3, 1) a{0.25w}]#(3,1)
    xlims --> (0.0,tsimmax*Δt-Δt)
    size --> (700,600)
    #linewidth --> 2
    xguidefontsize --> 15
    yguidefontsize --> 15
    legend --> false#:outertop
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
    hist_colors = [:cyan3, :maroon3, :blue, :gold2]
    # hist_colors = Dict(
    #     "nr" => :red,
    #     "first" => :maroon3,
    #     "full" => :gold2,
    #     "ours" => :cyan3
    # )

    num_plots = 3
    # Subplot 1
    for i in eachindex(hists)
        @series begin
            subplot := 1
            ylims --> [-side,side] + [-1,1]
            linewidth --> 4
            color --> hist_colors[i]
            tvec[i], x[i]
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :vline
            linestyle --> :dashdot
            color --> :black 
            [tmax[failtime]]
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :vline
            linestyle --> :dashdotdot
            color --> :blue
            [tmax[failtime + delaytime]]
        end
    end
    for i in 1:num_plots
        @series begin
            subplot := i
            seriestype:= :hline
            linestyle --> :dash
            color --> :firebrick
            [xref[i]]
        end
    end
    @series begin
        subplot := 1
        #opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,-side)
    end
    @series begin
        subplot := 1
        #opacity --> .5
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
            ylims --> [ylower,yupper] + [-1,1]
            linewidth --> 4
            color --> hist_colors[i]
            tvec[i], y[i]
        end
    end
    @series begin
        subplot := 2
        #opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,ylower)
    end
    @series begin
        subplot := 2
        #opacity --> .5
        color --> barrcolor
        ylabel --> "y(m)"
        ylims --> (ylower-1,yupper+1)
        legend --> false
        rectangle(recend,2,0,yupper)
    end

    # Subplot 3
    for i in eachindex(hists)
        @series begin
            subplot := 3
            ylims --> [-ground,ground] + [-1,1]
            xlabel --> "Time(s)"
            linewidth --> 4
            color --> hist_colors[i]
            tvec[i], z[i]
        end
    end
    @series begin
        subplot := 3
        #opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,-ground)
    end
    @series begin
        subplot := 3
        #opacity --> .5
        color --> barrcolor
        xlabel --> "Time(s)"
        ylabel --> "z(m)"
        ylims --> (-ground-2,ground+2)
        legend --> false
        rectangle(recend,2,0,ground)
    end
end

function rendezvous_crashfilter!(hists, tvec, ylower, yupper)
    for i ∈ eachindex(hists)
        crashind = findfirst(x -> any((x<ylower) || (x>yupper)), hists[i])
        if !isnothing(crashind)
            hists[i] = hists[i][1:crashind]
            tvec[i] = tvec[i][1:crashind]
        end
    end
    return hists, tvec
end

@recipe function plot(hists::Vector{ModeChangeSimHist}, params)
    Δt, barriers, xref, failtime, delaytime = params[1], params[2], params[3], params[4], params[5]
    x = Vector{Float64}[]
    y = Vector{Float64}[]
    z = Vector{Float64}[]
    tvec = Vector{Float64}[]

    ylower, yupper = barriers[1], barriers[2]

    for hist ∈ hists
        xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
        #xyz = crashfilter!(xyz, ground, side)
        xh=xyz[:,1];yh=xyz[:,2];zh=xyz[:,3]
        push!(x,xh) ; push!(y,yh) ; push!(z,zh)
        tsim = length(xh)
        push!(tvec, 0.0:Δt:tsim*Δt-Δt)
    end
    tmax_ind = findmax([length(x[i]) for i in eachindex(hists)])[2]
    tmax = tvec[tmax_ind]
    tsimmax = length(x[tmax_ind])
    recend = tmax[end] + 1
    barrcolor = colorant"#FFE6E6"#RGB(255.0, 230.0, 230.0)#"lightsalmon1"#"lightgray"

    labels = ["" "" "" "Rotor Fail" "IMM Delay" "Reference" "" ""]
    #layout := (3,1) #@layout [grid(3, 1) a{0.25w}]#(3,1)
    xlims --> (0.0,tsimmax*Δt-Δt)
    size --> (600,400)
    #linewidth --> 2
    xguidefontsize --> 15
    yguidefontsize --> 15
    legend --> false#:outertop
    #legend_columns --> -1 =#
    fontfamily --> "Computer Modern"
    label --> labels
    xlabel --> "Time(s)"
    
    hist_colors = [:cyan3, :maroon3, :blue, :gold2]
    #hist_colors = [:cyan3, :maroon3, :blue, :brown1]

    @series begin
        seriestype:= :vline
        linestyle --> :dashdot
        color --> :black 
        [tmax[failtime]]
    end
    @series begin
        seriestype:= :vline
        linestyle --> :dashdotdot
        color --> :blue
        [tmax[failtime + delaytime]]
    end
    @series begin
        seriestype:= :hline
        linestyle --> :dash
        color --> :firebrick
        [xref[2]]
    end

    y, tvec = rendezvous_crashfilter!(y, tvec, ylower, yupper)

    for i ∈ eachindex(hists)
        feasvec = [info.feas for info ∈ hists[i].info]
        infeas_ind = findfirst(x -> x==false, feasvec)
        if !isnothing(infeas_ind)
            @series begin
                ylims --> [ylower,yupper] + [-1,1]
                linewidth --> 4
                opacity --> .9
                color --> hist_colors[i]
                tvec[i][1:infeas_ind], y[i][1:infeas_ind]
            end
            @series begin
                ylims --> [ylower,yupper] + [-1,1]
                linewidth --> 4
                color --> hist_colors[i]
                #linestyle --> :dash
                opacity --> .2
                tvec[i][infeas_ind:end], y[i][infeas_ind:end]
            end
        else
            @series begin
                ylims --> [ylower,yupper] + [-1,1]
                linewidth --> 4
                color --> hist_colors[i]
                tvec[i], y[i]
            end
        end
    end
    @series begin
        #opacity --> .5
        color --> barrcolor
        rectangle(recend,-2,0,ylower)
    end
    @series begin
        #opacity --> .5
        color --> barrcolor
        ylabel --> "y(km)"
        ylims --> (ylower-1,yupper+1)
        legend --> false
        rectangle(recend,2,0,yupper)
    end

end