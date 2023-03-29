xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
#plot([hist.x'[:,[1,2]] -hist.x'[:,3]], label = ["x" "y" "z"])
x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
zlow = extrema(z)[1]

# Define the trajectory

pos = Observable(Point3f[(0,0,0)])

fig = Figure(); display(fig)
ax = Axis(fig[1,1])
limits!(ax, -1, 8, -1, 1, -1, 1)

cm.scatter!(ax, xyz)
frames = 1:simT

record(fig, "anim.mp4", frames; framerate = 10) do frame
    newp = Point3f(xyz[frame,:])
    pos[] = push!(pos[], newp)
end

