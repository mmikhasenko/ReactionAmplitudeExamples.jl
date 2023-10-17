using Plots
pgfplotsx()
plotly()
# pyplot()
gr()


using LaTeXStrings

const sth = 0.2
const m², g² = 1.8, 0.4

# f(s) = -sqrt(-(s-sth))/2
f(s) = 0.3 / (m²-s+g²*sqrt(-(s-sth)))
f_II(s) = 0.3/(m²-s-g²*sqrt(-(s-sth)))
# 
p(x,y,f) = imag(f(x+1im*y))


xmin, xmax = -0.5,3.5
ymin, ymax = -2,2
# 
x    = range(xmin, xmax, length=60)
y    = range(0.001, ymax, length=30)
y = [(-y[end:-1:1])..., y...]
# 
z    = p.(x',y,f)
z_II = p.(x',y,f_II)
# 

let
    plot(leg=:topright)
    ind = y .> 0
    nind = iszero.(ind)
    surface!(x,y[ind],z[ind,:], lab="", zlims=(-2,2), c=:RdBu_8, clim=(-2,2))
    surface!(x,y[nind],z[nind,:], lab="", zlims=(-2,2), c=:RdBu_8, clim=(-2,2))
    plot!(camera=(40,60))
end
let
    plot(leg=:topright)
    ind = y .> 0
    surface!(x,y[ind],z_II[ind,:], lab="", zlims=(-2,2), c=:RdBu_8, clim=(-2,2))
    nind = iszero.(ind)
    surface!(x,y[nind],z_II[nind,:], lab="", zlims=(-2,2), c=:RdBu_8, clim=(-2,2))
    plot!(camera=(40,60))
end



const sth2=2.2
g(s)     = - sqrt(-(s-sth))/2 - sqrt(-(s-sth2))/2
g_II(s)  = + sqrt(-(s-sth))/2 - sqrt(-(s-sth2))/2
g_III(s) = + sqrt(-(s-sth))/2 + sqrt(-(s-sth2))/2
g_IV(s)  = - sqrt(-(s-sth))/2 + sqrt(-(s-sth2))/2


t     = p.(x',y,g)
t_II  = p.(x',y,g_II)
t_III = p.(x',y,g_III)
t_IV  = p.(x',y,g_IV)



let
    plot(leg=:topright)
    ind = y .> 0
    nind = iszero.(ind)
    surface!(x,y[ind],t[ind,:], lab="", zlims=(-5,5), c=:green, clim=(-2,2))
    surface!(x,y[nind],t[nind,:], lab="", zlims=(-5,5), c=:green, clim=(-2,2))
    plot!(camera=(40,10))
end
let
    plot(leg=:topright)
    ind = y .> 0
    surface!(x,y[ind],t_II[ind,:], lab="", zlims=(-5,5), c=:blue, clim=(-2,2))
    nind = iszero.(ind)
    surface!(x,y[nind],t_II[nind,:], lab="", zlims=(-5,5), c=:blue, clim=(-2,2))
    plot!(camera=(40,60))
end
let
    plot(leg=:topright)
    ind = y .> 0
    surface!(x,y[ind],t_III[ind,:], lab="", zlims=(-5,5), c=:yellow, clim=(-2,2))
    nind = iszero.(ind)
    surface!(x,y[nind],t_III[nind,:], lab="", zlims=(-5,5), c=:yellow, clim=(-2,2))
    plot!(camera=(40,60))
end

let
    plot(leg=:topright)
    ind = y .> 0
    surface!(x,y[ind],t_IV[ind,:], lab="", zlims=(-5,5), c=:red, clim=(-2,2))
    nind = iszero.(ind)
    surface!(x,y[nind],t_IV[nind,:], lab="", zlims=(-5,5), c=:red, clim=(-2,2))
    plot!(camera=(40,60))
end

let
    plot(leg=:topright)
    ind = y .> 0
    nind = iszero.(ind)
    surface!(x,y[nind],t[nind,:], lab="", zlims=(-2.5,2.5), c=:green, clim=(-2,2))
    surface!(x,y[ind],t_II[ind,:], lab="", zlims=(-2.5,2.5), c=:blue, clim=(-2,2))
    surface!(x,y[ind],t_III[ind,:], lab="", zlims=(-2.5,2.5), c=:yellow, clim=(-2,2))
    surface!(x,y[ind],t_IV[ind,:], lab="", zlims=(-2.5,2.5), c=:red, clim=(-2,2))
    surface!(x,y[ind],t[ind,:], lab="", zlims=(-2.5,2.5), c=:green, clim=(-2,2))
    plot!(camera=(40,30))
end

let
    xmin, xmax = -0.5,3.5
    ymin, ymax = -2,2
    cmap = :deep
    # 
    x = range(xmin, xmax, length=30)
    y = range(ymin, ymax, length=30)
    z = p.(x',y)
    # 
    zmin, zmax = minimum(z), maximum(z)
    # 
    fig = Figure(resolution = (900,900))
    ax = Axis3(fig, aspect = :data, perspectiveness = 0.5, elevation = π/5,
        xlabel = L"Re $s$", ylabel = L"Im $s$"
        # viewmode = :fit
        # zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey
        )
    # 
    y1 = range(0.001, ymax, length=30)
    z1 = p.(x,y1')
    surface!(ax, x, y1, z1, colormap = :reds, colorrange = (zmin, zmax),
        # lightposition = Vec3f0(0, 0, 0),
        # ambient = Vec3f0(0.65, 0.65, 0.65),
        # backlight = 5f0
        )
    # 
    xm, ym, zm = minimum(ax.finallimits[])
    #
    y2 = range(ymin, -0.001, length=30)
    z2 = p.(x,y2')
    surface!(ax, x, y2, z2, colormap = :greens, colorrange = (zmin, zmax))
    #
    # ax.azimuth[] = -2π/3
    # contour!(ax, x, y1, z1, levels = 20, colormap = cmap, linewidth = 2,
        # colorrange=(zmin, zmax), transformation = (:xy, 2zm))
    # 
    wireframe!(ax, x, y1, z1, overdraw = true, transparency = true,
        color = (:black, 0.1))
    wireframe!(ax, x, y2, z2, overdraw = true, transparency = true,
        color = (:black, 0.1))
    # 
    fig[1,1] = ax
    fig
end


f_II(s) = 1/(m²-s-g²*sqrt(-(s-sth)))
p_II(x,y) = (v=f_II(x+1im*y); abs(imag(v)) > 1.1 ? NaN : imag(v))

let
    xmin, xmax = -0.5,3.5
    ymin, ymax = -2,2
    cmap = :deep
    # 
    x = range(xmin, xmax, length=60)
    y = range(ymin, ymax, length=60)
    z = p_II.(x',y)
    # 
    zmin, zmax = minimum(z), maximum(z)
    # 
    fig = Figure(resolution = (900,900))
    ax = Axis3(fig, aspect = :data, perspectiveness = 0.5, elevation = π/5,
        xlabel = L"Re $s$", ylabel = L"Im $s$"
        )
    # 
    y1 = range(0.001, ymax, length=60)
    z1 = p_II.(x,y1')
    surface!(ax, x, y1, z1, colormap = :reds, colorrange = (zmin, zmax),
        )
    # 
    xm, ym, zm = minimum(ax.finallimits[])
    #
    y2 = range(ymin, -0.001, length=60)
    z2 = p_II.(x,y2')
    surface!(ax, x, y2, z2, colormap = :greens, colorrange = (zmin, zmax))
    # 
    wireframe!(ax, x, y1, z1, overdraw = true, transparency = true,
        color = (:black, 0.1))
    wireframe!(ax, x, y2, z2, overdraw = true, transparency = true,
        color = (:black, 0.1))
    # 
    fig[1,1] = ax
    fig
end



myf(x,y) = (v=1/(x^2+y^2); v>2.1 ? NaN : v)

let
    xmin, xmax = -0.5,3.5
    ymin, ymax = -2,2
    # 
    x = range(xmin, xmax, length=60)
    y = range(ymin, ymax, length=60)
    z = myf.(x,y')
    # 
    zmin, zmax = minimum(z), maximum(z)
    # 
    fig = Figure(resolution = (900,900))
    ax = Axis3(fig, aspect = :data, perspectiveness = 0.5, elevation = π/5,
        xlabel = L"Re $s$", ylabel = L"Im $s$"
        )
    # 
    surface!(ax, x, y, z, colormap = :reds, colorrange = (zmin, zmax))
    fig[1,1] = ax
    fig
end

