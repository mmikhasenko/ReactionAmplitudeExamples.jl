
E(r,α,I) = I / (r^2) * cos(α) # 4π

E_of_r(r) = E(r,0,1.0)

using Plots
pyplot()

plot(E_of_r, 1, 3, xlabel="расстояние (м)", ylabel = "сила света (кд)")
E_of_α(α) = E(1,α/180*π,1.0)
plot(E_of_α, -90, 90, xlabel="угол (град.)", ylabel = "сила света (кд)",label="")


#  _|
#  _|    _|_|_|  _|_|_|  _|_|    _|_|_|
#  _|  _|    _|  _|    _|    _|  _|    _|
#  _|  _|    _|  _|    _|    _|  _|    _|
#  _|    _|_|_|  _|    _|    _|  _|_|_|
#                                _|
#                                _|
using Test
using Plots
using Interact

#
r1(x,y,x0,y0,sp) = sqrt(sp.h^2+(x-x0)^2+(y-y0)^2)
r2(x,y,i,sp) = r1(x,y,
    0.5*sp.a+sp.R*sin((i-1)*(2π/sp.N)),
    0.5*sp.b-sp.R*cos((i-1)*(2π/sp.N)),sp)
cosα(x,y,i,sp) =sp.h/r2(x,y,i,sp)
#
O(x,y,i,sp) = sp.I/r2(x,y,i,sp)^2*cosα(x,y,i,sp)
Otot(x,y,sp) = sum(O(x,y,i,sp) for i in 1:sp.N)
#
#
#
setup = (a = 4.0, b = 2.0, h = 0.5, R = 1.0, # meters
         N = 5, # amount
         I = 1.0) # kd
#
@test r1(0,0,0,0,setup) ≈ setup.h
@test r1(0,0,setup.a,setup.b,setup) ≈ sqrt(setup.h^2+setup.a^2+setup.b^2)
#

# static
let
    xv = range(0, setup.a, length=100)
    yv = range(0, setup.b, length=100)
    Ov = [Otot(x,y,(a = 4.0, b = 2.0, h = 1.0,
                    R = 1.0, N = 5, I = 1.0)) for y in yv, x in xv]
    heatmap(xv, yv, Ov, c=:viridis)
end

# dynamic
@manipulate for h=range(0,3, length=100)
    xv = range(0, setup.a, length=100)
    yv = range(0, setup.b, length=100)
    Ov = [Otot(x,y,(a = 4.0, b = 2.0, h = h,
                    R = 1.0, N = 8, I = 1.0)) for y in yv, x in xv]
    heatmap(xv, yv, Ov, clim=(0.05,2.0), c=:viridis)
end
