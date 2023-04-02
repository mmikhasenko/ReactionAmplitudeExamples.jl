using JSON

function readjson(path)
    f = read(path, String)
    return JSON.parse(f)
end

fin = readjson(joinpath("data", "complexplane.json"))
xi = Vector{Float64}(fin["x"])
yi = Vector{Float64}(fin["y"])
zi = Matrix{Float64}(hcat(fin["z"]...))
zr = Vector{Float64}(fin["z_on_real"])

using PlotlyJS, LaTeXStrings

layout = Layout(
    title="√s complex plane",
    # autosize=true,
    scene_camera_eye=attr(x=0.8, y=-2.0, z=0.6),
    width=600, height=600,
    scene=attr(xaxis_title="Re √s [MeV]",
        yaxis_title="Im √s [MeV]",
        zaxis_title="|A|²",
        #  
        zaxis_range=(2, 15),
        #  xaxis_showbackground=true,
        #  yaxis_showbackground=true,
        # #  
        #  xaxis_gridcolor="white",
        #  yaxis_gridcolor="white",
        #  zaxis_gridcolor="white",
        #  
        yaxis_zerolinewidth=4,
        yaxis_zerolinecolor="red",
        annotations=[
            attr(
                x=0.0,
                y=-0.04,
                z=4.0,
                text="Dˣ⁺D⁰ branch point",
                textangle=0,
                ax=0,
                ay=-75,
                font=attr(
                    color="black",
                    size=12
                ),
                arrowcolor="yellow",
                arrowsize=2,
                arrowwidth=1,
                arrowhead=2),
            attr(
                x=-0.36,
                y=-0.027,
                z=14.0,
                text="Tcc pole",
                textangle=0,
                ax=60,
                ay=0,
                font=attr(
                    color="black",
                    size=14
                ),
                arrowcolor="yellow",
                arrowsize=2,
                arrowwidth=1,
                arrowhead=2)],
    )
)
f = plot(surface(
        x=xi,
        y=yi,
        z=zi',
        showscale=false,
        contours=attr(
            y=attr(show=true, start=0.0, size=0.1, color="red", project_y=true),
            y_end=0.01,
            z=attr(show=true, start=3.0, size=0.5, usecolormap=true, project_z=true),
            z_end=12.5,
        ),
    ), layout)
# savefig(f, "complexplane.html")
