module DH32Utils

export drawSVG,
    drawSVGandPNG,
    screenareaPNG,
    loess_fit

using ImageMagick, FileIO, Loess

function drawSVG(figstem::String; edit::Bool=false)
    
    if (stat("$figstem.svg").size == 0) #generate new figure
        cp(joinpath(Pkg.dir("DH32Utils","data"),"template-drawing.svg"),"$figstem.svg")
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg`)
    elseif (edit == true) #edit existing figure
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg`)
    end

    #show figure
    graph=open(readall,"$figstem.svg") 
    display("image/svg+xml",graph)
    
end

function drawSVGandPNG(figstem::String; edit::Bool=false)
    
    if (stat("$figstem.svg").size == 0) #generate new figure
        cp(joinpath(Pkg.dir("DH32Utils","data"),"template-drawing.svg"),"$figstem.svg")
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg`)
    elseif (edit == true) #edit existing figure
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg -e $figstem.png`)
    end

    run(`inkscape $figstem.svg -e $figstem.png`)

    #show figure
    load("$figstem.png")
    
end

function screenareaPNG(figstem::String; edit::Bool=false)

    if (stat("$figstem.png").size == 0 || edit == true)
        run(`scrot "$figstem.png" -s`) #generate new screenshot
    end

    #show image
    load("$figstem.png")
    
end

"""
Produces a smooth fit (Loess) for a set of x and y values.

Input:

    - x : input x values

    - y : input y values

    - nsteps (default=100) : number of steps between min. and max. of x

    - normalize (default: true) [see Loess package]

    - span (default: 0.75) [see Loess package]

    - degree (default: 2) [see Loess package]

Output:

    - x steps

    - fitted y values at stepped x positions
"""
function loess_fit{T <: AbstractFloat}(
    x::AbstractVector{T},
    y::AbstractVector{T};
    nsteps::Int64=100,
    normalize::Bool=true, span::T=0.75, degree::Int=2
)
    if (nsteps < 2)
        error("loess_fit: need at least 2 steps")
    end
    if (length(x) != length(y))
        error("loess_fit: x and y input arrays must have same lengths")
    end
    maxx = maximum(x)
    minx = minimum(x)
    if (minx >= maxx)
        error("loess_fit: maximum x value not greater than minimum x value")
    end
    delta = (maxx-minx)/nsteps
    xsteps = collect(minx:delta:maxx)

    model = loess(x, y; normalize=normalize, span=span, degree=degree)
    ymodel = Loess.predict(model, xsteps)

    xsteps, ymodel
end

end # module
