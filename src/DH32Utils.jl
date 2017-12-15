module DH32Utils

export drawSVG,
    drawSVGandPNG,
    screenareaPNG,
    loess_fit,
    block_average_spread

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

"""
Computes block averages of a float vector and determines quantiles (default 0.05 and 0.95 quantilees) of the distribution of these averages for each block size.

The vector elements can be weighted with an additional vector of the same length. In this case the weights in each block are normalized to sum up to one.

Input:

- x: Float64 array

- optionally the following named parameters can be given:

    - dimin (default 10): the minimum block interval (= min. block size)

    - dimax (default 0, will be automatically calculated): the maximum block interval (= max. block size)

    - ddi (default 10): increment in block size

    - weighted (default false): should we weight each vector element

    - w: vector of weights, if not given will be set to [1., 1., ...., 1.]

    - qmin (default 0.05): lower quantile

    - qmax (defaults 0.95): upper quantile

"""
function block_average_spread(
            x::Array{Float64,1}; #quantity to be averaged
            dimin::Int64=10, #minimum width of blocks
            dimax::Int64=0, #maximum width of blocks
            ddi::Int64=10, #block increment 
            weighted::Bool=false, #shall we weight
            w::Array{Float64,1}=[1.], #weights
            qmin::Float64=0.05, #minimum quantile to be reported
            qmax::Float64=0.95 #maximum quantile to be reported
        )
    nx = length(x)
    if weighted == true
        if length(w) != nx
            error("error block_averages: length(w)!=length(x)")
        end
        if sum(w .> 0.) != nx
            error("error block_averages: some weights not positive")
        end
    else
        w = ones(nx)
    end
    if dimax == 0 #the last block average is between first and second half of data
        dimax = div(nx,2)
    end
    if dimax < dimin
        error("error block_averages: dimax < dimin")
    end
    
    ba = zeros(div(nx,dimin)) #array for analysis of block averages
    
    nout = length(collect(dimin:ddi:dimax))
    qmins = zeros(nout) #output of lower quantiles of block averages 
    qmaxs = zeros(nout) #output of upper quantiles of bl. av.
    
    jout=1
    for di in dimin:ddi:dimax #width of blocks
        ni = div(nx,di) # number of blocks
        bstart = 1
        bend = di
        for i in 1:ni 
            ba[i]=dot(x[bstart:bend],w[bstart:bend])/sum(w[bstart:bend])
            bstart += di
            bend += di
        end
        qmins[jout] = quantile(ba[1:ni],qmin)
        qmaxs[jout] = quantile(ba[1:ni],qmax)
        jout += 1
    end
    qmins, qmaxs
end

end # module
