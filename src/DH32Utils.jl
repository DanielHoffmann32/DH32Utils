module DH32Utils

export     block_average_spread,
    draw_mamba_model_graph,
    drawSVG,
    drawSVGandPNG,
    hclustplot,
    loess_fit,
    midpoints,
    save_script_and_fig_to_unique_filenames,
    screenareaPNG

using ImageMagick, FileIO, Loess, DataFrames, Clustering, StatPlots

"""
Takes a Mamba model (see Mamba.jl) and a string (stem name of files), and
draws an acyclic graph corresponding to the model.

NOTE: needs the GraphViz dot program!

Symbols in the graph:

- ellipses are stochastic variables

- squares are logic variables

- white symbols = monitored variables

- gray symbols = variables not monitored

- arrows = "generates"

"""

function draw_mamba_model_graph(model::Any, filestem::String)
    filename=filestem*".dot"
    draw(model,filename=filename)
    
    #we replace the "." in the model name because it causes
    #an error:
    m = replace(readstring(filename),"Mamba.Model","MambaModel")
    write(filename,m)

    #we use the original dot program because it plots with a better
    #layout:
    display("image/svg+xml",open(readstring,`dot $filename -Tsvg`))
end

"""
Plot recipe for cluster dendrogram, from
https://gist.github.com/kescobo/ -> hclustrecipe.jl
"""
function hclustplot(hc::Hclust, useheight::Bool)
    o = indexmap(hc.order)
    n = [x for x in 1:length(o)]

    pos = treepositions(hc, useheight)


    xs = []
    ys = []
    for i in 1: size(hc.merge, 1)
        x1 = pos[hc.merge[i,1]][1]
        x2 = pos[hc.merge[i,2]][1]
        append!(xs, [x1,x1,x2,x2])

        y1 = pos[hc.merge[i,1]][2]
        y2 = pos[hc.merge[i,2]][2]
        useheight ? h = hc.height[i] : h = 1
        newy = maximum([y1,y2]) + h
        append!(ys, [y1,newy,newy,y2])
    end
    return (reshape(xs, 4, size(hc.merge, 1)), reshape(ys, 4, size(hc.merge, 1)))
end

"""
Belongs to plot recipe for cluster dendrogram, from
https://gist.github.com/kescobo/ -> hclustrecipe.jl
"""
function treepositions(hc::Hclust, useheight::Bool)
    order = indexmap(hc.order)
    positions = Dict{}()
    for (k,v) in order
        positions[-k] = (v, 0)
    end
    for i in 1:size(hc.merge,1)
        xpos = mean([positions[hc.merge[i,1]][1], positions[hc.merge[i,2]][1]])
        if hc.merge[i,1] < 0 && hc.merge[i,2] < 0
            useheight ? ypos = hc.height[i] : ypos = 1
        else
            useheight ? h = hc.height[i] : h = 1
            ypos = maximum([positions[hc.merge[i,1]][2], positions[hc.merge[i,2]][2]]) + h
        end

        positions[i] = (xpos, ypos)
    end
    return positions
end

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
Computes block averages of a float vector and determines quantiles (default 0.05 and 0.95 quantiles) of the distribution of these averages for each block size.

The vector elements can be weighted with an additional vector of the same length. In this case the weights in each block are normalized to sum up to one.

Input:

- x: Float64 array

- optionally the following named parameters can be given:

    - dimin (default 10): the minimum block interval (= min. block size)

    - dimax (default 0, will be automatically calculated): the maximum block interval (= max. block size)

    - ddi (default 10): increment in block size

    - w (default [1.] = no weights used): vector of weights, if not given will be set to [1., 1., ...., 1.]

    - qmin (default 0.05): lower quantile

    - qmax (defaults 0.95): upper quantile

Output:

    - qmins: array of lower quantile of block averages

    - qmaxs: array of upper quantile of block averages
"""
function block_average_spread(
            x::Array{Float64,1}; #quantity to be averaged
            dimin::Int64=10, #minimum width of blocks
            dimax::Int64=0, #maximum width of blocks
            ddi::Int64=10, #block increment 
            w::Array{Float64,1}=[1.], #weights
            qmin::Float64=0.05, #minimum quantile to be reported
            qmax::Float64=0.95 #maximum quantile to be reported
        )
    nx = length(x)
    if w != [1.] #weighted
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
        bstart = 1 #index of block start
        bend = di #index of block end
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

"""

Compute midpoints of given range or vector.

"""
function midpoints(r::Range)
    r[1:length(r)-1] + 0.5*step(r)
end

function midpoints(v::AbstractVector)
    [0.5*(v[i] + v[i+1]) for i in 1:length(v)-1]
end

"""
Example:

    For given file "bla.jl" script the call

    save_script_and_fig_to_unique_filenames("bla")

will generate a copy of that file in "bla_RsTzzu.jl" and save the last figure to "bla_RsTzzu.pdf".

Input:

    - source_stem: stem name of source file, string; this should be the name of the existing script file (without extension_script) that will be written in first output

    - extension_script (optional, default ".jl"): ending of script file

    - extension_fig (optional, default ".pdf"): ending and type of figure

    - separator (optional, defaul "_"): string that separates source stem and unique component of file names

Output:

    - writes script file with name source*unique*extension_script in current directory

    - saves last figure to file named source*unique*extension_fig in current directory

"unique" in both outputs is the same unique filename string.

"""

function save_script_and_fig_to_unique_filenames(source_stem::String="";
                                                 extension_script::String=".jl",
                                                 extension_fig::String=".pdf",
                                                 separator::String="_")

    if extension_script == extension_fig
        error("extension_script and extension_fig must be different")
    end

    if source_stem==""
        error("source file stem name required as input")
    end

    if !isfile(source_stem*extension_script)
        error("source file "*source_stem*extension_script*" not present")
    end
    
    written=false
    trials=0
    while (!written && trials<1000)
        #propose unique file names
        proposed_unique=randstring(10)
        outfile_script = source_stem*separator*proposed_unique*extension_script
        outfile_fig = source_stem*separator*proposed_unique*extension_fig

        if !isfile(outfile_script) && !isfile(outfile_fig)
            #no files with the proposed names present
            cp(source_stem*extension_script, outfile_script)
            savefig(outfile_fig)
            written=true
            break
        end
        trials += 1
    end
end

end # module
