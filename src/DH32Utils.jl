module DH32Utils

export drawSVG,
       drawSVGandPNG,
       screenareaPNG

using Images

function drawSVG(figstem::ASCIIString; edit::Bool=false)
    
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

function drawSVGandPNG(figstem::ASCIIString; edit::Bool=false)
    
    if (stat("$figstem.svg").size == 0) #generate new figure
        cp(joinpath(Pkg.dir("DH32Utils","data"),"template-drawing.svg"),"$figstem.svg")
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg`)
    elseif (edit == true) #edit existing figure
        run(`inkscape $figstem.svg --with-gui -D -l $figstem.svg -e $figstem.png`)
    end

    run(`inkscape $figstem.svg -e $figstem.png`)

    #show figure
    imread("$figstem.png")
    
end


function screenareaPNG(figstem::ASCIIString; edit::Bool=false)

    if (stat("$figstem.png").size == 0 || edit == true)
        run(`scrot "$figstem.png" -s`) #generate new screenshot
    end

    #show image
    Images.imread("$figstem.png")
    
end
    
end # module
