library(shiny)

## R Helper Functions (for load before plotting) ## 6/5/17
## Load sample. e.g.: source("C:\\RData\\ALLRHFs.R")

# NOTE: 1. All generating/plotting functions (but not pIfsFractal()) are using
#          matrix mat to fill it with 0/nonzero integer and plot it.
#       2. plotv2() plotting helper functions is using 2 vectors X,Y
#          from the dump file.
#       3. The file names used are without extension (which will be added
#          as ".png", ".dmp" and ".dat" when needed).
#       4. Requesting dump file is useful if the generating/plotting time
#          is big. Having a dump file makes it easy and fast to repeat plotting with
#          different colors, titles and picture sizes.

## HFR#1 plotmat(): Simple plotting using matrix mat (filled with 0/nonzero int).
#  Where: mat - matrix; fn - file name (no extension); clr - color;
#         ttl - plot title; dflg - writing dump file flag (0-no/1-yes):
#         psz - picture size; cx - cex (scale).
#
plotmat <- function(mat, fn, clr, ttl, dflg=0, psz=600, cx=1.0) {
    m = nrow(mat); d = 0; X=NULL; Y=NULL;
    pf = paste0(fn, ".png"); df = paste0(fn, ".dmp");
    
    # Building X and Y arrays for plotting from not equal to zero values in mat.
    #
    for (i in 1:m) {
        for (j in 1:m) {if(mat[i,j]==0){next} else {d=d+1; X[d]=i; Y[d]=j} }
    };

    # Dumping if requested (dflg=1).
    #
    if (dflg==1) {
        dump(c("X","Y"), df); # cat(" *** Dump file:", df, "\n")
    };
    
    # Plotting
    #
    if (ttl!="") { 
        mp <- plot(X,Y, main=ttl, axes=FALSE, xlab="", ylab="", col=clr, pch=20, cex=cx)}
    else {
        par(mar=c(0,0,0,0));
        mp <- plot(X,Y, axes=FALSE, xlab=NULL, ylab=NULL, col=clr, pch=20, cex=cx)
    };
    return (mp);
}

## HFR#2 plotv2(): Simple plotting using 2 vectors (dumped into ".dmp" file).
# Where: fn - file name; clr - color; ttl - plot title; psz - picture size;
#         cx - cex or scale.
#
plotv2 <- function(fn, clr, ttl, psz=600, cx=1.0) {
    pf = paste0(fn, ".png"); df = paste0(fn, ".dmp");
    source(df); d = length(X);

    # Plotting
    #
    if (ttl!="") {
        mp <- plot(X,Y, main=ttl, axes=FALSE, xlab="", ylab="", col=clr, pch=20, cex=cx)}
    else {
        par(mar=c(0,0,0,0));
        mp <- plot(X,Y, axes=FALSE, xlab=NULL, ylab=NULL, col=clr, pch=20, cex=cx)
    };
    return (mp);
}
#=================================================

## Kronecker power of a matrix. 
## Where: m - initial matrix, n - power.
#
matkronpow <- function(m, n) {
    if (n<2) {return (m)};
    r = m; n = n-1; 
    for(i in 1:n) {r = r%x%m};
    return (r);
}

## Generate and plot Kronecker product based fractals. 8/12/16
## gpKronFractal(m, n, pf, clr, ttl, dflg, psz, cx):
## Where: m - initial matrix (filled with 0/int); n - order of the fractal;
## fn - plot file name (without extension); clr - color; ttl - plot title;
## dflg - writing dump file flag (0/1); psz - picture size; cx - cex
#
gpKronFractal <- function(m, n, fn, clr, ttl, dflg=0, psz=640, cx=1.0) {
    fign="Kpbf";
    if(fn=="") {fn=paste0(fign,"o", n)} else {fn=paste0(fn)};
    if(ttl!="") {ttl=paste0(ttl,", order ", n)};
    r = matkronpow(m, n);
    return (plotmat(r, fn, clr, ttl, dflg, psz, cx));
}

#=================================================

## Plotting fractals using IFS style   7/27/16
## Based on already calculated M x 7 table of coefficients in the input file.
## Note: 1. Input ifs-file should be dat-file; output is set as png-file.
##       2. Swap 2nd and 3rd column if you've got data used in Java, JavaScript, etc.

## pIfsFractal(fn,n,clr,ttl,cx): Plot fractals using IFS style.
## Where: fn - file name; n - number of dots; clr - color; ttl - plot title,
##        psz - plot size, cx - cex.

pIfsFractal <- function(mat, n, clr, ttl, psz=600, cx=0.5) {
    # Reading a complete data table from the file: space delimited, no header.
    # Table has any number of rows, but always 7 columns is a must.
    #
    Tb = mat
    tr = nrow(Tb)
    
    # Creating matrix M1 from 1st 4 columns of each row.
    #
    M1 = vector("list",tr);
    for (i in 1:tr) {M1[[i]] = matrix(c(Tb[i,1:4]),nrow=2)}
    
    # Creating matrix M2 from columns 5,6 of each row.
    #
    M2 = vector("list",tr);
    for (i in 1:tr) {M2[[i]] = matrix(c(Tb[i,5:6]),nrow=2)}
    
    ## Creating matrix M3 (actualy a vector) from column 7 of each row.
    #
    M3 = c(Tb[1:tr,7])
    
    x = numeric(n); y = numeric(n);
    x[1] = y[1] = 0;
    
    # Main loop
    #
    for (i in 1:(n-1)) {
        k = sample(1:tr, prob=M3, size=1);
        M = as.matrix(M1[[k]]);
        z = M%*%c(x[i],y[i]) + M2[[k]];
        x[i+1] = z[1];
        y[i+1] = z[2];
    }
    
    # Plotting
    #
    if (ttl!="") {
        mp <- plot(x,y, main=ttl, axes=FALSE, xlab="", ylab="", col=clr, pch=20, cex=cx)}
    else {
        par(mar=c(0,0,0,0));
        mp <- plot(x,y, axes=FALSE, xlab=NULL, ylab=NULL, col=clr, pch=20, cex=cx)
    };
    return (mp);
}
#=================================================

# BTALL.R  ## ALL 4 versions   7/27/16
# translation of PARI/GP: http://rosettacode.org/wiki/Brownian_tree#PARI.2FGP

# Generate and plot Brownian tree. Version #1.
# gpBrownianTree1(m, n, clr, fn, ttl, dflg)
# Where: m - defines matrix m x m; n - limit of the number of moves; 
#   fn - file name (.ext will be added); ttl - plot title; dflg - 0-no dump,
#   1-dump.
#
gpBrownianTree1 <- function(m, n, clr, fn, ttl, dflg=0)
{
    M <- matrix(c(0),ncol=m,nrow=m,byrow=T);
    
    # Seed in center
    #
    x <- m%/%2; y <- m%/%2;
    M[x,y]=1;
    pf=paste0(fn,".png");

    # Main loops: Generating matrix M
    #
    for (i in 1:n) {
        if(i>1) {
            x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)}
        while(1) {
            ox=x; oy=y;
            x <- x + sample(-1:1, 1, replace=F);
            y <- y + sample(-1:1, 1, replace=F);
            if(x<=m && y<=m && x>0 && y>0 && M[x,y]) 
            {if(ox<=m && oy<=m && ox>0 && oy>0){M[ox,oy]=1; break}}
            if(!(x<=m && y<=m && x>0 && y>0)) {break}
        }
    }
    return (plotmat(M, fn, clr, ttl, dflg)); 
}

# Generate and plot Brownian tree. Version #2.
# gpBrownianTree2(m, n, clr, fn, ttl, dflg)
# Where: m - defines matrix m x m; n - limit of the number of moves; 
#   fn - file name (.ext will be added); ttl - plot title; dflg - 0-no dump,
#   1-dump; seed - 0-center, 1-random.
#
gpBrownianTree2 <- function(m, n, clr, fn, ttl, dflg=0)
{
    M <- matrix(c(0),ncol=m,nrow=m,byrow=T);
    
    # Random seed
    #
    x <- sample(1:m, 1, replace=F); y <- sample(1:m, 1, replace=F); 
    M[x,y]=1;
    pf=paste0(fn,".png");

    # Main loops: Generating matrix M
    #
    for (i in 1:n) {
        if(i>1) {
            x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)}
        while(1) {
            dx <- sample(-1:1, 1, replace=F);
            dy <- sample(-1:1, 1, replace=F);
            nx=x+dx; ny=y+dy;
            if(!(nx<=m && ny<=m && nx>0 && ny>0)) {
                x <- sample(1:m, 1, replace=F); y <- sample(1:m, 1, replace=F)}
            else {if(M[nx,ny]) {M[x,y]=1; break}
                else{x=nx; y=ny;}}
        }
    }
    return (plotmat(M, fn, clr, ttl, dflg));
}

# Generate and plot Brownian tree. Version #3. 
# gpBrownianTree3(m, n, clr, fn, ttl, dflg, seed):
# Where: m - defines matrix m x m; n - limit of the number of moves; 
#   fn - file name (.ext will be added); ttl - plot title; dflg - 0-no dump,
#   1-dump; seed - 0-center, 1-random.
#
gpBrownianTree3 <- function(m, n, clr, fn, ttl, dflg=0, seed=0)
{
    M <- matrix(c(0),ncol=m,nrow=m,byrow=T);
    
    # Random seed
    #
    if(seed==1) {
        x <- sample(1:m, 1, replace=F);y <- sample(1:m, 1, replace=F)
    } 
    
    # Seed in center
    #
    else {
        x <- m%/%2; y <- m%/%2
    }
    M[x,y]=1;
    pf=paste0(fn,".png");

    # Main loops: Generating matrix M
    #
    for (i in 1:n) {
        if(i>1) {
            x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)}
        b <- 0
        while(b==0) {
            dx <- sample(-1:1, 1, replace=F)
            dy <- sample(-1:1, 1, replace=F)
            if(!(x+dx<=m && y+dy<=m && x+dx>0 && y+dy>0))
            { x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)
            }
            else{if(M[x+dx,y+dy]==1) {M[x,y]=1; b=1}
                else {x=x+dx; y=y+dy;} } 
        }
    }
    return (plotmat(M, fn, clr, ttl, dflg)); 
}

# Generate and plot Brownian tree. Version #4.
# gpBrownianTree4(m, n, clr, fn, ttl, dflg, seed)
# Where: m - defines matrix m x m; n - limit of the number of moves; 
#   fn - file name (.ext will be added); ttl - plot title; dflg - 0-no dump,
#   1-dump; seed - 0-center, 1-random.
#
gpBrownianTree4 <- function(m, n, clr, fn, ttl, dflg=0, seed=0)
{
    M <- matrix(c(0),ncol=m,nrow=m,byrow=T);
    
    # Random seed
    #
    if(seed==1) {
        x <- sample(1:m, 1, replace=F);y <- sample(1:m, 1, replace=F)
    }
    
    # Seed in center
    #
    else {
        x <- m%/%2; y <- m%/%2
    }
    M[x,y]=1;
    pf=paste0(fn,".png");

    # Main loops: Generating matrix M
    #
    for (i in 1:n) {
        if(i>1) {
            x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)}
        while((x<=m && y<=m && x>0 && y>0)) {
            if(!(x+1<=m && y+1<=m && x-1>0 && y-1>0)) {break;}
            b=M[x+1,y+1]+M[x,y+1]+M[x-1,y+1]+M[x+1,y];
            b=b+M[x-1,y-1]+M[x-1,y]+M[x,y-1]+M[x+1,y-1];
            if(b!=0) {break;}
            x <- x + sample(-1:1, 1, replace=F)
            y <- y + sample(-1:1, 1, replace=F)
            if(!(x<=m && y<=m && x>0 && y>0))
            { x <- sample(1:m, 1, replace=F)
            y <- sample(1:m, 1, replace=F)
            }
        }
        M[x,y]=1;
    }
    return (plotmat(M, fn, clr, ttl, dflg)); 
}

fractals <- c("Vicsek", "Sierpinski", "Rug (square pattern)", 
              "Landing at LaGuardia", "2 crosses based fractal",
              "Checkerboard", "Chessboard", "Triangular sibling",
              "Hexagon (destorted)", "Rings", "Barnsley Fern",
              "Crystal", "Leaf", "Tree 1", "Tree 2", "Tower", "Dragon",
              "Forest", "Batman", "Coral")

shinyServer(
    function(input, output) {
        output$fractalInfo <- renderPrint({
            cat(input$fractal, 
                (if (input$fractal %in% fractals[1:10]) c("- order", input$order) else ""), 
                "-", input$color)
        })
        output$fractal <- renderPlot({
            if (input$fractal == fractals[1]) {
                gpKronFractal(
                    matrix(c(0,1,0,1,1,1,0,1,0), ncol=3, nrow=3, byrow=TRUE), 
                    input$order, "VicsekFractal1", input$color, "")
            }
            if (input$fractal == fractals[2]) {
                gpKronFractal(
                    matrix(c(1,1,1,1,0,1,1,1,1), ncol=3, nrow=3, byrow=TRUE), 
                    input$order, "SierpinskiCF1", input$color, "")
            }
            if (input$fractal == fractals[3]) {
                gpKronFractal(
                    matrix(c(1,1,1,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1), 
                           ncol=5, nrow=5, byrow=TRUE), 
                    input$order, "RugF1", input$color, "", 1)
            }
            if (input$fractal == fractals[4]) {
                gpKronFractal(
                    matrix(c(1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1), ncol=4, nrow=4, byrow=TRUE), 
                    input$order, "LaGuardiaFractal", input$color, "", 1, 1024)
            }
            if (input$fractal == fractals[5]) {
                gpKronFractal(
                    matrix(c(0,1,0,1,1,1,0,1,0), ncol=3, nrow=3, byrow=TRUE) %x% 
                    matrix(c(1,0,1,0,1,0,1,0,1), ncol=3, nrow=3, byrow=TRUE), 
                    (if (input$order > 3) 3 else input$order), "CrossesF1", input$color, "", 1)
            }
            if (input$fractal == fractals[6]) {
                gpKronFractal(
                    matrix(c(0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0), ncol=4, nrow=4, byrow=TRUE), 
                    (if (input$order > 3) 3 else input$order), "CkBF", input$color, "")
            }
            if (input$fractal == fractals[7]) {
                gpKronFractal(
                    matrix(c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,
                             1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1),
                           ncol=8, nrow=8, byrow=TRUE), 
                    (if (input$order > 3) 3 else input$order), "CsBF", input$color, "");
            }
            if (input$fractal == fractals[8]) {
                gpKronFractal(
                    matrix(c(1,1,1,0,1,1,0,0,1), ncol=3, nrow=3, byrow=TRUE), 
                    input$order, "TSF1z", input$color, "", 1 , 0.5)
            }
            if (input$fractal == fractals[8]) {
                gpKronFractal(
                    matrix(c(1,1,1,0,1,1,0,0,1), ncol=3, nrow=3, byrow=TRUE), 
                    input$order, "TSF1z", input$color, "", 1 , 0.5)
            }
            if (input$fractal == fractals[9]) {
                gpKronFractal(
                    matrix(c(1,1,0,1,1,1,0,1,1), ncol=3, nrow=3, byrow=TRUE), 
                    input$order, "HGF", input$color, "")
            }
            if (input$fractal == fractals[10]) {
                gpKronFractal(
                    matrix(c(1,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1), ncol=4, nrow=4, byrow=TRUE), 
                    input$order, "RGSF", input$color, "")
            }
            if (input$fractal == fractals[11]) {
                pIfsFractal(
                    matrix(c(0.00, 0.00, 0.00, 0.16, 0.0, 0.00, 0.01,
                             0.85, 0.04, -0.04, 0.85, 0.0, 1.60, 0.85,
                             0.20, -0.26, 0.23, 0.22, 0.0, 1.60, 0.07,
                             -0.15, 0.28, 0.26, 0.24, 0.0, 0.44, 0.07), ncol=7, nrow=4, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[12]) {
                pIfsFractal(
                    matrix(c(0.696970, -0.481061, -0.393939, -0.662879, 2.147003, 10.310288, 0.747826,
                             0.090909, -0.443182, 0.515152, -0.094697, 4.286558, 2.925762, 0.252174), ncol=7, nrow=2, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[13]) {
                pIfsFractal(
                    matrix(c(0.555, 0.000, 0.000, 0.555, 0.000, 0.000, 0.20,
                             0.550, 0.000, 0.000, 0.550, 0.000, 0.185, 0.30,
                             0.353, 0.281, -0.295, 0.336, 0.068, 0.112, 0.25,
                             0.353, -0.281, 0.295, 0.336, -0.068, 0.112, 0.25), ncol=7, nrow=4, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[14]) {
                pIfsFractal(
                    matrix(c(.195, -.488, .344, .443, .4431, .2452, .2,
                             .462, .414, -.252, .361, .2511, .5692, .25,
                             -.058, -.070, .453, -.111, .5976, .0969, .2,
                             -.035, .070, -.469, -.022, .4884, .5069, .25,
                             -.637, 0, 0, .501, .8562, .2513, .1), ncol=7, nrow=5, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[15]) {
                pIfsFractal(
                    matrix(c(0.000, 0.000, 0.000, 0.600, 0.00, -0.065, 0.100,
                             0.440, 0.000, 0.000, 0.550, 0.00, 0.200, 0.180,
                             0.343, -0.248, 0.199, 0.429, -0.03, 0.100, 0.180,
                             0.343, 0.248, -0.199, 0.429, 0.03, 0.100, 0.180,
                             0.280, -0.350, 0.280, 0.350, -0.05, 0.000, 0.180,
                             0.280, 0.350, -0.280, 0.350, 0.05, 0.000, 0.180), ncol=7, nrow=6, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[16]) {
                pIfsFractal(
                    matrix(c(0.75, 0.00, 0.00, 0.30, -0.20, 0.00, 0.23,
                             0.75, 0.00, 0.00, 0.30, 0.20, 0.00, 0.23,
                             0.50, 0.00, 0.00, 0.80, 0.00, 0.20, 0.54), ncol=7, nrow=3, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[17]) {
                pIfsFractal(
                    matrix(c(0.824074, 0.581482, -0.212346, 0.864198, 1.882290, 0.110607, 0.787473,
                             0.088272, 0.420988, -0.463889, -0.377778, 0.785360, 8.095795,  0.212528), ncol=7, nrow=2, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[18]) {
                pIfsFractal(
                    matrix(c(-0.632407, -0.614815, -0.545370, 0.659259, 3.840822, 1.282321, 0.888128,
                             -0.036111, 0.444444, 0.210185, 0.037037, 2.071081, 8.330552, 0.111872), ncol=7, nrow=2, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[19]) {
                pIfsFractal(
                    matrix(c(0.5, 0, 0, 0.5, 0, 0, 0.25,
                             0.5, 0.5, 0, 0.5, 0.5, 0, 0.25,
                             0.5, 0, 0.5, 0.5, 0, 0.5, 0.25,
                             0.5, 0, 0, 0.5, 0.5, 0.5, 0.25), ncol=7, nrow=4, byrow=TRUE), 
                    100000, input$color, "")
            }
            if (input$fractal == fractals[20]) {
                pIfsFractal(
                    matrix(c(-0.16666667, -0.1666667, 0.16666667, -0.1666667, 0.0000000, 0.000000, 0.163,
                             0.83333333, 0.2500000, -0.25000000, 0.8333333, -0.1666667, -0.166667, 0.600,
                             0.33333333, -0.0833333, 0.08333333, 0.3333333, 0.0833333, 0.666667, 0.237), ncol=7, nrow=3, byrow=TRUE), 
                    100000, input$color, "")
            }
        })
    }
)