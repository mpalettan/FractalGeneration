Fractal Generation
========================================================
author: Mauricio Paletta
date: March 25, 2018
autosize: true

What are Fractals?
========================================================

Infinitely complex patterns that are self-similar across different scales. They are created by repeating a simple process over and over in an ongoing feedback loop.
For more information please visit:

<https://fractalfoundation.org/resources/what-are-fractals/>

<https://en.wikipedia.org/wiki/Fractal>

R package
========================================================

The R code used in this project was written by Voevudko A. E., Ph.D.
For more information please visit:

<https://www.codeproject.com/Articles/1195034/A-Few-Approaches-to-Generating-Fractal-Images-in-R>

How to use the Fractal Generator?
========================================================

The fractal is generated automatically using the parameters selected in the panel on the left. There are three parameters:

- The name of a known fractal. 
- The order or dimension of the fractal. The higher the order, the longer it takes the fractal to generate.
- The color with which you want to appear the fractal. Some fractals look better in certain colors.

Example
========================================================




```r
gpKronFractal(matrix(c(0,1,0,1,1,1,0,1,0), 
       ncol=3, nrow=3, byrow=TRUE), 4, 
       "VicsekFractal1", "red", "Vicsek Fractal")
```

<img src="FractalGeneration-figure/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
NULL
```
