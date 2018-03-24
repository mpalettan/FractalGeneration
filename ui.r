library(shiny)

fractals <- c("Vicsek", "Sierpinski", "Rug (square pattern)", 
              "Landing at LaGuardia", "2 crosses based fractal",
              "Checkerboard", "Chessboard", "Triangular sibling",
              "Hexagon (destorted)", "Rings", "Barnsley Fern",
              "Crystal", "Leaf", "Tree 1", "Tree 2", "Tower", "Dragon",
              "Forest", "Batman", "Coral")

shinyUI(pageWithSidebar(
    headerPanel("Fractal generation"),
    sidebarPanel(
        h4("1) Select the kind of fractal you want to generate"),
        selectInput("fractal", NULL, choices = fractals),
        h4("2) Indicates the order of the fractal (from 1 to 5)"),
        numericInput("order", NULL, value = 4, min = 1, max = 5, step = 1),
        h4("3) Select the desired color for the fractal"),
        radioButtons("color", NULL, 
                     c("Red 1" = "red",
                       "Red 2" = "darkred",
                       "Blue" = "blue",
                       "Navy" = "navy",
                       "Green 1" = "green",
                       "Green 2" = "dark green",
                       "Cyan" = "cyan",
                       "Magenta" = "darkmagenta",
                       "Orange 1" = "orange1",
                       "Orange 2" = "orange2",
                       "Maroon" = "maroon",
                       "Brown" = "brown"))
    ),
    mainPanel(
        h4('Fractal generated:'),
        verbatimTextOutput("fractalInfo"),
        plotOutput('fractal'),
        h5('Notes:'),
        h6('   1) Some fractals require more time than others to be generated.'),
        h6('   2) Some fractals look better in certain colors. For example, the leaf in green'),
        h6('   3) The order does not apply to all fractals. The higher the order, the longer it takes the fractal to generate.')
    )
))
