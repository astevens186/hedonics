
## Preliminaries
rm(list=ls())

# Change working directory to where you've stored ZTRAX
path<- "P:/Peter/Hedonics/Groundwater/"
#install.packages("dplyr", repos = "http://mran.revolutionanalytics.com")
## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c("DiagrammeR")
lapply(packages, pkgTest)

#library(statar)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors

memory.limit(10000000000000)


library(DiagrammeR)


grViz("
	digraph causal {
      
      # Nodes
      node [shape = rectangle,
style = filled,penwidth= 2, fillcolor=gainsboro,
fontcolor = black,fontsize=10]
      S [label = 'Spatial\n Variables']
      P [label = 'Price']
      T [label = 'Treatment']
      X [label = 'Housing\n Characteristics']
      N [label = 'Pre-treatment\n Neighboorhood\n Quality']
      L [label = 'Pre-treatment\n Lagged\n Prices']
      
      # Edges
      edge [color = black,
      arrowhead = vee,penwidth=2,
headclip=true,tailclip=true]
      rankdir = LR
      S->P
      T->P [color = blue]
      X->P
      L->P
      N->T [ color = red]

      edge [ color= black, arrowhead = vee ,penwidth=2,
headclip=true,tailclip=true,style=dashed]
      rankdir = LR
      N->P 
      N->S 
      S->T 
      L->T
      N->L
     


      # Graph
      graph [overlap = true, fontsize = 12]
      }")


#####################

grViz("
      digraph causal {
      
      # Nodes
      node [shape = oval,
      style = filled,penwidth= 1, fillcolor=gainsboro,
      fontcolor = black,fontsize=10]
      P [label = 'Price']
      T [label = 'Treatment']
      X [label = 'Housing\n Characteristics']
      N [label = 'Pre-treatment\n Neighboorhood\n Quality']
      L [label = 'Pre-treatment\n Lagged\n Prices']
      
      # Edges
      edge [color = black,
      arrowhead = vee,penwidth=1,
      headclip=true,tailclip=true]
      rankdir = LR
      T->P [color = blue]
      X->P
      L->P
      
      
      edge [ color= black, arrowhead = vee ,penwidth=1,
      headclip=true,tailclip=true,style=dashed]
      rankdir = LR
      N->P 
      L->T
      N->L
      X->L
      N->X
      X->T
      N->T [ color = red]
      
      # Graph
      graph [overlap = false, fontsize = 10]
      }")


export_svg(gr, file=paste0(path,'dag.svg'))

jpeg(paste0(path,'dag.png'))

mermaid("
        graph LR
        S(Spatial <br> Variables)-->P(Price)
        S-->T(Treatment)
        T-->P
        X(House <br> Characteristics)-->P
        N(Neighboorhood <br> Quality)-->L(Lagged Prices)
        L-->T
        L-->P
        N-->P
        N-->S
        ")
dev.off()
export_graph(m1, file_name = paste0(path,'dag.png'), file_type = NULL, title = NULL,
             width = 6, height = 6)

