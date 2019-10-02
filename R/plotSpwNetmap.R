##' plotSpwNetmap
##'
##'
##' @title Polt a subpathway network map
##' @description Visualize a subpathway network map.
##' @param spwid The subpathway id which the user wants to plot.
##' @param margin A numeric. The value is usually between -0.5 and 0.5, which is able to zoom in or out a subpathway graph. The
##' default is 0.
##' @param vertex.label.cex A numeric vector of node label size.
##' @param vertex.label.font A numeric vector of label font.
##' @param vertex.size A numeric vector of Node size. See \code{\link[igraph]{plot.igraph}}.
##' @param vertex.size2 A numeric vector of Node size.
##' @param edge.arrow.size Edge arrow size.The default is 0.2.
##' @param edge.arrow.width Edge arrow width. The default is 3.
##' @param edge.label.cex Edge label size.
##' @param layout A matrix of x-y coordinates with two dims. Determine the placement
##' of the nodes for drawing a graph.
##' @param vertex.label.color A vector of node label colors. The default is black.
##' @param vertex.color A vector of node colors. The default is the KEGG node color.
##' @param vertex.frame.color A vector of node frame color. The default is dimgray.
##' @param edge.color A vector of edge color. The default is dimgray.
##' @param edge.label.color A vector of edge label color. The default is dimgray.
##' @param sub A character string of subtitle.
##' @param main A character string of main title.
##' @details
##'   The function plotSpwNetmap is able to display a subpathway graph.
##' The argument layout is used to determine the placement of the nodes
##' for drawing a graph.The layouts provided in igraph include `layout_as_star`, `layout_as_tree`,
##' `layout_in_circle`, `layout_nicely`,`layout_on_grid`, `layout_on_sphere`, `layout_randomly`, `layout_with_dh`, `layout_with_fr`,
##' `layout_with_gem`, `layout_with_graphopt`, `layout_with_kk`, `layout_with_lgl`, `layout_with_mds`.
##' The `layout_as_tree` generates a tree-like layout, so it is mainly for
##' trees. The `layout_randomly` places the nodes randomly. The `layout_in_circle` places
##' the nodes on a unit circle. Detailed information on the parameters can be found in \code{\link[igraph]{layout_}}
##' @return a plot
##' @examples
##' # load depend package.
##' library(igraph)
##' # plot network graph of the subpathway 00982_2
##' plotSpwNetmap(spwid="00982_2",layout=layout_nicely)
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
##' @importFrom graphics plot
##' @importFrom igraph  plot.igraph
##' @export
plotSpwNetmap<-function(spwid,layout=NULL,margin = 0, vertex.label.cex = 0.6, vertex.label.font = 1,
                         vertex.size = 8, vertex.size2 = 6, edge.arrow.size = 0.2,
                         edge.arrow.width = 3, edge.label.cex = 0.6,
                         vertex.label.color = "black",
                         vertex.color ="#BFFFBF",vertex.frame.color = "dimgray",
                         edge.color = "dimgray", edge.label.color = "dimgray",sub = NULL,
                         main = NULL){
  spwnetworkdata<-get("spwnetworkdata")
  plot(spwnetworkdata[[spwid]],layout=layout,margin=margin,vertex.label.cex =vertex.label.cex ,vertex.label.font=vertex.label.font,
       vertex.size=vertex.size,vertex.size2=vertex.size2,edge.arrow.size=edge.arrow.size,edge.arrow.width=edge.arrow.width,
       edge.label.cex=edge.label.cex,vertex.label.color=vertex.label.color,vertex.color=vertex.color,
       vertex.frame.color=vertex.frame.color,edge.color=edge.color,edge.label.color=edge.label.color,sub=sub,main=main)
}
