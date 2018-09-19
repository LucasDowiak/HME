plot_network <- function(tree, title)
{
  expnms <- tree[unlist(is_terminal(tree, tree))]
  gatnms <- setdiff(tree, expnms)
  ng <- length(gatnms)
  ne <- length(expnms)
  
  plot_labels <- function(x)
  {
    f_ <-function(x)
    {
      paste(rev(x), collapse="|")
    }
    unlist(lapply(strsplit(x, "\\."), f_))
  }
  
  match_to_adult <- function(d)
  {
    childs <- unlist(parent(d, tree))
    if (is(childs, "NULL"))
      return(data.frame())
    return(data.frame(to=childs, from=d))
  }
  
  p_ <- function(x, y) c(rep(x, ng), rep(y, ne))
  nodes <- data.frame(
    id=c(gatnms, expnms),
    vertex.type   =p_("gate", "expert"),
    vertex.col    =p_("light blue", "orange"),
    vertex.shape  =p_("circle", "square"),
    stringsAsFactors=FALSE
  )
  
  links <- do.call(rbind, lapply(c(gatnms, expnms), match_to_adult))
  net <- graph_from_data_frame(links, vertices=nodes, directed=TRUE)
  plot(net, layout     =layout_(net, as_tree()),
       vertex.label    =plot_labels(c(gatnms, expnms)),
       vertex.shape    =p_("circle", "square"),
       vertex.color    =p_("light blue", "orange"),
       vertex.size     =25,
       edge.arrow.size =0.3,
       main=title
  )
}
debugonce(plot_network)
plot_network(tree2, "D")

