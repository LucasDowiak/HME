plot_network <- function(tree, hme_type=c("hme", "hmre"), title)
{
  require(igraph)
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
      return(data.frame(stringsAsFactors=FALSE))
    return(data.frame(from=childs, to=d, stringsAsFactors=FALSE))
  }
  
  p_ <- function(x, y) c(rep(x, ng), rep(y, ne))
  
  nodes <- data.frame(id=c(gatnms, expnms), vertex.type=p_("gate", "expert"),
                      stringsAsFactors=FALSE)
  
  links <- do.call(rbind, lapply(c(gatnms, expnms), match_to_adult))
  
  if (hme_type == "hmre") {
    expidx <- expert_index("hmre", tree)
    for (ee in unique(expidx[expidx > 0])) {
      nms <- names(expidx[expidx==ee])
      nodes[nodes$id %in% nms, "id"] <- paste("*", ee, sep=".")
      links[links$to %in% nms, "to"] <- paste("*", ee, sep=".")
    }
    nodes <- unique(nodes)
    expnms <- nodes[nodes$vertex.type == "expert", "id"]
    ne <- length(expnms)
  }
  net <- graph_from_data_frame(links, vertices=nodes, directed=TRUE)
  plot(net, layout     =layout_(net, as_tree()),
       vertex.label    =plot_labels(c(gatnms, expnms)),
       vertex.shape    =p_("circle", "square"),
       vertex.color    =p_("light blue", "orange"),
       vertex.size     =45,
       edge.arrow.size =0.2,
       main=title
  )
}


plot.hme <- function(obj)
{
  plot_network(obj[["tree"]], obj[["hme.type"]], "Needs a Title")
}

if (FALSE) {
  debugonce(plot_network)
  
  
  treeA <- c("0",
             "0.1", "0.2", "0.3", "0.4")
  
  
  treeB <- c("0",
             "0.1", "0.2",
             "0.1.1", "0.1.2", "0.2.1", "0.2.2")
  
  
  treeC <- c("0",
             "0.1", "0.2",
             "0.1.1", "0.1.2",
             "0.1.1.1", "0.1.1.2")
  
  
  treeD <- c("0",
             "0.1", "0.2", "0.3", "0.4",
             "0.1.1", "0.1.2", "0.2.1", "0.2.2",
             "0.3.1", "0.3.2", "0.4.1", "0.4.2")
  
  
  par(mfrow=c(2,2), mai=c(0.5, 0.3, 0.5, 0.3))
  plot_network(treeA, "hme", "A")
  plot_network(treeB, "hme", "B")
  plot_network(treeC, "hme", "C")
  plot_network(treeD, "hmre", "D")
}


