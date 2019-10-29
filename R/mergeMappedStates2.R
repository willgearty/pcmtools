#' Rename, merge, or split mapped states
#'
#' This function renames, merges, or splits mapped states on a tree.
#'
#' If \code{node} is specified, only mapped states of descendants edges of that node are modified. This function is based on the similarly named utility function by Liam Revell in \code{phytools}.
#'
#' @param tree an object of class "simmap" or "multiSimmap" containing one or more phylogenetic trees with a mapped discrete character
#' @param old.states state(s) to rename or merge
#' @param new.state name for new state
#' @param node a node index, used to specify a subclade of interest
#' @return An object of class "simmap" or "multiSimmap".
#' @export
#' @examples
#' library(phytools)
#' # simulate a mapped tree
#' set.seed(4)
#' Q <- matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
#' rownames(Q) <- colnames(Q) <- letters[1:3]
#' tree <- sim.history(pbtree(n=100,scale=1),Q)
#' cols <- setNames(c("blue","red","green","orange"),letters[1:4])
#'
#' # plot the mapping
#' plot(tree, cols, ftype="i", fsize=0.7)
#'
#' # split state c to state d within subclade
#' tree2 <- mergeMappedStates2(tree, "c", "d", 173)
#'
#' # plot the new mapping
#' plot(tree2, cols, ftype="i", fsize=0.7)
mergeMappedStates2<-function(tree,old.states,new.state,node=NULL){
  if(inherits(tree,"multiSimmap")){
    tree<-unclass(tree)
    tree<-lapply(tree,mergeMappedStates2,old.states=old.states,new.state=new.state,node=node)
    class(tree)<-c("multiSimmap","multiPhylo")
  } else if(inherits(tree,"simmap")) {
    maps <- tree$maps

    #if node is supplied, we only want to modify descendants of that node
    if(!is.null(node)){
      repl <- which(tree$edge[,2] %in% getDescendants(tree, node))
    } else {
      repl <- 1:Nedge(tree)
    }

    rr<-function(map,oo,nn){
      for(i in 1:length(map)) if(names(map)[i] %in% oo) names(map)[i]<-nn
      map
    }

    mm<-function(map){
      if(length(map)>1){
        new.map<-vector()
        j<-1
        new.map[j]<-map[1]
        names(new.map)[j]<-names(map)[1]
        for(i in 2:length(map)){
          if(names(map)[i]==names(map)[i-1]){
            new.map[j]<-map[i]+new.map[j]
            names(new.map)[j]<-names(map)[i]
          } else {
            j<-j+1
            new.map[j]<-map[i]
            names(new.map)[j]<-names(map)[i]
          }
        }
        map<-new.map
      }
      map
    }

    #rename old states to new state
    maps[repl]<-lapply(maps[repl],rr,oo=old.states,nn=new.state)
    #join any adjacent mapped elements in the same state
    maps[repl]<-lapply(maps[repl],mm)
    mapped.edge<-tree$mapped.edge
    #move old state values to new state in mapped.edge
    if (!(new.state %in% colnames(mapped.edge))) {
      mapped.edge<-cbind(mapped.edge,0)
      colnames(mapped.edge)[ncol(mapped.edge)]<-new.state
    }
    if (length(old.states) > 1) old.val<-mapped.edge[repl,old.states] else old.val<-matrix(mapped.edge[repl,old.states])
    mapped.edge[repl,new.state]<-rowSums(old.val)
    #set old state values to 0
    mapped.edge[repl,old.states]<-0
    #remove unused states
    mapped.edge <- mapped.edge[,!(colSums(mapped.edge) == 0)]

    tree$maps<-maps
    tree$mapped.edge<-mapped.edge
  } else stop("tree should be an object of class \"simmap\" or \"multiSimmap\".")
  return(tree)
}
