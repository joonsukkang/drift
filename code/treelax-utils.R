library(fastTopics)


################################################################################
# functions for preparing data  
################################################################################

# G (S x N): genotype data
# P (S x K): population allele frequencies
# Q (N x K): individuals' population memberships
# X (S x K): minor allele counts
prep_GPQX <- function(G, P, Q){
  
  popsize <- round(colSums(Q)) # effective population size
  X <- round(P %*% diag(2*popsize))
  colnames(X) <- colnames(P)
  
  out.list <- list(G=G, P=P, Q=Q, X=X, popsize=popsize)
  return(out.list)
}





# G (S x N): genotype data
# P (S x K): population allele frequencies
# Q (N x K): individuals' population memberships
# X (S x K): minor allele counts
prep_GPQX_sim <- function(filename, K, pop){
  
  # load genotyope matrix
  G <- read_plink(filename, verbose=FALSE)
  G <- G$X
  
  # load ADMIXTURE fit
  P <- read_table(paste0(filename,'.', K, '.P'),
                  col_names = FALSE, show_col_types = FALSE)
  P <- as.matrix(P)
  dim(P)
  
  Q <- read_table(paste0(filename,'.', K, '.Q'),
                  col_names = FALSE, show_col_types = FALSE)
  Q <- as.matrix(Q)
  
  # order factors
  data.frame(pop = pop,
             Q) %>%
    pivot_longer(cols=2:(ncol(Q)+1)) %>%
    group_by(pop, name) %>%
    summarize(mean = mean(value)) %>%
    group_by(name) %>%
    arrange(name) %>%
    filter(mean == max(mean)) %>%
    ungroup() %>%
    arrange(pop) -> df.cols
  
  Q <- Q[, df.cols$name]
  P <- P[, df.cols$name]
  
  rownames(Q) <- pop
  colnames(Q) <- paste0("P", 1:ncol(Q))
  colnames(P) <- colnames(Q)
  
  
  # change P coding
  # P <- 1-P to represent minor allele frequency (consistent with G)
  # plot(P[1:100,1], rowMeans(G[1:100, pop=='YRI'])/2)
  P <- 1-P
  
  # X: a matrix of minor allele counts by each estimated latent population
  popsize <- round(colSums(Q)) # effective population size
  X <- round(P %*% diag(2*popsize))
  colnames(X) <- colnames(P)

  out.list <- list(G=G, P=P, Q=Q, X=X, popsize=popsize)
  
  return(out.list)
}




################################################################################
# functions for TreeMix 
################################################################################


treemix_output_to_L <- function(out_treemix){
  
  out <- out_treemix
  
  V <- data.frame(
    id = out$d[,1],
    label = ifelse(is.na(out$d[,2]), out$d[,3], out$d[,2])
  )
  V[V[,2]=='NOT_ROOT',2] <- ''
  
  E <- data.frame(
    from = out$e[,1],
    to = out$e[,2],
    driftsize = out$e[,3],
    weight = out$e[,4]
  )
  
  
  L <- matrix(0, nrow=nrow(V), ncol = nrow(E))
  rownames(L) <- V$label
  
  v_idx_Outgroup <- which(V[,2]=='Outgroup')
  
  # change (from=ROOT, to=Outgroup) --> (from=Outgroup, to=ROOT)
  e_idx <- which(E$to==V$id[v_idx_Outgroup])
  if ((V$id[which(V[,2]=='ROOT')]) != E[e_idx, 'from']){
    warning('ROOT is not directly connected to Outgroup' )
    return(NULL)
  }
  
  E[e_idx, 'to'] <- E[e_idx, 'from']
  E[e_idx, 'from'] <- V$id[v_idx_Outgroup]
  
  queue <- which(E$from == V$id[v_idx_Outgroup]) # queue is for edges
  processed <- c()
  
  while(length(queue)>0){
    e <- queue[1]
    queue <- queue[-1]
    
    from_id <- E[e, 'from']
    to_id <- E[e, 'to']
    driftsize <- E[e, 'driftsize']
    weight <- E[e, 'weight']
    
    L[which(V$id==to_id),] <- L[which(V$id==to_id),] + L[which(V$id==from_id),] * weight
    L[which(V$id==to_id), e] <- L[which(V$id==to_id), e] + driftsize * weight
    
    processed <- c(processed, e)
    if (all(which(E$to==to_id) %in% processed)){
      queue <- c(queue, which(E$from==to_id))
    }
  }
  
  L <- L[,processed] # reorder drifts, using the processed order
  E <- E[processed,]
  
  # combine the (Outrgoup --> Root) and (Root --> mu0) branches
  L[,2] <- L[,1] + L[,2]
  L <- L[,-1]
  
  L <- L[rownames(L)!="" & rownames(L)!="Outgroup" & rownames(L)!="ROOT" ,]
  L <- L[order(rownames(L)),]
  L <- L[,apply(L, 2, max)>0]
  colnames(L) <- paste0("D",0:(ncol(L)-1))
  
  L <- sqrt(L) # convert drift sizes into drift memberships
  
  return(L)
}


estimate_L_from_X_treemix <- function(X, popsize, include.outgroup=TRUE,
                                      m1=TRUE, m2=TRUE){
  
  # prepare temporary treemix input file
  treemix_input <- c()
  for (k in 1:ncol(X)){
    
    ra <- 2*popsize[k] - X[,k]# reference (ancestral) allele count
    ma <- X[,k] # minor (derived) allele count
    treemix_input <- cbind(treemix_input,
                           c(colnames(X)[k], paste(ra, ',', ma, sep='')))
  }
  
  if(include.outgroup==TRUE){
    ra <- rep(2, nrow(X))
    ma <- rep(0, nrow(X))
    treemix_input <- cbind(treemix_input,
                           c(paste0('Outgroup'), paste(ra, ',', ma, sep='')))
  }
  
  write.table(treemix_input, here('output', 'treemix', 'temp_treemix_input'),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  R.utils::gzip(here('output', 'treemix', 'temp_treemix_input'), overwrite=TRUE)
  
  
  # run treemix for m=0,1,2
  if(include.outgroup==TRUE){
    rooting <- '-root Outgroup'
  }else{
    rooting <- ''
    stop('This function is not yet implemented')
  }
  
  
  paste0("cd '/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix';", 
         'treemix -i temp_treemix_input.gz -seed 1 ', rooting, ' -o temp_treemix_m0') %>%
    system(ignore.stdout=TRUE, wait=TRUE)
  # use 'plotting_funcs.R' of treemix package to import results
  out_m0 <- plot_tree('/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix/temp_treemix_m0')
  L_m0 <- treemix_output_to_L(out_m0)
  
  
  
  if(m1==TRUE){
    paste0("cd '/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix';", 
           'treemix -i temp_treemix_input.gz -seed 1 ', rooting, ' -m 1 -o temp_treemix_m1') %>%
      system(ignore.stdout=TRUE, wait=TRUE)
    
    out_m1 <- plot_tree('/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix/temp_treemix_m1')
    L_m1 <- treemix_output_to_L(out_m1)
  }else{
    L_m1 <- NULL
  }
  if(m2==TRUE){
    paste0("cd '/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix';", 
           'treemix -i temp_treemix_input.gz -seed 1 ', rooting, ' -m 2 -o temp_treemix_m2') %>%
      system(ignore.stdout=TRUE, wait=TRUE)
    
    out_m2 <- plot_tree('/Users/jkang/Library/CloudStorage/Box-Box/research/drift/output/treemix/temp_treemix_m2')
    L_m2 <- treemix_output_to_L(out_m2)
  }else{
    L_m2 <- NULL
  }
  
  R.utils::gunzip(here('output', 'treemix', 'temp_treemix_m0.treeout.gz'),
                  overwrite=TRUE, remove=FALSE)
  tree <- read.tree(here('output', 'treemix', 'temp_treemix_m0.treeout'))
  tree_wo_Outgroup <- drop.tip(tree, 'Outgroup')
  
  
  out.list <- list(L_m0 = L_m0, L_m1 = L_m1, L_m2 = L_m2, 
                   tree = tree, tree_wo_Outgroup = tree_wo_Outgroup)
  
  return(out.list)
}












################################################################################
# functions for plotting 
################################################################################

plot_colors <- c(
  
  '#008080', '#FFA500', '#008000','#C0C0C0',
  '#00FFFF', '#0000FF', '#FF00FF', '#800000', '#808000',
  '#4682B4', '#00FF00', '#FF0000', '#808080', '#FFC0CB',
  '#FFD700', '#B8860B', '#800080', '#DC143C', '#00FF7F',
  '#FF8C00', '#9400D3', '#00BFFF', '#FF69B4', '#D8BFD8',
  '#00FA9A', '#FF1493', '#7FFF00', '#4B0082', '#FFB6C1',
  '#9370DB', '#00CED1', '#BA55D3', '#FFA07A', '#7CFC00',
  '#8B008B', '#6B8E23', '#FFFF00', '#FFDAB9', '#ADFF2F',
  '#FFE4C4', '#FF00FF', '#32CD32', '#DB7093', '#87CEFA',
  '#48D1CC', '#FF6347', '#AFEEEE', '#DEB887', '#9932CC'
)


set_loadings_order <- function(L){
  
  set.seed(0)
  templist         <- list(L = L,F = L)
  class(templist)  <- c("multinom_topic_model_fit","list")
  y <- drop(umap_from_topics(templist, dims=1))
  loadings_order = order(y)
  
  return(loadings_order)
}


plot_Q <- function(Q, pop, loadings_order){
  
  templist         <- list(L = Q, F = Q)
  class(templist)  <- c("multinom_topic_model_fit","list")
  
  fig <- fastTopics::structure_plot(templist, 
                                    topics=ncol(Q):1,
                                    grouping = pop,
                                    colors = rev(plot_colors),
                                    loadings_order = loadings_order,
                                    gap = 16,
                                    verbose = FALSE) +
    labs(y = "Population membership", x = 'Individuals')
  
  fig$labels$fill <- 'Population'
  fig$labels$colour <- 'Population'
  return(fig)
}

reorder_L <- function(L, L.ref){
  
  if (any(rownames(L) != rownames(L.ref))){
    stop("L and L.ref must have the same row names.")
  }
  
  Lnew <- matrix(0, nrow=nrow(L), ncol=ncol(L))
  remaining_col <- 1:ncol(L)
  for (i in 1:min(ncol(L.ref), ncol(L))){
    diffs <- sapply(remaining_col, function(j){ sum((L.ref[,i] - L[,j])^2) })
    Lnew[,i] <- L[, remaining_col[which.min(diffs)]]
    remaining_col <- remaining_col[-which.min(diffs)]
  }
  
  if(ncol(L) > ncol(L.ref)){
    Lnew[,(ncol(L.ref)+1):ncol(L)] <- L[,remaining_col]
  }
  rownames(Lnew) <- rownames(L)
  colnames(Lnew) <- paste0('D',0:(ncol(Lnew)-1))
  
  return(Lnew)
}



plot_L <- function(L, legend = FALSE, remove.D0=TRUE){
  
  if(remove.D0 == TRUE){
    L <- L[,colnames(L)!='D0']
  }
  if(remove.D0 == FALSE){
    L <- cbind(L[, colnames(L)=='D0', drop=FALSE], L[, colnames(L)!='D0'])
    plot_colors <- c('black', plot_colors)
  }
  
  
  templist         <- list(L = L, F = L)
  class(templist)  <- c("multinom_topic_model_fit","list")
  
  fig <- fastTopics::structure_plot(templist, 
                                    topics=ncol(L):1,
                                    grouping = factor(rownames(L), levels=rownames(L)),
                                    colors = plot_colors,
                                    loadings_order = rep(1:nrow(L), each=5),
                                    gap = 1,
                                    verbose = FALSE) +
    theme(legend.position='left')+
    labs(y = "Drift membership", x = 'Populations')
  
  if(legend==FALSE){
    fig <- fig + theme(legend.position="none")
  }
  
  
  fig$labels$fill <- 'D'
  fig$labels$colour <- 'D'
  
  return(fig)
}


plot_M <- function(M, pop, loadings_order, remove.D0=TRUE, legend=FALSE){
  
  if(remove.D0 == TRUE){
    M <- M[,colnames(M)!='D0']
  }
  if(remove.D0 == FALSE){
    M <- cbind(M[, colnames(M)=='D0', drop=FALSE], M[, colnames(M)!='D0'])
    plot_colors <- c('black', plot_colors)
  }
  
  
  templist         <- list(L = M, F = M)
  class(templist)  <- c("multinom_topic_model_fit","list")
  
  fig <- fastTopics::structure_plot(templist, 
                                    topics=ncol(M):1,
                                    grouping = pop,
                                    colors = plot_colors,
                                    loadings_order = loadings_order,
                                    gap = 16,
                                    verbose = FALSE) +
    labs(y = "Drift membership", x = 'Individuals')
  if(legend==FALSE){
    fig <- fig + theme(legend.position="none")
  }
  
  
  fig$labels$fill <- 'D'
  fig$labels$colour <- 'D'
  return(fig)
}


plot_M2 <- function(L,
                    M,
                    P,
                    G,
                    pop,
                    loadings_order,
                    legend = FALSE
){
  
  # estimate D0
  if ('D0' %in% colnames(L)){
    M <- M[,colnames(M)!='D0']
  }
  D0 <- sqrt(max(0, mean(crossprod(P)/nrow(P) - tcrossprod(L))))
  M <- cbind(D0, M)
  colnames(M)[1] <- 'D0'
  
  M2.scaled <- 4*M^2
  gi2.scaled <- colSums(G^2)/nrow(G)
  residuals <- gi2.scaled - rowSums(M2.scaled)
  
  M2.scaled <- cbind(M2.scaled, residuals)
  colnames(M2.scaled) <- c(colnames(M), 'DX')
  
  plot_colors <- c('grey90', plot_colors) # color of the common drift D0
  plot_colors[ncol(M2.scaled)] <- 'beige' # color of the residuals
  
  templist         <- list(L = M2.scaled, F = M2.scaled)
  class(templist)  <- c("multinom_topic_model_fit","list")
  
  fig <- fastTopics::structure_plot(templist, 
                                    topics=ncol(M2.scaled):1,
                                    grouping = pop,
                                    colors = plot_colors,
                                    loadings_order = loadings_order,
                                    gap = 16,
                                    verbose = FALSE) +
    labs(y = "Drift membership squared", x = 'Individuals')
  if(legend==FALSE){
    fig <- fig + theme(legend.position="none")
  }
  
  fig$labels$fill <- 'D'
  fig$labels$colour <- 'D'
  return(fig)
}



plot_T <- function(out_treemix, edgelabels_adj = 1){
  tr <- out_treemix$tree_wo_Outgroup
  
  # adjust the edge length to be consistent with L plot
  tr$edge.length <- sqrt(tr$edge.length)
  
  # reorder edges to be consistent with L plot
  remaining_idx <- 1:nrow(tr$edge)
  sorted_edge_idx <- c()
  for (i in 1:(ncol(out_treemix$L_m0)-1)){
    driftsize_in_L <- max(out_treemix$L_m0[,i+1])
    closest_idx <- remaining_idx[which.min((tr$edge.length[remaining_idx] - driftsize_in_L)^2)]
    sorted_edge_idx <- c(sorted_edge_idx, closest_idx)
    remaining_idx <- remaining_idx[-which(remaining_idx == closest_idx)]
  }
  tr$edge <- tr$edge[sorted_edge_idx,]
  tr$edge.length <- tr$edge.length[sorted_edge_idx]
  
  plot(tr, edge.color = plot_colors, edge.width = 3, cex=1.5); 
  edgelabels(col='white', bg=plot_colors, cex=1.5, adj = edgelabels_adj)
}





