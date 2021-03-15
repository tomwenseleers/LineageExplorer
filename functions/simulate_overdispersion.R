# Functions to model an epidemic using a negative binomial
# distribution to take overdispersion into account.
#---------------------------------------------------------
library(tidygraph)
library(tibble)
library(ggraph)
library(patchwork)
library(scales)
library(dplyr)
library(tidyr)
#----------------------------------------------------------
#' Simulate total number of infections after n steps
#' @param nsteps the number of steps in the infection chain
#' @param nstart the number of infected patients
#' @param k the overdispersion parameter
#' @param R the reproductive number (seen as mean of the neg binom dist)
simul_cases <- function(nsteps = 5, nstart = 1,
                        k = 0.1, R = 1.3, lim = NULL){
  n <- nstart
  i <- 1
  while(i <= nsteps){
    sim <- rnbinom(n,size = k, mu = R)
    n <- if(is.null(lim)) sum(sim) else sum(pmin(sim,lim))
    i <- i + 1
  }
  n
}

#---------------------------------------------
#' Make distribution for a specific simulation
#' @param nsimul the number of simulations
#' @inheritParams simul_cases
dist_cases <- function(nsimul = 1000,
                       ...){
  args <- list(...)
  replicate(nsimul,
            do.call(simul_cases,args))
}

#--------------------------------------------------
#' Make probability distribution of dying epidemic
#' @param nstart a sequence of values for nstart (see simul_cases)
#' @inheritParams dist_cases

simul_prob_containment <- function(nstart = seq.int(50),
                                   fold = 1,
                       ...){
  nstart <- rep(nstart, fold)
  nout <- length(nstart)
  p <- numeric(nout)
  
  for(i in seq.int(nout)){
    tmp <- dist_cases(nstart = nstart[i],
                      ...)
      p[i] <- mean(tmp == 0)
  }
  tibble(nstart = nstart,
         p = p)
}

#-------------------------------
#' Simulate an infection chain
#' @inheritParams simul_cases
#' @return a tidygraph with a random infection chaing

simul_chain <- function(k=0.1, R = 1.3, nstart = 10, nsteps = 5){
  # Make the start cases
  ninfected <- nstart
  start <- paste0("r1case0-",seq.int(ninfected))
  
  # initialize the vectors to store results
  infected_by <- rep("r0case", length(start))
  id <- c("r0case",start)
  current <- start
  nmax <- nsteps
  
  n <- 2
  # Do the algorithm
  while(n <= nmax){
    # re-initialize the new id's
    newid <- character()
    
    for(i in seq_along(current)){
      # Draw from a neg binomial
      nnew <- rnbinom(1,size = k, mu = R)
      ninfected <- c(ninfected,nnew)
      if(nnew > 0){
        tmp <- paste0("r",n,"case",i,"-",seq.int(nnew))
        id <- c(id,tmp)
        infected_by <- c(infected_by,rep(current[i],nnew))
        newid <- c(newid,tmp)
      }
    }
    if(!length(newid)) break
    current <- newid
    n <- n + 1
  }
  if(n == nmax + 1){
    ninfected <- c(ninfected, rep(0,length(current)))
  }
  # Create vertices and edges
  vert <- tibble(
    name = unique(id),
    ninfected = ninfected,
    level = gsub("r(\\d+)case.*","\\1",name)
  )
  edge <- data.frame(from = match(infected_by, vert$name),
                     to = match(id[-1],vert$name),
                     length = 1)
  
  # Return
  tbl_graph(nodes = vert, 
            edges = edge,
            directed = TRUE)
  
}

#-------------------------------------------
#' Plot an infection chain
#' @param x a tidygraph object, coming from simul_chain
#' @param title title for the graph
#' @param maxcol the maximum number of infections getting their own
#' color. If maxcol is 2, you get the categories 0,1,2,>2. 
plot_chain <- function(x, title = "Spread",maxcol = 10){
  x <- x %>% activate(nodes) %>%
    mutate(ninfected = ifelse(level == 0,0,ninfected),
           dinfected = cut(ninfected,
                           breaks = c((-1):10,+Inf),
                           labels = c(as.character(0:10), ">10")))
  p1 <- ggraph(x, layout = "partition",
               direction = "out") +
    geom_edge_link() +
    geom_node_point(aes(color = dinfected,
                        size = sqrt(ninfected + 1) )) +
    theme_void() +
    scale_color_viridis_d(direction = -1, option = "B", end = 0.8, limits = c(as.character(0:10), ">10")) +
    guides(size = FALSE) +
    labs(color = "number\ninfected",
         title = title) 
  
  vert <- as_tibble(x) 
  p2 <- ggplot(vert, aes(y=level)) +
    geom_bar() +theme_minimal() +
    labs(y = "") +
    theme(axis.text.y = element_blank())
  p1 + p2 + plot_layout(widths = c(5,1)) 
}
