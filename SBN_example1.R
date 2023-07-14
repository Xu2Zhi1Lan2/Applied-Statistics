###################################################################
## Example of Static Bayesian Net Modeling
###################################################################

library(bnlearn)
data(marks) # Examination marks of 88 students on five different topics
str(marks) # 5 variables, each with 88 observations

# mech: mechanics
# vect: vectors
# alg: algebra
# anl: analysis
# stat: statistics

### example: create an undirected graph

ug = empty.graph(names(marks)) # create an empty network

# object ug belongs to bn class, it contains
# learning: results of structure learning algorithm and its tuning parameters
# nodes: markov blanket, neighbor, parents and children
# arcs: two-column matrix, means arcs from the first column to the second column

arcs(ug) = matrix(
   c("MECH", "VECT", "MECH", "ALG", "VECT", "MECH",
   "VECT", "ALG", "ALG", "MECH", "ALG", "VECT",
   "ALG", "ANL", "ALG", "STAT", "ANL", "ALG",
   "ANL", "STAT", "STAT", "ALG", "STAT", "ANL"),
   ncol = 2, byrow = TRUE,
   dimnames = list(c(), c("from", "to"))) # structure of arcs and why to create arcs 

ug$arcs # print arcs

ug # print ug

### example: create a directed acyclic graph

dag = empty.graph(names(marks))

arcs(dag) = matrix(
  c("VECT", "MECH", "ALG", "MECH", "ALG", "VECT",
  "ANL", "ALG", "STAT", "ALG", "STAT", "ANL"),
  ncol = 2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))

dag$arcs

dag

mat = matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
               0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
               nrow = 5,
               dimnames = list(nodes(dag), nodes(dag))) # create adjacent matrix

mat

dag2 = empty.graph(nodes(dag))
amat(dag2) = mat # create DAG by adjacent matrix
all.equal(dag, dag2) # check whether DAGs created by arcs and adjacent matrix are equal

### concepts related to dag

hl1 = list(arcs = vstructs(dag, arcs = TRUE),
           lwd = 4, col = "black")

graphviz.plot(dag, highlight = hl1, layout = "fdp",
              main = "dag") # plot dag

node.ordering(dag) # topological order

nbr(dag, "ANL") # return the neighbor of a node

mb(dag, "ANL") # return the Markov blanket of a node

"ALG" %in% mb(dag, "ANL") # determine whether a node is the MB of another node

score(dag, data = marks, type = "loglik-g") # log-likelihood score, 
                                            # larger score means more conditional dependence

### Structural Learning

bn.hc = hc(marks) # object bn.hc belongs to bn class

hlhc = list(arcs = vstructs(bn.hc, arcs = TRUE),
           lwd = 4, col = "black")

graphviz.plot(bn.hc, highlight = hlhc, layout = "fdp",
              main = "bn.hc") # plot bn.hc

fitted = bn.fit(bn.hc, data = marks) # each node names a list of 7 components

fitted$STAT







