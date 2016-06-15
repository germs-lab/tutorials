#Tutorial on co-occurrence analysis using Spearman's rank correlation and Hoeffding dependence analysis

##This analysis is all done in R

##Packages required:
```
library(Hmisc)
library(plyr)
library(reshape2)
```

###1. Let's simulate a simple set of data:
```
## a vector of 100 randomly normally distributed values with a mean of 1, standard deviation of 9.  
x<-rnorm(100, mean=1, sd=9)
## a vector of 100 randomly normally distributed values with a mean of 14, standard deviation of 9.  
y<-rnorm(100, mean=14, sd=9)
## a quick plot to check the relationship between x and y
plot(y~x)
```
We can easily see that these two vectors are very much independent of each other.
But what would the statisics say?
```
## let's run a quick Spearman's correlation analysis:
sp<-rcorr(x, y, "spearman")
sp
```
Here is the results:
a. the x to y relationship (spearman's coefficient rho)  matrix:
```
     x    y
x 1.00 0.17
y 0.17 1.00
```
b. the number of replicates:
```
n= 100 
```
c. the p value:
```
P
  x     y    
x       0.082
y 0.082   
```
d. also try these and what do you get?
```
sp$r
sp$n
sp$P
```
#####Note:
a. rho value ranges from -1 to 1, with 0 as no correlation at all, -1 as greatly negatively correlated, and 1 as strongly positively correlated.

b. Spearman's correlation says that x and y do not have a significant relationship.

###2. Let's simulate a linear relationship:
```
x<-rnorm(100, mean=1, sd=9)
y = rnorm(100, mean=17, sd=5) +0.8*x
plot(y~x)
```
#####Quick Exercise: 
a. What's the Spearman correlation coefficient of x and y? Are they significant (p value cutoff at 0.05)

b. What if you change the method "spearman" to "pearson"?

###3. Let's simulate a non-linear but monotonic relationship:
```
x <- rnorm(100, mean=1, sd=9)
y <- rnorm(100, mean=800, sd=1000) + x^3
plot (y ~ x)
```
#####Quick Exercise: 
a. What's the Spearman correlation coefficient of x and y? Are they significant (p value cutoff at 0.05)

b. What if you change the method "spearman" to "pearson"?

###4. Let's simulate a non-monotonic relationship:
```
x <- rnorm(100, mean=1, sd=9)
y <-  rnorm(100, mean=100, sd=100) - x^2
plot (y ~ x)
```
##### Quick Exercise: 
a. What's the Spearman correlation coefficient of x and y? Are they significant (p value cutoff at 0.05)

b. What if you change the method "spearman" to "pearson"?

c. now try: `hoeffd(x, y)`, what do you get?

#####What is happening?
1. Spearman's rank correlation is very robust, but it assumes monotonic relationships  
2. Pearson correlation detects linear relationship, which is a special case of monotonic relationship   
3. Hoeffding's D evaluates the dependence of two datasets. It ranges from -1 to 1, the smaller the value (e.g, -1), the less the two data were dependent of each other.    
4. Real world data (e.g., microbial populations) are frequently monotonic and non-monotonic. Evaluate relationships with both Spearman's correlation and Hoeffding dependence test. Pay attention to relationships with high Hoeffding's D and low Spearman's rho. 


###5. Let's put simulated data together:
```
x <- rnorm(100, mean=1, sd=9)
independent<-0+rnorm(100, mean=14, sd=9)
linear <- rnorm(100, mean=17, sd=5) +0.8*x
monotonic <- rnorm(100, mean=800, sd=1000) + x^3
non.mono <-  rnorm(100, mean=100, sd=100) - x^2
df<-data.frame(x, independent, linear, monotonic, non.mono)

results_sp<-rcorr(as.matrix(df),type="spearman")
results_hd<-hoeffd(as.matrix(df))
        
#make two seperate objects for p-value and correlation coefficients
rhos<-results_sp$r
sp_ps<-results_sp$P
ds<-results_hd$D
ds_ps<-results_hd$P

# going to melt these objects to 'long form' where the first two columns make up the pairs of relationships. also removing NA's as they are self-comparisons.
# let's first get all of the P values:
sp_melt<-na.omit(melt(sp_ps))
ds_melt<-na.omit(melt(ds_ps))

#creating a qvalue (adjusted pvalue) based on FDR
sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")

#making column names more relevant
names(sp_melt)[3]<-"spearman_pval"
names(ds_melt)[3]<-"hoeffding_pval"

# now let's get rho and d into the same format
rho.melt<-na.omit(melt(rhos))
d.melt<-na.omit(melt(ds))

names(rho.melt)[3]<-"rho"
names(d.melt)[3]<-"D"

#merging together 
sp_merged<-merge(sp_melt,rho.melt,by=c("Var1","Var2"))
ds_merged<-merge(ds_melt, d.melt,by=c("Var1","Var2"))

# final results
results<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
# subsetting for p value < 0.05, either spearman adj.p or hoeffding adj.p
final_results<-subset(results, spearman_qval < 0.05 | hoeffding_qval < 0.05)
```
##### What does final_results look like? Which relationship is the strongest?

###6. From relationships to network statistics. See `http://jfaganuk.github.io/2015/01/24/basic-network-analysis/` for detailed explanations.
```
library(igraph)
library(sna)

# turn relationships into a network graph
# the relationships here are undirected
temp.graph<-graph.edgelist(as.matrix(final_results[,c(1,2)]),directed=FALSE)
# E() assigns edge atributes
E(temp.graph)$weight<-abs(final_results$rho)
# because our network is undirected, every relationship is duplicated. "simplify" removes duplicated nodes and edges
temp.graph<-simplify(temp.graph)

# number of nodes
N_nodes<-vcount(temp.graph)
# number of edges
N_edges<-ecount(temp.graph)

# how well do our network cluster?
g.components <- clusters(temp.graph)
# number of clusters
N_clusters<-g.components$no
# size of the largest cluster in the network
Max_csize<-max(g.components$csize)
# and the number of edges associated with
Max_c_edges<-ecount(induced.subgraph(temp.graph, which(g.components$membership == which.max(g.components$csize))))

# the number of edges in a network divided by the total number of possible edges
g.density<-graph.density(temp.graph)

# the paths around the network on average
g.pathlength.avg<-average.path.length(temp.graph)

# centrality: how connected is the network?
# by betweenness: the number of shortest paths from all nodes to all others that pass through that node.
betcent<-centralization.betweenness(temp.graph)$centralization
# by degree: the number of edges a node has.
degcent<-centralization.degree(temp.graph)$centralization

# modularity: are nodes connected within their clusters (modules) more than between clusters than they would by chance?
g.mod<-modularity(edge.betweenness.community(temp.graph))

# effective size: the number of relathips between a node and other unrelated nodes
## effective size function ##
# definite a function like this
ego.effective.size <- function(g, ego, ...) {
  egonet <- induced.subgraph(g, neighbors(g, ego, ...))
  n <- vcount(egonet)
  t <- ecount(egonet)
  return(n - (2 * t) / n)
}
effective.size <- function(g, ego=NULL, ...) {
  if(!is.null(ego)) {
    return(ego.effective.size(g, ego, ...))
  }
  return (sapply(V(g), function(x) {ego.effective.size(g,x, ...)}))
}
# calculate average effective size
efsize_avg<-mean(effective.size(temp.graph, mode = "all"))

# combine them together
test<-cbind(N_nodes, N_edges, N_clusters, Max_csize, Max_c_edges, g.density, g.pathlength.avg, betcent, degcent, g.mod, efsize_avg)
```

###7. Plot it!
```
plot(temp.graph)
```

