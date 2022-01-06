library(WGCNA)
##load data
otu = read.table("otutab.txt", header=T, row.names=1, sep="\t")###

cor.mat <- corAndPvalue(otu,method = "spearman")#
#
library(neten)#
#
cor.neten <- Network_Enhancement(cor.mat$cor)#
#
write.table(cor.neten,file="7cor.neten.csv",sep = ",",row.names=FALSE)



column <- colnames(cor.mat$cor)

row <- rownames(cor.mat$cor)

colnames(cor.neten) <- column

rownames(cor.neten) <- row

plot(cor.neten,cor.mat$cor,xlim = c(.7,1))



require(RMThreshold)
rm.get.threshold(cor.neten)##
cor.p <- matrix(p.adjust(cor.mat$p,method = "fdr"),nrow = dim(cor.mat$p)[1])#p value adjusted

write.csv(cor.p,"8cor_p.csv")

cor.neten1 <- cor.neten
write.csv(cor.neten1,"9cor_neten1.csv")

cor.neten1[abs(cor.neten)< threshold|cor.p>0.05] <- 0

#Network Generation
library(igraph)
cor.net <- graph_from_adjacency_matrix(cor.neten1,diag = F,weighted = TRUE,mode = "undirected")
cor.net.clean <- induced_subgraph(cor.net, degree(cor.net)>0)

ly <- layout_with_kk(cor.net.clean)
plot.igraph(cor.net.clean,
            layout = ly,
            vertex.label = NA)

#Add properties for vertices

otutab_tax <- read.table("11cor-neten-sample.txt",header=T,sep = "\t",row.names = 1)
taxa <- as.matrix(otutab_tax)
taxa[1:10,1:6]

V(cor.net.clean)$Phylum <- taxa[V(cor.net.clean)$name,1]
V(cor.net.clean)$Order <- taxa[V(cor.net.clean)$name,2]
V(cor.net.clean)$Genus <- taxa[V(cor.net.clean)$name,3]

V(cor.net.clean)$degree <- degree(cor.net.clean)
V(cor.net.clean)$modularity <- membership(cluster_fast_greedy(cor.net.clean))
V(cor.net.clean)$betweenness = betweenness(cor.net.clean)


#Add properties for edges
edge.value <- cor.mat$cor[lower.tri(cor.mat$cor,diag = FALSE) & cor.neten1 != 0]
edge.dir <- rep("positive", length(E(cor.net.clean)))
edge.dir[edge.value < 0] <- "negative"
E(cor.net.clean)$edge.dir <- edge.dir


#Add properties for edges
cor.cut <- cor.mat$cor
positive.cut <- cor.cut[upper.tri(cor.cut)]
positive <- cor.mat$cor[upper.tri(cor.mat$cor)][positive.cut!=0]
positive[positive>0] <- "Positive"
positive[positive<0] <- "Negative"
E(cor.net.clean)$positve <- positive


#Export Network file
write.graph(cor.net.clean,file = "12bulk-0.02%1-strength.graphml",format = "graphml")#

source("netpro.R")
netpropoty <- netpro(cor.net.clean)
write.csv(netpropoty,"13netpropoty.csv")

#zp value
igraph <- cor.net.clean
V(igraph)$degree <- degree(igraph)

set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list)    

write.table(nodes_list, '14nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

source('zi_pi.r')

adjacency <- get.adjacency(cor.net.clean,sparse=FALSE)

nodes_list <- read.delim('14nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

nodes_list <- nodes_list[rownames(adjacency), ]


zi_pi <- zi.pi(nodes_list, adjacency, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, '15zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(zi_pi, '16zi_pi_result_type.txt', sep = '\t', row.names = FALSE, quote = FALSE)

p<-ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)
p
ggsave("17zp-reslut.pdf",p,width=7,height=6)
                                  


