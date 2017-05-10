graph_data<-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/h1_egor_site_entrenchment_graphs.csv")
graph_data<-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/smart_egor_n1_wr.csv")
site_data  <- graph_data[graph_data$site == 250, ]
node_data  <- site_data[site_data$node =="INTNODE2456",]
site_data  <- graph_data[graph_data$site == 261, ]
node_data  <- site_data[site_data$node =="INTNODE2435",]
plot(node_data$radius, node_data$density, type="l")

graph_data<-read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h3_neva_entrenchment_20koel.csv")

graph_data<-read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_neva_site_entrenchment_graphs_density_10.csv")
graph_syn_data<-read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/depth_groups_entremchment_test")


graph_data <-graph_data[graph_data$expected != 0,]
graph_syn_data <-graph_syn_data[graph_syn_data$expected != 0,]

plot(graph_data$observed, graph_data$expected,  xlim=c(0,0.18),ylim=c(0,0.18),asp =1)
plot(graph_syn_data$observed, graph_syn_data$expected,  xlim=c(0,0.18),ylim=c(0,0.18),asp =1)

rad = 2;
radius_data <- graph_data[graph_data$radius == rad, ]
radius_syn_data <- graph_syn_data[graph_syn_data$radius == rad, ]

par(mfrow=c(2,1)) 
plot(radius_data$observed, radius_data$expected, xlim=c(0,0.2),ylim=c(0,0.2),asp =1, main = paste("", rad*10, sep = " "))
plot(radius_syn_data$observed, radius_syn_data$expected, xlim=c(0,0.2),ylim=c(0,0.2),asp =1, main = paste("lifespan>150", rad*10, sep = " "))


graph_data_2 <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)

uniq <- unique(graph_data_2$site_node)
uniq <- as.character(uniq)
rich <- lapply(uniq, function(e) {if (length(graph_data_2$site_node[graph_data_2$site_node == e]) >3) {
                                      return (e)
}
else {
  return (NA);
}
})
rich <- rich[!is.na(rich)]
rich <-as.character(rich)
coeffs <- lapply(rich, function(e) {node_data <- graph_data_2[graph_data_2$site_node == e,]
                                  x <- node_data$radius
                                   y <- node_data$density
                                   ts <- mblm(y~x)
                                  print (ts[["coefficients"]]["x"])
                                   ts[["coefficients"]]["x"]}
                 )
coeffs <- sapply(coeffs, function(e) e)
hist(coeffs, seq(from = -0.02, to = 0.02, by = 0.0001))

mk <- lapply (rich, function(e) {node_data <- graph_data_2[graph_data_2$site_node == e,]
              MannKendall(node_data$density)}
              )
mk_vector <- sapply(mk, function(e) e[2]$sl[1])
results <-data.frame(rich, coeffs, mk_vector)
best <- results[results$mk_vector < 0.05,]