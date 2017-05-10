#probit
library(mblm)
library(Kendall)
list_substs <- function(graph_data, min_reversals){
    uniq <- unique(graph_data$site_node)
    uniq <- as.character(uniq)
    list <- lapply(uniq, function(e) {
            if (length(graph_data$site_node[graph_data$site_node == e]) >= min_reversals) {
                return (e)
            }
            else {
                return (NA);
            }
            })
    list <- list[!is.na(list)]
    list <-as.character(list)
}


kendall_stats <-function (graph_data, substs_list){
    coeffs <- lapply(substs_list, function(e) {
              node_data <- graph_data[graph_data$site_node == e,]
              x <- node_data$radius
              y <- node_data$density
              ts <- mblm(y~x)
              ts[["coefficients"]]["x"]}
              )
    ts_coeffs <- sapply(coeffs, function(e) e)
    mk <- lapply (substs_list, function(e) {node_data <- graph_data[graph_data$site_node == e,]
          MannKendall(node_data$density)}
          )
    mk_pvalue <- sapply(mk, function(e) e[2]$sl[1])
    mk_tau <-sapply(mk, function(e) e[1]$tau[1])
    results <-data.frame(substs_list, ts_coeffs, mk_tau, mk_pvalue)
}

plot_substs <- function(substs_list, graph_data, prefix){
  lapply(substs_list, function(e){
    data  <- graph_data[graph_data$site_node == e,]
    file <- file.path(dirname(getwd()), "output", "plots", paste(prefix, e, ".png", sep=""), fsep = .Platform$file.sep)
    #file <-paste(prefix, e, ".png", sep="")
    print (file)
    png(filename = file)
    print(data$radius)
    plot(data$radius, data$density, type="l")
    dev.off()
  })
}

#graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_diff_rings_h3.csv")
graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_n2.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)

list <- list_substs(graph_data, 3)
results <- kendall_stats(graph_data, list)
best_sites <- results[results$mk_pvalue < 0.05,]

hist(results$ts_coeffs, seq(from = min(results$ts_coeffs), to = max(results$ts_coeffs)+0.0001, by = 0.00005))
hist(results$mk_tau, seq(from = min(results$mk_tau), to = max(results$mk_tau)+0.0001, by = 0.00005))

hist(best_sites$ts_coeffs, seq(from = min(best_sites$ts_coeffs), to = max(best_sites$ts_coeffs)+0.0001, by = 0.00005))
hist(best_sites$mk_tau, seq(from = min(best_sites$mk_tau), to = max(best_sites$mk_tau)+0.0001, by = 0.00005))

plots <-plot_substs(as.character(best_sites$substs_list), graph_data, "n2_diff_rings_plot_")

#####
h1_longliving <- c("4_INTNODE2429","60_INTNODE2435","73_INTNODE2425","78_INTNODE2434","85_INTNODE2434","90_INTNODE2429","111_INTNODE2428","113_INTNODE2429","138_INTNODE2416","145_INTNODE2426","151_INTNODE2425","169_INTNODE2434","171_INTNODE2429","176_INTNODE2434","178_INTNODE2426","179_INTNODE2429","184_INTNODE2429","186_INTNODE2429","201_INTNODE2434","206_INTNODE2435","238_INTNODE2422","240_INTNODE2421","261_INTNODE2435","269_INTNODE2423","292_INTNODE2433","293_INTNODE2415","311_INTNODE2423","415_INTNODE2429","470_INTNODE2428","488_INTNODE2434")
plots <-plot_substs(as.character(h1_longliving), graph_data, "h1_longliving_diff_rings_plot_")

h3_longliving <- c("495_INTNODE4258","469_INTNODE4224","315_INTNODE4207","291_INTNODE4258","276_INTNODE4232","258_INTNODE4258","233_INTNODE4230","229_INTNODE4224","206_INTNODE4200","204_INTNODE4256","189_INTNODE4224","179_INTNODE4222","176_INTNODE4231","175_INTNODE4215","173_INTNODE4200","162_INTNODE4232","159_INTNODE4230","142_INTNODE4238","138_INTNODE4257","110_INTNODE4209","98_INTNODE4207","94_INTNODE4262","79_INTNODE4239","70_INTNODE4232","69_INTNODE4232","66_INTNODE4232","47_INTNODE4266","19_INTNODE4238","18_INTNODE4224","16_INTNODE4236","15_INTNODE4236","14_INTNODE4262")
plots <-plot_substs(as.character(h3_longliving), graph_data, "h3_longliving_diff_rings_plot_")

## plot sites with significant LR for exponential and weibull model + signif kendall trends
###
graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_h3.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)
h3_best_likelihood<- c("19_INTNODE4238", "69_INTNODE4232", "66_INTNODE4232", "110_INTNODE4209", "15_INTNODE4236", "175_INTNODE4215", "19_INTNODE4238", "69_INTNODE4232", "241_INTNODE3838", "466_INTNODE4200", "66_INTNODE4232", "16_INTNODE4236", "110_INTNODE4209", "15_INTNODE4236", "278_INTNODE4195", "175_INTNODE4215", "202_INTNODE3930", "546_INTNODE3932", "70_INTNODE4232", "99_INTNODE4207", "233_INTNODE4230", "140_INTNODE4152", "174_INTNODE4127", "138_INTNODE4257", "209_INTNODE4241", "209_INTNODE4201", "153_INTNODE4102", "", "153_INTNODE4232", "", "218_INTNODE3933", "149_INTNODE4200", "149_INTNODE4230", "159_INTNODE4230", "19_INTNODE4238", "291_INTNODE4258")
#h3_big_p<- c("149_INTNODE4230", "160_INTNODE4257", "236_INTNODE2473")
plots <-plot_substs(as.character(h3_best_likelihood), graph_data, "h3_plot_")

graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_h1.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)
h1_best_likelihood<- c("238_INTNODE2422", "205_INTNODE2406", "169_INTNODE2434", "169_INTNODE2394", "73_INTNODE2425", "151_INTNODE2425", "171_INTNODE2429", "202_INTNODE2406", "144_INTNODE2406", "155_INTNODE2428", "85_INTNODE2434", "111_INTNODE2428", "289_INTNODE2399", "113_INTNODE2429" )
plots <-plot_substs(as.character(h1_best_likelihood), graph_data, "h1_plot_")

graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_n1.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)
#n1_best_likelihood<- c("250_INTNODE2456", "16_INTNODE2390", "222_INTNODE2391", "84_INTNODE2390", "86_INTNODE2384", "382_INTNODE2390", "248_INTNODE2479", "366_INTNODE2455", "382_INTNODE2390", "200_INTNODE2384", "332_INTNODE2384", "336_INTNODE2442")
#big aposteriori for big p
n1_wood_likelihood<- c("336_INTNODE2442","369_INTNODE2384","388_INTNODE2384","34_INTNODE2360","250_INTNODE2390","367_INTNODE2455","369_INTNODE2383","222_INTNODE2064","222_INTNODE2390","339_INTNODE2390","434_INTNODE2391","249_INTNODE3063","248_INTNODE2473","151_INTNODE1891")
plots <-plot_substs(as.character(n1_wood_likelihood), graph_data, "n1_wood_plot_")

graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_n2.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)
#n2_best_likelihood<- c("466_INTNODE4636", "435_INTNODE4571", "197_INTNODE4596", "248_INTNODE4608", "329_INTNODE4636", "221_INTNODE4476", "26_INTNODE4675", "434_INTNODE4675", "199_INTNODE4587", "466_INTNODE4636", "367_INTNODE4582", "313_INTNODE4608", "331_INTNODE4582", "127_INTNODE4476", "435_INTNODE4571", "339_INTNODE3984", "197_INTNODE4596", "248_INTNODE4608", "329_INTNODE4636", "221_INTNODE4476", "249_INTNODE4367", "126_INTNODE4646", "347_INTNODE4637")
#big aposteriori for big p
n2_wood_likelihood <- c("199_INTNODE4587","248_INTNODE4608","197_INTNODE4596")
plots <-plot_substs(as.character(n2_wood_likelihood), graph_data, "n2_wood_plot_")


##

graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_h3.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)
h3_chosen <- c("161_INTNODE3854")
plots <-plot_substs(as.character(h3_chosen), graph_data, "h3_chosen_plot_")

##

graph_data <-read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/Output/egor_diff_rings_n2.csv")
graph_data <-data.frame("radius"=graph_data$radius, "site_node"=paste(graph_data$site, graph_data$node, sep="_"), "density"=graph_data$density)

list <- list_substs(graph_data, 3)
results <- kendall_stats(graph_data, list)
best_sites <- results[results$mk_pvalue < 0.05,]

hist(results$ts_coeffs, seq(from = min(results$ts_coeffs), to = max(results$ts_coeffs)+0.0001, by = 0.00005))
hist(results$mk_tau, seq(from = min(results$mk_tau), to = max(results$mk_tau)+0.0001, by = 0.00005))

hist(best_sites$ts_coeffs, seq(from = min(best_sites$ts_coeffs), to = max(best_sites$ts_coeffs)+0.0001, by = 0.00005))
hist(best_sites$mk_tau, seq(from = min(best_sites$mk_tau), to = max(best_sites$mk_tau)+0.0001, by = 0.00005))

plots <-plot_substs(as.character(best_sites$substs_list), graph_data, "n2_diff_rings_plot_")
