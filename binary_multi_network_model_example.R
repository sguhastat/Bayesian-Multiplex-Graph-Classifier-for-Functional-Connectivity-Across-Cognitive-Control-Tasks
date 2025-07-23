

source("binary_multi_network_sim_data_gen.R")
source("binary_multi_network_function_with_results.R")

# Generate Data

V    <- 20      
N    <- 250       
R.star    <- 4  
pi.nz <- 0.2

data = binary_multi_network_sim_data_gen(V = V, N = N, R.star = R.star, pi.nz = pi.nz, mu.tr = 0.5, seed = 1235)


# For each value of R fit the model

model_outputs = binary_multilayer_network_model(y = data$y, Xmat1 = data$Xmat1, Xmat2 = data$Xmat2, b1.tr = data$b1.tr,
                                                b2.tr = data$b2.tr,
                                                R = 6, niter = 1000, n_burnin = 500)


model_outputs

which(model_outputs$node_probabilities > 0.5)
sort(data$nz.index, decreasing = FALSE)
