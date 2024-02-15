
max_patch_lifetime <- 0.1

##FF16 0.044
p0 <- scm_base_parameters("FF16")
p0$max_patch_lifetime <- max_patch_lifetime
env <- make_environment("FF16")
ctrl <- scm_base_control()

# one species
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar, 
                        birth_rate_list = list(20))

system.time(out <- run_scm(p1, env, ctrl))

##TF24 - single depths 0.046

p0 <- scm_base_parameters("TF24")
p0$max_patch_lifetime <- max_patch_lifetime
env <- make_environment("TF24", soil_number_of_depths = 1)
ctrl <- scm_base_control()

# one species
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, TF24_hyperpar, 
                        birth_rate_list = list(20))

system.time(out <- run_scm(p1, env, ctrl))

##TF24 - two depths 10.457

p0 <- scm_base_parameters("TF24")
p0$max_patch_lifetime <- max_patch_lifetime
env <- make_environment("TF24", soil_number_of_depths = 2)
ctrl <- scm_base_control()

# one species
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, TF24_hyperpar, 
                        birth_rate_list = list(20))

system.time(out <- run_scm(p1, env, ctrl))

##TF24 - 10 depths 18.068

p0 <- scm_base_parameters("TF24")
p0$max_patch_lifetime <- max_patch_lifetime
env <- make_environment("TF24", soil_number_of_depths = 10)
ctrl <- scm_base_control()

# one species
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, TF24_hyperpar, 
                        birth_rate_list = list(20))

system.time(out <- run_scm_collect(p1, env, ctrl))

##TF24 - 100 depths - small time decrease? 14.987

p0 <- scm_base_parameters("TF24")
p0$max_patch_lifetime <- max_patch_lifetime
env <- make_environment("TF24", soil_number_of_depths = 100)
ctrl <- scm_base_control()

# one species
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, TF24_hyperpar, 
                        birth_rate_list = list(20))

system.time(out <- run_scm_collect(p1, env, ctrl))
