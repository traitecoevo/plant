library(purrr)
library(tidyverse)
library(gganimate)
devtools::load_all()

extr <- function(m) {
  t = m$time
  h = m$species[[1]]["height", , ]
  d = m$species[[1]]["log_density", , ]
  n = ncol(d)
  
  df = tibble(cohort = rep(1:n, each = length(t)),
                   time = rep(t, times = n),
                   height = as.vector(h),
                   log_density = as.vector(d)) %>%
    mutate(density = exp(log_density)) %>%
    filter(!is.na(density))
}

mulga <- function(sla = 15.5, 
                  wood_density = 1100,
                  leaf_n = 3.6,
                  seed_mass = 9.26e-06, 
                  mat_height = 5,
                  rep_alloc = 6,
                  relative_belowground_importance = 1, 
                  recruitment_constraint = 0, 
                  canopy_dimension = 4, 
                  seed_rain_density = 1) {
  
  lma = sla^-1
  mol_n = 14.0067
  leaf_n_area = leaf_n * (mol_n * 1e-3) * lma
  
  p <- scm_base_parameters("FF16", "FF16_Env")
  
  m_h = make_FF16_hyperpar(lma_0 = sla^-1,
                            rho_0 = wood_density,
                            narea = leaf_n_area,
                            narea_0 = leaf_n_area,
                            latitude = coords$lat)
  
  m_t <- trait_matrix(c(mat_height, rep_alloc,
                        seed_mass, canopy_dimension,
                        relative_belowground_importance,
                        recruitment_constraint),
                      c("hmat", "a_f2", "omega", "eta", "k_2", "a_f4"))
  
  
  m <- expand_parameters(m_t, p, m_h, FALSE)
  m$seed_rain <- seed_rain_density
  
  return(m)
}

pars <- data.frame(label = factor(1:5, labels = c("Light", "low", "Light+Soil", "high", "very-high BG")),
                   k_2 = c(0, 0.05, 0.2, 2, 5)) %>%
  crossing(., a_f4 = c(0, 1)) %>%
  crossing(., seed_rain = c(1, 10, 100, 1000, 10000)) %>%
  rowid_to_column("run")


res <- pars %>% 
  mutate(res = pmap(list(k_2, a_f4, seed_rain),
                ~ {p = mulga(relative_belowground_importance = ..1,
                             recruitment_constraint = ..2,
                             seed_rain_density = ..3)
                   m = run_scm_collect(p)
                   return(m)}),
         out = map(res, extr))

long <- bind_rows(res$out, .id = "run") %>%
  mutate(across(run, as.numeric)) %>%
  left_join(pars)


## Animate one individual

p_ind <- mulga(relative_belowground_importance = 0,
               recruitment_constraint = 1,
               seed_rain_density = 1)

a <- FF16_Individual(p_ind$strategies[[1]])

tt <- seq(0, 110, length.out=101)
z <- grow_plant_to_time(a, tt, FF16_fixed_environment())

ind <- data.frame(time = tt,
           height = z$state[, 1],
           density = 1, 
           cohort = 1,
           timestep = tt * exp(-seq(0, 1, len = 101)) * 2)


gif <- ggplot(ind, aes(x = time, y = height)) +
  geom_line(size = 1.2, color = "forestgreen") +
  geom_point(aes(size = height), color = "forestgreen") + 
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep)

animate(gif, duration = 3.5, fps = 22, detail = 40, start_pause = 30,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/ind.gif",
                                   loop = F))


cohort <- filter(long, run == 1, cohort %in% c(1, 40, 60, 80, 
                                          seq(80, 140, by = 2))) %>%
  arrange(cohort, time) %>%
  group_by(cohort) %>%
  mutate(above_min = cumall(density > 0.3),
         age = time - min(time),
         timestep = age * exp(-seq(0, 1, len = n())) * 2) %>%
  filter(above_min) %>%
  ungroup()

gif <- ggplot(cohort, aes(x = time, y = height,
                group = cohort)) +
  geom_line(data = select(ind, -timestep), aes(group = NULL),
            size = 1.2, color = "forestgreen") +
  geom_line(size = 1, alpha = .1) +
  geom_point(aes(size = height)) +
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep) 

animate(gif, duration = 5, fps = 22, detail = 40, start_pause = 10,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/cohort_by_age.gif",
                                   loop = F))


low_cohort <- filter(long, run %in% c(6),
                 cohort %in% c(1, 40, 60, 62, 66, seq(68, 86, by = 2))) %>%
  arrange(cohort, time) %>%
  group_by(cohort) %>%
  mutate(above_min = cumall(density > 0.5),
         age = time - min(time),
         timestep = age * exp(-seq(0, 1, len = n())) * 2,
         across(seed_rain, ~ "Low density"),
         across(label, ~ paste(., " competition"))) %>%
  filter(above_min) %>%
  ungroup()

gif <- ggplot(low_cohort, aes(x = time, y = height,
                          group = cohort)) +
  geom_line(data = select(ind, -timestep), aes(group = NULL),
            size = 1.2, color = "forestgreen") +
  geom_line(size = 1, alpha = .1) +
  geom_point(aes(size = height)) +
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep) 

animate(gif, duration = 5, fps = 22, detail = 40, start_pause = 10,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/cohort2.gif",
                                   loop = F))


high_cohort <- filter(long, run %in% c(8),
                 cohort %in% c(1, 40, 60, 62, 66, seq(68, 86, by = 2))) %>%
  arrange(cohort, time) %>%
  group_by(cohort) %>%
  mutate(above_min = cumall(density > 0.5),
         age = time - min(time),
         timestep = age * exp(-seq(0, 1, len = n())) * 2,
         across(seed_rain, ~ "High density")) %>%
  filter(above_min) %>%
  ungroup()

gif <- ggplot(high_cohort, aes(x = time, y = height,
                          group = cohort)) +
  geom_line(data = select(ind, -timestep), aes(group = NULL),
            size = 1.2, color = "forestgreen") +
  geom_line(data = select(low_cohort, -timestep),
            size = 1, alpha = 0.1) +
  geom_line(size = 1, alpha = .4, colour = "steelblue2") +
  geom_point(aes(size = height), fill = "steelblue2") +
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep) 

animate(gif, duration = 5, fps = 22, detail = 40, start_pause = 10,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/cohort3.gif",
                                   loop = F))



bg_cohort <- filter(long, run %in% c(48),
                      cohort %in% c(1, 40, 60, 62, 66, seq(68, 86, by = 2))) %>%
  arrange(cohort, time) %>%
  group_by(cohort) %>%
  mutate(above_min = cumall(density > 0.1),
         age = time - min(time),
         timestep = age * exp(-seq(0, 1, len = n())) * 2,
         across(seed_rain, ~ "High density")) %>%
  filter(above_min) %>%
  ungroup()

gif <- ggplot(bg_cohort, aes(x = time, y = height,
                               group = cohort)) +
  geom_line(data = select(ind, -timestep), aes(group = NULL),
            size = 1.2, color = "forestgreen") +
  geom_line(data = select(low_cohort, -timestep),
            size = 1, alpha = 0.1) +
  geom_line(size = 1, alpha = .4, colour = "steelblue2") +
  geom_point(aes(size = height), fill = "steelblue2") +
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep) 

animate(gif, duration = 5, fps = 22, detail = 40, start_pause = 10,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/cohort_bg.gif",
                                   loop = F))


bg_cohort_all <- filter(long, run %in% c(33),
                   cohort %in% c(1, 40, 60, 80, 
                                 seq(80, 140, by = 2))) %>%
  arrange(cohort, time) %>%
  mutate(timestep = time * exp(-seq(0, 1, len = n())) * 2) %>%
  group_by(cohort) %>%
  mutate(above_min = cumall(density > 0.1),
         age = time - min(time),
         across(seed_rain, ~ "High density")) %>%
  filter(above_min) %>%
  ungroup()

gif <- ggplot(bg_cohort_all, aes(x = time, y = height,
                             group = cohort)) +
  geom_line(data = select(ind, -timestep), aes(group = NULL),
            size = 1.2, color = "forestgreen") +
  geom_line(data = select(low_cohort, -timestep),
            size = 1, alpha = 0.1) +
  geom_line(size = 1, alpha = .4, colour = "steelblue2") +
  geom_point(aes(size = height), fill = "steelblue2") +
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.2, 9))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(timestep) 

animate(gif, duration = 3, fps = 22, detail = 40, start_pause = 10,
        width = 7, height = 5, units = "in", res = 200,
        renderer = gifski_renderer(file = "vignettes/cohort_bg2.gif",
                                   loop = F))


slice(x, 100)
  group_by(cohort) %>%
  mutate(min = min(density),
         above_min = cumall(density > min)) %>%
  ungroup() %>%
  filter(above_min)

ggplot(x, aes(time, height, alpha = density)) +
  geom_line() + facet_wrap(~ cohort)
  
  gplot(., aes(x = time, y = height, group = cohort)) +
  geom_line(size = 1, alpha = .1) +
  geom_point(aes(size = height)) + 
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  scale_x_continuous(expand = c(0, 1), limits = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0.1))  +
  scale_size(range = c(1, 3)) +
  guides(alpha = F, size = F) +
  theme_classic(base_size = 18) +
  transition_reveal(time) +
  shadow_wake(wake_length = 0.2)

animate(x, duration = 2.5, detail = 5, 
        width = 7, height = 5, 
        units = "in", res = 150,
        renderer = gifski_renderer(file = "vignettes/cohort.gif",
                                   loop = F))

           
filter(long, run == 3) %>%
  ggplot(., aes(x = time, y = height, color = density,
                   group = cohort, alpha = log_density)) +
  geom_line() + 
  labs(x = "Time (years)",
       y = "Height (m)",
       color = "Density (ind.m-2)") +
  facet_grid(seed_rain ~ label) +
  scale_color_viridis_c(direction = -1) +
  guides(alpha = F) +
  theme_bw()



p1 <- expand_parameters(trait_matrix(4, "eta"), p0, FF16bg_hyperpar,FALSE)

p1$seed_rain <- 20
p1_out <- run_scm_collect(p1)

g <- patch(p1_out)
ggsave(g, filename = "vignettes/base.png")


# Increasing below-ground importance
p2 <- expand_parameters(trait_matrix(c(4, 0.05), c("eta", "k_2")),
                        p0, mutant = FALSE)
p2$seed_rain <- 20
p2_out <- run_scm_collect(p2)

g <- patch(p2_out)
ggsave(g, filename = "vignettes/small_bg.png")


# Increasing below-ground importance
p3 <- expand_parameters(trait_matrix(c(4, 0.2), c("eta", "k_2")),
                        p0, mutant = FALSE)
p3$seed_rain <- 20
p3_out <- run_scm_collect(p3)

g <- patch(p3_out)
ggsave(g, filename = "vignettes/med_bg.png")



# Increasing below-ground importance
p3 <- expand_parameters(trait_matrix(c(4, 2), c("eta", "k_2")),
                        p0, mutant = FALSE)
p3$seed_rain <- 20
p3_out <- run_scm_collect(p3)

g <- patch(p3_out)
ggsave(g, filename = "vignettes/large_bg.png", width = 7, height = 7)


# Increasing below-ground importance
p3 <- expand_parameters(trait_matrix(c(4, 1, 1), c("eta", "a_dG1", "k_2")),
                        p0, mutant = FALSE)
p3$seed_rain <- 20
p3_out <- run_scm_collect(p3)

g <- patch(p3_out)
ggsave(last_plot(), filename = "vignettes/large_bg.png", width = 7, height = 7)


