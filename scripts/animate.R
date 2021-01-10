library(ggplot2)
library(gganimate)


x <- filter(long, label %in% c("Light", "Light+Soil"), seed_rain == 100) %>%
  group_by(cohort) %>%
  mutate(age = time - min(time),
         across(time, floor),
         label2 = if_else(a_f4 == 1, "One cohort", "Constant\n recruitment")) %>%
  filter(density < 2, log_density > -14)

y <- ggplot(x, aes(age, height, size = height, alpha = density, color = age)) +
  geom_point(shape = 20) +
  scale_color_viridis_c(direction = -1) +
  scale_size(range = c(0, 5), guide = F) +
  scale_alpha(range = c(0.4, 1), guide = F) +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = .5) +
  facet_grid(label2 ~ label) +
  # Here comes the gganimate specific bits
  labs(title = 'Comparison of stand structure over time', 
       subtitle = 'Patch age: {round(frame_time, 0)} years', 
       x = "Cohort age (yrs)", y = 'Height (m)', 
       color = "Cohort\n age (yrs)") +
  transition_time(time) +
  ease_aes('linear')

animate(y, height = 6, width = 8, units = "in", res = 300)

anim_save("vignettes/plant.gif")


z <- ggplot(x, aes(log_density, height, size = height, alpha = density, color = age)) +
  geom_point(shape = 20) +
  scale_color_viridis_c(direction = -1) +
  scale_size(range = c(0, 5), guide = F) +
  scale_alpha(range = c(0.4, 1), guide = F) +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = .5) +
  facet_grid(label2 ~ label) +
  # Here comes the gganimate specific bits
  labs(title = 'Comparison of stand structure over time', 
       subtitle = 'Patch age: {round(frame_time, 0)} years', 
       x = "log(density)", y = 'Height (m)', 
       color = "Cohort\n age (yrs)") +
  transition_time(time) +
  ease_aes('linear')

animate(z, height = 6, width = 8, units = "in", res = 300)

anim_save("vignettes/plant2.gif")
