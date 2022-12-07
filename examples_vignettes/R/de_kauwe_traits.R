library(tidyverse)
austraits <- readRDS("examples_vignettes/data/austraits_develop_recent.rds")
austraits$traits %>%
  filter(grepl("Eucalyptus", taxon_name)) %>%
  drop_na(value, trait_name, taxon_name) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(taxon_name, trait_name) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  ungroup() %>% 
  filter(taxon_name %in% c("Eucalyptus blakelyi",
                           "Eucalyptus camaldulensis",
                           "Eucalyptus crebra",
                           "Eucalyptus dunnii",
                           "Eucalyptus globulus",
                           "Eucalyptus grandis",
                           "Eucalyptus largiflorens",
                           "Eucalyptus macrorhyncha",
                           "Eucalyptus melliodora",
                           "Eucalyptus obliqua",
                           "Eucalyptus populnea",
                           "Eucalyptus saligna",
                           "Eucalyptus sideroxylon",
                           "Eucalyptus tereticornis",
                           "Eucalyptus viminalis")) %>% 
  filter(trait_name %in% c("huber_value",
                           "leaf_N_per_dry_mass",
                           "leaf_N_per_area",
                           "seed_mass",
                           "specific_leaf_area",
                           "wood_density"
  )) %>%
  pivot_wider(names_from = trait_name, values_from = value) %>%
  mutate(leaf_N_per_area = if_else(is.na(leaf_N_per_area), leaf_N_per_dry_mass/specific_leaf_area, leaf_N_per_area)) %>%
  select(-leaf_N_per_dry_mass) %>%
  drop_na()


         