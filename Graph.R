Results <- read_rds("Results.rds")

Results2 <- read_rds("Results2.rds")

Results3 <- Results %>% mutate(Type = case_when(Meassurement %in% c("Ecos", "EcosTemp", "EcosTrop") ~ "Ecological entities",Meassurement %in% c("Richness", "RichnessTemp", "RichnessTrop") ~ "Richness")) %>% mutate(Meassurement = case_when(Meassurement %in% c("Ecos", "Richness") ~ "Total", Meassurement %in% c("EcosTrop", "RichnessTrop") ~ "Tropical", Meassurement %in% c("EcosTemp", "RichnessTemp") ~ "Temperate")) %>% mutate(Type = relevel(factor(Type), "Richness")) %>% mutate(Meassurement = relevel(factor(Meassurement), "Total"))


Results3 <- bind_rows(Results3)
saveRDS(Results3, "Results3.rds")


ggplot(Results3, aes(y = Mean, x = Time, group = Meassurement)) + geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = Meassurement), alpha = 0.3) + geom_path(aes(color = Meassurement)) + geom_point(aes(color = Meassurement)) + theme_bw() + facet_grid(Type ~ ., scales = "free_y") + theme(legend.position = "bottom") + scale_x_reverse() + xlab("Time [Ma]") + ylab("") + geom_vline(xintercept = 250, lty = 2, color = "red") + geom_vline(xintercept = 200, lty = 2, color = "red")

ggplot(Results2, aes(y = Mean, x = Time, group = Meassurement)) + geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = Meassurement), alpha = 0.3) + geom_path(aes(color = Meassurement)) + geom_point(aes(color = Meassurement)) + theme_bw() + facet_grid(Type ~ ., scales = "free_y") + theme(legend.position = "bottom") + scale_x_reverse() + xlab("Time [Ma]") + ylab("") + geom_vline(xintercept = 250, lty = 2, color = "red") + geom_vline(xintercept = 200, lty = 2, color = "red")


View(Results3)
