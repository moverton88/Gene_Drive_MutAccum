# Make figures for poster

###############################################################################
# Global LOH counts per clone ####
LOHcounts_in$dot_color <- txPal[as.numeric(LOHcounts_in$Tx_ID)]
LOHcounts_in$fill_color <- txPal[as.numeric(LOHcounts_in$Tx_ID)]
LOHcounts_in$fill_color[LOHcounts_in$ID %in% DA_clones] <- NA
LOHcounts_in <- LOHcounts_in %>% arrange(Tx_ID, n_LOH)
LOH_final_counts_stats$fill_color <- rev(txPal[as.numeric(LOH_final_counts_stats$Tx_ID)])

lab_y <- max(LOHcounts_in$n_LOH)
nudge <- c(2, 6, 10)
max_break <- ceiling(lab_y/5)*5
sr <- 1.1
x_diff <- (sr - 1)/50

LOHcount_dotplot_mean <- LOHcounts_in %>%
  # mutate(n_LOH = n_Inter) %>%
  # filter(Tx_ID == "C") %>%
  ggplot(aes(x = x_diff)) + 
  geom_dotplot(aes(y = n_LOH, group = Tx_ID), 
               fill = LOHcounts_in$fill_color,
               color = LOHcounts_in$dot_color,
               stroke = 2,
               binwidth = 1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize = 0.5, stackratio = sr) +
  stat_summary(aes(y = n_LOH), color = "grey30", fun = "mean", fun.min = "mean", fun.max= "mean",
               size = 0.2, width = 0.9, geom = "crossbar") +
  scale_fill_manual(values = txPal) +
  # scale_color_manual(values = txPal) +
  scale_y_continuous(breaks = seq(0, max_break, 5), limits = c(0, NA),
                     name = "Number of LOH events") +
  scale_x_continuous(breaks = 0, name = "Number of Clones") +
  theme(axis.ticks.y = element_blank(), 
        # axis.ticks.x = element_line(size = 0.25, color = "grey80"),
        # axis.ticks.length = unit(4, "mm"),
        plot.margin = unit(c(t = 0, r = 5, b = 5, l = 5), "mm"),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.title.y = element_text(hjust = 0.17, vjust = 1.5),
        # axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 30), 
        strip.text = element_text(size = 32),
        # axis.text = element_text(size = 12),
        legend.position = "none") +
  # facet_grid(.~Line %in% c("H_F", "H_H"))
  facet_grid(.~Tx_ID)

LOHcount_dotplot_mean

ggsave(file.path(outIntDir, "LOHcount_all_dotplot_2022_05.png"),
       plot = LOHcount_dotplot_mean,
       device = "png",
       width = 16, height = 9,
       units = "in",
       dpi = 600)


###############################################################################
# Chr V D strain LOH hotspot ####

BPrate_SW_ChrV <- ggplot() + 
  geom_rect(data = chrom_arms[5, ],
            aes(xmin = 0, xmax = (abs(arm_1_cent) + arm_2_cent)/1000,
                ymin = -Inf, ymax = Inf),
            alpha = 0.2, fill = "white", color = "black", size = 0.5) +
  geom_vline(data = chrom_arms[5, ],
             aes(xintercept = abs(arm_1_cent)/1000), size = 0.35, color = "blue4", alpha = 0.5) +
  # geom_segment(data = scale_bar[5, ], 
  #              aes(x = .x, xend = .xend, y = .y, yend = .yend),
  #              size = 1, color = "grey30") +
  geom_line(data = LOH_BPrate_SW %>%
              filter(rom_CHROM == "V") %>%
              filter(length > 40000),
            aes(x = mid_POS/1000, y = bp_rate/2, color = Tx_ID), size = 1) +
  scale_x_continuous(name = "Window Position (kb)", 
                     limits = c(0, NA),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = format_sci_10,
                     breaks = c(0, 1E-9, 2E-9),
                     limits = c(0, 2E-9),
                     expand = expansion(mult = c(0.05, 0)),
                     name = "Breakpoint rate /bp/gen") +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 0),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 24, vjust = 2.5),
        axis.title.x = element_text(size = 24, vjust = -1.7),
        axis.text = element_text(size = 18),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = "none")

BPrate_SW_ChrV

ggsave(file.path(outDir, "2022_PEQG_poster/BPrate_SW_SW50x5_ChrV.png"),
       plot = BPrate_SW_ChrV,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)
