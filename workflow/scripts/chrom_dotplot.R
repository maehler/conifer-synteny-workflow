library(tidyverse)

source_species <- snakemake@wildcards[["species1"]]
target_species <- snakemake@wildcards[["species2"]]
source_seq <- snakemake@wildcards[["seq"]]

anchor_filename <- snakemake@input[["anchors"]]
source_bed_filename <- snakemake@input[["source_bed"]]
target_bed_filename <- snakemake@input[["target_bed"]]
source_fai_filename <- snakemake@input[["source_fai"]]
target_fai_filename <- snakemake@input[["target_fai"]]

source_column <- which(str_split(basename(anchor_filename), "\\.")[[1]] == source_species)
target_column <- ifelse(source_column == 1, 2, 1)

source_bed <- read_tsv(source_bed_filename,
    col_names = c("seq", "start", "end", "transcript"),
    col_types = "cddc") %>%
    filter(seq == source_seq)
target_bed <- read_tsv(target_bed_filename,
    col_names = c("seq", "start", "end", "transcript"),
    col_types = "cddc")

anchor_lines <- read_lines(anchor_filename)
block_rows <- which(str_detect(anchor_lines, "^#"))
anchor_block_lines <- split(anchor_lines,
                            cumsum(seq_along(anchor_lines) %in% block_rows))

source_blocks <- map_df(anchor_block_lines, ~ {
        d <- read_tsv(str_c(.[2:length(.)], collapse = "\n"),
            col_types = "ccd",
            col_names = c(ifelse(source_column == 1, "source_transcript", "target_transcript"),
                          ifelse(source_column == 2, "source_transcript", "target_transcript"),
                          "score")) %>%
            inner_join(source_bed, by = c("source_transcript" = "transcript"))
        if (nrow(d) > 0) {
            d
        } else {
            NULL
        }
    }, .id = "block")

summarised_blocks <- source_blocks %>%
    inner_join(target_bed %>% rename(t_seq = seq,
                                     t_start = start,
                                     t_end = end,
                                     target_transcript = transcript),
        by = "target_transcript") %>%
    group_by(block, seq, t_seq) %>%
    summarise(s_start = min(start),
              s_end = max(end),
              t_start = min(t_start),
              t_end = max(t_end),
              mean_score = mean(score), .groups = "drop")

target_lengths <- read_tsv(target_fai_filename,
                           col_types = "cd---",
                           col_names = c("t_seq", "end")) %>%
  filter(t_seq %in% summarised_blocks$t_seq)

source_length <- read_tsv(source_fai_filename,
    col_types = "cd---",
    col_names = c("seq", "end")) %>%
    filter(seq == source_seq) %>%
    pull(end)

n_panels <- nrow(target_lengths)

p <- ggplot(summarised_blocks) +
  facet_grid(rows = vars(t_seq), space = "free_y", scales = "free_y") +
  # Make sure the limits are set for each panel
  geom_blank(aes(x = source_length, y = 0)) +
  geom_blank(data = target_lengths,
             mapping = aes(x = 0, y = end)) +
  geom_segment(aes(x = s_start, xend = s_end,
                   y = t_start, yend = t_end),
               size = 0.5, colour = "firebrick2") +
  labs(x = str_c(source_species, source_seq, "position", sep = " "),
       y = str_c(target_species, "sequence position", sep = " "),
       title = str_c(source_species, source_seq,
                     "vs", target_species, "synteny", sep = " ")) +
  theme_void() +
  theme(axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(angle = 90, margin = margin(r = 5)),
        axis.line.x = element_line(),
        axis.ticks.length = unit(3, "pt"),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 10),
        panel.border = element_rect(fill = NA, size = 0.1),
        strip.text.y = element_text(hjust = 0, margin = margin(l = 5)),
        plot.title = element_text(margin = margin(b = 5)),
        plot.margin = margin(10, 10, 10, 7),
        plot.background = element_rect(fill = "white"))

ggsave(filename = snakemake@output[["png"]], plot = p,
       width = 7, height = 2 + 0.2 * n_panels)

