library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Overview table + figure for the 41 isolates retained after RevTrans filtering
#
# The retained set is defined by the taxa actually present in the final
# alignments (SML_*_noout_align_rtex.fasta), not by the metadata, so the table
# can never drift from what was fed to IQ-TREE.
# ------------------------------------------------------------------------------

# Run from the repository root:  Rscript amanda/script/07_isolate_overview.R
stopifnot(dir.exists("helene/out"), file.exists("complete_isolates_SML2.csv"))

dir.create("amanda/out", showWarnings = FALSE, recursive = TRUE)
dir.create("amanda/plot", showWarnings = FALSE, recursive = TRUE)

read_taxa <- function(path) {
  ln <- readLines(path, warn = FALSE)
  str_trim(str_remove(ln[str_starts(ln, ">")], "^>"))
}

taxa <- list(
  S = read_taxa("helene/out/SML_S_noout_align_rtex.fasta"),
  M = read_taxa("helene/out/SML_M_noout_align_rtex.fasta"),
  L = read_taxa("helene/out/SML_L_noout_align_rtex.fasta")
)
stopifnot(lengths(taxa) == 41)

meta <- read_csv("complete_isolates_SML2.csv", show_col_types = FALSE)

# Isolates retained = those whose S accession survived into the final alignment
kept_isolates <- meta %>%
  filter(Segment == "S", Accession %in% taxa$S) %>%
  pull(Isolate)
stopifnot(length(kept_isolates) == 41,
          !any(duplicated(kept_isolates)))

# One row per isolate, the three segment accessions side by side
tbl <- meta %>%
  filter(Isolate %in% kept_isolates) %>%
  select(Isolate, Organism_Name, Country, Collection_Date, Collection_Year,
         Segment, Accession) %>%
  pivot_wider(names_from = Segment, values_from = Accession,
              names_prefix = "acc_")

# Every retained isolate must carry all three segments, and each of those
# accessions must be the one present in that segment's alignment.
stopifnot(!any(is.na(c(tbl$acc_S, tbl$acc_M, tbl$acc_L))),
          setequal(tbl$acc_S, taxa$S),
          setequal(tbl$acc_M, taxa$M),
          setequal(tbl$acc_L, taxa$L))

virus_levels <- c("Schmallenberg virus", "Shamonda virus",
                  "Sathuperi virus", "Douglas virus")
virus_abbr <- c("Schmallenberg virus" = "SBV", "Shamonda virus" = "SHAV",
                "Sathuperi virus" = "SATV", "Douglas virus" = "DOUV")

tbl <- tbl %>%
  mutate(
    Virus   = factor(Organism_Name, levels = virus_levels),
    Abbr    = virus_abbr[as.character(Virus)],
    # read_csv types these columns differently depending on how many "NA"
    # strings they hold, so normalise to character before blanking them
    Country = na_if(as.character(Country), "NA"),
    Year    = suppressWarnings(as.integer(na_if(as.character(Collection_Year), "NA"))),
    Date    = na_if(as.character(Collection_Date), "NA")
  ) %>%
  arrange(Virus, Year, Country, Isolate)

stopifnot(nrow(tbl) == 41)
print(count(tbl, Virus, Abbr))

# ------------------------------------------------------------------------------
# Table 1 - csv, markdown and LaTeX
# ------------------------------------------------------------------------------

out_tbl <- tbl %>%
  transmute(
    Virus = Abbr,
    Isolate,
    Country = replace_na(Country, "n.d."),
    Year = ifelse(is.na(Year), "n.d.", as.character(Year)),
    `S segment` = acc_S, `M segment` = acc_M, `L segment` = acc_L
  )

write_csv(out_tbl, "amanda/out/table1_isolates.csv")

md <- c(
  paste0("| ", paste(names(out_tbl), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(out_tbl)), collapse = "|"), "|"),
  apply(out_tbl, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
)
writeLines(md, "amanda/out/table1_isolates.md")

tex <- c(
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{The 41 isolates retained for phylogenetic analysis, each represented",
  "by all three genome segments. Accession numbers are GenBank/ENA.",
  "n.d., no date or country recorded in the source record.}",
  "\\label{tab:isolates}",
  "\\begin{tabular}{lllllll}", "\\toprule",
  paste0(paste(names(out_tbl), collapse = " & "), " \\\\"), "\\midrule",
  apply(out_tbl, 1, function(r)
    paste0(paste(str_replace_all(r, "_", "\\\\_"), collapse = " & "), " \\\\")),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"
)
writeLines(tex, "amanda/out/table1_isolates.tex")

# ------------------------------------------------------------------------------
# Figure 1 - collection year x country, one panel per virus
#
# Species identity is carried by the panel, not by hue, so a single validated
# categorical step is used for every mark (four hues cannot clear the
# all-pairs colour-vision floors).
# ------------------------------------------------------------------------------

dot_col <- "#2a78d6"
ink     <- "#0b0b0b"
ink2    <- "#52514e"

plot_dat <- tbl %>%
  filter(!is.na(Year), !is.na(Country)) %>%
  count(Virus, Abbr, Country, Year, name = "n")

n_dropped <- nrow(tbl) - sum(plot_dat$n)
message("isolates not placed on the figure (no date/country): ", n_dropped)

panel_lab <- tbl %>%
  count(Virus, Abbr, name = "n_tot") %>%
  mutate(lab = glue("{Abbr} - {Virus}  (n = {n_tot})")) %>%
  select(Virus, lab) %>%
  deframe()

plot_dat <- plot_dat %>%
  # draw the largest marks first so smaller neighbours keep their ring and label
  arrange(desc(n)) %>%
  mutate(
    panel   = factor(panel_lab[as.character(Virus)], levels = unname(panel_lab)),
    # levels pinned alphabetically, so the draw order above cannot move the axis
    Country = factor(Country, levels = rev(sort(unique(Country))))
  )

p <- ggplot(plot_dat, aes(x = Year, y = Country)) +
  # shape 21 + white stroke = the surface ring that keeps marks from adjacent
  # years reading as one blob where they overlap
  geom_point(aes(size = n), shape = 21, fill = dot_col, colour = "white",
             stroke = 0.8) +
  geom_text(
    data = ~ filter(.x, n > 1), aes(label = n),
    colour = "white", size = 2.6, fontface = "bold"
  ) +
  facet_grid(rows = vars(panel), scales = "free_y", space = "free_y",
             switch = "y") +
  scale_size_area(max_size = 8.5, breaks = c(1, 2, 4, 6), name = "Isolates") +
  scale_x_continuous(breaks = seq(1980, 2026, 5), limits = c(1978, 2028),
                     expand = expansion(mult = 0.01)) +
  labs(
    title = "Sampling of the 41 isolates used for phylogenetic analysis",
    subtitle = "One point per country and year of collection; point area and label give the number of isolates",
    x = "Year of collection", y = NULL,
    caption = glue(
      "{n_dropped} isolates without a recorded date or country are omitted here ",
      "and listed as n.d. in Table 1."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 0, face = "bold",
                                      colour = ink, size = 9.5),
    panel.spacing.y    = unit(10, "pt"),
    plot.title         = element_text(face = "bold", colour = ink),
    plot.subtitle      = element_text(colour = ink2, size = 9),
    axis.text          = element_text(colour = ink2),
    legend.position    = "bottom",
    legend.justification = "left",
    plot.caption       = element_text(colour = ink2, size = 8, hjust = 0)
  )

ggsave("amanda/plot/fig1_isolate_overview.png", p, width = 9, height = 6.2,
       dpi = 400, bg = "white")
ggsave("amanda/plot/fig1_isolate_overview.pdf", p, width = 9, height = 6.2,
       bg = "white")

message("wrote amanda/out/table1_isolates.{csv,md,tex} and amanda/plot/fig1_isolate_overview.{png,pdf}")
