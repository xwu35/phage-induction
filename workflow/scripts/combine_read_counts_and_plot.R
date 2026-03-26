#!/usr/bin/env Rscript

# this script is used to combine read counts from different steps and plot the composition of reads

suppressPackageStartupMessages(library("tidyverse"))

# define input
trimmed <- snakemake@input$trimmed
mapped <- snakemake@input$mapped

# trimmed reads
counts_trimmed <- read.table(trimmed, 
                                   sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         trimmed_reads = num_seqs * 2) %>%
  select(sample, trimmed_reads)

# mapped to the genome
counts_mapped <- read.table(mapped, sep = " ", header = T) %>%
  rename(mapped = "read_counts") 
 
# merge all tables
df_lists = list(counts_trimmed, counts_mapped)

combined_table <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), df_lists)

# calculate number of removed reads at each step
combined_table <- combined_table %>% 
  mutate(unmapped = trimmed_reads - mapped)

# barplot of removed reads at each step
# define colors
colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999", "#CC79A7")

# plot
p_composition_reads <- combined_table %>%
  select(sample, mapped, unmapped) %>%
  pivot_longer(!sample, names_to = "type", values_to = "count") %>% 
  ggplot(aes(x=sample, y=count, fill=type)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  theme(axis.title = element_text(size = 11),
        axis.line = element_line(colour = "black"), 
        axis.text.x=element_text(size=8, colour = "black", angle=90, vjust=0.5),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x=element_blank(), 
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.position = "bottom") + 
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Proportion of reads") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors) # + coord_flip()  

# select few columns to write out
output_table <- combined_table %>%
  mutate(pct_mapped = mapped / trimmed_reads,
         pct_unmapped = unmapped / trimmed_reads)
  
# save table
write.table(output_table, snakemake@output$table,
            sep="\t", row.names = F, col.names=T, quote = F)

# save figure
ggsave(file=snakemake@output$figure,
       plot=p_composition_reads, 
       width=9, height=5)
