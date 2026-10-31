library(tidyverse)
#filtrelenmis olanlari gorsellestirmek uzere anotasyonlarini gerceklestir
filtered <- read_tsv("results/variants/filtered_calls.tsv", col_names = c("pos", "qual", "dp", "filter"),
                     col_types = "dddc")
filtered |> count(filter)

#Nokta grafiginin olusuturulmasi
ggplot(filtered, aes(dp, qual, colour = filter)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_y_log10() +
  geom_hline(yintercept = 30, linetype = "dashed") +
  geom_vline(xintercept = c(10, 100), linetype = "dashed") +
  labs(x = "INFO/DP", y = "QUAL (log scale)", title = "Which calls the chosen filter removes") +
  theme_bw()

#yapisal elementleri yazdır
annot <- read_tsv("results/annot/annotated.tsv", comment = "#", na = ".",
                  col_names = c("pos", "rsid", "ref", "alt", "gt", "bcsq", "gnomad_af", "gnomad_faf95",
                                "clnsig", "clnrevstat", "clnvid"),
                  col_types = "dcccccddccc") |>
  mutate(first_csq = str_split_i(bcsq, ",", 1)) |>
  separate_wider_delim(first_csq, delim = "|", too_few = "align_start", too_many = "drop",
                       names = c("consequence", "gene", "transcript", "biotype", "strand",
                                 "protein_change", "dna_change"))
annot |> count(consequence, sort = TRUE)

#kodlayan bölgelerdeki mutasyonlar
coding <- annot |>
  filter(str_detect(consequence, "missense|synonymous|stop|frameshift|inframe|splice|start_lost")) |>
  select(pos, rsid, gt, gene, consequence, protein_change, gnomad_af, gnomad_faf95, clnsig)
coding

#Bu kisinin varyantlari genel olarak database'de ne kadar yaygin
annot <- annot |>
  mutate(frequency = case_when(is.na(gnomad_af) ~ "not in gnomAD",
                               gnomad_af < 0.001 ~ "below 0.1 %",
                               gnomad_af < 0.01 ~ "0.1-1 %",
                               gnomad_af < 0.05 ~ "1-5 %",
                               .default = "5 % or more"),
         frequency = factor(frequency, levels = c("not in gnomAD", "below 0.1 %", "0.1-1 %", "1-5 %", "5 % or more")))
annot |> count(frequency)


#grafikle de gorelim:
ggplot(filter(annot, !is.na(gnomad_af)), aes(gnomad_af)) +
  geom_histogram(bins = 50) +
  scale_x_log10(labels = scales::label_percent()) +
  geom_vline(xintercept = c(0.001, 0.01, 0.05), linetype = "dashed") +
  labs(x = "allele frequency in gnomAD v4.1 genomes (log scale)", y = "variants",
       title = "How common are this person's variants?") +
  theme_bw()




