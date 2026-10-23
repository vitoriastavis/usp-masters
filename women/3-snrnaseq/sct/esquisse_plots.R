# Script para visualizar dados de DEGs no esquisse

# Instalar pacotes se necessário
if (!require("esquisse")) {
  install.packages("esquisse")
}
if (!require("readr")) {
  install.packages("readr")
}

# Carregar bibliotecas
library(esquisse)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Ler o arquivo CSV
dados_deg <- read_csv("degs_data.csv")

# Visualizar os dados
print(dados_deg)

dados_long <- dados_deg %>%
  pivot_longer(
    cols = -`Analysis`,
    names_to = "DEGs",
    values_to = "Number of DEGs"
  )

library(ggplot2)

pp <- ggplot(dados_long) +
 aes(x = Analysis, y = `Number of DEGs`, fill = DEGs) +
 geom_col(position = "dodge") +
 scale_fill_manual(values = c(Shared = "#87C909", GSE144136 = "#E69F00", GSE213982 = "#0072B2")) +
 theme_classic() +
 theme(legend.position = "top", axis.title.y = element_text(size = 14L),
       axis.title.x = element_text(size = 14L),
       axis.text.y = element_text(size = 14L),
       axis.text.x = element_text(size = 14L),
       legend.text = element_text(size = 14L),
      legend.title = element_blank())
ggsave("plots/Figure_2.tiff", width = 10)


load("intersection_all_degs.RData")
all_degs_df <- table(lengths(all_degs))
all_degs_df_long <- data.frame(
  N_alvos = names(all_degs_df),
  Contagem = as.integer(all_degs_df),
  row.names = NULL
)

# esquisse::esquisser(all_degs_df_long)

p_dist_degs <- ggplot(all_degs_df_long) +
  aes(x = N_alvos, y = Contagem) +
  geom_col(fill = "#6FAF9C") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 14L),
    axis.title.x = element_text(size = 14L),
    axis.text.y = element_text(size = 14L),
    axis.text.x = element_text(size = 14L)
  )

sizes_degs <- lengths(all_degs)
df_sizes_degs <- data.frame(N_alvos = sizes_degs)

p_dist_degs <- ggplot(df_sizes_degs, aes(x = N_alvos)) +
  geom_histogram(binwidth = 1, fill = "#6FAF9C", color = "white") +
  theme_classic() +
  labs(
    x = "Número de alvos",
    y = "Frequência"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 14L),
    plot.subtitle = element_text(size = 14L),
    axis.title.y = element_text(size = 14L),
    axis.title.x = element_text(size = 14L),
    axis.text.y = element_text(size = 14L),
    axis.text.x = element_text(size = 14L),
    legend.text = element_text(size = 14L),
    legend.title = element_text(size = 14L)
    )

ggsave("dist.png", width = 3, height = 2)

