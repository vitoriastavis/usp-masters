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
    cols = -`Análise`,
    names_to = "DEGs",
    values_to = "Número_DEGs"
  )

# Abrir no esquisse para criar plots interativos
esquisse::esquisser(dados_long)

p <- ggplot(dados_long) +
  aes(x = Análise, y = Número_DEGs, fill = DEGs) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = c(DEGs_comum = "#A665B5",
               DEGs_GSE144136 = "#367B85",
               DEGs_GSE213982 = "#87CEE0")
  ) +
  labs(
    title = "Tabela 1: DEGs encontrados nas análises de expressão diferencial para tipos celulares e pseudo-bulk",
    subtitle = " Valores para normalização SCTransform,
    p_adj <= 0.05,
    logFC >= 0.6",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 14L),
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    axis.text.y = element_text(size = 12L),
    axis.text.x = element_text(size = 12L),
    legend.text = element_text(size = 12L),
    legend.title = element_text(size = 12L)
  )


library(ggplot2)

pp <- ggplot(dados_long) +
 aes(x = Análise, y = Número_DEGs, fill = DEGs) +
 geom_col(position = "dodge") +
 scale_fill_manual(values = c(Comum = "#A665B5", GSE144136 = "#6FAF9C", GSE213982 = "#87CEE0")) +
 theme_classic() +
 theme(legend.position = "top", axis.title.y = element_text(size = 14L), axis.title.x = element_text(size = 14L), 
 axis.text.y = element_text(size = 14L), axis.text.x = element_text(size = 14L), legend.text = element_text(size = 14L))
ggsave("degs.png")


load("intersection_all_degs.RData")
all_degs_df <- table(lengths(all_degs))
all_degs_df_long <- data.frame(
  N_alvos = names(all_degs_df),
  Contagem = as.integer(all_degs_df),
  row.names = NULL
)

esquisse::esquisser(all_degs_df_long)

p_dist <- ggplot(all_degs_df_long) +
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
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    axis.text.y = element_text(size = 12L),
    axis.text.x = element_text(size = 12L),
    legend.text = element_text(size = 12L),
    legend.title = element_text(size = 12L)
    )

ggsave("dist.png", width = 3, height = 2)

