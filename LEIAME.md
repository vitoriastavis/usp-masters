# Scripts e dados para processar metabólitos da base GutMGene,
# prever seus alvos moleculares e realizar enriquecimento funcional

Primeiro, foi rodado o script target_prediction.Rmd até gerar os SMILES. Como alguns metabólitos não possuem
CID para gerar diretamente o SMILES pelo PubChem, foi obtida uma lista em data/names_ids.txt. Essa lista foi
usada como entrada no script get_smiles.ipynb, que gera os SMILES do OPSIN. Os SMILES foram concatenados na
lista data/smiles_ids.txt e esses dados foram usados para predição nas plataformas SEA e STP. Os resultados
são processados no script target_prediction.Rmd, que gera uma lista da interseção desses alvos. Essa lista é
usada no script enrichment_ora.Rmd, que realiza o enriquecimento por Over-representation analysis com os gene
sets GO-BP, GO-MF, GO-CC, KEGG, reactome e DisGeNET. Esse script gera grafos e gráficos de barra, tanto para
as análises sem filtro, quanto para selecionando apenas categorias que contém termos relacionados ao
Transtorno Depressivo Maior (TDM).

projeto-vitoria
│
├── data/
│   ├── gutmgene/
│   │   ├── gutmgene.csv                   # Tabela da base de dados GutMGene
│   │   ├── gutmgene_na.RData              # GutMGene com apenas NAs
│   │   └── gutmgene_nodup.RData           # GutMGene sem duplicados e sem NAs
│   ├── sea/
│   │   ├── cid.xls                        # Predições de alvos com CID
│   │   ├── nocid.xls                      # Predições de alvos sem CID
│   │   └── sea-results.xls                # Predições consolidadas (CID + noCID)
│   ├── stp/
│   │   ├── *.csv                          # Um arquivo por metabólito (246 no total)
│   │   └── stp-results.csv                # Resultados concatenados do STP (sem processamento)
│   ├── smiles_ids.txt                     # Lista de SMILES e CIDs ou IDs
│   ├── names_ids.txt                      # Nomes e identificadores dos metabólitos sem CID
│   └── targets_intersection.RData         # Interseção dos alvos preditos por SEA e STP
│
├── plots/
│   ├── gutmgene/
│   │   └── microbe_distribution.png       # Histograma de distribuição microbiana
│   └── enrichment-predicted-genes/
│       ├── enrichments.RData              # Objeto em R que contém os enriquecimentos 
│       ├── filtered/                 
│           └── *.png / *.pdf              # Grafos e gráficos de barra (filtrando termos TDM)
│       └── no-filter/                 
│           └── *.png / *.pdf              # Grafos e gráficos de barra sem filtro dos termos
│
├── enrichment_ora.Rmd                     # Enriquecimento funcional (ORA)
├── target_prediction.Rmd                  # Processamento de GutMGene e processamento da predição de alvos
└── get_smiles.ipynb                       # Obtenção de SMILES a partir de nomes/IDs dos metabólitos em CID
