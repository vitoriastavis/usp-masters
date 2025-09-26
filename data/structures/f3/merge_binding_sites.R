library(stringr)
library(bio3d)

############################## BINDING SITES

# PRANK WEB

prankweb_result <- read.csv("binding_sites/prankweb/structure.cif_predictions.csv")
prankweb_residues <- list()
prankweb_residues[[1]] <- unname(as.numeric(
                            vapply(tail(
                            strsplit(prankweb_result$residue_ids, " ")[[1]], -1),
                            function(x) strsplit(x, "_")[[1]][2],
                            FUN.VALUE = character(1))))

# DOGSITE3

dogsite_result <- read.csv("binding_sites/dogsite3/1boy_desc.txt", sep="\t")

directory <- "binding_sites/dogsite3/residues"
files <- list.files(path = directory, pattern = "\\.pdb$", full.names = TRUE)
dogsite_residues <- list()
for (i in 1:length(files)) {
  file_path <- files[i]
  pdb <- read.pdb(file_path)
  
  residues <- unique(pdb$atom[, c("resno")])
  
  dogsite_residues[[i]] <- residues
}

# DOGSITESCORER

dogsitescorer_result <- read.csv("binding_sites/dogsitescorer/1boy_desc.txt", sep="\t")

directory <- "binding_sites/dogsitescorer/residues"
files <- list.files(path = directory, pattern = "\\.pdb$", full.names = TRUE)
dogsitescorer_residues <- list()
for (i in 1:length(files)) {
  file_path <- files[i]
  pdb <- read.pdb(file_path)
  
  residues <- unique(pdb$atom[, c("resno")])
  
  dogsitescorer_residues[[i]] <- residues
}

# CAVITY PLUS

directory <- "binding_sites/cavityplus"
files <- list.files(path = directory, pattern = "^this_cavity", full.names = TRUE)
cavityplus_residues <- list()
for (i in 1:length(files)) {
  file_path <- files[i]
  pdb <- read.pdb(file_path)
  
  residues <- unique(pdb$atom[, c("resno")])
  
  cavityplus_residues[[i]] <- residues
}

# FPOCKET

directory <- "binding_sites/fpocket/pockets"
files <- list.files(path = directory, pattern = "atm\\.cif$", full.names = TRUE)
fpocket_residues <- list()

for (i in 1:length(files)) {
  
  file_path <- files[i]

  atoms <- read.table(file_path, skip = 41, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  fpocket_residues[[i]] <- atoms$V15
}

############################## ALLOSTERIC SITES

# ASD

asd_results <- read.table("allosteric_sites/asd/ASD_Release_202306_AP.txt", sep="\t", header=T)
asd_f3 <- asd_results[asd_results$pdb_id=="1BOY", ]

asd_residues <- lapply(asd_f3$site_residue, function(site) {
  residues <- strsplit(site, split = ";")[[1]]
  positions <- sapply(residues, function(res) {
    strsplit(res, split = ",")[[1]][3]
  })
  as.integer(positions)
})

# ALLO - use results from dog

allo_results <- read.csv("allosteric_sites/allo/1boy_desc_nb_out.txt", header=T)
allo_results_a <- allo_results[allo_results$prediction == "A", 2]

directory <- "binding_sites/dogsitescorer/residues"
files <- list.files(path = directory, pattern = "\\.pdb$", full.names = TRUE)
allo_residues <- list()
for (i in 1:length(files)) {
  file_path <- files[i]
  
  if (any(sapply(allo_results_a, function(key) grepl(key, file_path)))){
    pdb <- read.pdb(file_path)
    residues <- unique(pdb$atom[, c("resno")])
    allo_residues[[i]] <- residues
  }
}


# PASSER

directory <- "allosteric_sites/passer/pockets"
files <- list.files(path = directory, pattern = "atm\\.pdb$", full.names = TRUE)
passer_residues <- list()
for (i in 1:length(files)) {
  file_path <- files[i]
  pdb <- read.pdb(file_path)
  
  residues <- unique(pdb$atom[, c("resno")])
  
  passer_residues[[i]] <- residues
}

# Analyze binding sites

ortho_predictions <- list(prankweb = prankweb_residues,
                        dogsite = dogsite_residues,
                        dogsitescorer = dogsitescorer_residues,
                        cavityplus = cavityplus_residues,
                        fpocket = fpocket_residues)

allo_predictions <- list(asd = asd_residues,
                     allo = allo_residues,
                     passer = passer_residues)

build_consensus_sites <- function(predictions, jaccard_threshold = 0.3){
  
  # Gather info of all sites 
  sites <- data.frame(
    tool = rep(names(predictions), sapply(predictions, length)),
    site_id = unlist(lapply(seq_along(predictions), function(i){
      paste0(names(predictions)[i], "_", seq_along(predictions[[i]]))
    })),
    stringsAsFactors = FALSE
  )
  sites$residues <- unlist(lapply(predictions, function(x) lapply(x, unique)),
                           recursive = FALSE)
  # print(sites)
  # Compare all sites pairwise using Jaccard
  n <- nrow(sites)
  overlaps <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      a <- sites$residues[[i]]
      b <- sites$residues[[j]]
      jac <- length(intersect(a, b)) / length(union(a, b))
      overlaps[i, j] <- jac
      overlaps[j, i] <- jac
    }
  }
  
  # print(overlaps)
  
  # Cluster sites based on threshold
  visited <- rep(FALSE, n)
  clusters <- list()
  
  for (i in 1:n) {
    if (!visited[i]) {
      cluster <- i
      visited[i] <- TRUE
      
      # Find all sites similar to current one
      similar <- which(overlaps[i, ] >= jaccard_threshold & !visited)
      while (length(similar) > 0) {
        cluster <- c(cluster, similar)
        visited[similar] <- TRUE
        new_similar <- unique(unlist(lapply(similar, function(k) {
          which(overlaps[k, ] >= jaccard_threshold & !visited)
        })))
        similar <- new_similar
      }
      
      clusters[[length(clusters) + 1]] <- cluster
    }
  }
  
  # print(clusters)
  
  # Build consensus sites
  consensus <- lapply(clusters, function(idx) {
    list(
      tools = sites$tool[idx],
      site_ids = sites$site_id[idx],
      residues = sort(unique(unlist(sites$residues[idx])))
    )
  })
  
  return(consensus)
}

allosteric_sites <- build_consensus_sites(allo_predictions, jaccard_threshold = 0.2)

allosteric_sites_filtered <- lapply(allosteric_sites, function(residue){
  if (length(residue$tools) > 1){
    residue
} else {NULL}
})