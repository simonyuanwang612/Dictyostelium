library("sleuth")

sample_id <- dir(file.path("results"))
sample_id
kal_dirs <- file.path("results", sample_id)
kal_dirs
s2c <- read.table(file.path("metadata", "Dicty_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
kal_dirs <- file.path("results", sample_id, "kallisto")
s2c <- read.table(file.path("metadata", "Dicty_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
