install.packages("fs")
library(fs)
library(Giotto)

args <- commandArgs(trailingOnly = TRUE)
CPUs <- as.integer(args[1])

set_up_giotto <- function(work_dir, instrs) {

  exp_path <- fs::path(work_dir, "exp.txt")
  loc_path <- fs::path(work_dir, "loc.txt")
  meta_path <- fs::path(work_dir, "meta.txt")

  obj <- createGiottoObject(raw_exprs = exp_path, spatial_locs = loc_path, instructions = instrs)
  metadata <- data.table::fread(file = meta_path)
  obj <- addCellMetadata(obj, new_metadata = metadata, by_column = T, column_cell_ID = 'cell_id')

  return(obj)
}

REPEAT <- 3
K <- 10
TIMES <- 1000
heat <- FALSE
dataset <- c('Stereo-seq-MouseEmbryo', 'IMC-BreastCancer')

nns_tb <- list(software = c(), dataset = c(), cpu_count = c(), exec_time = c(), exec_mem = c())
cci_tb <- list(software = c(), dataset = c(), cpu_count = c(), exec_time = c(), exec_mem = c())

python_path <- "/usr/bin/python3"
results_folder <- 'result/'

instructions <- createGiottoInstructions(save_plot = FALSE,
                                         show_plot = FALSE,
                                         save_dir = results_folder,
                                         python_path = python_path)

# Testing Stero-seq dataset
print("Running for Stereo-seq MouseEmbryo dataset", quote = FALSE)
data_dir <- "data/E16-5"
obj <- set_up_giotto(data_dir, instructions)

if (!heat) {
  temp_obj <- createSpatialKNNnetwork(gobject = obj, k = 1, name = "knn_network",)
  temp_obj <- cellProximityEnrichment(gobject = temp_obj,
                                      cluster_column = 'cell_type',
                                      spatial_network_name = "knn_network",
                                      number_of_simulations = 10)
  heat <- TRUE
  print("Finish HEATING!")
}

for (i in 1:REPEAT) {
  sprintf("Running %s Time", i)
  # Neighbor search
  Rprof(tf <- "rprof.log", memory.profiling = TRUE, interval = 0.01)
  obj <- createSpatialKNNnetwork(gobject = obj, k = K, name = "knn_network",)
  Rprof(NULL)
  p <- summaryRprof(tf, memory = "both")
  nns_mem <- max(p$by.total$mem.total)
  nns_time <- max(p$by.total$total.time)


  # Cell Cell interaction
  Rprof(tf <- "rprof.log", memory.profiling = TRUE, interval = 0.01)
  obj <- cellProximityEnrichment(gobject = obj,
                                 cluster_column = 'cell_type',
                                 spatial_network_name = "knn_network",
                                 number_of_simulations = TIMES)
  Rprof(NULL)
  p <- summaryRprof(tf, memory = "both")
  cci_mem <- max(p$by.total$mem.total)
  cci_time <- max(p$by.total$total.time)

  nns_tb$software <- append(nns_tb$software, "Giotto")
  nns_tb$dataset <- append(nns_tb$dataset, dataset[1])
  nns_tb$cpu_count <- append(nns_tb$cpu_count, CPUs)
  nns_tb$exec_time <- append(nns_tb$exec_time, nns_time)
  nns_tb$exec_mem <- append(nns_tb$exec_mem, nns_mem)

  cci_tb$software <- append(cci_tb$software, "Giotto")
  cci_tb$dataset <- append(cci_tb$dataset, dataset[1])
  cci_tb$cpu_count <- append(cci_tb$cpu_count, CPUs)
  cci_tb$exec_time <- append(cci_tb$exec_time, cci_time)
  cci_tb$exec_mem <- append(cci_tb$exec_mem, cci_mem)

}

write.csv(data.frame(nns_tb), sprintf("result/giotto_embryo_nns_%score.csv", CPUs), quote = FALSE, row.names = FALSE)
write.csv(data.frame(cci_tb), sprintf("result/giotto_embryo_cci_%score.csv", CPUs), quote = FALSE, row.names = FALSE)

# Testing IMC Breast Cancer datset, Only run when CPUS >= 8

# if (CPUs >= 8) {
#   roi_list <- list.files(path = "data/IMC-scp-basel", full.names = T)
#   print("Running for IMC Breast Cancer dataset", quote = FALSE)
#   Rprof(tf <- "rprof.log", memory.profiling = TRUE, interval = 0.01)
#   for (roi in roi_list) {
#     obj <- set_up_giotto(roi, instructions)
#     obj <- createSpatialKNNnetwork(gobject = obj, k = K, name = "knn_network")
#     obj <- cellProximityEnrichment(gobject = obj,
#                                    cluster_column = 'cell_type',
#                                    spatial_network_name = "knn_network",
#                                    number_of_simulations = TIMES)
#   }
#   Rprof(NULL)
#   p <- summaryRprof(tf, memory = "both")
#   multi_mem <- max(p$by.total$mem.total)
#   multi_time <- max(p$by.total$total.time)
#   multi_tb <- list(software = c("Giotto"),
#                    dataset = c(dataset[2]),
#                    cpu_count = c(CPUs),
#                    exec_time = c(multi_time),
#                    exec_mem = c(multi_mem))
#
#   write.csv(data.frame(multi_tb), sprintf("result/giotto_multi_all_%score.csv", CPUs), quote = FALSE, row.names = FALSE)
#
# }
