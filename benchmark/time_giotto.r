library(Giotto)

set_up_giotto <- function(work_dir, instrs) {
    
    exp_path = fs::path(work_dir, "exp.txt")
    loc_path = fs::path(work_dir, "loc.txt")
    meta_path = fs::path(work_dir, "meta.txt")
    
    obj <- createGiottoObject(raw_exprs = exp_path,spatial_locs = loc_path,instructions = instrs)
    metadata = data.table::fread(file = meta_path)
    obj <- addCellMetadata(obj, new_metadata = metadata, by_column = T, column_cell_ID = 'cell_id')
    
    return(obj)
}

cell_num = c('1000', '5000', '10000', '50000', '100000')

nns_time = list(software=c(),cell_num=c(),exec_time=c())
cci_time = list(software=c(),cell_num=c(),exec_time=c())

python_path = "/usr/bin/python3"
results_folder = 'result/'

instrs = createGiottoInstructions(save_plot = FALSE,
                                  show_plot = FALSE,
                                  save_dir = results_folder,
                                  python_path = python_path)

for (c in cell_num) {
    write(sprintf("Running for %s \n", c), stdout())
    work_dir = sprintf('data_%s', c)

    obj = set_up_giotto(work_dir, instrs)

    for (i in 1:3) {
        start_nns = as.numeric(as.POSIXct(Sys.time()))
        temp_obj <- createSpatialKNNnetwork(gobject = obj, k=5, name = "knn_network",)
        end_nns = as.numeric(as.POSIXct(Sys.time()))
        temp_obj <- cellProximityEnrichment(gobject = temp_obj,
                                           cluster_column = 'cell_type',
                                            spatial_network_name = "knn_network",
                                           number_of_simulations = 500)
        end_cci = as.numeric(as.POSIXct(Sys.time()))


        write(sprintf("nns time is %s \n", end_nns - start_nns), stdout())
        nns_time$software <- append(nns_time$software, 'Giotto')
        nns_time$cell_num <- append(nns_time$cell_num, c)
        nns_time$exec_time <- append(nns_time$exec_time, end_nns - start_nns)

        write(sprintf("cci time is %s \n", end_cci - end_nns), stdout())

        cci_time$software <- append(cci_time$software, 'Giotto')
        cci_time$cell_num <- append(cci_time$cell_num, c)
        cci_time$exec_time <- append(cci_time$exec_time, end_cci - end_nns)
    }

}

nns_time = data.frame(nns_time)
cci_time = data.frame(cci_time)
write.csv(nns_time, "giotto_nns_time.csv", quote=FALSE, row.names = FALSE)
write.csv(cci_time, "giotto_cci_time.csv", quote=FALSE, row.names = FALSE)