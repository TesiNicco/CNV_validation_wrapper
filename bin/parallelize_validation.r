# Function to guide parallelization of the validation_CNV.R script

# Libraries
    library(data.table)
    library(stringr)

# Set working directory
    setwd('/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/')

# Data paths
    input_path <- '/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/input_data'
    snps_path <- paste0(input_path, '/SNPs_input.txt')
    cnvs_path <- paste0(input_path, '/CNVs_info_pennCNV.txt')
    samples_path <- paste0(input_path, '/IDs_pacbio_passed_sampleQC_updated.txt')
    model_path <- str_replace_all(input_path, 'input_data', 'trained_models/joint.rds')
    output_path <- '/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/results/'

# Read data
    snps <- fread(snps_path, header = TRUE, stringsAsFactors = FALSE)
    cnvs <- fread(cnvs_path, header = TRUE, stringsAsFactors = FALSE)
    samples <- fread(samples_path, header = TRUE, stringsAsFactors = FALSE)

# Define parallelization parameters -- number of batches = number of cpus
    n_cpu = 6
    # Create batches of same number of samples per batch
    samples_per_batch <- ceiling(nrow(samples) / n_cpu)
    # create n_cpu lists of samples
    samples_list <- split(samples, rep(1:n_cpu, each = samples_per_batch, length.out = nrow(samples)))
    # Create folders, one for each batch in the output path and place the samples in the corresponding batch folder
    for (i in 1:n_cpu) {
        batch_path <- paste0(output_path, 'batch_', i)
        dir.create(batch_path, recursive = TRUE)
        write.table(samples_list[[i]], file = paste0(batch_path, '/samples.txt'), row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

# Finally run the validation_CNV.R script in parallel
    # iterate over the batches and submit the jobs
    for (i in 1:n_cpu) {
        cmd = paste0('Rscript bin/validation_CNV.R --batch ', i, ' --snps ', snps_path, ' --cnvs ', cnvs_path, ' --samples ', paste0(output_path, 'batch_', i, '/samples.txt'), ' --model ', model_path, ' --outdir ', paste0(output_path, 'batch_', i))
        message(paste("Running batch", i))
        system(cmd, wait = FALSE)
    }
    # Define expected output flag files
    expected_files <- paste0(output_path, 'batch_', 1:n_cpu, '/predictions.txt')
    # Wait for completion
    while (TRUE) {
        done <- file.exists(expected_files)
        message(sum(done), " of ", n_cpu, " batches complete")
        if (all(done)) break
        Sys.sleep(5)
    }
    # Combine results from all batches
    all_preds <- lapply(expected_files, fread)
    final_preds <- rbindlist(all_preds)
    # Write final predictions to a single file
    final_output_path <- paste0(output_path, 'final_predictions.txt')
    write.table(final_preds, final_output_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    message("All batches completed and results combined into ", final_output_path)

