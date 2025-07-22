# Script to prepare data for CNV validation

# Make sure tabix and bgzip are installed and accessible in the system -- can be through a conda environment, or globally installed
    # Load conda environment before starting the script -- required as bgzip and tabix are already installed in the conda environment
    conda activate newtreat
    # Alternative is to install tabix and bgzip through conda, in a new or existing environment
    # conda env list
    # conda activate {env_name}
    # conda install -c bioconda tabix
    # conda install -c bioconda bgzip
#####################################

# Make sure the required packages are installed
    required_packages <- c("BiocParallel", "torch", "luz", 'BiocManager', 'devtools', 'argparse', 'data.table', 'stringr')
    new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
    # Check if BiocManager is installed
    if ("BiocManager" %in% new_packages) { install.packages("BiocManager") }
    # Then load BiocManager package
    library(BiocManager)
    # Install the remaining packages through BiocManager
    if (length(new_packages) >0) { BiocManager::install(new_packages) }
    # Load the required packages to check everything is OK
    library(BiocParallel)
    library(torch)
    library(luz)
    library(data.table)
    library(stringr)
    library(devtools)
    # Then install the tool for the validation if it's not installed
    if ('CNValidatron' %in% rownames(installed.packages()) == FALSE){
        devtools::install_github("sinomem/CNValidatron_fl")
    }
###################################

# Make sure the data is in the right format
    # Define original samples path
    samples_path <- '/project/holstegelab/Share/joan/CNVs/data/IDs_pacbio_passed_sampleQC.txt'
    snps_path <- '/project/holstegelab/Share/gwas_array/Raw_arraydata_CNVs/IQ_dementia/data/Input_PennCNV_GSAV1_customcontent/CNVinput.Index'
    cnvs_path <- '/project/holstegelab/Share/gwas_array/Raw_arraydata_CNVs/IQ_dementia/data/Output_PennCNV_GSAV1_customcontent/GSA1_customcontent_20210824.rawcnv'
    # Define new data path where to place reformatted data and the tabix-indexed files
    new_data_path_raw <- "/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/input_data/raw_intensity_files"
    new_data_path <- "/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/input_data/tabix_indexed_files"
    # Define path to the input data for SNPs, CNVs and samples
    outdir = '/project/holstegelab/Share/nicco/collaborations/Joan/CNV_validation/input_data'
    # Create the directory if it does not exist
    if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
    if (!dir.exists(new_data_path)) { dir.create(new_data_path, recursive = TRUE)}
    if (!dir.exists(new_data_path_raw)) { dir.create(new_data_path_raw, recursive = TRUE)}
    # Read the samples
    samples = read.table(samples_path, header = TRUE, stringsAsFactors = FALSE)
    # Add tabix column to samples
    samples$file_path_tabix = NA
    # Iterate through the samples, and copy the data path to the tabix_indexed_files in the same directory
    for (i in 1:nrow(samples)){
        cat(sprintf("Processing sample: %s\r", i))
        flush.console()
        # Get the path to the sample
        tmp_path = sub("'", "", samples$file_path[i])
        # Get the base name of the file
        base_name = basename(tmp_path)
        # Check if the new file exists
        if (file.exists(paste0(new_data_path, '/', base_name, '.gz'))) {
            cat(sprintf("File %s already exists in %s. Skipping...\n", base_name, new_data_path))
            samples$file_path_tabix[i] = paste0(new_data_path, '/', base_name, '.gz')
            samples$file_path[i] = paste0(new_data_path_raw, '/', base_name)
            next
        }
        # Read the file
        tmp = fread(tmp_path, h=T, stringsAsFactors=F)
        # Update the column names
        colnames(tmp) <- c("Name", "Chr", "Start", "Type", "B Allele Freq", "Log R Ratio")
        tmp$End = tmp$Start
        # Reorder
        tmp = tmp[, c('Chr', 'Start', 'End', 'Log R Ratio', 'B Allele Freq')]
        # Add a LRR adjusted
        tmp$"Log R Ratio Adjusted" = tmp$"Log R Ratio"
        # Make sure chromosome and position are numeric
        tmp$Chr = as.numeric(tmp$Chr)
        tmp$Start = as.numeric(tmp$Start)
        tmp$End = as.numeric(tmp$End)
        # Exclude NAs
        tmp = tmp[!is.na(tmp$Start),]
        # Sort by chromosome, then position
        tmp = tmp[order(tmp$Chr, tmp$Start), ]
        # Write the file back
        write.table(tmp, paste0(new_data_path_raw, '/', base_name), row.names = FALSE, col.names = TRUE, quote = FALSE, sep= "\t")
        # Compress with bgzip
        system(paste0('bgzip -c -f ', paste0(new_data_path_raw, '/', base_name), ' > ', paste0(new_data_path, '/', base_name, '.gz')))
        # Index with tabix
        system(paste0('tabix -S 1 -s 1 -b 2 -e 3 ', paste0(new_data_path, '/', base_name, '.gz')))
        # Udate the file path in the samples file
        samples$file_path_tabix[i] = paste0(new_data_path, '/', base_name, '.gz')
        samples$file_path[i] = paste0(new_data_path_raw, '/', base_name)
    }
    # Write the updated samples file
    outpath = paste0(outdir, '/IDs_pacbio_passed_sampleQC_updated.txt')
    write.table(samples, outpath, row.names = FALSE, col.names = TRUE, quote = FALSE, sep= "\t")

    # Read the SNP file
    snps = fread(snps_path, h=T, stringsAsFactors=F)
    head(snps)
    snps$Index = NULL
    write.table(snps, paste0(outdir, '/SNPs_input.txt'), quote=F, row.names=F, sep="\t")

    # Finally the CNV file
    cnvs = fread(cnvs_path, h=F, stringsAsFactors=F)
    # they need to be: sample_ID, chr, start, end, GT, CN, numsnp
    cnvs_parsed = data.table(sample_ID = str_replace_all(str_replace_all(cnvs$V5, 'Input_PennCNV_GSAV1_customcontent/CNVinput.', ''), '.gz', ''), chr = str_replace_all(str_split_fixed(cnvs$V1, ':', 2)[, 1], 'chr', ''), start = as.integer(sub(".*:(\\d+)-.*", "\\1", cnvs$V1)), end = as.integer(sub(".*-(\\d+)", "\\1", cnvs$V1)), CN = str_split_fixed(str_split_fixed(cnvs$V4, ',', 2)[, 2], '=', 2)[, 2], numsnp = str_split_fixed(cnvs$V2, '=', 2)[, 2])
    # add GT: 1 when CN is 0 or 1, 2 otherwise
    cnvs_parsed$GT = ifelse(cnvs_parsed$CN %in% c(0, 1), 1, 2)
    # reorder columns
    cnvs_parsed = cnvs_parsed[, c('sample_ID', 'chr', 'start', 'end', 'GT', 'CN', 'numsnp')]
    # remove CN and nums
    cnvs_parsed$chrom = as.integer(cnvs_parsed$chr)
    write.table(cnvs_parsed, paste0(outdir, '/CNVs_info_pennCNV.txt'), quote=F, row.names=F, sep="\t")
###################################