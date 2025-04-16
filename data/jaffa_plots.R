# JAFFA fusion gene visualization functions

# Load necessary packages
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(RCircos)
    library(gtools)
    library(dplyr)
})

# ---------------------- JAFFA Data Import Function ----------------------

import_jaffa_data <- function(filename, genome_version, limit) {
    # Validate genome version
    valid_genomes <- c("hg19", "hg38", "mm10")
    if (is.na(match(tolower(genome_version), tolower(valid_genomes)))) {
        stop("Invalid genome version. Please use: hg19, hg38 or mm10")
    }
    
    # Validate limit parameter
    if (!missing(limit) && (!is.numeric(limit) || limit <= 0)) {
        stop("limit must be a numeric value greater than 0")
    }
    
    # Read JAFFA result file
    col_types <- c(
        "sample" = "character",
        "fusion genes" = "character",
        "chrom1" = "character",
        "base1" = "integer",
        "chrom2" = "character",
        "base2" = "integer",
        "spanning pairs" = "character",
        "spanning reads" = "integer",
        "inframe" = "character",
        "classification" = "character",
        "nb_tag" = "character"  # Add nb_tag column type
    )
    
    # First read one row to check columns
    header <- read.csv(filename, nrows = 1)
    # Filter out non-existent columns
    col_types <- col_types[names(col_types) %in% names(header)]
    
    report <- if (missing(limit)) {
        data.table::fread(
            input = filename,
            colClasses = col_types,
            showProgress = FALSE
        )
    } else {
        data.table::fread(
            input = filename,
            colClasses = col_types,
            showProgress = FALSE,
            nrows = limit
        )
    }
    
    # Create fusion gene list
    fusion_list <- vector("list", dim(report)[1])
    
    # Process each row of data, create fusion gene entry
    for (i in 1:dim(report)[1]) {
        # Gene names
        gene_names <- unlist(strsplit(report[[i, "fusion genes"]], split = ":"))
        
        # Determine if reading frame is maintained
        if (report[[i, "inframe"]]=="True") {
            inframe <- TRUE
        } else if (report[[i, "inframe"]]=="False") {
            inframe <- FALSE
        } else {
            inframe <- NA
        }
        
        # Supporting reads count
        split_reads_count <- as.numeric(report[[i, "spanning reads"]])
        spanning_reads_count <- tryCatch(
            as.numeric(report[[i, "spanning pairs"]]),
            warning = function(w) { 0 }
        )
        
        # Get nb_tag if exists
        nb_tag <- NULL
        if ("nb_tag" %in% colnames(report)) {
            nb_tag <- report[[i, "nb_tag"]]
        }
        
        # Create simplified fusion event list
        fusion_list[[i]] <- list(
            id = as.character(i),
            fusion_tool = "jaffa",
            genome_version = genome_version,
            
            # Upstream gene information
            gene_upstream = list(
                name = gene_names[1],
                chromosome = report[[i, "chrom1"]],
                breakpoint = report[[i, "base1"]],
                strand = "*"  # JAFFA does not provide strand direction
            ),
            
            # Downstream gene information
            gene_downstream = list(
                name = gene_names[2],
                chromosome = report[[i, "chrom2"]],
                breakpoint = report[[i, "base2"]],
                strand = "*"  # JAFFA does not provide strand direction
            ),
            
            # Supporting evidence for fusion event
            spanning_reads_count = spanning_reads_count,
            split_reads_count = split_reads_count,
            inframe = inframe,
            
            # JAFFA specific data
            classification = report[[i, "classification"]],
            nb_tag = nb_tag  # Add nb_tag field
        )
    }
    
    return(fusion_list)
}

# ---------------------- Numeric Scaling Helper Function ----------------------

scale_list_to_interval <- function(the_list, new_min, new_max) {
    if (length(the_list) <= 1) {
        stop("Invalid list. This function requires at least two values.")
    }
    (new_max - new_min) * (the_list - min(the_list)) / (max(the_list) - min(the_list)) + new_min
}

# ---------------------- Create Link Data ----------------------

create_link_data <- function(fusion_list, min_link_width = 1, max_link_width = 10) {
    chromosome <- vector(mode = "character", length = length(fusion_list))
    chrom_start <- vector(mode = "numeric", length = length(fusion_list))
    chrom_end <- vector(mode = "numeric", length = length(fusion_list))
    
    chromosome_1 <- vector(mode = "character", length = length(fusion_list))
    chrom_start_1 <- vector(mode = "numeric", length = length(fusion_list))
    chrom_end_1 <- vector(mode = "numeric", length = length(fusion_list))
    
    link_width <- vector(mode = "numeric", length = length(fusion_list))
    
    for (i in seq_along(fusion_list)) {
        fusion <- fusion_list[[i]]
        
        # Upstream gene data
        chromosome[[i]] <- fusion$gene_upstream$chromosome
        chrom_start[[i]] <- fusion$gene_upstream$breakpoint
        chrom_end[[i]] <- fusion$gene_upstream$breakpoint + 1
        
        # Downstream gene data
        chromosome_1[[i]] <- fusion$gene_downstream$chromosome
        chrom_start_1[[i]] <- fusion$gene_downstream$breakpoint
        chrom_end_1[[i]] <- fusion$gene_downstream$breakpoint + 1
        
        # Link width = total supporting reads for fusion
        link_width[[i]] <- fusion$spanning_reads_count + fusion$split_reads_count
    }
    
    # Normalize link widths
    if (length(link_width) > 1) {
        link_width <- scale_list_to_interval(link_width, min_link_width, max_link_width)
    } else {
        link_width[[1]] <- max_link_width
    }
    
    data.frame(chromosome,chrom_start,chrom_end,chromosome_1,chrom_start_1,chrom_end_1,link_width)
}

# ---------------------- Create Gene Label Data ----------------------

create_gene_label_data <- function(fusion_list) {
    original_length <- length(fusion_list)
    new_length <- original_length * 2
    
    chromosome <- vector(mode = "character", length = new_length)
    chrom_start <- vector(mode = "numeric", length = new_length)
    chrom_end <- vector(mode = "numeric", length = new_length)
    gene <- vector(mode = "character", length = new_length)
    
    for (i in seq_along(fusion_list)) {
        fusion <- fusion_list[[i]]
        
        # Upstream gene label data
        chromosome[[i]] <- fusion$gene_upstream$chromosome
        chrom_start[[i]] <- fusion$gene_upstream$breakpoint
        chrom_end[[i]] <- fusion$gene_upstream$breakpoint + 1
        gene[[i]] <- fusion$gene_upstream$name
        
        # Downstream gene label data
        chromosome[[i + original_length]] <- fusion$gene_downstream$chromosome
        chrom_start[[i + original_length]] <- fusion$gene_downstream$breakpoint
        chrom_end[[i + original_length]] <- fusion$gene_downstream$breakpoint + 1
        gene[[i + original_length]] <- fusion$gene_downstream$name
    }
    
    data.frame(chromosome,chrom_start,chrom_end,gene)
}

# ---------------------- Create Simple Ideogram ----------------------
create_simple_ideogram <- function(chromosomes) {
    n_chr <- length(chromosomes)
    # Create simple chromosome dataframe
    data.frame(
        Chromosome = chromosomes,
        ChromStart = rep(0, n_chr),
        ChromEnd = rep(100000000, n_chr),  # Simple chromosome length
        Band = paste0(chromosomes, "p"), # Simple band label
        Stain = rep("gpos50", n_chr)  # Default staining
    )
}

# ---------------------- Circos Plot Function ----------------------
plot_fusion_circos <- function(csv_file, genome_version="hg38", output_file=NULL, cytoband_path=NULL) {
    # Import JAFFA data
    fusion_list <- import_jaffa_data(csv_file, genome_version)
    
    # Determine cytogenetic data file path
    if (!is.null(cytoband_path) && dir.exists(cytoband_path)) {
        # Use provided path
        if (genome_version == "hg19") {
            cytoband_file <- file.path(cytoband_path, "UCSC.HG19.Human.CytoBandIdeogram.txt")
        } else if (genome_version == "hg38") {
            cytoband_file <- file.path(cytoband_path, "UCSC.HG38.Human.CytoBandIdeogram.txt")
        } else if (genome_version == "mm10") {
            cytoband_file <- file.path(cytoband_path, "UCSC.MM10.Mus.musculus.CytoBandIdeogram.txt")
        }
    }
    
    # Check if file exists
    if (!file.exists(cytoband_file)) {
        warning(paste("Cannot find", cytoband_file, ", will use simplified chromosome structure"))
        # Create simplified chromosome structure data
        if (genome_version == "hg19" || genome_version == "hg38") {
            chromosomes <- paste0("chr", c(1:22, "X", "Y"))
            cytoband <- create_simple_ideogram(chromosomes)
        } else if (genome_version == "mm10") {
            chromosomes <- paste0("chr", c(1:19, "X", "Y"))
            cytoband <- create_simple_ideogram(chromosomes)
        }
    } else {
        # Read cytogenetic data
        cytoband <- utils::read.table(cytoband_file)
    }
    
    # Set column names expected by RCircos
    names(cytoband) <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")
    
    # Add RCircos.Env object to global environment
    assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)
    
    # Sort cytogenetic data to make chromosomes appear in order
    cytoband <- RCircos::RCircos.Sort.Genomic.Data(
        genomic.data = cytoband,
        is.ideo = TRUE
    )
    
    # Initialize plotting device
    if (!is.null(output_file)) {
        pdf(output_file, width = 10, height = 10)
    }
    
    # Initialize components
    cyto_info <- cytoband
    chr_exclude <- NULL
    tracks_inside <- 3
    tracks_outside <- 0
    RCircos::RCircos.Set.Core.Components(
        cyto_info,
        chr_exclude,
        tracks_inside,
        tracks_outside
    )
    
    # Set plot area
    RCircos::RCircos.Set.Plot.Area()
    
    # Draw chromosome ideogram
    RCircos::RCircos.Chromosome.Ideogram.Plot()
    
    # Create gene label data
    gene_label_data <- create_gene_label_data(fusion_list)
    
    # Draw gene name connectors
    track_num <- 1
    side <- "in"
    RCircos::RCircos.Gene.Connector.Plot(
        gene_label_data,
        track_num,
        side
    )
    
    # Draw gene names
    track_num <- 2
    name_col <- 4
    RCircos::RCircos.Gene.Name.Plot(
        gene_label_data,
        name_col,
        track_num,
        side
    )
    
    # Create link data
    link_data <- create_link_data(fusion_list)
    
    # Ensure correct sorting
    multi_mixedorder <- function(..., na.last = TRUE, decreasing = FALSE) {
        do.call(
            order,
            c(
                lapply(
                    list(...),
                    function(l) {
                        if (is.character(l)) {
                            factor(
                                l,
                                levels = gtools::mixedsort(unique(l))
                            )
                        } else {
                            l
                        }
                    }
                ),
                list(na.last = na.last, decreasing = decreasing)
            )
        )
    }
    
    ordered_link_width <- link_data[
        multi_mixedorder(
            as.character(link_data$chromosome),
            link_data$chrom_start
        ),
    ]$link_width
    
    # Draw link lines
    track_num <- 3
    RCircos::RCircos.Link.Plot(
        link.data = link_data,
        track.num = track_num,
        by.chromosome = TRUE,
        start.pos = NULL,
        genomic.columns = 3,
        is.sorted = FALSE,
        lineWidth = ordered_link_width
    )
    
    # If output to file, close plotting device
    if (!is.null(output_file)) {
        dev.off()
        cat(paste("Circos plot saved to:", output_file, "\n"))
    }
    
    # Remove RCircos.Env object from global environment
    remove("RCircos.Env", envir = .GlobalEnv)
}

# ---------------------- Cell Fusion Gene Count Statistics Plot ----------------------

plot_fusion_per_cell <- function(csv_file, output_file=NULL) {
    # Read JAFFA result file
    dat <- read.csv(csv_file)
    
    # Check if nb_tag column exists
    if (!"nb_tag" %in% colnames(dat)) {
        stop("nb_tag column not found in file")
    }
    
    # Count fusion genes per cell
    cell_counts <- dat %>%   
        dplyr::count(nb_tag) %>%   
        dplyr::rename(Cell_Type = nb_tag, Count = n)
    
    # Set group intervals
    breaks <- c(0,1,2,3,4,5,6,7,9,13,Inf)
    labels <- c("(0,1]", "(1,2]", "(2,3]", "(3,4]", "(4,5]", "(5,6]","(6,7]","(7,9]","(9,13]","13+")
    cell_counts$count_group <- cut(cell_counts$Count, breaks = breaks, labels = labels, include.lowest = FALSE, right = TRUE)
    
    # Create plot
    p <- ggplot(cell_counts, aes(x = count_group, fill = count_group)) +   
        geom_bar(stat = "count") +   
        theme_minimal() +   
        labs(title = "The number of fusion genes detected per cell", x = "", y = "Count") +   
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.4, color = "black", size = 15),
            axis.text.y = element_text(color = "black", size = 15),
            axis.title.y = element_text(color = "black", size = 18)) +   
        scale_fill_brewer(palette = "Set3") +
        guides(fill = FALSE)
    
    # Save chart
    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 10, height = 6)
        cat(paste("Fusion gene cell distribution plot saved to:", output_file, "\n"))
    } else {
        print(p)
    }
}

# ---------------------- Fusion Gene Classification Statistics Plot ----------------------

plot_fusion_classification <- function(csv_file, output_file=NULL) {
    # Read JAFFA result file
    dat <- read.csv(csv_file)
    
    # Create plot
    p <- ggplot(dat, aes(x = classification, fill = classification)) +   
        geom_bar(stat = "count") +   
        theme_minimal() +   
        labs(title = "Classification of fusion genes", x = "", y = "Count") +   
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 15),
            axis.text.y = element_text(color = "black", size = 15),
            axis.title.y = element_text(color = "black", size = 18)) +   
        scale_fill_brewer(palette = "Set3") +
        guides(fill = FALSE)
    
    # Save chart
    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 10, height = 6)
        cat(paste("Fusion gene classification plot saved to:", output_file, "\n"))
    } else {
        print(p)
    }
}

# ---------------------- Command Line Entry Function ----------------------
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Usage: Rscript jaffa_plots.R <function_name> <arguments>")
}

function_name <- args[1]
args <- args[-1]  # Remove function name, leave parameters

if (function_name == "circos") {
    if (length(args) < 1) {
        stop("Usage: Rscript jaffa_plots.R circos <csv_file> [genome_version] [output_file] [cytoband_path]")
    }
    
    csv_file <- args[1]
    genome_version <- if (length(args) >= 2) args[2] else "hg38"
    output_file <- if (length(args) >= 3) args[3] else NULL
    cytoband_path <- if (length(args) >= 4) args[4] else NULL
    
    plot_fusion_circos(csv_file, genome_version, output_file, cytoband_path)
} else if (function_name == "cell_stats") {
    if (length(args) < 1) {
        stop("Usage: Rscript jaffa_plots.R cell_stats <csv_file> [output_file]")
    }
    
    csv_file <- args[1]
    output_file <- if (length(args) >= 2) args[2] else NULL
    
    plot_fusion_per_cell(csv_file, output_file)
} else if (function_name == "class_stats") {
    if (length(args) < 1) {
        stop("Usage: Rscript jaffa_plots.R class_stats <csv_file> [output_file]")
    }
    
    csv_file <- args[1]
    output_file <- if (length(args) >= 2) args[2] else NULL
    
    plot_fusion_classification(csv_file, output_file)
} else {
    stop(paste("Unknown function:", function_name))
}

