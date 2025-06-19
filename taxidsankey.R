library(argparse)
library(networkD3)

# Get the directory of this script, even when run via Rscript
get_script_dir <- function() {
  # Works for Rscript execution
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(dirname(sub("--file=", "", cmdArgs[match]))))
  }
  # Fallback for interactive use
  return(getwd())
}

script_dir <- get_script_dir()
taxdump_dir <- file.path(script_dir, "taxdump")
rankedlineage_path <- file.path(taxdump_dir, "rankedlineage.dmp")

# Ensure rankedlineage.dmp exists, otherwise download and extract it
if (!file.exists(rankedlineage_path)) {
  dir.create(taxdump_dir, showWarnings = FALSE)
  owd <- getwd()
  setwd(taxdump_dir)
  if (!file.exists("new_taxdump.tar.gz")) {
    system("wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz")
  }
  system("tar zxvf new_taxdump.tar.gz")
  setwd(owd)
}

taxidsankey <- function(
  taxids = NULL,
  taxid_file = NULL,
  rankedlineage_file = rankedlineage_path,
  output_html = "sankey_diagram.html",
  ranks = c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
  fill = FALSE
) {
  # Define all possible ranks in order from domain to species
  all_ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # Allow ranks to be specified as numbers (1 = domain, 9 = species)
  if (is.numeric(ranks)) {
    ranks <- all_ranks[ranks]
  }
  
  # Load the rankedlineage.dmp file
  print("Loading rankedlineage.dmp file...")
  rankedlineage <- read.table(
    rankedlineage_file,
    sep = "|",
    quote = "",
    fill = TRUE,
    strip.white = TRUE,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  print("Ranked lineage file loaded.")
  rankedlineage <- rankedlineage[, -11]
  colnames(rankedlineage) <- c(
    "taxid", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "domain"
  )
  
  # Get taxids from file if not provided as vector
  if (is.null(taxids) && !is.null(taxid_file)) {
    taxids <- scan(taxid_file, what = character(), quiet = TRUE)
  }
  if (is.null(taxids)) stop("No taxids provided.")
  
  # Filter
  filtered_rows <- rankedlineage[rankedlineage$taxid %in% taxids, ]
  
  # Build links for the selected ranks
  links_list <- list()
  for (i in 1:nrow(filtered_rows)) {
    row <- filtered_rows[i, ]
    # print(paste("row",row,sep=""))
    if (fill) {
      blank_ranks <- row==""
      newrow <- c(row[blank_ranks], row[!blank_ranks])
      # print(paste("newrow",row,sep=""))
      names(newrow) <- names(row)
      row <- newrow
    }
    lineage <- row[ranks]
    # lineage <- as.character(filtered_rows[i, ranks])
    lineage <- lineage[lineage != ""]
    if (length(lineage) < 2) next
    for (j in 1:(length(lineage) - 1)) {
      links_list[[length(links_list) + 1]] <- data.frame(
        source = lineage[j],
        target = lineage[j + 1],
        value = 1
      )
    }
  }
  if (length(links_list) == 0) stop("No links to plot.")
  links <- do.call(rbind, links_list)
  links_agg <- aggregate(value ~ source + target, data = links, FUN = sum)
  nodes_names <- unique(c(links_agg$source, links_agg$target))
  nodes <- data.frame(name = nodes_names, stringsAsFactors = FALSE)
  links_agg$source <- match(links_agg$source, nodes$name) - 1
  links_agg$target <- match(links_agg$target, nodes$name) - 1
  
  sankey <- sankeyNetwork(
    Links = links_agg,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    fontSize = 30,      # Reduced font size for less overlap
    nodeWidth = 40,     # Slightly wider nodes
    nodePadding = 30,
    # width = 2400,       # Wider diagram
    # height = 900,       # Taller diagram
    iterations = 1000    # More layout iterations for better node placement
    # sinksRight = FALSE # Uncomment if you want to try left/right layout
  )
  saveNetwork(sankey, file = output_html)
}

# argparse CLI
if (interactive() == FALSE) {
  parser <- ArgumentParser(description = "Create a taxonomic Sankey diagram from a list of NCBI taxids.")
  parser$add_argument("-i", "--input", required=TRUE, help="Input file with taxids (one per line)")
  parser$add_argument("-o", "--output", default="sankey_diagram.html", help="Output HTML file [default: sankey_diagram.html]")
  parser$add_argument("-r", "--ranks", default="domain,superkingdom,kingdom,phylum,class,order,family,genus,species",
                      help="Comma-separated list of ranks (by name or number, e.g. 1,5,9 or domain,phylum,species)")
  parser$add_argument("-l", "--lineage", default=rankedlineage_path, help=paste("Ranked lineage file [default:", rankedlineage_path, "]"))
  parser$add_argument("-f", "--fill", action="store_true", help="Remove blank values in lineage (skip missing ranks)")
  args <- parser$parse_args()
  
  # Parse ranks argument
  ranks_arg <- strsplit(args$ranks, ",")[[1]]
  all_ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (all(grepl("^[0-9]+$", ranks_arg))) {
    ranks <- as.numeric(ranks_arg)
  } else {
    ranks <- ranks_arg
  }
  
  taxidsankey(
    taxid_file = args$input,
    output_html = args$output,
    rankedlineage_file = args$lineage,
    ranks = ranks,
    fill = args$fill
  )
}