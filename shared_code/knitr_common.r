library(yaml)

options(knitr.project_name = "ChIP-nexus methods paper")

figure_path <- function(filename="") {
  file.path(getOption("knitr.figure_dir"), filename)
}

data_path <- function(filename="") {
  file.path(yaml.load_file("_output.yaml")$data$base_path, filename)
}

# Default xtable output
options(xtable.type = 'html')

# Output HTML table (requires xtable package)
html_table <- function(df, digits=2, row.names=FALSE, col.names=TRUE) {
  if(class(df)[1] != "xtable") df <- xtable(df)
  digits(df) <- digits
  print(df, include.rownames=row.names, include.colnames=col.names)
}

# Format number with commas
pn <- function(i, ...) {
  prettyNum(i, big.mark=",", ...)
}

# Wrap output and code
options(width=80)

# Force knitr to stop evaluation when an error is encountered
knitr::opts_chunk$set(error=FALSE)

# Don't show code blocks by default
knitr::opts_chunk$set(echo=FALSE)

# Don't reformat R code
knitr::opts_chunk$set(tidy=FALSE)

# Set up figure defaults 
knitr::opts_chunk$set(fig.width=7, fig.height=5, fig.path=figure_path())

# Create output directory if it doesn't exist
if(!file.exists(getOption("knitr.figure_dir"))) dir.create(getOption("knitr.figure_dir"))

source("shared_code/caching.r")
source("shared_code/parallel.r")
