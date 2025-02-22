% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectFeatures.r
\name{selectFeatures}
\alias{selectFeatures}
\title{Select genes with large absolute fold change for each spatial domain}
\usage{
selectFeatures(spadesign, logfc_cutoff, mean_in_cutoff, max_num_gene)
}
\arguments{
\item{spadesign}{A \code{spaDesign} object}

\item{logfc_cutoff}{Numeric value specifying the threshold for log fold change, used to select genes with significant differential expression within vs outside domain}

\item{mean_in_cutoff}{Numeric value specifying the cutoff for within-domain mean expression (log-transformed), used to select domain-informative genes}

\item{max_num_gene}{Integer specifying the maximum number of genes that can be selected for each spatial domain.}
}
\value{
A object with selected features for each domain
}
\description{
Select genes with large absolute fold change for each spatial domain
}
\examples{
# Assuming you have a spaDesign object named `my_spadesign`
my_spadesign <- selectFeatures(my_spadesign, logfc_cutoff = .7, mean_in_cutoff = 2, max_num_gene = 10)
# You can now access the top genes for each domain using:
top_genes <- my_spadesign@topGenes
}
