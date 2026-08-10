\name{differential_region_analysis_wilcoxon}
\alias{differential_region_analysis_wilcoxon}
\title{Subject-level Wilcoxon sensitivity analysis across spatial regions}
\description{
Aggregates spatial samples at an explicit inference unit and performs fixed-design
paired signed-rank or unpaired rank-sum sensitivity analyses with effect sizes,
multiple-testing correction, and auditable missingness diagnostics.
}
\usage{
differential_region_analysis_wilcoxon(
  sample_matrix,
  group_column = "roi_id",
  subject_column = NULL,
  section_column = NULL,
  reference_group = NULL,
  aggregate_fun = c("mean", "median"),
  exact_method = c("auto", "exact", "asymptotic"),
  min_replicates = 5L,
  confidence_level = 0.95,
  p_adjust_method = "BH"
)
}
\arguments{
  \item{sample_matrix}{A data frame containing \code{sample_id}, a region
  column, and numeric \code{mz_} feature columns.}
  \item{group_column}{Name of the region or condition column.}
  \item{subject_column}{Optional biological-subject identifier. When supplied,
  each contrast is paired on the fixed intersection of subjects.}
  \item{section_column}{Optional section identifier used only when no biological
  subject identifier is available. Section-level inference is flagged as
  potentially pseudoreplicated.}
  \item{reference_group}{Optional reference group. If omitted, all pairwise
  contrasts are evaluated.}
  \item{aggregate_fun}{Function used to summarize multiple rows within each
  group and replicate: \code{"mean"} or \code{"median"}.}
  \item{exact_method}{Whether to select exact tests automatically, require an
  exact test, or use the asymptotic approximation. Exact tests are considered
  eligible only for small data without relevant ties or zero paired differences.}
  \item{min_replicates}{Minimum complete pairs, or minimum observations in each
  unpaired group. The default is five.}
  \item{confidence_level}{Confidence level for the location-shift interval.}
  \item{p_adjust_method}{Method passed to \code{p.adjust}. Adjustment is
  performed within each group contrast across features.}
}
\value{
A data frame with descriptive medians, the Hodges--Lehmann location-shift
estimate and interval when available, rank-biserial correlation, p-value and
within-contrast FDR, replication counts, tie and zero-difference diagnostics,
and explicit fit or skip status. Positive effects mean that \code{group_b}
tends to have larger values than \code{group_a}.
}
\details{
Rows are first aggregated within group and inference unit. Pairing is fixed at
the contrast level and is never changed to an unpaired test because a feature
has missing values. When no subject identifier is supplied, results are marked
as potentially pseudoreplicated and should be treated as exploratory.

The rank-sum test is a test of ranks or distributions; the reported difference
in medians is descriptive. A simple location-shift interpretation of the
Hodges--Lehmann estimate needs the corresponding distributional assumptions.
The paired signed-rank test assumes symmetry of paired differences. Agreement
with a t-test or spatial mixed model is a sensitivity analysis, not independent
validation.
}
