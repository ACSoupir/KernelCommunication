#' Clear Cell Renal Cell Carcinoma CosMx Data
#'
#' A core from a TMA of clear cell renal cell carcinoma samples subjected to single cell spatial transcriptomics with CosMx SMI.
#'
#' Coordinates are in \strong{millimetres}, spanning roughly 1.02 x 1.02 mm with
#' a median nearest-neighbour distance of 0.0073 mm. Bandwidths for this dataset
#' therefore belong in the range 0.01 to 0.25, which corresponds to the 10-250
#' micrometre range used in the companion manuscript. Passing bandwidths of 10
#' or more drives every pairwise weight to approximately 1, making the kernel a
#' constant matrix and collapsing Moran's I toward zero.
#'
#' @name ccRCCExample
#' @docType data
#' @keywords datasets
#' @format ## `ccRCCExample`
#' A data frame with 7,122 rows and 12 columns:
#' \describe{
#'   \item{Slide}{Our slide identifier}
#'   \item{x}{centroid location in the x direction - in mm}
#'   \item{y}{centroid location in the y direction - in mm}
#'   \item{cell_id}{unique identifier for each cell}
#'   \item{CD274}{SCTransform scale.data (residuals) of CD274 (slide normalized)}
#'   \item{PDCD1}{SCTransform scale.data (residuals) of PDCD1 (slide normalized)}
#'   \item{slide}{run slide number}
#'   \item{fov}{slide FOV number}
#'   \item{cell}{FOV cell number}
#'   \item{unique_fov}{unique identifier for CosMx SMI imaging FOV}
#'   \item{tma_core}{combined FOV number relating back to the TMA core}
#'   \item{unique_core}{unique identifier for the TMA core given the slide and core}
#' }
#' @source Mitchell Hayes
#' @usage data(ccRCCExample)
NULL
