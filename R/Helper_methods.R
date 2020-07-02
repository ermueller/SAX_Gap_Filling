#' @include Helper_methods.R
NULL

# Internal function to get the ppm limits given a mass and delta ppm
.get_ppmlimits <- function (mz, dppm) 
{
    return(c(mz * (1 - dppm * 1e-06), mz * (1 + dppm * 1e-06)))
}