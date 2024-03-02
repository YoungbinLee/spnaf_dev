#' @importFrom magrittr %>%
#' @importFrom rlang .data

SpatialWeight <- function(df, shape, snap, queen){

    ## this function will include other types of inputs (e.g. points),
    ## therefore it remains concise for now.

    cat("(1) Creating Spatial Weights... ")

    nb <- unclass(spdep::poly2nb(shape, row.names = NULL, snap = snap, queen = snap)) # create nb data from input shape
    names(nb) <- shape$id

    cat("Done !\n")
    return(nb)
}
