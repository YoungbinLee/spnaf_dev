#' @importFrom magrittr %>%

SpatialWeight <- function(df, shape, snap, queen){

    w <- NULL
    nb <- unclass(spdep::poly2nb(shape, snap, queen))
    original_id <- shape$id
    nb <- lapply(nb, function(data){
        return(original_id[data])
    })
    names(nb) <- original_id

    did <- stack(nb)
    did <- setNames(did$values, did$ind)
    oid <- names(did)
    nb_frame <- data.frame(oid = oid)
    nb_frame$did <- did
    nb_frame$w <- 1

    result <- nb_frame %>%
        dplyr::left_join(df, by = c("oid" = "oid", "did" = "did")) %>%
        dplyr::mutate(w = ifelse(oid == did, 1, w)) %>%
        dplyr::arrange(oid, did)
    ## Union is created by input polygons^2
    U <- merge(shape$id, shape$id)
    names(U) <- c("oid", "did")
    ### Join OD data to Union
    result <- U %>%
        dplyr::left_join(result, by = c("oid" = "oid", "did" = "did")) %>%
        dplyr::mutate(n = ifelse(is.na(n), 0, n),
                      w = ifelse(oid == did, 1, w),
                      w = ifelse(is.na(w), 0, 1))

    return(result)

}
