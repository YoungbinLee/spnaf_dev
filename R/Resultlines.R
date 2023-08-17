#' @importFrom magrittr %>%

Resultlines <- function(shape, result_frame){

    cat("(3) Generating result in lines ... ")

    x <- y <- geometry <- coords <- id <- oid <- did <- n <- Gij <- pval_positive <- pval_negative <- NULL
    x.x <- y.x <- x.y <- y.y <- oid_x <- oid_y <- did_x <- did_y <- origin_neighbor <- destination_neighbor <- NULL
    centroids <- suppressWarnings(sf::st_centroid(shape, of_largest_polygon = T)) %>%
        dplyr::mutate(coords = as.character(geometry)) %>%
        as.data.frame() %>%
        dplyr::select(-geometry) %>%
        dplyr::mutate(coords = gsub("c|\\(|\\)", "", coords)) %>%
        tidyr::separate(col = coords, into = c("x", "y"), sep = ",") %>%
        dplyr::select(id, x, y)
    centroids$x <- trimws(centroids$x)
    centroids$y <- trimws(centroids$y)
    result_lines <- result_frame %>%
        dplyr::left_join(centroids, by = c("oid" = "id")) %>%
        dplyr::left_join(centroids, by = c("did" = "id")) %>%
        dplyr::select(oid, did, n, origin_neighbor, destination_neighbor, Gij, pval_positive, pval_negative, oid_x = x.x, oid_y = y.x, did_x = x.y, did_y = y.y) %>%
        dplyr::mutate(sfc = paste0("LINESTRING (", oid_x, " ", oid_y, ", ", did_x, " ", did_y, ")")) %>%
        dplyr::select(-oid_x, -oid_y, -did_x, -did_y) %>%
        sf::st_as_sf(wkt = "sfc")
    sf::st_crs(result_lines) <- sf::st_crs(shape)

    cat("Done !\n")

    return(result_lines)
}
