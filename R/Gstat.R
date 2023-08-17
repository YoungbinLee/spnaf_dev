#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom rlang .data

Gstat <- function(od, shape, spatialweights, method, n.cores, R){

    cat("(2) Calculating Gij ... \n")

    ## Union is created by input polygons^2
    U <- merge(shape$id, shape$id)
    names(U) <- c("oid", "did")
    U <- U %>% dplyr::filter(.data$oid != .data$did) # remove flow from i to i

    ### Join OD data to Union
    U <- U %>%
        dplyr::left_join(
            od, by = c("oid" = "oid", "did" = "did")) %>%
        dplyr::mutate(n = ifelse(is.na(.data$n), 0, .data$n)) # fill NA flows with zeros

    ## Alarm the size of the union set
    cat(paste0("Note: Total ",
               nrow(U),
               " network combinations are ready to be analyzed which is calculated by the input polygon data\n"))

    ## create list to iteration
    union_list <- split(U, f = paste0(U$oid, "-", U$did))

    ## calculate fixed values
    # t: n^2 = rows of Union
    t <- nrow(U)
    r_bar <- sum(U$n)/t
    s_sq <- sum((U$n - r_bar)**2)/(t-1)
    s <- sqrt(s_sq)

    # iterate for every flow
    subframe <- function(l,  shp = shape, ref = spatialweights, m = method, union = U){

        o <- l$oid
        d <- l$did

        if(m == 't'){
            ## origins
            l$origin_neighbor <- paste(shp$id[ref[[o]]], collapse = " ")
            origins <- c(o, shp$id[ref[[o]]])
            ## destinations
            l$destination_neighbor <- paste(shp$id[ref[[d]]], collapse = " ")
            destinations <- c(d, shp$id[ref[[d]]])
        }else if(m == 'o'){
            ## origins
            l$origin_neighbor <- o
            origins <- o
            ## destinations
            l$destination_neighbor <- paste(shp$id[ref[[d]]], collpase = " ")
            destinations <- c(d, shp$id[ref[[d]]])
        }else if(m == 'd'){
            ## origins
            l$origin_neighbor <- paste(shp$id[ref[[o]]], collpase = " ")
            origins <- c(o, shp$id[ref[[o]]])
            ## destinations
            l$destination_neighbor <- d
            destinations <- d
        }

        ## Merge valid networks
        set <- union %>%
            dplyr::filter(.data$oid %in% origins & .data$did == d |
                          .data$oid == o, .data$did %in% destinations)

        ## Calculate Wij*
        Wij_star <-  nrow(set) # = length(origins) + length(destinations) -1
        ## Wij*^2
        Wij_star_sq <- Wij_star**2
        ## S1
        S1 <- Wij_star
        ## Gij calculation
        sigma <- sum(set$n, na.rm = T)
        ## numerator
        numerator <- sigma - Wij_star*r_bar
        ## denominator
        denominator <- s * sqrt((t*S1 - Wij_star_sq)/(t-1))

        l$Gij <- numerator/denominator # calculated Gij


        ## create conditional permutation
        pseudo_Gij <- vector(mode = "numeric", R)
        set.seed(23)
        for(i in 1:R){
            keep <- set[set$oid == o & set$did ==d,]
            shuffle <- set[!(set$oid == o & set$did ==d),]
            shuffle$n <- sample(U[U$oid != o & U$did !=d, "n"],
                                replace = FALSE, size = nrow(shuffle))
            set <- rbind(keep, shuffle)
            sigma <- sum(set$n, na.rm = T) # recalculation of sum
            ## numerator
            numerator <- sigma - Wij_star*r_bar # recalculation of numerator
            ## denominator
            denominator <- s * sqrt((t*S1 - Wij_star_sq)/(t-1)) # recalculation of denominator

            ## add value to pseudo_Gij
            pseudo_Gij[i] <- numerator/denominator
        }

        ## create pseudo p-value by conditional permutation
        Gij_set <- c(l$Gij, pseudo_Gij)

        ## positive p-value
        l$pval_positive <- (sum(Gij_set[-1] >= Gij_set[1]) + 1)/(R+1)

        ## negative p-value
        l$pval_negative <- (sum(Gij_set[-1] < Gij_set[1]) + 1)/(R+1)

        return(l)
    }

    # use multi-core computation
    if(n.cores > 1){
        cl <- parallel::makeCluster(n.cores)
        cat(paste0("NOTE: using ", n.cores, " cores of the machine which has total ",
                   parallel::detectCores(), " cores \n"))

        parallel::clusterExport(cl,
                                varlist = c("union_list", "shape", "spatialweights",
                                            "U", "method", "n.cores", "R"),
                                envir = environment())
        parallel::clusterEvalQ(cl, library("dplyr"))

        G <- parallel::parLapply(cl = cl, X = union_list, fun = subframe)

        parallel::stopCluster(cl)
    }else{
        G <- lapply(union_list, FUN = subframe)
    }

    ## create result in data.frame type
    G <- do.call("bind_rows", G)

    cat("Done! \n")
    # cat(paste0("note: Total ", nrow(SpatialWeightsl), " networks are considered (zeros are excluded)\n"))

    return(G)
}

