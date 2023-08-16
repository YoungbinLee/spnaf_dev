#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows

Gstat <- function(SpatialWeights, method, n.cores){

    cat("(2) Calculating Gij ... \n")

    # t: n^2 = rows of Union
    t <- nrow(SpatialWeights)
    r_bar <- sum(SpatialWeights$n)/t
    s_sq <- sum((SpatialWeights$n - r_bar)**2)/(t-1)
    s <- sqrt(s_sq)

    SpatialWeightsl <- SpatialWeights %>%
        dplyr::mutate(seq = paste(oid, did, sep="-"))

    SWL <- split(SpatialWeightsl, f = SpatialWeightsl$seq)

    oid <- did <- w <- NULL

    subframe <- function(l, ref = SpatialWeights, m = method){
        o <- l$oid
        d <- l$did

        if(m == 't'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            origins <- unique(c(o, origins)) # include o
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- unique(c(d, destinations)) # include d
        }else if(m == 'o'){
            origins <- o
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- unique(c(d, destinations)) # include d
        }else if(m == 'd'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            origins <- unique(c(o, origins)) # include o
            destinations <- d
        }

        ## Merge valid networks
        set1 <- ref %>%
            dplyr::filter(oid %in% origins, did == d)
        set2 <- ref %>%
            dplyr::filter(oid == o, did %in% destinations)
        set <- rbind(set1, set2) %>%
            dplyr::distinct()

        ## Calculate Wij* for each case
        if(m == 't'){
            Wij_star <- length(origins) + length(destinations) -1 # remove duplicated i --> j
        }else if(m == 'o'){
            Wij_star <- length(origins)
        }else{
            Wij_star <- length(destinations)
        }

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

        l$Gij <- numerator/denominator

        return(l)
    }

    # use multi-core computation
    if(n.cores > 1){
        cl <- parallel::makeCluster(n.cores)
        cat(paste0("using ", n.cores, " cores of the machine which has total ", parallel::detectCores(), " cores \n"))

        parallel::clusterExport(cl, varlist = c("t", "s_sq", "s", "r_bar", "SWL", "oid", "did","w", "method", "SpatialWeights", "subframe"), envir = environment())
        parallel::clusterEvalQ(cl, library("dplyr"))

        G <- parallel::parLapply(cl = cl, SWL, fun = subframe)
        parallel::stopCluster(cl)
    }else{
        G <- lapply(SWL, FUN = subframe)
    }

    ## create result in data.frame type
    G <- do.call("bind_rows", G)

    cat("Done! \n")
    # cat(paste0("note: Total ", nrow(SpatialWeightsl), " networks are considered (zeros are excluded)\n"))

    return(G)
}
