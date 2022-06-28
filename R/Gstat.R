#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows

Gstat <- function(SpatialWeights, method, star = TRUE){

    cat("(2) Calculating Gij ... ")

    # t: n^2 = rows of Union
    t <- nrow(SpatialWeights)
    r_bar <- sum(SpatialWeights$n)/t
    s_sq <- sum((SpatialWeights$n - r_bar)**2)/(t-1)
    s <- sqrt(s_sq)

    SpatialWeightsl <- SpatialWeights %>%
        dplyr::mutate(seq = paste(oid, did, sep="-"))
    SWL <- split(SpatialWeightsl,
                 f = SpatialWeightsl$seq)

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
        }else{
            stop("method must be one of c('t', 'o', 'd') \n")
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
    G <- do.call("bind_rows", lapply(SWL, subframe))

    cat("Done! \n")

    return(G)
}
