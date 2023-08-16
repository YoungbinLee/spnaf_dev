Boot <- function(rf, R, n.cores){

    cat("(3) Do bootstrapping ... \n")

    G <- rf$Gij
    Gset <- data.frame(no = as.numeric(row.names(rf)), G = G)
    Gset_list <- split(Gset, Gset$no)

    boot <- function(l){
        sample <- c(l$G, sample(x = G[G != l$G], R-1, replace = T))
        sample[is.na(sample)] <- -100
        sorted <- sort(sample, decreasing = T)
        cr <- which(sorted == sample[1])/R
        return(mean(cr))
    }


    # use multi-core computation
    if(n.cores > 1){
        cl <- parallel::makeCluster(n.cores)
        cat(paste0("using ", n.cores, " cores of the machine which has total ", parallel::detectCores(), " cores \n"))

        parallel::clusterExport(cl, varlist = c("Gset_list", "boot", "R", "G"), envir = environment())
        pval <- parallel::parLapply(cl = cl, Gset_list, fun = boot)
        parallel::stopCluster(cl)
    }else{
        pval <- lapply(Gset_list, FUN = boot)
    }

    ## create results
    pval <- do.call("c", pval)

    ## append a pval column to the result data.frame
    result <- cbind(rf, pval)

    cat("Done ! \n")

    return(result)
}
