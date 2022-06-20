Boot <- function(rf, R){
    G <- rf$Gij
    Gset <- data.frame(no = as.numeric(row.names(rf)), G = G)
    Gset_list <- split(Gset, Gset$no)

    boot <- function(l){
        sample <- c(l$G, sample(x = G[G != l$G], R-1, replace = T))
        sample[is.na(sample)] <- -100
        sorted <- sort(sample, decreasing = T)
        cr <- which(sorted == sample[1])/R
        return(cr)
    }

    pval <- do.call("c", lapply(Gset_list, boot))

    result <- cbind(rf, pval)

    return(result)
}
