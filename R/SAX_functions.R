#' @include SAX_functions.R
NULL

# Get the SAX sequence of a given curve
.Func.SAX2 <- function (x, w, a, eps, norm) 
{
    sd_x <- sd(x)
    if (sd_x <= eps) {
        sym <- rep(letters[round((1 + a)/2, digits = 0)], w)
    }
    else {
        if (norm == TRUE) {
            data.nor <- as.vector(scale(x))
        }
        else {
            data.nor <- x
        }
        
        ind <- round(seq(from = 1, to = length(data.nor), length.out = w + 1), digits = 0)
        
        pieces <- rep(FALSE,length(ind)-1)
        
        for(i in 1:(length(ind) - 1)){
            if (i != (length(ind) - 1)) {
                piece <- data.nor[ind[i]:(ind[i + 1] - 1)]
            }
            else {
                piece <- data.nor[ind[i]:ind[i + 1]]
            }
            pieces[i] <- mean(piece, na.rm = T)
        }
        let <- letters[1:a]
        bks <- round(qnorm(p = seq(from = 0, to = 1, length.out = a + 
                                       1)), digits = 2)
        sym <- sapply(pieces, function(obs){
            let[max(which(bks < obs))]
        })
    }
    return(sym)
}

# Calculates the distance between 2 SAX representations
.Func.dist2 <- function (x, y, mat, n){

    w <- length(x)
    d <- sapply(1:length(x), function(i){mat[x[i], y[i]]})
    return(sqrt(sum(d^2)) * sqrt(n/w))
}

.Func.matrix2 <- function (a){
    # Get segments in quantile normal distribution
    alphabet <- letters[1:a]
    segments_pos <- seq(from = 0, to = 1, length.out = (a +  1))
    
    # Get block sizes
    blocks <- round(qnorm(p = segments_pos), digits = 2)[-c(1, a+1)]
    
    # Create empty matrix
    dist_mat <- matrix(data = NA, nrow = a, ncol = a, dimnames = list(alphabet, alphabet))
    
    # Fill matrix with corresponding values
    for (i in 1:a){
        for (j in 1:a){
            dist_mat[i, j] <- ifelse(abs(i - j) <= 1, 0, blocks[max(c(i, j)-1)] - blocks[min(c(i, j))])
        }
    }
    return(dist_mat)
}

# Generates a consensus sequence between multiple SAX strings of the same length,
# handling the letters as if they were numbers
pseudo_consensus <- function(SAX_mat, ...){
    alph <- letters
    num <- 1:26
    names(num) <- alph
    sax_nums <- apply(SAX_mat, 2, function(x) num[x])
    return(alph[round(colMeans(sax_nums))])
}

# Generates a consensus sequence, based on majority vote. 
# Tie handler is either minimum, maximum or mean of the present SAX sequence
consensus <- function(SAX_mat, tie_handler=c("mean","min","max")){
    
    # Variables for SAX to numbers and vice versa
    tie_handler <- tie_handler[1]
    alph <- letters
    num <- 1:26
    names(num) <- alph
    
    # Go over SAX matrix and get consensus letter for each position
    apply(SAX_mat,2,function(x){
        tab_x <- table(x)
        if(tie_handler == "mean"){
            maxima_val  <- num[names(tab_x[which(tab_x == max(tab_x))])]
            mean_maxima <- num[round(mean(maxima_val))]
            return(names(mean_maxima))
        } else{
            return(names(tab_x[getMethod(tie_handler)(which(tab_x == max(tab_x)))]))
        }
    })
}

# Function that scores SAX by mean distance from each sequence to the other
mean_distance_sax <- function(SAXtable, consensus_fun, EIClength, training_index, dist_mat){
    score_mat <- matrix(0,nrow(SAXtable),nrow(SAXtable))
    
    for(s1 in 1:(nrow(score_mat) - 1)){
        s2_s <- s1 + 1
        for(s2 in s2_s:nrow(score_mat)){
            score_mat[s1,s2] <- .Func.dist2(
                x = SAXtable[s1,], y = SAXtable[s2,], 
                mat = dist_mat, n = EIClength)
            score_mat[s2,s1]  <- score_mat[s1,s2] 
        }
    }
    
    mean_p  <-  mean(as.dist(score_mat[training_index, training_index]))
    sd_p    <-    sd(as.dist(score_mat[training_index, training_index]))
    
    return(list(score_mat, mean_p, sd_p))
}

# Function that scores SAX by mean distance from each sequence to a consensus sequence
consensus_distance_sax <- function(SAXtable, consensus_fun, EIClength, training_index, dist_mat){
    
    consensus_sequence <- consensus_fun(SAXtable[training_index,,drop=FALSE], "mean")
    
    score_vec <- vector()
    
    for(s1 in 1:(nrow(SAXtable))){
            score_vec[s1] <- .Func.dist2(
                x = SAXtable[s1,], y = consensus_sequence, 
                mat = dist_mat, n = EIClength)
    }
    
    mean_p  <-  mean(score_vec[training_index])
    sd_p    <-    sd(score_vec[training_index])
    
    return(list(as.matrix(score_vec), mean_p, sd_p))
}

classify_peaks <- function(SAX, a = 4, w = 7, single_matrix=FALSE){
    
    
}

align_EIC <- function(EIClist){
    
    # Go through all peaks and align the EICs
    for(i in 1:length(EIClist[[1]])){
        n_list <- lapply(EIClist, function(x) x[[i]])
        fftXcor()
        A <- sapply(n_list,ncol)
        lowestIndex    <- which(A == min(A))
        notLowestIndex <- which(A != min(A))
        mean_RT <- colMeans(do.call(rbind,lapply(n_list[lowestIndex], function(x) x[1,])))

    }
}

# Fast fourier transform cross correlation 
fftXcor <- function(x, y) {
    
    # Length of vector x
    n <- length(x)
    
    # Enlarge with 0's to size 2*n-1 to account for periodicity of convolution
    x <- c(x, rep(0, length(x) - 1))
    y <- c(y, rep(0, length(y) - 1))
    
    # Cross-correlation via convolution
    crosscor <- fft(Conj(fft(x)) * fft(y), inverse=T) / length(x)
    crosscor <- Re(crosscor) / (n - 1)
    
    # Slice it up to make it for lags -n:n not 0:(2n-1)
    crosscor <- c(crosscor[(n+1):length(crosscor)], crosscor[1:n])
    
    # Store lag as names attribute of vector
    names(crosscor) <- (1-n):(n-1)
    return(crosscor)
}


# Alternative: findShiftStepFFT from speaq
# library(speaq)
# 
# system.time(for(i in 1:1000){
#     x <- findShiftStepFFT(allEICs[[16]][[313]][2,],allEICs[[31]][[313]][2,-29])
# })
