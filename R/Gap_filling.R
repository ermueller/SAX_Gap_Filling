#' @include Gap_filling.R
NULL

#' @import progress
#' @import data.table
#' @importFrom tools file_path_sans_ext
#' @importFrom xcms xcmsRaw
#' @importFrom xcms rawEIC
NULL

#' @title Fill MS-gaps using SAX
#' 
#' @description Function that fills gaps utilizing SAX (Symbol Aggregate approXimation)
#' 
#' @param peak_table A mandatory character vector giving the path to an aligned peak table in csv format or a data.table that is the aligned peak table. See the Details section for more information on how this parameter and the filenames parameter interact.
#' @param mzML_filenames The names of the .mzML-mzML_filenames that are references in peak table. See the Details section for more information
#' @param SAXGAP_params A list containing the main parameters for the SAXGAP function. See the Details section for more information
#' @param output_file An optional charcter vector giving the filepath of where the output table should be stored.
#' 
#' @details
#' \subsection{peak_table and mzML_filename:}{
#' 
#' The peak_table and mzML_filenames parameter need to have the same underlying names for the storage of the relationship of the intensities and retention times to the files. 
#' If, for example, one supplies the file \code{Sample_1.mzML} then the peak_table MUST contain a column named either \code{Intensity_Sample_1} or \code{Intensity_Sample_1.mzML}. The holds true for
#' retention time, i.e. the table needs to have a column named either \code{RT_Sample_1.mzML} or \code{RT_Sample_1}. The peak table also needs general \code{mz} and \code{RT} columns. All retention times must provided in seconds!
#' 
#' }
#' 
#' \subsection{SAXGAP_params:}{
#' The SAXGAP_params paarmeter must be a named list containing the parameters:
#' \describe{
#'  \item{a}{The alphabet size of the SAX conversion and comparison}
#'  \item{w}{The word length for peaks}
#'  \item{t}{The length of the extracted EICs}
#'  \item{pe}{The percentage cutoff for peak removal}
#'  \item{dppm}{The instrument-specific allowed ppm error of the peaks}
#'  \item{intensity_threshold}{The minimum intensity for peaks}
#' }
#' }
#' 
#' @return The function returns the gap filled table with the gap filled peaks stored as a negative value and an added \code{Removal_Candidate} column that signifies whether a peak was likely picked in error. 
#' If output_file is set it writes the same table to the path specified by output_file.
#' 
#' @export
SAXGAP <- function(peak_table, mzML_filenames, SAXGAP_params = list(a=6,w=6,t=14,pe=0.9,intensity_threshold = 1e4,dppm=5), output_file){
    
    # Parameter check and read section

        # Check the peak_table parameter
        # Is it a filepath? If it is, read the file
        if(is.character(peak_table)){
            if(file.exists(peak_table)){
                peak_table <- fread(peak_table)
            } else{
                stop(paste("File",peak_table,"does not exist"))
            }
        } else{
            if(!is.data.table(peak_table)){
                stop("peak_table must be an object of class data.table")
            }
        }
        
        # Check the mzML_filenames parameter
        # Do all filenames exist?
        if(is.character(mzML_filenames)){
            mzML_filenames_exist <- file.exists(mzML_filenames)
            
            if(!all(mzML_filenames_exist)){
                not_existing_mzML_files <- mzML_filenames[!mzML_filenames_exist]
                if(length(not_existing_mzML_files) == 1){
                    stop(paste0("File does not exist:\n", not_existing_mzML_files))
                }
                if(length(not_existing_mzML_files) > 1){
                    stop(paste0("Files do not exist:\n", paste(not_existing_mzML_files,collapse = "\n")))
                }
            } 
            
        } else{
            stop("mzML_filenames must be a character vector containg filenames")
        }
        
        # Check the SAXGAP_params parameter
        # Is SAXGAP_params a list with the proper named elements?
        
        if(!is.list(SAXGAP_params) || !all(c("a","w","t","pe","intensity_threshold","dppm") %in% names(SAXGAP_params))){
            stop("SAXGAP_params must be a list containing named elements \"a\",\"w\",\"t\",\"pe\",\"intensity_threshold\",\"dppm\"")
        }
        
        # Check the output_file parameter (if it is supplied)
        # Does the folder exist where the table should go?
        if(!missing(output_file)){
            dirname_output_file <- dirname(output_file)
            if(!dir.exists(dirname_output_file)){
                stop(paste("The target directory of the output file", dirname_output_file, "does not exist"))
            }
        }
    
    # Get intensity and rt column names
    INT_column_names  <- paste0("Intensity_",file_path_sans_ext(basename(mzML_filenames)))
    RT_column_names   <- paste0("RT_",file_path_sans_ext(basename(mzML_filenames)))
    
    # If the column names are not found in the table check if they are found if we don't remove the file extension from the name
    # It's easy to check this, so we allow this with file extensions and without, for user convenience.
    
    if(!any(INT_column_names %in% colnames(peak_table))){
        INT_column_names  <- paste0("Intensity_",basename(mzML_filenames))
    }
    if(!any(RT_column_names %in% colnames(peak_table))){
        RT_column_names  <- paste0("RT_",basename(mzML_filenames))
    }
    
    # Get intensity table and retention time table
    INT_Table  <- peak_table[,  INT_column_names,with=FALSE]
    RT_Table   <- peak_table[,   RT_column_names,with=FALSE]
    nPeaks <- nrow(INT_Table)
    
    # Check if the intensity table has a column for every file
    if(ncol(INT_Table) != length(INT_column_names)){
        stop("The peak_table does not contain a properly named intensity column for every provided mzML_filename\nThe columns must have the same name as the mzML-files with \"Intensity_\" added in front")
    }
    
    # Check if the retention time table has a column for every file
    if(ncol(RT_Table) != length(RT_column_names)){
        stop("The peak_table does not contain a properly named retention time column for every provided mzML_filename\nThe columns must have the same name as the mzML-files with \"RT_\" added in front")
    }
    
    if(!all(c("mz","RT") %in% colnames(peak_table))){
        stop("The peak_table does not contain mandatory columns with names \"mz\",\"RT\"")
    }
    
    if(all(peak_table$RT < 60)){
        warning("All retention times are smaller than 60, please note that retention times must be supplied in seconds")
    }
    
    # Set up variables for EIC extraction
    SAX <- list()
    detected_peaks      <- list() # Which peaks have already been found?
    detected_peaks_l    <- matrix(FALSE, nPeaks, length(mzML_filenames)) # Which peaks have already been found (but as a logical vector)
    
    # Much faster with a matrix than with a data.table, so a one time conversion just for this purpose
    INT_matrix <- as.matrix(INT_Table)
    
    for(i in 1:nPeaks){
        SAX[[i]]        <- vector()
        detected_peaks_l[i,] <- INT_matrix[i,] > 0
        detected_peaks[[i]] <- which(detected_peaks_l[i,])
    }
    
    ### Set up the EIC extraction and SAX conversion loop
    
    # progress_bar stuff
    length_progressbar <- length(mzML_filenames)*nPeaks
    pb <- progress_bar$new(
        format = " [NonTplus] Generating SAX sequences    [:bar] :percent ETA: :eta Elapsed: :elapsed ",
        total = length_progressbar, clear = FALSE, width = 80)
    
    # Variable setup, create new data.tables and EIClengths
    i <- 0
    EIClengths <- matrix(0,nPeaks,length(mzML_filenames))
    INT_dt_new <- copy(INT_Table)
    RT_dt_new  <- data.table(matrix(0,nrow(RT_Table),ncol(RT_Table)))
    colnames(RT_dt_new) <- RT_column_names
    
    # Calculate the ppmlimits of the mz of each peak and get the SAX string for empty EICs
    mzranges <- t(sapply(peak_table$mz,.get_ppmlimits,SAXGAP_params$dppm))
    empty_SAX <- paste0(rep("a",SAXGAP_params$w),collapse="")
    
    # Go through each file
    for(f in mzML_filenames){
        
        # Read the file
        suppressMessages(XR <- xcmsRaw(f))
        
        # Go through each peak
        for(p in 1:nPeaks){
            
            # Only update the progress bar every 1000 peaks, or else the updating takes a noticable amount of time. 
            if(p %% 1000 == 0){
                pb$update((i*nPeaks+p)/length_progressbar)
            }

            
            # Get the mzrange for each peak
            current_mzrange <- mzranges[p,]
            
            # Get the RT - either it is peak-specific or we take the "general" supplied RT
            if(detected_peaks_l[p,i+1]){
                RT <- RT_Table[[i+1]][p]
            } else{
                RT <- peak_table$RT[p]
            }
            
            # Retrieve the currently relevant EIC from the file
            tempEIC <- rawEIC(XR,mzrange=current_mzrange,rtrange = c(RT-2*SAXGAP_params$t,RT+2*SAXGAP_params$t))
            
            # Move the retention time window to center around the weighted means of intensities - i.e.
            # "Center the middle of the peak if there is one"
            intensities <- tempEIC$intensity
            retention_t <- XR@scantime[tempEIC$scan]
            ind <- .Internal(which(retention_t > RT - SAXGAP_params$t & retention_t < RT + SAXGAP_params$t))
            new_RT <- weighted.mean(retention_t[ind],intensities[ind])
            
            if(!is.nan(new_RT)){
                ind2 <- which(retention_t > new_RT - SAXGAP_params$t & retention_t < new_RT + SAXGAP_params$t)
                intensities <- intensities[ind2]
                retention_t <- retention_t[ind2]
            } else{
                new_RT <- 0
                ind <- .Internal(which(retention_t > RT - SAXGAP_params$t & retention_t < RT + SAXGAP_params$t))
                intensities <- intensities[ind]
            }
            
            sumint <- sum(intensities)
            
            # Save the new Intensities and RT in their tables if we assume a peak to be here
            # If there is a peak, get the new RT with the weighted mean strategy
            if(detected_peaks_l[p,i+1]){
                set(RT_dt_new,p,RT_column_names[i+1], new_RT)
            } else{
                set(INT_dt_new,p,INT_column_names[i+1], sumint)
                set(RT_dt_new,p,RT_column_names[i+1], new_RT)
            }
            
            # If the intensity is below the threshold, pretend that it is empty, otherwise create a SAX string
            if(sumint < SAXGAP_params$intensity_threshold){
                SAX[[p]][i+1] <- empty_SAX
            } else{
                SAX[[p]][i+1] <- paste0(.Func.SAX2(intensities, SAXGAP_params$w, SAXGAP_params$a, 0, T), collapse="")
            }
            
            # The EIC length is technically constant, but SAX is dependent on the number of data points and 
            # in MS experiments the same timerange can have a different number of datapoints at different RTs
            # because scans are not exactly taken at a constant rate but at an only *somewhat* constant rate
            EIClengths[p,i+1] <- length(intensities)
        }
        
        # Clear the memory
        gc(reset = T,full = T)
        i <- i + 1
    }
    
    # Finish up the progress bar in case it is not, because we only update it every 1000 steps
    if(!pb$finished){
        pb$update(1)
    }
    
    # Get a list of all SAX strings that belong to a "detected peak"
    valid_SAX <- list()
    for(p in 1:length(SAX)){
        valid_SAX[[p]] <- SAX[[p]][detected_peaks[[p]]]
    }
    
    # Reduce the "valid" SAX strings by the accepted percentage
    valid_SAX_table <- sort(table(do.call(c,valid_SAX)),decreasing = T)
    valid_SAX_table_reduced <- valid_SAX_table[1:which(cumsum(valid_SAX_table)/sum(valid_SAX_table) > SAXGAP_params$p)[1]]
    valid_SAX_sequences <- names(valid_SAX_table_reduced)
    
    # Setup for finding the consensus sequences
    consensus_mean <- list()
    consensus_majority <- list()
    length_progressbar <- nPeaks
    pb <- progress_bar$new(
        format = " [NonTplus] Finding consensus sequences [:bar] :percent ETA: :eta Elapsed: :elapsed",
        total = length_progressbar, clear = FALSE, width = 80)
    
    # Go through each peak
    for(p in 1:nPeaks){
        
        # If more than one peak has been detected, get the consensus
        if(length(detected_peaks[[p]]) > 1){
            SAX_correct <- do.call(rbind,strsplit(SAX[[p]][detected_peaks[[p]]],split = ""))
            consensus_mean[[p]] <- pseudo_consensus(SAX_correct)
            consensus_majority[[p]] <- consensus(SAX_correct)
        } else{
        # If only one was detected, it IS the consensus
            consensus_mean[[p]]     <- strsplit(SAX[[p]][detected_peaks[[p]]],split = "")
            consensus_majority[[p]] <- strsplit(SAX[[p]][detected_peaks[[p]]],split = "")
        }
        
        if(p %% 100 == 0){
            pb$update(p/length_progressbar)
        }
    }
    
    # Finish progress bar
    if(!pb$finished){
        pb$update(1)
    }
    
    
    
    # Get distance matrix
    dist_mat <- .Func.matrix2(SAXGAP_params$a)
    
    # Set up variables for classification
    classification_function <- consensus_distance_sax
    consensus_fun <- consensus
    mean_s <- vector()
    sd_s <- vector()
    cutoff <- vector()
    A <- vector()
    count <- 0
    classFunc <- median
    res <- list()
    classes <- list()
    single_score_mats <- list()
    length_progressbar <- nrow(INT_Table)
    pb <- progress_bar$new(
        format = " [NonTplus] Adding new peaks            [:bar] :percent ETA: :eta Elapsed: :elapsed",
        total = length_progressbar, clear = FALSE, width = 80)
    
    # Classification loop
    
    rempeaks <- rep(FALSE,times=nPeaks)
    Removed_Peaks <- logical(length = nPeaks)
    for(p in 1:nPeaks){
        
        # Get the index of the filex where the peak has been detected
        training_index <- detected_peaks[[p]]
        
        # In case something with the peak picking has gone horribly wrong (Sometimes the RT doesn't fit)
        
        empty_index <- which(SAX[[p]][training_index] == empty_SAX)
        if(length(empty_index) != 0){
            training_index <- training_index[-empty_index]
        }
        
        test_index <- which(!detected_peaks_l[p,])
        if(length(test_index) == 0){
            next
        }
        
        # Convert the SAX sequences to a table
        SAXtable <- do.call(rbind,strsplit(SAX[[p]], split=""))
        
        # Differentiate method based on whether the current peak is low-coverage or not
        
        if(length(training_index) < 0.05 * length(mzML_filenames)){

            # Low-coverage: See if the low-coverage SAX sequences are plausibly peaks
            
            current_SAX <- SAX[[p]][training_index]
            current_SAX_test <- SAX[[p]][test_index]
            
            if(!any(current_SAX %in% valid_SAX_sequences)){
                
                # If not, remove the peak
                rempeaks[p] <- TRUE
                count <- count + 1
                Removed_Peaks[p] <- TRUE
                classes[[p]] <- rep(FALSE,length(mzML_filenames))
                
            }else{
                
                # If they are plausibly peaks, just look which other peaks match any peak of the valid ones (using the SAX distance function)
                # See if the not yet classified peaks have the same sequence as one of the found peaks
                
                SAX_matches <- which(sapply(current_SAX_test, function(SAX_test_sequence){
                    any(sapply(current_SAX, function(SAX_training_sequence){
                        .Func.dist2(strsplit(SAX_test_sequence, split="")[[1]], strsplit(SAX_training_sequence, split="")[[1]], dist_mat,n = mean(EIClengths[p,]))
                    }) == 0)
                }))
                
                classes[[p]] <- rep(FALSE,length(mzML_filenames))
                if(length(SAX_matches) != 0){
                    classes[[p]][test_index][SAX_matches] <- TRUE
                }
                
            }
        } else{
            
            # High-coverage: Mintest with a consensus sequence
            res[[p]] <- classification_function(SAXtable, consensus_fun, mean(EIClengths[p,]), training_index, dist_mat)
            
            # Get mean distances and standard deviations
            single_score_mats[[p]] <- res[[p]][[1]]
            mean_s <- res[[p]][[2]]
            sd_s <- res[[p]][[3]]
            
            # If the mean is 0
            if(is.nan(mean_s)){
                mean_s <- 1
                sd_s <- 1
            }
            
            if(is.na(sd_s)){
                sd_s <- 0
            }
            
            # The cutoff is the maximum distance of one of training peaks
            cutoff[p] <- max(single_score_mats[[p]][detected_peaks[[p]]])
            A[p] <- all(single_score_mats[[p]][detected_peaks[[p]]] <= cutoff[p])
            classes[[p]] <- (apply(single_score_mats[[p]],1, classFunc) <= cutoff[p]) & !detected_peaks_l[p,]
        }
        
        # Set the new intensities and retention times in the new tables
        for(s in which(classes[[p]])){
            set(peak_table,p,INT_column_names[s],-INT_dt_new[[s]][p])
            set(peak_table,p,RT_column_names[s],RT_dt_new[[s]][p])
        }
        
        if(p %% 100){
            pb$update(p/nPeaks)
        }
    }
    
    
    if(!pb$finished){
        pb$update(1)
    }

    # Clean up
    peak_table[, Removal_Candidate := rempeaks]
    
    if(!missing(output_file)){
        fwrite(peak_table,file = output_file)
    }
    return(peak_table)
}