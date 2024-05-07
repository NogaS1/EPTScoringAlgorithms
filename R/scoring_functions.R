library(dplyr)
library(abind)

# Pre-processing ----------------------------------------------------------

# Check data's shape
check_shape <- function(data, blockNum){
  if(!"block" %in% colnames(data)){
    data <- data %>%
      dplyr::group_by(sid) %>%
      dplyr::mutate(block = rep(1:blockNum, length.out = n()))
  }
  data <- data %>% dplyr::select(c("sid", "error", "rt", "prime_cat", "target_cat", "block"))
  
  return(data)
}


# General calculations for algorithms -------------------------------------

# add 600ms penalty to errors
add_penalty <- function(correct, errors){
  # Calculate mean rt of block for each participant's block (after trimming)
  ep_mean_rt <- subset(correct, select = c("sid", "block", "rt"))
  mean_rt <- ep_mean_rt %>% dplyr::group_by(sid, block, .group = TRUE) %>% dplyr::mutate(mean_rt = mean(rt))
  mean_rt <- unique(mean_rt[,c(1,2,5)])
  errors_penalty <- merge(x=errors, y=mean_rt, by.x=c("sid", "block"), by.y = c("sid", "block"), all.x = TRUE)
  errors_penalty <- errors_penalty %>% dplyr::mutate(rt = mean_rt+600)
  errors_penalty <- dplyr::select(errors_penalty, -"mean_rt")
  
  correct <- dplyr::select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  res <- rbind(correct, errors_penalty)          # Add error trials (after penalty) to correct trials trimmed df
  
  return(res)
}

# remove problematic participants
remove_participants <- function(data){
  # Remove participants with less than 2 correct trials per prime-target condition
  correct_trials <- subset(data, error == 0) %>% dplyr::group_by(sid, prime_cat, target_cat)
  no_correct_trials <- unique(as.matrix(data[,1]))[!unique(as.matrix(data[,1])) %in% unique(as.matrix(correct_trials[,1]))]           # remove participants who has only error trials after treatments
  c <- aggregate(cbind(count = rt) ~ sid+prime_cat+target_cat, data = correct_trials, FUN = function(x){NROW(x)})
  
  # less than 2 trials of prime-target
  not_ok <- unique(select(subset(c, count < 2), "sid"))                                                                          
  not_ok_sid <- append(not_ok$sid, no_correct_trials)
  
  # missing prime-target combinations
  missing <- c %>% dplyr::group_by(sid, prime_cat) %>% dplyr::mutate(has_neg = any(target_cat == "neg"), has_pos = any(target_cat == "pos")) %>% dplyr::filter(!(has_neg & has_pos))
  not_ok_sid <- append(missing$sid, not_ok_sid)
  data <- data[!data$sid %in% not_ok_sid,]
  
  # missing blocks
  max_block <- max(data$block)
  all_parts <- list()
  
  # Iterate over each block to get unique sids
  for (i in 1:max_block) {
    all_parts[[i]] <- unique(data[data$block == i, "sid"])
  }
  
  # Find missing parts
  missing_parts <- vector("list", length = max_block)
  for (i in 1:max_block) {
    other_blocks <- setdiff(1:max_block, i)
    missing_parts[[i]] <- lapply(other_blocks, function(j) all_parts[[j]][!all_parts[[j]] %in% all_parts[[i]]])
  }
  
  all_missing_parts <- unique(unlist(unlist(missing_parts)))
  
  # Remove missing
  data <- data[!data$sid %in% all_missing_parts, ]
  
  return(data)
}

# Calculate basic G score
my.Gscores <- function(inID, inDV, inPair, inCond, include){
  dt <- data.frame(sessionId=inID, dv=inDV, pair=inPair, cond=inCond, include=include)
  dt <- dt[which(dt$include),]
  # Comments in "quotes" from Nosek, Bar-Anan, Sriram, Greenwald (2013),
  # "Understanding and using the brief implicit association test: I.
  # recommended scoring procedures". Table 9. http://ssrn.com/abstract=2196002
  
  Gaussianranks <- function(x)
  {
    # It handles NA values by leaving them in the same place as found
    y <- x[!is.na(x)]    
    
    # "1. Assign fractional ranks to N latencies. The longest latency will be
    # assigned a value of 1.0 and the shortest will be assigned a value of 1/N.
    # In the case of ties, ranks are averaged across tied values"
    N <- length(y)
    Fr <- rank(y)/N
    
    # "2. Subtract 1/2N from each fractional rank. Assuming untied values, the
    # largest latency will now have a value of 1-1/2N or (2N-1)/2N.
    # The 1/2N downward adjustment applies generally, even when tied values
    # exist".
    Fr <- Fr - 1/(2*N)
    
    # "3. For each of the N observations, compute the standard normal deviate
    # (mean = 0 and standard deviation = 1) corresponding to the adjusted
    # fractinal rank latency".
    Gr <- scale(Fr)
    
    # Remove the attributes given by the scale() function
    attr(Gr, "scaled:center") <- NULL
    attr(Gr, "scaled:scale") <- NULL
    x[!is.na(x)] <- Gr
    return(x)
  }
  
  # "4. G1 is the mean of the normal deviates in condition 1. G2 is the mean
  # of the normal deviates in condition 2"
  Mranks <- dt %>% 
    dplyr::group_by(sessionId) %>% # here do NOT group by blockcode
    dplyr::mutate(Gr = Gaussianranks(dv)) %>% # compute gaussian ranks
    dplyr::group_by(sessionId, pair, cond) %>% # here do group by blockcode
    dplyr::summarize(Mean = mean(Gr, na.rm = TRUE)) # mean rank in each block
  
  Mranks <- reshape2::dcast(Mranks, sessionId + pair ~ cond, value.var = "Mean")
  
  # "5. G = G2-G1"
  Mranks <- mutate(Mranks, Gscore = incong - cong)
  
  # Put data in wide format
  Gsc <- reshape2::dcast(Mranks, sessionId ~ pair, value.var = "Gscore")
  return(Gsc)
}

# Overall
calc_g_overall <- function(data){
  # Include all trials
  data$include <- TRUE
  
  # G score function computes the G score by comparing the conditions "incong" & "cong"
  data$cond <- ifelse(data$target_cat == "neg", "incong", ifelse(data$target_cat == "pos", "cong", NA))
  
  # Calc G score & fix col names
  overall <- my.Gscores(inID = data$sid, inDV = data$rt, inPair = data$prime_cat, inCond = data$cond, include = data$include)
  
  tmp <- data.frame(sessionId=data$sid, dv=data$rt, pair=data$prime_cat, cond=data$cond, include=data$include)
  # Comments in "quotes" from Nosek, Bar-Anan, Sriram, Greenwald (2013),
  # "Understanding and using the brief implicit association test: I.
  # recommended scoring procedures". Table 9. http://ssrn.com/abstract=2196002
  
  Gaussianranks <- function(x){
    # It handles NA values by leaving them in the same place as found
    y <- x[!is.na(x)]    
    
    # "1. Assign fractional ranks to N latencies. The longest latency will be
    # assigned a value of 1.0 and the shortest will be assigned a value of 1/N.
    # In the case of ties, ranks are averaged across tied values"
    N <- length(y)
    Fr <- rank(y)/N
    
    # "2. Subtract 1/2N from each fractional rank. Assuming untied values, the
    # largest latency will now have a value of 1-1/2N or (2N-1)/2N.
    # The 1/2N downward adjustment applies generally, even when tied values
    # exist".
    Fr <- Fr - 1/(2*N)
    
    # "3. For each of the N observations, compute the standard normal deviate
    # (mean = 0 and standard deviation = 1) corresponding to the adjusted
    # fractinal rank latency".
    Gr <- scale(Fr)
    
    # Remove the attributes given by the scale() function
    attr(Gr, "scaled:center") <- NULL
    attr(Gr, "scaled:scale") <- NULL
    x[!is.na(x)] <- Gr
    return(x)
  }
  
  # "4. G1 is the mean of the normal deviates in condition 1. G2 is the mean
  # of the normal deviates in condition 2"
  Mranks <- tmp %>% 
    dplyr::group_by(sessionId) %>% # here do NOT group by blockcode
    dplyr::mutate(Gr = Gaussianranks(dv)) %>% # compute gaussian ranks
    dplyr::group_by(sessionId, pair, cond) %>% # here do group by blockcode
    dplyr::summarize(Mean = mean(Gr, na.rm = TRUE)) # mean rank in each block
  
  Mranks <- reshape2::dcast(Mranks, sessionId + pair ~ cond, value.var = "Mean")
  
  # "5. G = G2-G1"
  Mranks <- Mranks %>% dplyr::mutate(Gscore = incong - cong)
  
  # Put data in wide format
  overall <- reshape2::dcast(Mranks, sessionId ~ pair, value.var = "Gscore")
  
  # change col names
  colnames(overall)[1] <- "sid"
  names(overall)[names(overall) != "sid"] <- paste0("g_", names(overall)[names(overall) != "sid"])
  
  return(overall)
}

# Average 3 parts of G score
calc_g_parts <- function(data){
  max_block <- max(data$block)  # Get the maximum block number
  
  # Initialize an empty list to store scores for each block
  score_list <- list()
  
  # Calculate scores for each block
  for (i in 1:max_block) {
    scores <- calc_g_overall(data[data$block == i, ])
    score_list[[i]] <- scores
  }
  
  # Calculate the row means across blocks
  scores <- Reduce(`+`, score_list) / max_block
  
  # Initialize an empty data frame to store the blocks' scores
  all_block_scores <- data.frame(sid = unique(data$sid))
  
  # Iterate over each block and merge scores by 'sid'
  for (i in 1:max_block) {
    block_scores <- score_list[[i]]
    colnames(block_scores)[-1] <- paste0(names(block_scores)[-1], "_", i)  # Exclude the first column (sid)
    all_block_scores <- merge(all_block_scores, block_scores, by = "sid", all.x = TRUE)
  }
  
  # Merge all scores - final (average) and seperated by blocks
  scores <- merge(scores, all_block_scores, by = "sid", all.x = TRUE)
  
  
  return(scores)
}


# Algorithms' functions ----------------------------------------------------

# Algo 1 - e0.5_lowNone_2SD_noWinsorize_600msPenaltyErrors_noLog_g_overall
calc_algo_1 <- function(data){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- dplyr::select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- dplyr::select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- dplyr::select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_overall(final_data)
  
  # Add the blocks' scores
  max_block <- max(final_data$block)  # Get the maximum block number
  # Initialize an empty data frame to store the blocks' scores
  all_block_scores <- data.frame(sid = unique(final_data$sid))
  
  # Iterate over each block and merge scores by 'sid'
  for (i in 1:max_block) {
    block_scores <- calc_g_overall(data[final_data$block == i, ])
    colnames(block_scores)[-1] <- paste0(names(block_scores)[-1], "_", i)  # Exclude the first column (sid)
    all_block_scores <- merge(all_block_scores, block_scores, by = "sid", all.x = TRUE)
  }
  
  score <- merge(score, all_block_scores, by = "sid", all.x = TRUE)
  
  return(score)
}

# Algo 2 - e0.4_300_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_2 <- function(data){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- subset(correct, 300 <= rt)   # remove rt<300
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 3 - e0.5_300_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_3 <- function(data){
  # remove participant with error rate higher than 0.5
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  correct <- correct %>% mutate(rt = ifelse(rt < 300, 300, rt))   # rt<300 >> change rt to 300
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 4 - e0.4_lowNone_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_4 <- function(data){
  
  # remove participant with error rate higher than 0.4
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 5 - e0.4_300_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_5 <- function(data){
  
  # remove participant with error rate higher than 0.4
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  correct <- correct %>% mutate(rt = ifelse(rt < 300, 300, rt))   # rt<300 >> change rt to 300
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 6 - e0.5_lowNone_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_6 <- function(data){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 7 - e0.5_300_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_7 <- function(data){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- subset(correct, 300 <= rt)   # remove rt<300
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 8 - best treatments no winsorize - e0.4_lowNone_1000_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_8 <- function(data){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- subset(correct, rt <= 1000)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}

# Algo 9 - best tratments with winsorize - e0.4_lowNone_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_9 <- function(data){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)
  
  return(score)
}


# Users' functions --------------------------------------------------------

#' EPT's G score calculation
#'
#' Calculate EPT scores using one of the 9 algorithms recommended in Segal-Gordon et al. (2024).
#' Segal-Gordon et al. recommended Algorithms 1 and 7. Algorithm 1 is set as default.
#' 
#'
#' @param algoNum The chosen algorithm, with default value of 1.
#' @param blockNum The number of blocks, with default value of max(data$block).
#' @param data The processed data, with columns:
#' \itemize{
#'   \item \code{sid} (Participant's unique ID)
#'   \item \code{error} (Indicator - 1 is error, 0 is correct)
#'   \item \code{rt} (Reaction time)
#'   \item \code{prime_cat} (The primes category)
#'   \item \code{target_cat} (pos / neg)
#'   \item \code{block - optional} (Block number. If column is missing, it is calculated automatically by the "blockNum" argument.)
#' }
#' @return The results include EPT scores for each prime_cat, for each sid, overall and by block (for computing internal consistency).
#' 
#' @examples 
#' scores1 <- calc_ept_score(myProcessedData);
#' scores2 <- calc_ept_score(myProcessedData, algoNum = 7);
#' scores3 <- calc_ept_score(myProcessedDataWithoutBlockCol, blockNum = 5, algoNum = 3);
#' @export
#' 

calc_ept_score <- function(data, blockNum = max(data$block), algoNum = 1) {
  print(blockNum)
  data <- check_shape(data, blockNum)
  scores <- switch(as.character(algoNum),
                   "1" = calc_algo_1(data),
                   "2" = calc_algo_2(data),
                   "3" = calc_algo_3(data),
                   "4" = calc_algo_4(data),
                   "5" = calc_algo_5(data),
                   "6" = calc_algo_6(data),
                   "7" = calc_algo_7(data),
                   "8" = calc_algo_8(data),
                   "9" = calc_algo_9(data),
                   "Unsupported algorithm number")
  return(scores)
}




