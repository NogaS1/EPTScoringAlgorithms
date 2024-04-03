library(dplyr)
library(abind)

# Pre-processing ----------------------------------------------------------

# Data's shape:
# c("sid", "error", "rt", "prime_cat", "target_cat", "block", "trial", "part")
#  sid - subjects' unique id
# error - 1 is wrong, 0 is right
# rt - reaction time
# prime_cat - primes' categories
# target_cat - pos / neg
# block = part (need both)
# trial - trial number




# General calculations for algorithms -------------------------------------

# add 600ms penalty to errors
add_penalty <- function(correct, errors){
  # Calculate mean rt of block for each participant's block (after trimming)
  ep_mean_rt <- subset(correct, select = c("sid", "block", "rt"))
  mean_rt <- ep_mean_rt %>% dplyr::group_by(sid, block, .group = TRUE) %>% dplyr::mutate(mean_rt = mean(rt))
  mean_rt <- unique(mean_rt[,c(1,2,5)])
  errors_penalty <- merge(x=errors, y=mean_rt, by.x=c("sid", "block"), by.y = c("sid", "block"), all.x = TRUE)
  errors_penalty <- errors_penalty %>% mutate(rt = mean_rt+600)
  errors_penalty <- select(errors_penalty, -"mean_rt")
  
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
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
  
  
  # check sid have all 3 parts
  part_1 <- unique(data[data$part == 1, c("sid")])
  part_2 <- unique(data[data$part == 2, c("sid")])
  part_3 <- unique(data[data$part == 3, c("sid")])
  
  missing_parts_12 <- part_1[!part_1 %in% part_2]
  missing_parts_13 <- part_1[!part_1 %in% part_3]
  missing_parts_21 <- part_2[!part_2 %in% part_1]
  missing_parts_23 <- part_2[!part_2 %in% part_3]
  missing_parts_31 <- part_3[!part_3 %in% part_1]
  missing_parts_32 <- part_3[!part_3 %in% part_2]
  missing_parts <- unique(abind(missing_parts_12, missing_parts_13, missing_parts_21, missing_parts_23, missing_parts_31, missing_parts_32))
  data <- data[!data$sid %in% missing_parts,]
  
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
  Mranks <- mutate(Mranks, Gscore = incong - cong)

  # Put data in wide format
  overall <- reshape2::dcast(Mranks, sessionId ~ pair, value.var = "Gscore")

  
  # change col names
  colnames(overall)[1] <- "sid"
  names(overall)[names(overall) != "sid"] <- paste0("g_overall_", names(overall)[names(overall) != "sid"])
  
  return(overall)
}


# Average 3 parts of G score
calc_g_parts <- function(data){
  
  # Calculate scores for each part
  scores_1 <- calc_g_overall(data[data$part == 1,])
  scores_2 <- calc_g_overall(data[data$part == 2,])
  scores_3 <- calc_g_overall(data[data$part == 3,])
  
  # Average the scores
  arr <- abind(scores_1, scores_2, scores_3, along = 3)
  scores <- data.frame(rowMeans(arr, dims = 2))
  names(scores) <- gsub(x = names(scores), pattern = "overall", replacement = "parts")
  
  return(scores)
}




# Algorithms' functions ----------------------------------------------------

# Algo 1 - e0.5_lowNone_2SD_noWinsorize_600msPenaltyErrors_noLog_g_overall
calc_algo_1 <- function(data, outcome){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_overall(final_data)
  
  return(score)
}


# Algo 2 - e0.4_300_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_2 <- function(data, outcome){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- subset(correct, 300 <= rt)   # remove rt<300
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 3 - e0.5_300_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_3 <- function(data, outcome){
  # remove participant with error rate higher than 0.5
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  correct <- correct %>% mutate(rt = ifelse(rt < 300, 300, rt))   # rt<300 >> change rt to 300
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 4 - e0.4_lowNone_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_4 <- function(data, outcome){
  
  # remove participant with error rate higher than 0.4
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 5 - e0.4_300_1000_winsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_5 <- function(data, outcome){
  
  # remove participant with error rate higher than 0.4
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  correct <- correct %>% mutate(rt = ifelse(rt < 300, 300, rt))   # rt<300 >> change rt to 300
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 6 - e0.5_lowNone_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_6 <- function(data, outcome){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 7 - e0.5_300_2SD_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_7 <- function(data, outcome){
  
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.5)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- subset(correct, 300 <= rt)   # remove rt<300
  
  sd <- correct %>% dplyr::group_by(sid, prime_cat, target_cat, .group = TRUE) %>% dplyr::summarise(mean = mean(rt, na.rm = TRUE), sd = sd(rt, na.rm = TRUE))      # SD within participants, within each prime-target condition of that participant
  sd$high <- sd$mean + sd$sd*2
  correct <- select(correct, c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- select(merge(x = correct, y = sd, by = c("sid", "prime_cat", "target_cat"), all.x = TRUE), -c(".group"))
  correct <- subset(correct, rt <= correct$high)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 8 - best treatments no winsorize - e0.4_lowNone_1000_noWinsorize_600msPenaltyErrors_noLog_g_parts
calc_algo_8 <- function(data, outcome){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- subset(correct, rt <= 1000)
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Algo 9 - best tratments with winsorize - e0.4_lowNone_1000_winsorize_600msPenaltyErrors_noLog_g_parts

calc_algo_9 <- function(data, outcome){
  # add error rate to fit other functions
  data <- as.data.frame(data %>% dplyr::group_by(sid) %>% dplyr::mutate(error_rate = mean(error, na.rm = TRUE)))
  data <- subset(data, error_rate <= 0.4)
  
  # Keep errors aside
  errors <- subset(data, error == 1)
  
  # Outlier trial treatment - only for correct trials
  correct <- select(subset(data, error == 0), c("sid", "prime_cat", "target_cat", "error", "rt", "block", "trial", "part", "error_rate"))
  correct <- correct %>% mutate(rt = ifelse(1000 < rt, 1000, rt)) # 1000<rt >> change rt to 1000
  
  final_data <- remove_participants(add_penalty(correct, errors))
  
  score <- calc_g_parts(final_data)

  return(score)
}


# Users' functions --------------------------------------------------------

#' EPT's G score calculation
#'
#' Calculate scores according to the users' choice of algorithm, with algorithm 1 as default
#'
#' @param algoNum The chosen algorithm, with default value of 1.
#' @param data The processed data, with columns:
#' sid (subjects' unique id), 
#' error (1 is wrong, 0 is right), 
#' rt (reaction time), 
#' prime_cat (primes' categories), 
#' target_cat (pos / neg), 
#' block = part (EPT with 3 parts), 
#' trial - trial number.
#' The order should be: sid, error, rt, prime_cat, target_cat, block, trial, part.
#' @return Data frame of scores by subject (sid)
#' @examples 
#' scores1 <- calc_ept_score(myProcessedData);
#' scores2 <- calc_ept_score(7, myProcessedData);
#' @export
calc_ept_score <- function(data, algoNum = 1) {
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


