#Initialize

library('stringr')
library("phytools")
library("pROC")

###########################################################################################
###########################################################################################

#Function: convert_matrix
#The function converts the matrix into a long-format data frame ('lower.tri')
#where each row represents a unique pairwise comparison between two matrix 
#elements that are not on the diagonal. The pairwise comparisons are determined
#by the lower triangle of the matrix.
#Arguments: 'matrix' - a pairwise identity square matrix (output from geneious)
#           'groups' - a data frame that associates each row of the matrix 
#                      (column 1) with a group (column 2).
#           'group_column' - a column within 'groups' to define the 'comparison' column (see 
#                            below). The ccolumn values can be numeric or strings. Call can be
#                            done with column names (group_column = "column name") or numbers 
#                            (group_column = 2)
#The resulting data frame has two columns: 'comparison', which is a factor 
#variable that combines the names of the two groups being compared; and 
#'identity', which is the value of the matrix element being compared.
#The 'comparison' column is created by finding the minimum and maximum group 
#names for each pairwise comparison and concatenating them (needs 'str_c' function
#from package 'stringr'). The resulting 'comparison' column is then converted to a 
#factor variable.

convert_matrix <- function(matrix, groups, group_column, type_column) {
  x <- which(lower.tri(matrix, diag = FALSE), arr.ind = TRUE)
  matrix_long <- cbind(x, matrix[x])
  nrow_matrix_long <- nrow(matrix_long)
  
  groupid <- data.frame(group1 = rep(NA, nrow_matrix_long), 
                        group2 = rep(NA, nrow_matrix_long),
                        type1 = rep(NA, nrow_matrix_long),
                        type2 = rep(NA, nrow_matrix_long),
                        comparison = rep(NA, nrow_matrix_long),
                        comparison_type = rep(NA, nrow_matrix_long)) 
  
  for (i in 1:nrow_matrix_long) {
    value1 = matrix_long[i, 1]
    value2 = matrix_long[i, 2]
    groupid[i, 1] <- groups[value1, group_column]
    groupid[i, 2] <- groups[value2, group_column]
    groupid[i, 3] <- groups[value1, type_column]
    groupid[i, 4] <- groups[value2, type_column]
  }
  
  for (i in 1:nrow_matrix_long) {
    maxx = max(groupid[i, 1], groupid[i, 2])
    minn = min(groupid[i, 1], groupid[i, 2])
    groupid[i, 5] <- str_c(minn, maxx)
    maxx_type = max(groupid[i, 3], groupid[i, 4])
    minn_type = min(groupid[i, 3], groupid[i, 4])
    groupid[i, 6] <- str_c(minn_type, maxx_type)
    
  }
  
  groupid[, 5] <- as.factor(groupid[, 5])
  groupid[, 6] <- as.factor(groupid[, 6])
  
  matrix_long2 <- data.frame(comparison = factor(groupid$comparison), 
                             comparison_type = groupid$comparison_type,
                             identity = matrix_long[, 3])
  return(matrix_long2)
}



###########################################################################################
###########################################################################################

#Function: compare_allergen
#This function creates a vector 'x' meant to be attached to the long-form 
#data frame created with the 'convert_matrix' function created previously, 
#containing information about the kind of comparison being made (AA, ANA, NANA) .
#Arguments: matrix: a data frame created with the 'convert_matrix' function 
#           AA: a vector containing the comparison levels representing an 
#               allergenic vs allergenic comparison
#           ANA: a vector containing the comparison levels representing an 
#                allergenic vs non-allergenic comparison
#           NANA: a vector containing the comparison levels representing a 
#                 non-allergenic vs non-allergenic comparison
#The function loops over each row of the input matrix (output from 'convert_matrix'),
#checks if each element of the row is in the AA, ANA, or NANA vectors (%in%) 
#and assigns the corresponding string value ("AA", "ANA", or "NANA") to the 'x' vector. 
#If the value is not in any of the vectors, the function assigns "UNKNOWN" to 
#the Allergen_comparison vector.

compare_allergen <- function(matrix, AA, ANA, NANA) {
  
  x <- c()
  for (i in 1:nrow(matrix)) {
    value = matrix[i,1]
    if (value %in% AA) { x[i] <- "AA" }
    else if (value %in% ANA){ x[i] <- "ANA" }
    else if (value %in% NANA) { x[i] <- "NANA" }
    else { x[i] <- "UNKNOWN" }
  }
  x <- factor(x, levels = c("AA", "ANA", "NANA", "UNKNOWN"))
  return(x)
}

###########################################################################################
###########################################################################################

#Function: organism_split
#This function takes a long-format data frame 'x' and splits it into a list of data frames based on the provided 'org' list.
#Each data frame in the resulting list represents a subset of 'x' containing rows where the 'comparison' column matches the corresponding organism in 'org'.
#Arguments: x - a long-format data frame with 'comparison' and 'identity' columns.
#           org - a list of organism names used to split the 'x' data frame.
#Returns: A list of data frames, where each data frame corresponds to a specific organism from the 'org' list.

organism_split <- function(x, org) {
  dflistorg <- list()
  
  # Iterate over each organism in the 'org' list
  for (i in org) {
    # Create a unique data frame name based on the organism
    df_name <- paste("x_", i, sep = "")
    
    # Subset 'x' to include rows where the 'comparison' column matches the current organism
    x_sub <- x[grepl(i, x$comparison), ]
    
    # Assign the subset data frame to the unique name and add it to the list
    assign(df_name, x_sub)
    dflistorg[[i]] <- get(df_name)
  }
  
  # Return the list of data frames, each corresponding to an organism in 'org'
  return(dflistorg)
}

###########################################################################################
###########################################################################################

# Function: organism_row_plot
# This function generates row-wise plots for a list of data frames, where each data frame corresponds to the subset of data for a specific organism.
# It creates a grid of plots, with each row representing a different organism's data.
# Arguments: org_list - a list of data frames obtained using the organism_split function, each
#                       containing 'identity' and 'comparison' columns.
#            minvalue - the minimum value for the y-axis (identity) in the plots, defined       
#                       beforehand.
#            puntos - a logical value indicating whether to include points on the plots (default: #                     FALSE).
# Returns: None (plots are displayed).

organism_row_plot <- function(org_list, minvalue, org_names, puntos = c('TRUE','FALSE')) {
  
  rownum <- length(org_list)
  
  # Set up plot layout and margins
  par(mfrow = c(rownum, 1), mai = c(1, 0.2, 0.2, 0.2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(5,5,0.2,0.2))
  
  # Iterate over each data frame in org_list and create plots
  for (i in 1:rownum) {
    
    # Create a scatter plot
    plot(identity ~ comparison, org_list[[i]], 
         xlab = FALSE, xaxt = "n", ylab = org_names[i], 
         ylim = c(minvalue, 100), lwd = 0.75, las = 1)
    
    # Add points to the plot if puntos is TRUE
    if (puntos == "TRUE") {points(identity ~ comparison, org_list[[i]], pch = 20)}
    
    # Only add x-axis labels to bottom plot
    if(i == rownum) {axis(1, at=1:length(org_names), labels=org_names, las = 2)}
  }
}


