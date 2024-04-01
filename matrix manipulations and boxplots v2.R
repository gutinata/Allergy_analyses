
library('stringr')

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
###########################################################################################
###########################################################################################

#Tropomyosin

trop <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_identity.csv", header = T, row.names = 1)
trop <- as.matrix(trop)

trop_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_groups.csv", header = T)

tropomyosin <- convert_matrix(trop, trop_groups, group_column = "group_type", type_column = "allergenic")
head(tropomyosin)
levels(tropomyosin$comparison_type)
summary(tropomyosin$comparison)

trop_org_order <-c ("Mollusk","Crustacean","Insect","Mite","Worm","Fish","Bird","Mammal")
trop_org <-  organism_split(tropomyosin, trop_org_order)
summary(trop_org)

#reorder comparison column levels for boxplot
for (i in trop_org) {
  print(levels(droplevels(i$comparison)))
}

trop_org$Mollusk[,1] <- factor(trop_org$Mollusk[,1], levels = c("MolluskMollusk","CrustaceanMollusk","InsectMollusk","MiteMollusk","MolluskWorm","FishMollusk","BirdMollusk","MammalMollusk"))
trop_org$Crustacean[,1] <- factor(trop_org$Crustacean[,1], levels = c("CrustaceanMollusk","CrustaceanCrustacean","CrustaceanInsect","CrustaceanMite","CrustaceanWorm","CrustaceanFish","BirdCrustacean","CrustaceanMammal"))
trop_org$Insect[,1] <- factor(trop_org$Insect[,1], levels = c("InsectMollusk","CrustaceanInsect","InsectInsect","InsectMite","InsectWorm","FishInsect","BirdInsect","InsectMammal"))
trop_org$Mite[,1] <- factor(trop_org$Mite[,1], levels = c("MiteMollusk","CrustaceanMite","InsectMite","MiteMite","MiteWorm","FishMite","BirdMite","MammalMite"))
trop_org$Worm[,1] <- factor(trop_org$Worm[,1], levels = c("MolluskWorm","CrustaceanWorm","InsectWorm","MiteWorm","WormWorm","FishWorm","BirdWorm","MammalWorm"))
trop_org$Fish[,1] <- factor(trop_org$Fish[,1], levels = c("FishMollusk","CrustaceanFish","FishInsect","FishMite","FishWorm","FishFish","BirdFish","FishMammal"))
trop_org$Bird[,1] <- factor(trop_org$Bird[,1], levels = c("BirdMollusk","BirdCrustacean","BirdInsect","BirdMite","BirdWorm","BirdFish","BirdBird","BirdMammal"))
trop_org$Mammal[,1] <- factor(trop_org$Mammal[,1], levels = c("MammalMollusk","CrustaceanMammal","InsectMammal","MammalMite","MammalWorm","FishMammal","BirdMammal","MammalMammal"))

trop_min_iden <- min(tropomyosin$identity)

#plot pairwise boxplots
organism_row_plot(trop_org, trop_min_iden, trop_org_order, puntos = TRUE)

#plot comparison bloxplot
plot(identity ~ comparison_type, tropomyosin, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

TukeyHSD(aov(identity~allergen_comparison, tropomyosin2))

###########################################################################################
###########################################################################################

#Triosephosphate Isomerase

triose <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tim_identity.csv", header = T, row.names = 1)
triose <- as.matrix(triose)

triose_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tim_groups.csv", header = T)

tim <- convert_matrix(triose, triose_groups, group_column = "group_type", type_column = "allergenic", type_column = "allergenic")
head(tim)
levels(tim$comparison)
summary(tim$comparison)

tim_org_order <-c ("Fungi","Plant","Decapod","Mite","Fish","Mammal")
tim_org <-  organism_split(tim, tim_org_order)
summary(tim_org)

#reorder comparison column levels for boxplot

for (i in tim_org) {
  print(levels(droplevels(i$comparison)))
}

tim_org$Fungi[,1] <- factor(tim_org$Fungi[,1], levels = c("FungiFungi","FungiPlant","DecapodFungi","FungiMite","FishFungi","FungiMammal"))
tim_org$Plant[,1] <- factor(tim_org$Plant[,1], levels = c("FungiPlant","PlantPlant","DecapodPlant","MitePlant","FishPlant","MammalPlant"))
tim_org$Decapod[,1] <- factor(tim_org$Decapod[,1], levels = c("DecapodFungi","DecapodPlant","DecapodDecapod","DecapodMite","DecapodFish","DecapodMammal"))
tim_org$Mite[,1] <- factor(tim_org$Mite[,1], levels = c("FungiMite","MitePlant","DecapodMite","MiteMite", "FishMite","MammalMite"))
tim_org$Fish[,1] <- factor(tim_org$Fish[,1], levels = c("FishFungi","DecapodFish","DecapodMite","FishMite","FishFish" ,"FishMammal"))
tim_org$Mammal[,1] <- factor(tim_org$Mammal[,1], levels = c("FungiMammal","MammalPlant","DecapodMammal","MammalMite","FishMammal","MammalMammal"))

#check level order again

tim_min_iden <- min(tim$identity)

#plot pairwise boxplots
organism_row_plot(tim_org, tim_min_iden, tim_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, tim, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Beta Parvalbumin

bparv <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/bparvalbumin_identity.csv", header = T, row.names = 1)
bparv <- as.matrix(bparv)

bparv_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/bparvalbumin_groups.csv", header = T)

bparvalbumin <- convert_matrix(bparv, bparv_groups, group_column = "group_type", type_column = "allergenic")
head(bparvalbumin)
levels(bparvalbumin$comparison_type)
summary(bparvalbumin$comparison)

bparv_org_order <-c ("Otom.","Euacan.","Paracan.","Protacan.","Mammal")
bparv_org_order2 <-c ("Otomorpha","Euacanthomorphacea","Paracanthopterygii","Protacanthopterygii","Mammal")
bparv_org <-  organism_split(bparvalbumin, bparv_org_order)
summary(bparv_org)

#reorder comparison column levels for boxplot
for (i in bparv_org) {
  print(levels(droplevels(i$comparison)))
}

bparv_org$Otom.[,1] <- factor(bparv_org$Otom.[,1], levels = c("Otom.Otom.","Euacan.Otom.","Otom.Paracan.","Otom.Protacan.","MammalOtom."))
bparv_org$Euacan.[,1] <- factor(bparv_org$Euacan.[,1], levels = c("Euacan.Otom.","Euacan.Euacan.","Euacan.Paracan.","Euacan.Protacan.","Euacan.Mammal"))
bparv_org$Paracan.[,1] <- factor(bparv_org$Paracan.[,1], levels = c("Otom.Paracan.","Euacan.Paracan.","Paracan.Paracan.","Paracan.Protacan.","MammalParacan."))
bparv_org$Protacan.[,1] <- factor(bparv_org$Protacan.[,1], levels = c("Otom.Protacan.","Euacan.Protacan.","Paracan.Protacan.","Protacan.Protacan.","MammalProtacan."))
bparv_org$Mammal[,1] <- factor(bparv_org$Mammal[,1], levels = c("MammalOtom.","Euacan.Mammal","MammalParacan.","MammalProtacan.","MammalMammal"))

#check level order again

bparv_min_iden <- min(bparvalbumin$identity)

#plot pairwise boxplots
organism_row_plot(bparv_org, bparv_min_iden, bparv_org_order2, puntos = FALSE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, bparvalbumin, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#HSP

hsp <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/HSP_identity.csv", header = T, row.names = 1)
hsp <- as.matrix(hsp)

hsp_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/HSP_groups.csv", header = T)

heatshockp <- convert_matrix(hsp, hsp_groups, group_column = "group_type", type_column = "allergenic")
head(heatshockp)
levels(heatshockp$comparison_type)
summary(heatshockp$comparison)

hsp_org_order <-c ("AFungi","NAlFungi","Arthropod","Human")
hsp_org_order2 <-c ("Allergenic Fungi","Non Allergenic Fungi","Arthropod","Human")
hsp_org <-  organism_split(heatshockp, hsp_org_order)
summary(hsp_org)


#reorder comparison column levels for boxplot
for (i in hsp_org) {
  print(levels(droplevels(i$comparison)))
}

hsp_org$AFungi[,1] <- factor(hsp_org$AFungi[,1], levels = c("AFungiAFungi","AFungiNAlFungi","AFungiArthropod","AFungiHuman"))
hsp_org$NAlFungi[,1] <- factor(hsp_org$NAlFungi[,1], levels = c("AFungiNAlFungi","NAlFungiNAlFungi","ArthropodNAlFungi","HumanNAlFungi"))
hsp_org$Arthropod[,1] <- factor(hsp_org$Arthropod[,1], levels = c("AFungiArthropod","ArthropodNAlFungi","ArthropodArthropod","ArthropodHuman"))
hsp_org$Human[,1] <- factor(hsp_org$Human[,1], levels = c("AFungiHuman","HumanNAlFungi","ArthropodHuman","HumanHuman"))
  
#check level order again

hsp_min_iden <- min(heatshockp$identity)

#plot pairwise boxplots
organism_row_plot(hsp_org, hsp_min_iden, hsp_org_order2, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, heatshockp, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Casein alpha

c_a <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_alpha_identity.csv", header = T, row.names = 1)
c_a <- as.matrix(c_a)

c_a_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_alpha_groups.csv", header = T)

caseina <- convert_matrix(c_a, c_a_groups, group_column = "group_type", type_column = "allergenic")
head(caseina)
levels(caseina$comparison_type)
summary(caseina$comparison)

c_a_org_order <-c ("Bovinae","Caprinae","Camelidae","Equidae","Human")
c_a_org <-  organism_split(caseina, c_a_org_order)
summary(c_a_org)

#reorder comparison column levels for boxplot
for (i in c_a_org) {
  print(levels(droplevels(i$comparison)))
}

c_a_org$Bovinae[,1] <- factor(c_a_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeCamelidae","BovinaeEquidae","BovinaeHuman"))
c_a_org$Caprinae[,1] <- factor(c_a_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CamelidaeCaprinae","CaprinaeEquidae","CaprinaeHuman"))
c_a_org$Camelidae[,1] <- factor(c_a_org$Camelidae[,1], levels = c("BovinaeCamelidae","CamelidaeCaprinae","CamelidaeCamelidae","CamelidaeEquidae","CamelidaeHuman"))
c_a_org$Equidae[,1] <- factor(c_a_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","CamelidaeEquidae","EquidaeEquidae","EquidaeHuman"))
c_a_org$Human[,1] <- factor(c_a_org$Human[,1], levels = c("BovinaeHuman","CaprinaeHuman","CamelidaeHuman","EquidaeHuman","HumanHuman"))

#check level order again

c_a_min_iden <- min(caseina$identity)

#plot pairwise boxplots
organism_row_plot(c_a_org, c_a_min_iden, c_a_org_order, puntos = FALSE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, caseina, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Casein beta

c_b <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_beta_identity.csv", header = T, row.names = 1)
c_b <- as.matrix(c_b)

c_b_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_beta_groups.csv", header = T)

caseinb <- convert_matrix(c_b, c_b_groups, group_column = "group_type", type_column = "allergenic")
head(caseinb)
levels(caseinb$comparison_type)
summary(caseinb$comparison)

c_b_org_order <-c ("Bovinae","Caprinae","Camelidae","Equidae","Human")
c_b_org <-  organism_split(caseinb, c_b_org_order)
summary(c_b_org)

#reorder comparison column levels for boxplot
for (i in c_b_org) {
  print(levels(droplevels(i$comparison)))
}

c_b_org$Bovinae[,1] <- factor(c_b_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeCamelidae","BovinaeEquidae","BovinaeHuman"))
c_b_org$Caprinae[,1] <- factor(c_b_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CamelidaeCaprinae","CaprinaeEquidae","CaprinaeHuman"))
c_b_org$Camelidae[,1] <- factor(c_b_org$Camelidae[,1], levels = c("BovinaeCamelidae","CamelidaeCaprinae","CamelidaeCamelidae","CamelidaeEquidae","CamelidaeHuman"))
c_b_org$Equidae[,1] <- factor(c_b_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","CamelidaeEquidae","EquidaeEquidae","EquidaeHuman"))
c_b_org$Human[,1] <- factor(c_b_org$Human[,1], levels = c("BovinaeHuman","CaprinaeHuman","CamelidaeHuman","EquidaeHuman","HumanHuman"))

#check level order again

c_b_min_iden <- min(caseinb$identity)

#plot pairwise boxplots
organism_row_plot(c_b_org, c_b_min_iden, c_b_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, caseinb, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Arginine kinase

ak <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/argkin_identity.csv", header = T, row.names = 1)
ak <- as.matrix(ak)

ak_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/argkin_groups.csv", header = T)
levels(as.factor(ak_groups$group_type))

argkin <- convert_matrix(ak, ak_groups, group_column = "group_type", type_column = "allergenic")
head(argkin)
levels(argkin$comparison_type)
summary(argkin$comparison)

ak_org_order <-c ("Mollusk","Decapod","Insect","Mite","Bird","Mammal")
ak_org <-  organism_split(argkin, ak_org_order)
summary(ak_org)

#reorder comparison column levels for boxplot
for (i in ak_org) {
  print(levels(droplevels(i$comparison)))
}

ak_org$Mollusk[,1] <- factor(ak_org$Mollusk[,1], levels = c("MolluskMollusk","DecapodMollusk","InsectMollusk","MiteMollusk","BirdMollusk","MammalMollusk"))
ak_org$Decapod[,1] <- factor(ak_org$Decapod[,1], levels = c("DecapodMollusk","DecapodDecapod","DecapodInsect","DecapodMite","BirdDecapod","DecapodMammal"))
ak_org$Insect[,1] <- factor(ak_org$Insect[,1], levels = c("InsectMollusk","DecapodInsect","InsectInsect","InsectMite","BirdInsect","InsectMammal"))
ak_org$Mite[,1] <- factor(ak_org$Mite[,1], levels = c("MiteMollusk","DecapodMite","InsectMite","MiteMite","BirdMite","MammalMite"))
ak_org$Bird[,1] <- factor(ak_org$Bird[,1], levels = c("BirdMollusk","BirdDecapod","BirdInsect","BirdMite","BirdBird","BirdMammal"))
ak_org$Mammal[,1] <- factor(ak_org$Mammal[,1], levels = c("MammalMollusk","DecapodMammal","InsectMammal","MammalMite","BirdMammal","MammalMammal"))

#check level order again

ak_min_iden <- min(argkin$identity)

#plot pairwise boxplots
organism_row_plot(ak_org, ak_min_iden, ak_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, argkin, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Vicilin 7SC

vic_7SC <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SC_all_identity.csv", header = T, row.names = 1)

vic_7SC <- as.matrix(vic_7SC)

vic_7SC_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SC_groups.csv", header = T)


vicilin7SC <- convert_matrix(vic_7SC, vic_7SC_groups, group_column = "group_type", type_column = "allergenic")
head(vicilin7SC)

#subset to exclude non-allergenic groups
vicC <- vicilin7SC[vicilin7SC$comparison_type == "ALAL",]
levels(vicC$comparison_type) = droplevels(vicC$comparison_type)
summary(vicC)

vicC_org_order <- c("Fabales","Fagales","Rosales","Sapindales","Lamiales","Proteales","Pinopsida")
vicC_org <-  organism_split(vicC, vicC_org_order)
summary(vicC_org)

#reorder comparison column levels for boxplot
for (i in vicC_org) {
  print(levels(droplevels(i$comparison)))
}

vicC_org$Fabales[,1] <- factor(vicC_org$Fabales[,1], levels = c("FabalesFabales","FabalesFagales","FabalesRosales","FabalesSapindales","FabalesLamiales","FabalesProteales","FabalesPinopsida"))
vicC_org$Fagales[,1] <- factor(vicC_org$Fagales[,1], levels = c("FabalesFagales","FagalesFagales","FagalesRosales","FagalesSapindales","FagalesLamiales","FagalesProteales","FagalesPinopsida"))
vicC_org$Rosales[,1] <- factor(vicC_org$Rosales[,1], levels = c("FabalesRosales","FagalesRosales","RosalesRosales","RosalesSapindales","LamialesRosales","ProtealesRosales","PinopsidaRosales"))
vicC_org$Sapindales[,1] <- factor(vicC_org$Sapindales[,1], levels = c("FabalesSapindales","FagalesSapindales","RosalesSapindales","SapindalesSapindales","LamialesSapindales","ProtealesSapindales", "PinopsidaSapindales"))
vicC_org$Lamiales[,1] <- factor(vicC_org$Lamiales[,1], levels = c("FabalesLamiales","FagalesLamiales","LamialesRosales","LamialesSapindales","LamialesLamiales","LamialesProteales","LamialesPinopsida"))
vicC_org$Proteales[,1] <- factor(vicC_org$Proteales[,1], levels = c("FabalesProteales","FagalesProteales","ProtealesRosales","ProtealesSapindales","LamialesProteales","ProtealesProteales","PinopsidaProteales"))
vicC_org$Pinopsida[,1] <- factor(vicC_org$Pinopsida[,1], levels = c("FabalesPinopsida","FagalesPinopsida","PinopsidaRosales","PinopsidaSapindales","LamialesPinopsida","PinopsidaProteales","PinopsidaPinopsida"))

for (i in vicC_org) {
  print(levels(i$comparison))
}

vicC_min_iden <- min(vicC$identity)

#plot pairwise boxplots
organism_row_plot(vicC_org, vicC_min_iden, vicC_org_order, puntos = TRUE)

levels(vicC$comparison)

#replace values in "comparison_type" column to between/within

within <- c("FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales","LamialesLamiales","ProtealesProteales","PinopsidaPinopsida")
vicC$comparison_type <- NA
for (i in 1:nrow(vicC)){
  value = vicC$comparison[i]
  if (value %in% within) {vicC$comparison_type[i] <- "Within"}
  else {vicC$comparison_type[i] <- "Between"}
}
vicC$comparison_type <- factor(vicC$comparison_type, levels = c("Within","Between"))
levels(vicC$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, vicC, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Vicilin 7SN

vic_7SN <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SN_all_identity.csv", header = T, row.names = 1)

vic_7SN <- as.matrix(vic_7SN)

vic_7SN_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SN_groups.csv", header = T)


vicilin7SN <- convert_matrix(vic_7SN, vic_7SN_groups, group_column = "group_type", type_column = "allergenic")
head(vicilin7SN)
summary(vicilin7SN$comparison_type)

#subset to exclude non-allergenic groups
vicN <- vicilin7SN[vicilin7SN$comparison_type == "ALAL",]
levels(vicN$comparison_type) = droplevels(vicN$comparison_type)
summary(vicN)

vicN_org_order <- c("Fabales","Fagales","Rosales","Sapindales","Lamiales","Proteales","Pinopsida")
vicN_org <-  organism_split(vicN, vicN_org_order)
summary(vicN_org)

#reorder comparison column levels for boxplot
for (i in vicN_org) {
  print(levels(droplevels(i$comparison)))
}

vicN_org$Fabales[,1] <- factor(vicN_org$Fabales[,1], levels = c("FabalesFabales","FabalesFagales","FabalesRosales","FabalesSapindales","FabalesLamiales","FabalesProteales","FabalesPinopsida"))
vicN_org$Fagales[,1] <- factor(vicN_org$Fagales[,1], levels = c("FabalesFagales","FagalesFagales","FagalesRosales","FagalesSapindales","FagalesLamiales","FagalesProteales","FagalesPinopsida"))
vicN_org$Rosales[,1] <- factor(vicN_org$Rosales[,1], levels = c("FabalesRosales","FagalesRosales","RosalesRosales","RosalesSapindales","LamialesRosales","ProtealesRosales","PinopsidaRosales"))
vicN_org$Sapindales[,1] <- factor(vicN_org$Sapindales[,1], levels = c("FabalesSapindales","FagalesSapindales","RosalesSapindales","SapindalesSapindales","LamialesSapindales","ProtealesSapindales", "PinopsidaSapindales"))
vicN_org$Lamiales[,1] <- factor(vicN_org$Lamiales[,1], levels = c("FabalesLamiales","FagalesLamiales","LamialesRosales","LamialesSapindales","LamialesLamiales","LamialesProteales","LamialesPinopsida"))
vicN_org$Proteales[,1] <- factor(vicN_org$Proteales[,1], levels = c("FabalesProteales","FagalesProteales","ProtealesRosales","ProtealesSapindales","LamialesProteales","ProtealesProteales","PinopsidaProteales"))
vicN_org$Pinopsida[,1] <- factor(vicN_org$Pinopsida[,1], levels = c("FabalesPinopsida","FagalesPinopsida","PinopsidaRosales","PinopsidaSapindales","LamialesPinopsida","PinopsidaProteales","PinopsidaPinopsida"))

for (i in vicN_org) {
  print(levels(i$comparison))
}

vicN_min_iden <- min(vicN$identity)

#plot pairwise boxplots
organism_row_plot(vicN_org, vicN_min_iden, vicN_org_order, puntos = TRUE)

#replace values in "comparison_type" column to between/within

within <- c("FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales","LamialesLamiales","ProtealesProteales","PinopsidaPinopsida")
vicN$comparison_type <- NA
for (i in 1:nrow(vicN)){
  value = vicN$comparison[i]
  if (value %in% within) {vicN$comparison_type[i] <- "Within"}
  else {vicN$comparison_type[i] <- "Between"}
}
vicN$comparison_type <- factor(vicN$comparison_type, levels = c("Within","Between"))
levels(vicN$comparison_type)


#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, vicN, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#2S albumin

twoSalbumin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/2Salb_identity.csv", header = T, row.names = 1)

twoSalbumin <- as.matrix(twoSalbumin)

twoSalb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/2Salb_groups.csv", header = T)


twoSalbumin_all <- convert_matrix(twoSalbumin, twoSalb_groups, group_column = "group_type", type_column = "allergenic")
head(twoSalb)
summary(twoSalbumin_all$comparison_type)

#subset to exclude non-allergenic groups
twoSalb <- twoSalbumin_all[twoSalbumin_all$comparison_type == "ALAL",]
levels(twoSalb$comparison_type) = droplevels(twoSalb$comparison_type)
summary(twoSalb)

twoSalb_org_order <- c("Fabales","Fagales","Malpighiales","Cucurbitales","Sapindales","Brassicales","Lamiales","Ericales")
twoSalb_org <-  organism_split(twoSalb, twoSalb_org_order)
summary(twoSalb_org)

#reorder comparison column levels for boxplot
for (i in twoSalb_org) {
  print(levels(droplevels(i$comparison)))
}

twoSalb_org$Fabales[,1] <- factor(twoSalb_org$Fabales[,1], levels = c("FabalesFabales","FabalesFagales","FabalesMalpighiales","CucurbitalesFabales","FabalesSapindales","BrassicalesFabales","FabalesLamiales","EricalesFabales"))
twoSalb_org$Fagales[,1] <- factor(twoSalb_org$Fagales[,1], levels = c("FabalesFagales","FagalesFagales","FagalesMalpighiales","CucurbitalesFagales","FagalesSapindales","BrassicalesFagales","FagalesLamiales","EricalesFagales"))
twoSalb_org$Malpighiales[,1] <- factor(twoSalb_org$Malpighiales[,1], levels = c("FabalesMalpighiales","FagalesMalpighiales","MalpighialesMalpighiales","CucurbitalesMalpighiales","MalpighialesSapindales","BrassicalesMalpighiales","LamialesMalpighiales","EricalesMalpighiales"))
twoSalb_org$Cucurbitales[,1] <- factor(twoSalb_org$Cucurbitales[,1], levels = c("CucurbitalesFabales","CucurbitalesFagales","CucurbitalesMalpighiales","CucurbitalesCucurbitales","CucurbitalesSapindales","BrassicalesCucurbitales","CucurbitalesLamiales","CucurbitalesEricales"))
twoSalb_org$Sapindales[,1] <- factor(twoSalb_org$Sapindales[,1], levels = c("FabalesSapindales","FagalesSapindales","MalpighialesSapindales","CucurbitalesSapindales","SapindalesSapindales","BrassicalesSapindales","LamialesSapindales","EricalesSapindales"))
twoSalb_org$Brassicales[,1] <- factor(twoSalb_org$Brassicales[,1], levels = c("BrassicalesFabales","BrassicalesFagales","BrassicalesMalpighiales","BrassicalesCucurbitales","BrassicalesSapindales","BrassicalesBrassicales","BrassicalesLamiales","BrassicalesEricales"))
twoSalb_org$Lamiales[,1] <- factor(twoSalb_org$Lamiales[,1], levels = c("FabalesLamiales","FagalesLamiales","LamialesMalpighiales","CucurbitalesLamiales","LamialesSapindales","BrassicalesLamiales","LamialesLamiales","EricalesLamiales"))
twoSalb_org$Ericales[,1] <- factor(twoSalb_org$Ericales[,1], levels = c("EricalesFabales","EricalesFagales","EricalesMalpighiales","CucurbitalesEricales","EricalesSapindales","BrassicalesEricales","EricalesLamiales","LamialesLamiales"))

for (i in twoSalb_org) {
  print(levels(droplevels(i$comparison)))
}

twoSalb_min_iden <- min(twoSalb$identity)

#plot pairwise boxplots
organism_row_plot(twoSalb_org, twoSalb_min_iden, twoSalb_org_order, puntos = TRUE)

#replace values in "comparison_type" column to between/within

within <- c("FabalesFabales","FagalesFagales","MalpighialesMalpighiales","CucurbitalesCucurbitales","SapindalesSapindales","BrassicalesBrassicales","LamialesLamiales","EricalesEricales")
twoSalb$comparison_type <- NA
for (i in 1:nrow(twoSalb)){
  value = twoSalb$comparison[i]
  if (value %in% within) {twoSalb$comparison_type[i] <- "Within"}
  else {twoSalb$comparison_type[i] <- "Between"}
}
twoSalb$comparison_type <- factor(twoSalb$comparison_type, levels = c("Within","Between"))
levels(twoSalb$comparison_type)


#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, twoSalb, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#11S globulin

eleven <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/11S_identity.csv", header = T, row.names = 1)

eleven <- as.matrix(eleven)

eleven_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/11S_groups.csv", header = T)

eleven_all <- convert_matrix(eleven, eleven_groups, group_column = "group_type", type_column = "allergenic")
summary(eleven_all)

eleven_org_order <- c("Fabales","Fagales","Rosales","Cucurbitales","Sapindales","Brassicales","Lamiales","Ericales","Caryophyllales")
eleven_org <-  organism_split(eleven_all, eleven_org_order)
summary(eleven_org)

#reorder comparison column levels for boxplot
for (i in eleven_org) {
  print(levels(droplevels(i$comparison)))
}

eleven_org$Fabales[,1] <- factor(eleven_org$Fabales[,1], levels = c("FabalesFabales","FabalesFagales","FabalesRosales","CucurbitalesFabales","FabalesSapindales","BrassicalesFabales","FabalesLamiales","EricalesFabales","CaryophyllalesFabales"))
eleven_org$Fagales[,1] <- factor(eleven_org$Fagales[,1], levels = c("FabalesFagales","FagalesFagales","FagalesRosales","CucurbitalesFagales","FagalesSapindales","BrassicalesFagales","FagalesLamiales","EricalesFagales","CaryophyllalesFagales"))
eleven_org$Rosales[,1] <- factor(eleven_org$Rosales[,1], levels = c("FabalesRosales","FagalesRosales","RosalesRosales","CucurbitalesRosales","RosalesSapindales","BrassicalesRosales","LamialesRosales","EricalesRosales","CaryophyllalesRosales"))
eleven_org$Cucurbitales[,1] <- factor(eleven_org$Cucurbitales[,1], levels = c("CucurbitalesFabales","CucurbitalesFagales","CucurbitalesRosales","CucurbitalesCucurbitales","CucurbitalesSapindales","BrassicalesCucurbitales","CucurbitalesLamiales","CucurbitalesEricales","CaryophyllalesCucurbitales"))
eleven_org$Sapindales[,1] <- factor(eleven_org$Sapindales[,1], levels = c("FabalesSapindales","FagalesSapindales","RosalesSapindales","CucurbitalesSapindales","SapindalesSapindales","BrassicalesSapindales","LamialesSapindales","EricalesSapindales","CaryophyllalesSapindales"))
eleven_org$Brassicales[,1] <- factor(eleven_org$Brassicales[,1], levels = c("BrassicalesFabales","BrassicalesFagales","BrassicalesRosales","BrassicalesCucurbitales","BrassicalesSapindales","BrassicalesBrassicales","BrassicalesLamiales","BrassicalesEricales","BrassicalesCaryophyllales"))
eleven_org$Lamiales[,1] <- factor(eleven_org$Lamiales[,1], levels = c("FabalesLamiales","FagalesLamiales","LamialesRosales","CucurbitalesLamiales","LamialesSapindales","BrassicalesLamiales","LamialesLamiales","EricalesLamiales","CaryophyllalesLamiales"))
eleven_org$Ericales[,1] <- factor(eleven_org$Ericales[,1], levels = c("EricalesFabales","EricalesFagales","EricalesRosales","CucurbitalesEricales","EricalesSapindales","BrassicalesEricales","EricalesLamiales","EricalesEricales","CaryophyllalesEricales"))
eleven_org$Caryophyllales[,1] <- factor(eleven_org$Caryophyllales[,1], levels = c("CaryophyllalesFabales","CaryophyllalesFagales","CaryophyllalesRosales","CaryophyllalesCucurbitales","CaryophyllalesSapindales","BrassicalesCaryophyllales","CaryophyllalesLamiales","CaryophyllalesEricales","CaryophyllalesCaryophyllales"))

eleven_min_iden <- min(eleven_all$identity)

#plot pairwise boxplots
organism_row_plot(eleven_org, eleven_min_iden, eleven_org_order, puntos = TRUE)

#replace values in "comparison_type" column to between/within
within <- c("FabalesFabales","FagalesFagales","RosalesRosales","CucurbitalesCucurbitales","SapindalesSapindales","BrassicalesBrassicales","LamialesLamiales","EricalesEricales","CaryophyllalesCaryophyllales")
eleven_all$comparison_type <- NA
for (i in 1:nrow(eleven_all)){
  value = eleven_all$comparison[i]
  if (value %in% within) {eleven_all$comparison_type[i] <- "Within"}
  else {eleven_all$comparison_type[i] <- "Between"}
}
eleven_all$comparison_type <- factor(eleven_all$comparison_type, levels = c("Within","Between"))
levels(eleven_all$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, eleven_all, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#thaumatin

thaumatin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/thaumatin_identity.csv", header = T, row.names = 1)

thaumatin <- as.matrix(thaumatin)

thau_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/thaumatin_groups.csv", header = T)

thau <- convert_matrix(thaumatin, thau_groups, group_column = "group_type", type_column = "allergenic")
summary(thau)


thau_org_order <- c("Lamiales","Solanales","Ericales","Rosales","Poales","Zingiberales","Cupressales")
thau_org <-  organism_split(thau, thau_org_order)
summary(thau_org)

for (i in thau_org) {
  print(levels(droplevels(i$comparison)))
}

thau_org$Lamiales[,1] <- factor(thau_org$Lamiales[,1], levels = c("LamialesLamiales","LamialesSolanales","EricalesLamiales","LamialesRosales","LamialesPoales","LamialesZingiberales","CupressalesLamiales"))
thau_org$Solanales[,1] <- factor(thau_org$Solanales[,1], levels = c("LamialesSolanales","SolanalesSolanales","EricalesSolanales","RosalesSolanales","PoalesSolanales","SolanalesZingiberales","CupressalesSolanales"))
thau_org$Ericales[,1] <- factor(thau_org$Ericales[,1], levels = c("EricalesLamiales","EricalesSolanales","EricalesEricales","EricalesRosales","EricalesPoales","EricalesZingiberales","CupressalesEricales"))
thau_org$Rosales[,1] <- factor(thau_org$Rosales[,1], levels = c("LamialesRosales","RosalesSolanales","EricalesRosales","RosalesRosales","PoalesRosales","RosalesZingiberales","CupressalesRosales"))
thau_org$Poales[,1] <- factor(thau_org$Poales[,1], levels = c("LamialesPoales","PoalesSolanales","EricalesPoales","PoalesRosales","PoalesPoales","PoalesZingiberales","CupressalesPoales"))
thau_org$Zingiberales[,1] <- factor(thau_org$Zingiberales[,1], levels = c("LamialesZingiberales","SolanalesZingiberales","EricalesZingiberales","RosalesZingiberales","PoalesZingiberales","ZingiberalesZingiberales","CupressalesZingiberales"))
thau_org$Cupressales[,1] <- factor(thau_org$Cupressales[,1], levels = c("CupressalesLamiales","CupressalesSolanales","CupressalesEricales","CupressalesRosales","CupressalesPoales","CupressalesZingiberales","CupressalesCupressales"))

thau_min_iden <- min(thau$identity)

#plot pairwise boxplots
organism_row_plot(thau_org, thau_min_iden, thau_org_order, puntos = TRUE)

#replace values in "comparison_type" column to between/within
within <- c("LamialesLamiales","SolanalesSolanales","EricalesEricales","RosalesRosales","PoalesPoales","ZingiberalesZingiberales","CupressalesCupressales")

thau$comparison_type <- NA

for (i in 1:nrow(thau)){
  value = thau$comparison[i]
  if (value %in% within) {thau$comparison_type[i] <- "Within"}
  else {thau$comparison_type[i] <- "Between"}
}

thau$comparison_type <- factor(thau$comparison_type, levels = c("Within","Between"))
levels(thau$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, thau, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#gibberellin-like

gibberellin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/gibberellin_identity.csv", header = T, row.names = 1)
gibberellin <- as.matrix(gibberellin)

gib_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/gibberellin_groups.csv", header = T)

gib <- convert_matrix(gibberellin, gib_groups, group_column = "group_type", type_column = "allergenic")
summary(gib)

gib_org_order <- c("Rosales","Sapindales","Solanales","Cupressales")
gib_org <-  organism_split(gib, gib_org_order)
summary(gib_org)

for (i in gib_org) {
  print(levels(droplevels(i$comparison)))
}

gib_org$Rosales[,1] <- factor(gib_org$Rosales[,1], levels = c("RosalesRosales","RosalesSapindales",  "RosalesSolanales","CupressalesRosales"))
gib_org$Sapindales[,1] <- factor(gib_org$Sapindales[,1], levels = c("RosalesSapindales","SapindalesSapindales","SapindalesSolanales","CupressalesSapindales" ))
gib_org$Solanales[,1] <- factor(gib_org$Solanales[,1], levels = c("RosalesSolanales","SapindalesSolanales","SolanalesSolanales","CupressalesSolanales"))
gib_org$Cupressales[,1] <- factor(gib_org$Cupressales[,1], levels = c("CupressalesRosales","CupressalesSapindales","CupressalesSolanales","CupressalesCupressales"))

gib_min_iden <- min(gib$identity)

#plot pairwise boxplots
organism_row_plot(gib_org, gib_min_iden, gib_org_order, puntos = TRUE)

#replace values in "comparison_type" column to between/within
within <- c("RosalesRosales","SapindalesSapindales","SolanalesSolanales","CupressalesCupressales")

gib$comparison_type <- NA

for (i in 1:nrow(gib)){
  value = gib$comparison[i]
  if (value %in% within) {gib$comparison_type[i] <- "Within"}
  else {gib$comparison_type[i] <- "Between"}
}

gib$comparison_type <- factor(gib$comparison_type, levels = c("Within","Between"))
levels(gib$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, gib, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Lipocalin

lipoc <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lipocalin_identity.csv", header = T, row.names = 1)
lipoc <- as.matrix(lipoc)

lipoc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lipocalin_groups.csv", header = T)

lipocalin <- convert_matrix(lipoc, lipoc_groups, group_column = "group_type", type_column = "allergenic")
head(lipocalin)
levels(lipocalin$comparison_type)
summary(lipocalin$comparison)

lipoc_org_order <-c ("Carnivora","Perissodactyla","Rodentia","Lagomorpha","Primates")
lipoc_org <-  organism_split(lipocalin, lipoc_org_order)
summary(lipoc_org)

#reorder comparison column levels for boxplot
for (i in lipoc_org) {
  print(levels(droplevels(i$comparison)))
}

lipoc_org$Carnivora[,1] <- factor(lipoc_org$Carnivora[,1], levels = c("CarnivoraCarnivora","CarnivoraPerissodactyla","CarnivoraRodentia","CarnivoraLagomorpha","CarnivoraPrimates"))
lipoc_org$Perissodactyla[,1] <- factor(lipoc_org$Perissodactyla[,1], levels = c("CarnivoraPerissodactyla","PerissodactylaPerissodactyla","PerissodactylaRodentia","LagomorphaPerissodactyla","PerissodactylaPrimates"))
lipoc_org$Rodentia[,1] <- factor(lipoc_org$Rodentia[,1], levels = c("CarnivoraRodentia","PerissodactylaRodentia","RodentiaRodentia","LagomorphaRodentia","PrimatesRodentia"))
lipoc_org$Lagomorpha[,1] <- factor(lipoc_org$Lagomorpha[,1], levels = c("CarnivoraLagomorpha","LagomorphaPerissodactyla","LagomorphaRodentia","LagomorphaLagomorpha","LagomorphaPrimates"))
lipoc_org$Primates[,1] <- factor(lipoc_org$Primates[,1], levels = c("CarnivoraPrimates","PerissodactylaPrimates","PrimatesRodentia","LagomorphaPrimates","PrimatesPrimates"))

lipoc_min_iden <- min(lipocalin$identity)

#plot pairwise boxplots
organism_row_plot(lipoc_org, lipoc_min_iden, lipoc_org_order, puntos = TRUE)

#plot comparison bloxplot
plot(identity ~ comparison_type, lipocalin, 
     names = c("Allergenic", "Allergenic v. Non-allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#LTP

LTP <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/LTP_identity.csv", header = T, row.names = 1)

LTP <- as.matrix(LTP)

LTP_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/LTP_groups.csv", header = T)


LTPs <- convert_matrix(LTP, LTP_groups, group_column = "group_type", type_column = "allergenic")
head(LTPs)

LTP_org_order <- c("Poales","Proteales","Apiales","Asterales","Brassicales","Ericales","Fabales","Fagales","Malpighiales","Myrtales","Rosales","Sapindales","Solanales","Vitales")
LTP_org <-  organism_split(LTPs, LTP_org_order)
summary(LTP_org)

#reorder comparison column levels for boxplot
for (i in LTP_org) {
  print(levels(droplevels(i$comparison)))
}

LTP_org$Poales[,1] <- factor(LTP_org$Poales[,1], levels = c("PoalesPoales","PoalesProteales","ApialesPoales","AsteralesPoales","BrassicalesPoales","EricalesPoales","FabalesPoales","FagalesPoales","MalpighialesPoales","MyrtalesPoales","PoalesRosales","PoalesSapindales","PoalesSolanales","PoalesVitales"))
LTP_org$Proteales[,1] <- factor(LTP_org$Proteales[,1], levels = c("PoalesProteales","ProtealesProteales","ApialesProteales","AsteralesProteales","BrassicalesProteales","EricalesProteales","FabalesProteales","FagalesProteales","MalpighialesProteales","MyrtalesProteales","ProtealesRosales","ProtealesSapindales","ProtealesSolanales","ProtealesVitales"))
LTP_org$Apiales[,1] <- factor(LTP_org$Apiales[,1], levels = c("ApialesPoales","ApialesProteales","ApialesApiales","ApialesAsterales","ApialesBrassicales","ApialesEricales","ApialesFabales","ApialesFagales","ApialesMalpighiales","ApialesMyrtales","ApialesRosales","ApialesSapindales","ApialesSolanales","ApialesVitales"))
LTP_org$Asterales[,1] <- factor(LTP_org$Asterales[,1], levels = c("AsteralesPoales","AsteralesProteales","ApialesAsterales","AsteralesAsterales","AsteralesBrassicales","AsteralesEricales","AsteralesFabales","AsteralesFagales","AsteralesMalpighiales","AsteralesMyrtales","AsteralesRosales","AsteralesSapindales","AsteralesSolanales","AsteralesVitales"))
LTP_org$Brassicales[,1] <- factor(LTP_org$Brassicales[,1], levels = c("BrassicalesPoales","BrassicalesProteales","ApialesBrassicales","AsteralesBrassicales","BrassicalesBrassicales","BrassicalesEricales","BrassicalesFabales","BrassicalesFagales","BrassicalesMalpighiales","BrassicalesMyrtales","BrassicalesRosales","BrassicalesSapindales","BrassicalesSolanales","BrassicalesVitales"))
LTP_org$Ericales[,1] <- factor(LTP_org$Ericales[,1], levels = c("EricalesPoales","EricalesProteales","ApialesEricales","AsteralesEricales","BrassicalesEricales","EricalesEricales","EricalesFabales","EricalesFagales","EricalesMalpighiales","EricalesMyrtales","EricalesRosales","EricalesSapindales","EricalesSolanales","EricalesVitales"))
LTP_org$Fabales[,1] <- factor(LTP_org$Fabales[,1], levels = c("FabalesPoales","FabalesProteales","ApialesFabales","AsteralesFabales","BrassicalesFabales","EricalesFabales","FabalesFabales","FabalesFagales","FabalesMalpighiales","FabalesMyrtales","FabalesRosales","FabalesSapindales","FabalesSolanales","FabalesVitales"))
LTP_org$Fagales[,1] <- factor(LTP_org$Fagales[,1], levels = c("FagalesPoales","FagalesProteales","ApialesFagales","AsteralesFagales","BrassicalesFagales","EricalesFagales","FabalesFagales","FagalesFagales","FagalesMalpighiales","FagalesMyrtales","FagalesRosales","FagalesSapindales","FagalesSolanales","FagalesVitales"))
LTP_org$Malpighiales[,1] <- factor(LTP_org$Malpighiales[,1], levels = c("MalpighialesPoales","MalpighialesProteales","ApialesMalpighiales","AsteralesMalpighiales","BrassicalesMalpighiales","EricalesMalpighiales","FabalesMalpighiales","FagalesMalpighiales","MalpighialesMalpighiales","MalpighialesMyrtales","MalpighialesRosales","MalpighialesSapindales","MalpighialesSolanales","MalpighialesVitales"))
LTP_org$Myrtales[,1] <- factor(LTP_org$Myrtales[,1], levels = c("MyrtalesPoales","MyrtalesProteales","ApialesMyrtales","AsteralesMyrtales","BrassicalesMyrtales","EricalesMyrtales","FabalesMyrtales","FagalesMyrtales","MalpighialesMyrtales","MyrtalesMyrtales","MyrtalesRosales","MyrtalesSapindales","MyrtalesSolanales","MyrtalesVitales"))
LTP_org$Rosales[,1] <- factor(LTP_org$Rosales[,1], levels = c("PoalesRosales","ProtealesRosales","ApialesRosales","AsteralesRosales","BrassicalesRosales","EricalesRosales","FabalesRosales","FagalesRosales","MalpighialesRosales","MyrtalesRosales","RosalesRosales","RosalesSapindales","RosalesSolanales","RosalesVitales"))
LTP_org$Sapindales[,1] <- factor(LTP_org$Sapindales[,1], levels = c("PoalesSapindales","ProtealesSapindales","ApialesSapindales","AsteralesSapindales","BrassicalesSapindales","EricalesSapindales","FabalesSapindales","FagalesSapindales","MalpighialesSapindales","MyrtalesSapindales","RosalesSapindales","SapindalesSapindales","SapindalesSolanales","SapindalesVitales"))
LTP_org$Solanales[,1] <- factor(LTP_org$Solanales[,1], levels = c("PoalesSolanales","ProtealesSolanales","ApialesSolanales","AsteralesSolanales","BrassicalesSolanales","EricalesSolanales","FabalesSolanales","FagalesSolanales","MalpighialesSolanales","MyrtalesSolanales","RosalesSolanales","SapindalesSolanales","SolanalesSolanales","SolanalesVitales"))
LTP_org$Vitales[,1] <- factor(LTP_org$Vitales[,1], levels = c("PoalesVitales","ProtealesVitales","ApialesVitales","AsteralesVitales","BrassicalesVitales","EricalesVitales","FabalesVitales","FagalesVitales","MalpighialesVitales","MyrtalesVitales","RosalesVitales","SapindalesVitales","SolanalesVitales","VitalesVitales"))

for (i in LTP_org) {
  print(levels(i$comparison))
}

LTP_min_iden <- min(LTPs$identity)

#plot pairwise boxplots
organism_row_plot(LTP_org, LTP_min_iden, LTP_org_order, puntos = TRUE)

levels(LTPs$comparison)

#replace values in "comparison_type" column to between/within

within <- c("PoalesPoales","ProtealesProteales","ApialesApiales","AsteralesAsterales","BrassicalesBrassicales","EricalesEricales","FabalesFabales","FagalesFagales","MalpighialesMalpighiales","MyrtalesMyrtales","RosalesRosales","SapindalesSapindales","SolanalesSolanales","VitalesVitales")
LTPs$comparison_type <- NA
for (i in 1:nrow(LTPs)){
  value = LTPs$comparison[i]
  if (value %in% within) {LTPs$comparison_type[i] <- "Within"}
  else {LTPs$comparison_type[i] <- "Between"}
}
LTPs$comparison_type <- factor(LTPs$comparison_type, levels = c("Within","Between"))
levels(LTPs$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, LTPs, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Ovalbumin

ovalbumin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/ovalbumin_identity.csv", header = T, row.names = 1)
ovalbumin <- as.matrix(ovalbumin)

ovalb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/ovalbumin_groups.csv", header = T)

ovalb <- convert_matrix(ovalbumin, ovalb_groups, group_column = "group_type", type_column = "allergenic")
head(ovalb)
levels(ovalb$comparison_type)
summary(ovalb$comparison)

ovalb_org_order <-c ("Galliformes","Anseriformes","Charadriiformes","Cow","Human")
ovalb_org <-  organism_split(ovalb, ovalb_org_order)
summary(ovalb_org)

#reorder comparison column levels for boxplot
for (i in ovalb_org) {
  print(levels(droplevels(i$comparison)))
}

ovalb_org$Galliformes[,1] <- factor(ovalb_org$Galliformes[,1], levels = c("GalliformesGalliformes","AnseriformesGalliformes","CharadriiformesGalliformes","CowGalliformes","GalliformesHuman"))
ovalb_org$Anseriformes[,1] <- factor(ovalb_org$Anseriformes[,1], levels = c("AnseriformesGalliformes","AnseriformesAnseriformes","AnseriformesCharadriiformes","AnseriformesCow","AnseriformesHuman"))
ovalb_org$Charadriiformes[,1] <- factor(ovalb_org$Charadriiformes[,1], levels = c("CharadriiformesGalliformes","AnseriformesCharadriiformes","CharadriiformesCharadriiformes","CharadriiformesCow","CharadriiformesHuman"))
ovalb_org$Cow[,1] <- factor(ovalb_org$Cow[,1], levels = c("CowGalliformes","AnseriformesCow","CharadriiformesCow","CowCow","CowHuman"))
ovalb_org$Human[,1] <- factor(ovalb_org$Human[,1], levels = c("GalliformesHuman","AnseriformesHuman","CharadriiformesHuman","CowHuman","HumanHuman"))

ovalb_min_iden <- min(ovalb$identity)

#plot pairwise boxplots
organism_row_plot(ovalb_org, ovalb_min_iden, ovalb_org_order, puntos = TRUE)

#plot comparison bloxplot
plot(identity ~ comparison_type, ovalb, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Profilin

profilin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/profilin_identity.csv", header = T, row.names = 1)

profilin <- as.matrix(profilin)

prof_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/profilin_groups.csv", header = T)


prof <- convert_matrix(profilin, prof_groups, group_column = "group_type", type_column = "allergenic")
head(prof)

prof_org_order <- c("Apiales","Asterales","Solanales","Lamiales","Ericales","Cucurbitales","Fabales","Fagales","Malpighiales","Rosales","Brassicales","Sapindales","Caryophyllales","Asparagales","Arecales","Poales","Zingiberales","Proteales")
prof_org <-  organism_split(prof, prof_org_order)
summary(prof_org)

#reorder comparison column levels for boxplot
for (i in prof_org) {
  print(levels(droplevels(i$comparison)))
}

prof_org$Apiales[,1] <- factor(prof_org$Apiales[,1], levels = c("ApialesApiales","ApialesAsterales","ApialesSolanales","ApialesLamiales","ApialesEricales","ApialesCucurbitales","ApialesFabales","ApialesFagales","ApialesMalpighiales","ApialesRosales","ApialesBrassicales","ApialesSapindales","ApialesCaryophyllales","ApialesAsparagales","ApialesArecales","ApialesPoales","ApialesZingiberales","ApialesProteales"))
prof_org$Asterales[,1] <- factor(prof_org$Asterales[,1], levels = c("ApialesAsterales","AsteralesAsterales","AsteralesSolanales","AsteralesLamiales","AsteralesEricales","AsteralesCucurbitales","AsteralesFabales","AsteralesFagales","AsteralesMalpighiales","AsteralesRosales","AsteralesBrassicales","AsteralesSapindales","AsteralesCaryophyllales","AsparagalesAsterales","ArecalesAsterales","AsteralesPoales","AsteralesZingiberales","AsteralesProteales"))
prof_org$Solanales[,1] <- factor(prof_org$Solanales[,1], levels = c("ApialesSolanales","AsteralesSolanales","SolanalesSolanales","LamialesSolanales","EricalesSolanales","CucurbitalesSolanales","FabalesSolanales","FagalesSolanales","MalpighialesSolanales","RosalesSolanales","BrassicalesSolanales","SapindalesSolanales","CaryophyllalesSolanales","AsparagalesSolanales","ArecalesSolanales","PoalesSolanales","SolanalesZingiberales","ProtealesSolanales"))
prof_org$Lamiales[,1] <- factor(prof_org$Lamiales[,1], levels = c("ApialesLamiales","AsteralesLamiales","LamialesSolanales","LamialesLamiales","EricalesLamiales","CucurbitalesLamiales","FabalesLamiales","FagalesLamiales","LamialesMalpighiales","LamialesRosales","BrassicalesLamiales","LamialesSapindales","CaryophyllalesLamiales","AsparagalesLamiales","ArecalesLamiales","LamialesPoales","LamialesZingiberales","LamialesProteales"))
prof_org$Ericales[,1] <- factor(prof_org$Ericales[,1], levels = c("ApialesEricales","AsteralesEricales","EricalesSolanales","EricalesLamiales","EricalesEricales","CucurbitalesEricales","EricalesFabales","EricalesFagales","EricalesMalpighiales","EricalesRosales","BrassicalesEricales","EricalesSapindales","CaryophyllalesEricales","AsparagalesEricales","ArecalesEricales","EricalesPoales","EricalesZingiberales","EricalesProteales"))
prof_org$Cucurbitales[,1] <- factor(prof_org$Cucurbitales[,1], levels = c("ApialesCucurbitales","AsteralesCucurbitales","CucurbitalesSolanales","CucurbitalesLamiales","CucurbitalesEricales","CucurbitalesCucurbitales","CucurbitalesFabales","CucurbitalesFagales","CucurbitalesMalpighiales","CucurbitalesRosales","BrassicalesCucurbitales","CucurbitalesSapindales","CaryophyllalesCucurbitales","AsparagalesCucurbitales","ArecalesCucurbitales","CucurbitalesPoales","CucurbitalesZingiberales","CucurbitalesProteales"))
prof_org$Fabales[,1] <- factor(prof_org$Fabales[,1], levels = c("ApialesFabales","AsteralesFabales","FabalesSolanales","FabalesLamiales","EricalesFabales","CucurbitalesFabales","FabalesFabales","FabalesFagales","FabalesMalpighiales","FabalesRosales","BrassicalesFabales","FabalesSapindales","CaryophyllalesFabales","AsparagalesFabales","ArecalesFabales","FabalesPoales","FabalesZingiberales","FabalesProteales"))
prof_org$Fagales[,1] <- factor(prof_org$Fagales[,1], levels = c("ApialesFagales","AsteralesFagales","FagalesSolanales","FagalesLamiales","EricalesFagales","CucurbitalesFagales","FabalesFagales","FagalesFagales","FagalesMalpighiales","FagalesRosales","BrassicalesFagales","FagalesSapindales","CaryophyllalesFagales","AsparagalesFagales","ArecalesFagales","FagalesPoales","FagalesZingiberales","FagalesProteales"))
prof_org$Malpighiales[,1] <- factor(prof_org$Malpighiales[,1], levels = c("ApialesMalpighiales","AsteralesMalpighiales","MalpighialesSolanales","LamialesMalpighiales","EricalesMalpighiales","CucurbitalesMalpighiales","FabalesMalpighiales","FagalesMalpighiales","MalpighialesMalpighiales","MalpighialesRosales","BrassicalesMalpighiales","MalpighialesSapindales","CaryophyllalesMalpighiales","AsparagalesMalpighiales","ArecalesMalpighiales","MalpighialesPoales","MalpighialesZingiberales","MalpighialesProteales"))
prof_org$Rosales[,1] <- factor(prof_org$Rosales[,1], levels = c("ApialesRosales","AsteralesRosales","RosalesSolanales","LamialesRosales","EricalesRosales","CucurbitalesRosales","FabalesRosales","FagalesRosales","MalpighialesRosales","RosalesRosales","BrassicalesRosales","RosalesSapindales","CaryophyllalesRosales","AsparagalesRosales","ArecalesRosales","PoalesRosales","RosalesZingiberales","ProtealesRosales"))
prof_org$Brassicales[,1] <- factor(prof_org$Brassicales[,1], levels = c("ApialesBrassicales","AsteralesBrassicales","BrassicalesSolanales","BrassicalesLamiales","BrassicalesEricales","BrassicalesCucurbitales","BrassicalesFabales","BrassicalesFagales","BrassicalesMalpighiales","BrassicalesRosales","BrassicalesBrassicales","BrassicalesSapindales","BrassicalesCaryophyllales","AsparagalesBrassicales","ArecalesBrassicales","BrassicalesPoales","BrassicalesZingiberales","BrassicalesProteales"))
prof_org$Sapindales[,1] <- factor(prof_org$Sapindales[,1], levels = c("ApialesSapindales","AsteralesSapindales","SapindalesSolanales","LamialesSapindales","EricalesSapindales","CucurbitalesSapindales","FabalesSapindales","FagalesSapindales","MalpighialesSapindales","RosalesSapindales","BrassicalesSapindales","SapindalesSapindales","CaryophyllalesSapindales","AsparagalesSapindales","ArecalesSapindales","PoalesSapindales","SapindalesZingiberales","ProtealesSapindales"))
prof_org$Caryophyllales[,1] <- factor(prof_org$Caryophyllales[,1], levels = c("ApialesCaryophyllales","AsteralesCaryophyllales","CaryophyllalesSolanales","CaryophyllalesLamiales","CaryophyllalesEricales","CaryophyllalesCucurbitales","CaryophyllalesFabales","CaryophyllalesFagales","CaryophyllalesMalpighiales","CaryophyllalesRosales","BrassicalesCaryophyllales","CaryophyllalesSapindales","CaryophyllalesCaryophyllales","AsparagalesCaryophyllales","ArecalesCaryophyllales","CaryophyllalesPoales","CaryophyllalesZingiberales","CaryophyllalesProteales"))
prof_org$Asparagales[,1] <- factor(prof_org$Asparagales[,1], levels = c("ApialesAsparagales","AsparagalesAsterales","AsparagalesSolanales","AsparagalesLamiales","AsparagalesEricales","AsparagalesCucurbitales","AsparagalesFabales","AsparagalesFagales","AsparagalesMalpighiales","AsparagalesRosales","AsparagalesBrassicales","AsparagalesSapindales","AsparagalesCaryophyllales","AsparagalesAsparagales","ArecalesAsparagales","AsparagalesPoales","AsparagalesZingiberales","AsparagalesProteales"))
prof_org$Arecales[,1] <- factor(prof_org$Arecales[,1], levels = c("ApialesArecales","ArecalesAsterales","ArecalesSolanales","ArecalesLamiales","ArecalesEricales","ArecalesCucurbitales","ArecalesFabales","ArecalesFagales","ArecalesMalpighiales","ArecalesRosales","ArecalesBrassicales","ArecalesSapindales","ArecalesCaryophyllales","ArecalesAsparagales","ArecalesArecales","ArecalesPoales","ArecalesZingiberales","ArecalesProteales"))
prof_org$Poales[,1] <- factor(prof_org$Poales[,1], levels = c("ApialesPoales","AsteralesPoales","PoalesSolanales","LamialesPoales","EricalesPoales","CucurbitalesPoales","FabalesPoales","FagalesPoales","MalpighialesPoales","PoalesRosales","BrassicalesPoales","PoalesSapindales","CaryophyllalesPoales","AsparagalesPoales","ArecalesPoales","PoalesPoales","PoalesZingiberales","PoalesProteales"))
prof_org$Zingiberales[,1] <- factor(prof_org$Zingiberales[,1], levels = c("ApialesZingiberales","AsteralesZingiberales","SolanalesZingiberales","LamialesZingiberales","EricalesZingiberales","CucurbitalesZingiberales","FabalesZingiberales","FagalesZingiberales","MalpighialesZingiberales","RosalesZingiberales","BrassicalesZingiberales","SapindalesZingiberales","CaryophyllalesZingiberales","AsparagalesZingiberales","ArecalesZingiberales","PoalesZingiberales","ZingiberalesZingiberales","ProtealesZingiberales"))
prof_org$Proteales[,1] <- factor(prof_org$Proteales[,1], levels = c("ApialesProteales","AsteralesProteales","ProtealesSolanales","LamialesProteales","EricalesProteales","CucurbitalesProteales","FabalesProteales","FagalesProteales","MalpighialesProteales","ProtealesRosales","BrassicalesProteales","ProtealesSapindales","CaryophyllalesProteales","AsparagalesProteales","ArecalesProteales","PoalesProteales","ProtealesZingiberales","ProtealesProteales"))

for (i in prof_org) {
  print(levels(i$comparison))
}

prof_min_iden <- min(prof$identity)

#plot pairwise boxplots
organism_row_plot(prof_org, prof_min_iden, prof_org_order, puntos = TRUE)

levels(prof$comparison)

#replace values in "comparison_type" column to between/within

within <- c("ApialesApiales","AsteralesAsterales","SolanalesSolanales","LamialesLamiales","EricalesEricales","CucurbitalesCucurbitales","FabalesFabales","FagalesFagales","MalpighialesMalpighiales","RosalesRosales","BrassicalesBrassicales","SapindalesSapindales","CaryophyllalesCaryophyllales","AsparagalesAsparagales","ArecalesArecales","PoalesPoales","ZingiberalesZingiberales","ProtealesProteales")
prof$comparison_type <- NA
for (i in 1:nrow(prof)){
  value = prof$comparison[i]
  if (value %in% within) {prof$comparison_type[i] <- "Within"}
  else {prof$comparison_type[i] <- "Between"}
}
prof$comparison_type <- factor(prof$comparison_type, levels = c("Within","Between"))
levels(prof$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, prof, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#PR10 Bet v 1 like

pr10_m <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/PR10_identity.csv", header = T, row.names = 1)

pr10_m <- as.matrix(pr10_m)

pr10_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/PR10_groups.csv", header = T)

pr10 <- convert_matrix(pr10_m, pr10_groups, group_column = "group_type", type_column = "allergenic")
head(pr10)

pr10_org_order <- c("Apiales","Ericales","Solanales","Fabales","Fagales","Rosales","Sapindales")
pr10_org <-  organism_split(pr10, pr10_org_order)
summary(pr10_org)

#reorder comparison column levels for boxplot
for (i in pr10_org) {
  print(levels(droplevels(i$comparison)))
}

pr10_org$Apiales[,1] <- factor(pr10_org$Apiales[,1], levels = c("ApialesApiales","ApialesEricales","ApialesSolanales","ApialesFabales","ApialesFagales","ApialesRosales","ApialesSapindales"))
pr10_org$Ericales[,1] <- factor(pr10_org$Ericales[,1], levels = c("ApialesEricales","EricalesEricales","EricalesSolanales","EricalesFabales","EricalesFagales","EricalesRosales","EricalesSapindales"))
pr10_org$Solanales[,1] <- factor(pr10_org$Solanales[,1], levels = c("ApialesSolanales","EricalesSolanales","SolanalesSolanales","FabalesSolanales","FagalesSolanales","RosalesSolanales","SapindalesSolanales"))
pr10_org$Fabales[,1] <- factor(pr10_org$Fabales[,1], levels = c("ApialesFabales","EricalesFabales","FabalesSolanales","FabalesFabales","FabalesFagales","FabalesRosales","FabalesSapindales"))
pr10_org$Fagales[,1] <- factor(pr10_org$Fagales[,1], levels = c("ApialesFagales","EricalesFagales","FagalesSolanales","FabalesFagales","FagalesFagales","FagalesRosales","FagalesSapindales"))
pr10_org$Rosales[,1] <- factor(pr10_org$Rosales[,1], levels = c("ApialesRosales","EricalesRosales","RosalesSolanales","FabalesRosales","FagalesRosales","RosalesRosales","RosalesSapindales"))
pr10_org$Sapindales[,1] <- factor(pr10_org$Sapindales[,1], levels = c("ApialesSapindales","EricalesSapindales","SapindalesSolanales","FabalesSapindales","FagalesSapindales","RosalesSapindales"))

for (i in pr10_org) {
  print(levels(i$comparison))
}

pr10_min_iden <- min(pr10$identity)

#plot pairwise boxplots
organism_row_plot(pr10_org, pr10_min_iden, pr10_org_order, puntos = TRUE)

levels(pr10$comparison)

#replace values in "comparison_type" column to between/within
within <- c("ApialesApiales","EricalesEricales","SolanalesSolanales","FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales")

pr10$comparison_type <- NA
for (i in 1:nrow(pr10)){
  value = pr10$comparison[i]
  if (value %in% within) {pr10$comparison_type[i] <- "Within"}
  else {pr10$comparison_type[i] <- "Between"}
}
pr10$comparison_type <- factor(pr10$comparison_type, levels = c("Within","Between"))
levels(pr10$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, pr10, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Cyclophilin

cyclophilin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/cyclophilin_identity.csv", header = T, row.names = 1)
cyclophilin <- as.matrix(cyclophilin)

cyc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/cyclophilin_groups.csv", header = T)
levels(as.factor(cyc_groups$group_type))

cyc <- convert_matrix(cyclophilin, cyc_groups, group_column = "group_type", type_column = "allergenic")
head(cyc)
levels(cyc$comparison_type)
summary(cyc$comparison)

cyc_org_order <-c ("Plant","Fungi","Invertebrate","Human")
cyc_org <-  organism_split(cyc, cyc_org_order)
summary(cyc_org)

#reorder comparison column levels for boxplot
for (i in cyc_org) {
  print(levels(droplevels(i$comparison)))
}

cyc_org$Plant[,1] <- factor(cyc_org$Plant[,1], levels = c("PlantPlant","FungiPlant","InvertebratePlant","HumanPlant"))
cyc_org$Fungi[,1] <- factor(cyc_org$Fungi[,1], levels = c("FungiPlant","FungiFungi","FungiInvertebrate","FungiHuman"))
cyc_org$Invertebrate[,1] <- factor(cyc_org$Invertebrate[,1], levels = c("InvertebratePlant","FungiInvertebrate","InvertebrateInvertebrate","HumanInvertebrate"))
cyc_org$Human[,1] <- factor(cyc_org$Human[,1], levels = c("HumanPlant","FungiHuman","HumanInvertebrate","HumanHuman"))

#check level order again

cyc_min_iden <- min(cyc$identity)

#plot pairwise boxplots
organism_row_plot(cyc_org, cyc_min_iden, cyc_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, cyc, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)



###########################################################################################
###########################################################################################

#Polcalcin

polcalcin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/polcalcin_identity.csv", header = T, row.names = 1)

polcalcin <- as.matrix(polcalcin)

polc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/polcalcin_groups.csv", header = T)

polc <- convert_matrix(polcalcin, polc_groups, group_column = "group_type", type_column = "allergenic")
head(polc)

polc_org_order <- c("Asterales","Lamiales","Fagales","Rosales","Brassicales","Caryophyllales","Poales")
polc_org <-  organism_split(polc, polc_org_order)
summary(polc_org)

#reorder comparison column levels for boxplot
for (i in polc_org) {
  print(levels(droplevels(i$comparison)))
}

polc_org$Asterales[,1] <- factor(polc_org$Asterales[,1], levels = c("AsteralesAsterales","AsteralesLamiales","AsteralesFagales","AsteralesRosales","AsteralesBrassicales","AsteralesCaryophyllales","AsteralesPoales"))
polc_org$Lamiales[,1] <- factor(polc_org$Lamiales[,1], levels = c("AsteralesLamiales","LamialesLamiales","FagalesLamiales","LamialesRosales","BrassicalesLamiales","CaryophyllalesLamiales","LamialesPoales"))
polc_org$Fagales[,1] <- factor(polc_org$Fagales[,1], levels = c("AsteralesFagales","FagalesLamiales","FagalesFagales","FagalesRosales","BrassicalesFagales","CaryophyllalesFagales","FagalesPoales"))
polc_org$Rosales[,1] <- factor(polc_org$Rosales[,1], levels = c("AsteralesRosales","LamialesRosales","FagalesRosales","RosalesRosales","BrassicalesRosales","CaryophyllalesRosales","PoalesRosales"))
polc_org$Brassicales[,1] <- factor(polc_org$Brassicales[,1], levels = c("AsteralesBrassicales","BrassicalesLamiales","BrassicalesFagales","BrassicalesRosales","BrassicalesBrassicales","BrassicalesCaryophyllales","BrassicalesPoales"))
polc_org$Caryophyllales[,1] <- factor(polc_org$Caryophyllales[,1], levels = c("AsteralesCaryophyllales","CaryophyllalesLamiales","CaryophyllalesFagales","CaryophyllalesRosales","BrassicalesCaryophyllales","CaryophyllalesCaryophyllales","CaryophyllalesPoales"))
polc_org$Poales[,1] <- factor(polc_org$Poales[,1], levels = c("AsteralesPoales","LamialesPoales","FagalesPoales","PoalesRosales","BrassicalesPoales","CaryophyllalesPoales","PoalesPoales"))

for (i in polc_org) {
  print(levels(i$comparison))
}

polc_min_iden <- min(polc$identity)

#plot pairwise boxplots
organism_row_plot(polc_org, polc_min_iden, polc_org_order, puntos = TRUE)

levels(polc$comparison)

#replace values in "comparison_type" column to between/within
within <- c("AsteralesAsterales","LamialesLamiales","FagalesFagales","RosalesRosales","BrassicalesBrassicales","CaryophyllalesCaryophyllales","PoalesPoales")

polc$comparison_type <- NA
for (i in 1:nrow(polc)){
  value = polc$comparison[i]
  if (value %in% within) {polc$comparison_type[i] <- "Within"}
  else {polc$comparison_type[i] <- "Between"}
}
polc$comparison_type <- factor(polc$comparison_type, levels = c("Within","Between"))
levels(polc$comparison_type)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, polc, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL)



###########################################################################################
###########################################################################################

#Enolase

enolase <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/enolase_identity.csv", header = T, row.names = 1)
enolase <- as.matrix(enolase)

eno_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/enolase_groups.csv", header = T)
levels(as.factor(eno_groups$group_type))

eno <- convert_matrix(enolase, eno_groups, group_column = "group_type", type_column = "allergenic")
tail(eno)
levels(eno$comparison_type)
summary(eno$comparison)

eno_org_order <-c ("Plant","Fungi","Invertebrate","Fish","Chicken","Human")
eno_org <-  organism_split(eno, eno_org_order)
summary(eno_org)

#reorder comparison column levels for boxplot
for (i in eno_org) {
  print(levels(droplevels(i$comparison)))
}

eno_org$Plant[,1] <- factor(eno_org$Plant[,1], levels = c("PlantPlant","FungiPlant","InvertebratePlant","FishPlant","ChickenPlant","HumanPlant"))
eno_org$Fungi[,1] <- factor(eno_org$Fungi[,1], levels = c("FungiPlant","FungiFungi","FungiInvertebrate","FishFungi","ChickenFungi","FungiHuman"))
eno_org$Invertebrate[,1] <- factor(eno_org$Invertebrate[,1], levels = c("InvertebratePlant","FungiInvertebrate","InvertebrateInvertebrate","FishInvertebrate","ChickenInvertebrate","HumanInvertebrate"))
eno_org$Fish[,1] <- factor(eno_org$Fish[,1], levels = c("FishPlant","FishFungi","FishInvertebrate","FishFish","ChickenFish","FishHuman"))
eno_org$Chicken[,1] <- factor(eno_org$Chicken[,1], levels = c("ChickenPlant","ChickenFungi","ChickenInvertebrate","ChickenFish","ChickenChicken","ChickenHuman"))
eno_org$Human[,1] <- factor(eno_org$Human[,1], levels = c("HumanPlant","FungiHuman","HumanInvertebrate","FishHuman","ChickenHuman","HumanHuman"))

#check level order again

eno_min_iden <- min(eno$identity)

#plot pairwise boxplots
organism_row_plot(eno_org, eno_min_iden, eno_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, eno, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Alpha Lactalbumin

alphalac <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/alpha_lactalbumin_identity.csv", header = T, row.names = 1)
alphalac <- as.matrix(alphalac)

alac_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/alpha_lactalbumin_groups.csv", header = T)

alac <- convert_matrix(alphalac, alac_groups, group_column = "group_type", type_column = "allergenic")
head(alac)
levels(alac$comparison_type)
summary(alac$comparison)

alac_org_order <-c ("Bovinae","Caprinae","Camelidae","Equidae","Human")
alac_org <-  organism_split(alac, alac_org_order)
summary(alac_org)

#reorder comparison column levels for boxplot
for (i in alac_org) {
  print(levels(droplevels(i$comparison)))
}

alac_org$Bovinae[,1] <- factor(alac_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeCamelidae","BovinaeEquidae","BovinaeHuman"))
alac_org$Caprinae[,1] <- factor(alac_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CamelidaeCaprinae","CaprinaeEquidae","CaprinaeHuman"))
alac_org$Camelidae[,1] <- factor(alac_org$Camelidae[,1], levels = c("BovinaeCamelidae","CamelidaeCaprinae","CamelidaeCamelidae","CamelidaeEquidae","CamelidaeHuman"))
alac_org$Equidae[,1] <- factor(alac_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","CamelidaeEquidae","EquidaeEquidae","EquidaeHuman"))
alac_org$Human[,1] <- factor(alac_org$Human[,1], levels = c("BovinaeHuman","CaprinaeHuman","CamelidaeHuman","EquidaeHuman","HumanHuman"))

#check level order again

alac_min_iden <- min(alac$identity)

#plot pairwise boxplots
organism_row_plot(alac_org, alac_min_iden, alac_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, alac, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Beta Lactoglobulin

betalac <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/beta_lactoglobulin_identity.csv", header = T, row.names = 1)
betalac <- as.matrix(betalac)

blac_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/beta_lactoglobulin_groups.csv", header = T)

blac <- convert_matrix(betalac, blac_groups, group_column = "group_type", type_column = "allergenic")
head(blac)
levels(blac$comparison_type)
summary(blac$comparison)

blac_org_order <-c ("Bovinae","Caprinae","Equidae")
blac_org <-  organism_split(blac, blac_org_order)
summary(blac_org)

#reorder comparison column levels for boxplot
for (i in blac_org) {
  print(levels(droplevels(i$comparison)))
}

blac_org$Bovinae[,1] <- factor(blac_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeEquidae"))
blac_org$Caprinae[,1] <- factor(blac_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CaprinaeEquidae"))
blac_org$Equidae[,1] <- factor(blac_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","EquidaeEquidae"))

blac_min_iden <- min(blac$identity)

#plot pairwise boxplots
organism_row_plot(blac_org, blac_min_iden, blac_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, blac, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)


###########################################################################################
###########################################################################################

#Serum albumin

serumalb <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/serum_albumin_identity.csv", header = T, row.names = 1)
serumalb <- as.matrix(serumalb)

sealb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/serum_albumin_groups.csv", header = T)

sealb <- convert_matrix(serumalb, sealb_groups, group_column = "group_type", type_column = "allergenic")
head(sealb)
levels(sealb$comparison_type)
summary(sealb$comparison)

sealb_org_order <-c ("Bovinae","Caprinae","Camelidae","Equidae","Human")
sealb_org <-  organism_split(sealb, sealb_org_order)
summary(sealb_org)

#reorder comparison column levels for boxplot
for (i in sealb_org) {
  print(levels(droplevels(i$comparison)))
}

sealb_org$Bovinae[,1] <- factor(sealb_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeCamelidae","BovinaeEquidae","BovinaeHuman"))
sealb_org$Caprinae[,1] <- factor(sealb_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CamelidaeCaprinae","CaprinaeEquidae","CaprinaeHuman"))
sealb_org$Camelidae[,1] <- factor(sealb_org$Camelidae[,1], levels = c("BovinaeCamelidae","CamelidaeCaprinae","CamelidaeCamelidae","CamelidaeEquidae","CamelidaeHuman"))
sealb_org$Equidae[,1] <- factor(sealb_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","CamelidaeEquidae","EquidaeEquidae","EquidaeHuman"))
sealb_org$Human[,1] <- factor(sealb_org$Human[,1], levels = c("BovinaeHuman","CaprinaeHuman","CamelidaeHuman","EquidaeHuman","HumanHuman"))

#check level order again

sealb_min_iden <- min(sealb$identity)

#plot pairwise boxplots
organism_row_plot(sealb_org, sealb_min_iden, sealb_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, sealb, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)

###########################################################################################
###########################################################################################

#Lactotransferrin

lactran <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lactotransferrin_identity.csv", header = T, row.names = 1)
lactran <- as.matrix(lactran)

lactr_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lactotransferrin_groups.csv", header = T)

lactr <- convert_matrix(lactran, lactr_groups, group_column = "group_type", type_column = "allergenic")
head(lactr)
levels(lactr$comparison_type)
summary(lactr$comparison)

lactr_org_order <-c ("Bovinae","Caprinae","Camelidae","Equidae","Human")
lactr_org <-  organism_split(lactr, lactr_org_order)
summary(lactr_org)

#reorder comparison column levels for boxplot
for (i in lactr_org) {
  print(levels(droplevels(i$comparison)))
}

lactr_org$Bovinae[,1] <- factor(lactr_org$Bovinae[,1], levels = c("BovinaeBovinae","BovinaeCaprinae","BovinaeCamelidae","BovinaeEquidae","BovinaeHuman"))
lactr_org$Caprinae[,1] <- factor(lactr_org$Caprinae[,1], levels = c("BovinaeCaprinae","CaprinaeCaprinae","CamelidaeCaprinae","CaprinaeEquidae","CaprinaeHuman"))
lactr_org$Camelidae[,1] <- factor(lactr_org$Camelidae[,1], levels = c("BovinaeCamelidae","CamelidaeCaprinae","CamelidaeCamelidae","CamelidaeEquidae","CamelidaeHuman"))
lactr_org$Equidae[,1] <- factor(lactr_org$Equidae[,1], levels = c("BovinaeEquidae","CaprinaeEquidae","CamelidaeEquidae","EquidaeEquidae","EquidaeHuman"))
lactr_org$Human[,1] <- factor(lactr_org$Human[,1], levels = c("BovinaeHuman","CaprinaeHuman","CamelidaeHuman","EquidaeHuman","HumanHuman"))

#check level order again

lactr_min_iden <- min(lactr$identity)

#plot pairwise boxplots
organism_row_plot(lactr_org, lactr_min_iden, lactr_org_order, puntos = TRUE)

#plot comparison bloxplot
dev.off()
plot(identity ~ comparison_type, lactr, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)



