#plots for milk proteins

#Alpha Lactalbumin

lactalbumin_data <- read.csv("~/Documents/UNL/FARRP/Datasets/alpha lactalbumin/alpha_lactalbumin_blasttohuman.csv", header = T)

l_query <- as.factor(lactalbumin_data[,1])
l_hit <- as.factor(lactalbumin_data[,2])
seq_names <- c("NP_001371279.1_alpha-lactalbumin_precursor_[Homo_sapiens]", "pdb|4L41|A_Chain_A_Alpha-lactalbumin", "pdb|3B0I|A_Chain_A_Alpha-lactalbumin", "pdb|3B0O|A_Chain_A_Alpha-lactalbumin", "pdb|1A4V|A_Chain_A_ALPHA-LACTALBUMIN", "pdb|1B9O|A_Chain_A_Protein_(alpha-lactalbumin)", "KAI2565395.1_lactalbumin_alpha_[Homo_sapiens]", "pdb|1I22|A_Chain_A_LYSOZYME_C", "pdb|1DI4|A_Chain_A_LYSOZYME_C", "pdb|1DI5|A_Chain_A_LYSOZYME_C", "NP_000230.1_lysozyme_C_precursor_[Homo_sapiens]", "pdb|2LHM|A_Chain_A_HUMAN_LYSOZYME", "AAA36188.1_lysozyme_precursor_(EC_3.2.1.17)_[Homo_sapiens]", "pdb|1GDW|A_Chain_A_LYSOZYME_C", "AAC63078.1_lysozyme_precursor_[Homo_sapiens]", "pdb|1B7S|A_Chain_A_LYSOZYME", "pdb|1I20|A_Chain_A_LYSOZYME_C", "pdb|1W08|A_Chain_A_LYSOZYME", "pdb|1I1Z|A_Chain_A_LYSOZYME_C", "pdb|1LAA|A_Chain_A_HUMAN_LYSOZYME")
l_hitname <- factor(lactalbumin_data[,3], levels = seq_names)
l_iden <- lactalbumin_data[,4]
l_ev <- lactalbumin_data[,5]

l_data <- cbind(l_query,l_hit, l_hitname, l_iden, l_ev)
colnames(l_data) = colnames(lactalbumin_data)

#diferentes símbolos para cada proteina "query", diferentes colores para cada parámetro
par(mar = c(15, 4, 4, 4) + 0.3)
plot(Percent_identity~hit_name, l_data[1:20,], col = "black", pch = 1, cex = 0.8, axes = FALSE, ylab = "% identity", xlab = "")
axis(side = 1, at = 1:20 , labels = seq_names, las = 2, cex.axis = 0.5)
axis(side = 2, at = pretty(range(l_data[,4])))
points(Percent_identity~hit_name, l_data[21:40,], col = "black", pch = 5)
legend(1.5, 55, title = "Cow Lactalbumin:", legend = c("CAA29664.1", "AAA30615.1"), pch = c(1, 5), col = c("black", "black"), cex = 0.7, box.lty = 0, bg = "transparent")
legend(0.8, 48, legend = "Parameter:", cex = 0.7, box.lty = 0, bg = "transparent")
legend(1.45, 45.5, legend = c("% identity", "E-value"), fill = c("black", "dark grey"), cex = 0.7, box.lty = 0, bg = "transparent")
rect(1.1, 39.3, 5.6, 56)
par(new = TRUE)
plot(E_value~hit_name, l_data[21:40,], col = "dark grey", axes = FALSE, xlab = "", ylab = "", pch = 5)
points(E_value~hit_name, l_data[1:20,], col = "dark grey", pch = 1, cex = 0.8)
axis(side = 4, at = pretty(range(l_data[,5])))
mtext("E-value", side = 4, line = 3)  

###########################################################################################
###########################################################################################

#serum albumin

serum_albumin_data <- read.csv("~/Documents/UNL/FARRP/Datasets/serum albumin/serum_albumin_blasttohuman.csv", header = T)

sa_query <- as.factor(serum_albumin_data[,1])
sa_hit <- as.factor(serum_albumin_data[,2])
seq_names_sa <- c("KAI2534612.1_albumin_[Homo_sapiens]", "AAA98798.1_alloalbumin_Venezia_[Homo_sapiens]", "AAA64922.1_albumin-like,_partial_[Homo_sapiens]", "EAX05669.1_albumin,_isoform_CRA_k_[Homo_sapiens]", "EAX05660.1_albumin,_isoform_CRA_b_[Homo_sapiens]", "KAI2534615.1_albumin_[Homo_sapiens]", "NP_000468.1_albumin_preproprotein_[Homo_sapiens]", "AAH41789.1_ALB_protein_[Homo_sapiens]", "CAA23754.1_serum_albumin_[Homo_sapiens]", "BAG37325.1_unnamed_protein_product_[Homo_sapiens]", "pdb|1TF0|A_Chain_A,_Serum_albumin", "EAX05661.1_albumin,_isoform_CRA_c_[Homo_sapiens]", "AAN17825.1_serum_albumin_[Homo_sapiens]", "CAA23753.1_unnamed_protein_product_[Homo_sapiens]", "AAF01333.1_serum_albumin_precursor_[Homo_sapiens]", "pdb|6ZL1|A_Chain_A,_Albumin", "BAF85444.1_unnamed_protein_product_[Homo_sapiens]", "pdb|5GIX|A_Chain_A,_Serum_albumin", "BAG60674.1_unnamed_protein_product_[Homo_sapiens]", "AAX63425.1_serum_albumin_[Homo_sapiens]", "EAX05666.1_albumin,_isoform_CRA_h_[Homo_sapiens]", "pdb|3CX9|A_Chain_A,_Serum_albumin", "pdb|4N0U|D_Chain_D,_Serum_albumin", "KAI2534609.1_albumin_[Homo_sapiens]", "pdb|1BKE|A_Chain_A,_SERUM_ALBUMIN", "pdb|1AO6|A_Chain_A,_SERUM_ALBUMIN", "pdb|2BXI|A_Chain_A,_SERUM_ALBUMIN", "pdb|2I2Z|A_Chain_A,_Serum_albumin", "CAH18185.1_hypothetical_protein_[Homo_sapiens]", "pdb|6M5D|A_Chain_A,_Serum_albumin", "pdb|2VDB|A_Chain_A,_Serum_Albumin", "pdb|6JE7|A_Chain_A,_Serum_albumin", "pdb|7AAE|AAA_Chain_AAA,_Albumin", "pdb|5Z0B|A_Chain_A,_Serum_albumin", "pdb|1HK2|A_Chain_A,_SERUM_ALBUMIN", "AEE60908.1_albumin,_partial_[Homo_sapiens]", "pdb|1HK3|A_Chain_A,_SERUM_ALBUMIN", "pdb|4K71|A_Chain_A,_Serum_albumin", "BAG60658.1_unnamed_protein_product_[Homo_sapiens]", "AAG35503.1_PRO2619_[Homo_sapiens]", "EAX05668.1_albumin,_isoform_CRA_j_[Homo_sapiens]", "EAX05678.1_albumin,_isoform_CRA_t_[Homo_sapiens]", "EAX05663.1_albumin,_isoform_CRA_e_[Homo_sapiens]", "AAH35969.1_ALB_protein_[Homo_sapiens]", "NP_001341646.2_alpha-fetoprotein_isoform_2_[Homo_sapiens]", "BAG60036.1_unnamed_protein_product_[Homo_sapiens]", "BAG60031.1_unnamed_protein_product_[Homo_sapiens]")
sa_hitname <- factor(serum_albumin_data[,3], levels = seq_names_sa)
sa_iden <- serum_albumin_data[,4]
sa_ev <- serum_albumin_data[,5]

sa_data <- cbind(sa_query,sa_hit, sa_hitname, sa_iden, sa_ev)
colnames(sa_data) = colnames(serum_albumin_data)
summary(sa_query)

#diferentes símbolos para cada proteina "query", diferentes colores para cada parámetro
par(mar = c(15, 4, 4, 4) + 0.3)
plot(Percent_identity~hit_name, sa_data[1:47,], col = "black", pch = 1, cex = 0.8, axes = FALSE, ylab = "% identity", xlab = "")
axis(side = 1, at = 1:47 , labels = seq_names_sa, las = 2, cex.axis = 0.5)
axis(side = 2, at = pretty(range(sa_data[,4])))
points(Percent_identity~hit_name, sa_data[48:94,], col = "black", pch = 5)
legend(0.5, 52, title = "Cow Serum Albumin:", legend = c("AAA51411.1", "CAA76847.1"), pch = c(1, 5), col = c("black", "black"), cex = 0.7, box.lty = 0, bg = "transparent")
rect(0.2, 44, 7.3, 52.5)

#serum albumin E-values not plotted

###########################################################################################
###########################################################################################

#beta lactoglobulin

lactoglobulin_data <- read.csv("~/Documents/UNL/FARRP/Datasets/beta lactoglobulin/beta_lactoglobulin_blasttohuman.csv", header = T)

lg_query <- as.factor(lactoglobulin_data[,1])
lg_hit <- as.factor(lactoglobulin_data[,2])
seq_names_lg <- c("XP_016870272.1_glycodelin_isoform_X7_[åHomo_sapiens]", "KAI2554567.1_progestagen_associated_endometrial_protein_[Homo_sapiens]", "EAW88163.1_progestagen-associated_endometrial_protein_[Homo_sapiens]", "NP_001018059.1_glycodelin_isoform_1_precursor_[Homo_sapiens]", "CAB43305.1_hypothetical_protein_partial_[Homo_sapiens]", "EAW88165.1_progestagen-associated_endometrial_protein_[Homo_sapiens]", "pdb|4R0B|A_Chain_A_Glycodelin", "BAG34756.1_unnamed_protein_product_[Homo_sapiens]", "AAA60147.1_placental_protein_14_[Homo_sapiens]", "EAW88170.1_hCG1654757_partial_[Homo_sapiens]", "KAI2554568.1_progestagen_associated_endometrial_protein_partial_[Homo_sapiens]", "XP_011517047.1_glycodelin_isoform_X1_[Homo_sapiens]", "NP_001018058.1_glycodelin_isoform_2_precursor_[Homo_sapiens]", "BAG65432.1_unnamed_protein_product_[Homo_sapiens]", "XP_011517054.1_glycodelin_isoform_X5_[Homo_sapiens]", "XP_011517050.1_glycodelin_isoform_X2_[Homo_sapiens]", "EAW50553.1_hCG1795014_partial_[Homo_sapiens]")
lg_hitname <- factor(lactoglobulin_data[,3], levels = seq_names_lg)
lg_iden <- lactoglobulin_data[,4]
lg_ev <- lactoglobulin_data[,5]

lg_data <- cbind(lg_query,lg_hit, lg_hitname, lg_iden, lg_ev)
colnames(lg_data) = colnames(lactoglobulin_data)

#diferentes símbolos para cada proteina "query", diferentes colores para cada parámetro
par(mar = c(20, 4, 1, 4) + 0.3)
plot(Percent_identity~hit_name, lg_data[1:17,], col = "black", pch = 1, cex = 0.8, axes = FALSE, ylab = "% identity", xlab = "")
axis(side = 1, at = 1:17 , labels = seq_names_lg, las = 2, cex.axis = 0.5)
axis(side = 2, at = pretty(range(lg_data[,4])))
points(Percent_identity~hit_name, lg_data[18:34,], col = "black", pch = 5)
points(Percent_identity~hit_name, lg_data[35:51,], col = "black", pch = 4)
legend(2.2, 34.5, title = "Cow lactoglobulin:", legend = c("CAA32835.1", "P02754.3", "ACG59280.1"), pch = c(1, 5, 4), col = c("black", "black", "black"), cex = 0.7, box.lty = 0, bg = "transparent")
rect(2, 27.7, 5.5, 34.8)


###########################################################################################
###########################################################################################

#lactotransferrin

lactotransferrin_data <- read.csv("~/Documents/UNL/FARRP/Datasets/lactotransferrin/lactotransferrin_blasttohuman.csv", header = T)

lt_query <- as.factor(lactotransferrin_data[,1])
lt_hit <- as.factor(lactotransferrin_data[,2])
seq_names_lt <- c ("BAG64912.1_unnamed_protein_product_[Homo_sapiens]", 
                   "AAA86665.1_lactoferrin,_partial_[Homo_sapiens]", 
                   "AAB57795.1_lactoferrin,_partial_[Homo_sapiens]", 
                   "KAI4029334.1_lactotransferrin_[Homo_sapiens]", 
                   "NP_001186078.1_lactotransferrin_isoform_2_[Homo_sapiens]", 
                   "BAH14701.1_unnamed_protein_product_[Homo_sapiens]", 
                   "KAI2529273.1_lactotransferrin_[Homo_sapiens]", 
                   "pdb|1B0L|A_Chain_A,_PROTEIN_(LACTOFERRIN)", 
                   "CAA37116.1_unnamed_protein_product,_partial_[Homo_sapiens]", 
                   "KAI4029336.1_lactotransferrin,_partial_[Homo_sapiens]", 
                   "pdb|1BKA|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1N76|A_Chain_A,_Lactoferrin", 
                   "pdb|1SQY|A_Chain_A,_lactoferrin", 
                   "NP_001308051.1_lactotransferrin_isoform_4_[Homo_sapiens]", 
                   "NP_002334.2_lactotransferrin_isoform_1_preproprotein_[Homo_sapiens]", 
                   "BAG52774.1_unnamed_protein_product,_partial_[Homo_sapiens]", 
                   "KAI2529274.1_lactotransferrin,_partial_[Homo_sapiens]", 
                   "pdb|1FCK|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1LCF|A_Chain_A,_LACTOFERRIN", 
                   "pdb|2BJJ|X_Chain_X,_Lactotransferrin", 
                   "AAA36159.1_lactoferrin_[Homo_sapiens]", 
                   "AAB60324.1_lactoferrin_[Homo_sapiens]", 
                   "AAG48753.1_lactoferrin_precursor_[Homo_sapiens]", 
                   "AAW71443.1_lactoferrin_[Homo_sapiens]", 
                   "BAF83548.1_unnamed_protein_product_[Homo_sapiens]", 
                   "EAW64767.1 lactotransferrin [Homo sapiens]", 
                   "AAA58656.1_HLF2,_partial_[Homo_sapiens]", 
                   "pdb|1LFG|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1LFH|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1LFI|A_Chain_A,_LACTOFERRIN", 
                   "AAA59511.1_lactoferrin_[Homo_sapiens]", 
                   "AAH22347.1_Lactotransferrin_[Homo_sapiens]", 
                   "AAS72878.1_growth-inhibiting_protein_12_[Homo_sapiens]", 
                   "ACF19793.1_lactoferrin_[Homo_sapiens]", 
                   "CAA37914.1_precursor_(AA_-19_to_692)_[Homo_sapiens]", 
                   "KAI2529272.1_lactotransferrin_[Homo_sapiens]", 
                   "NP_001308050.1_lactotransferrin_isoform_3_preproprotein_[Homo_sapiens]", 
                   "AAH15822.1_Lactotransferrin_[Homo_sapiens]", 
                   "AAH15823.1_Lactotransferrin_[Homo_sapiens]", 
                   "KAI4029335.1_lactotransferrin_[Homo_sapiens]", 
                   "KAI2529271.1_lactotransferrin_[Homo_sapiens]", 
                   "AAR12276.1_lactoferrin_[Homo_sapiens]", 
                   "ACC95967.1_truncated_lactoferrin_[Homo_sapiens]", 
                   "pdb|1EH3|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1H44|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1L5T|A_Chain_A,_lactoferrin", 
                   "pdb|1H43|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1H45|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1LCT|A_Chain_A,_LACTOFERRIN", 
                   "pdb|2PMS|A_Chain_A,_Lactotransferrin", 
                   "ABF69105.1_lactoferrin,_partial_[Homo_sapiens]", 
                   "pdb|1DSN|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1VFE|A_Chain_A,_HUMAN_LACTOFERRIN", 
                   "pdb|1VFD|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1HSE|A_Chain_A,_LACTOFERRIN", 
                   "pdb|1D4N|A_Chain_A,_TRANSFERRIN", 
                   "pdb|1A8E|A_Chain_A,_SERUM_TRANSFERRIN", 
                   "pdb|1N7W|A_Chain_A,_Serotransferrin", 
                   "pdb|1BP5|A_Chain_A,_PROTEIN_(SERUM_TRANSFERRIN)", 
                   "pdb|1N7X|A_Chain_A,_Serotransferrin", 
                   "pdb|1N84|A_Chain_A,_Serotransferrin", 
                   "pdb|1RYO|A_Chain_A,_Serotransferrin", 
                   "pdb|1B3E|A_Chain_A,_PROTEIN_(SERUM_TRANSFERRIN)", 
                   "pdb|1D3K|A_Chain_A,_SERUM_TRANSFERRIN", 
                   "pdb|1FQE|A_Chain_A,_SEROTRANSFERRIN", 
                   "pdb|1FQF|A_Chain_A,_SEROTRANSFERRIN", 
                   "pdb|1JQF|A_Chain_A,_TRANSFERRIN", 
                   "pdb|1OQG|A Chain A, Serotransferrin", 
                   "pdb|1OQH|A_Chain_A,_Serotransferrin", 
                   "pdb|2O84|X Chain X, Serotransferrin", 
                   "AAF22007.1_PRO1400_[Homo_sapiens]", 
                   "NP_001341632.2_serotransferrin_isoform_2_[Homo_sapiens]", 
                   "BAG58369.1_unnamed_protein_product_[Homo_sapiens]", 
                   "NP_001341633.2_serotransferrin_isoform_3_[Homo_sapiens]", 
                   "AAA61141.1_transferrin,_partial_[Homo_sapiens]", 
                   "pdb|1DTG|A_Chain_A,_TRANSFERRIN", 
                   "pdb|2O7U|A_Chain_A,_Serotransferrin", 
                   "pdb|3FGS|A_Chain_A,_Serotransferrin", 
                   "pdb|2HAV|A_Chain_A,_Serotransferrin", 
                   "pdb|6SOY_C_Chain_C,_Serotransferrin", 
                   "pdb|3QYT|A_Chain_A,_Serotransferrin", 
                   "ABI97197.1_transferrin_[Homo_sapiens]", 
                   "sp|P02787.4|TRFE_HUMAN_Serotransferrin;_Short=Transferrin;_AltName:_Full=Siderophilin",
                   "pdb|7Q1L|A_Chain_A,_Serotransferrin",
                   "pdb|3VE1_B_Chain_B,_Serotransferrin",
                   "pdb|4H0W|A_Chain_A,_Serotransferrin",
                   "AAH59367.1_Transferrin_[Homo_sapiens]",
                   "BAD96475.1_transferrin_variant,_partial_[Homo_sapiens]",
                   "KAI4031630.1_transferrin_[Homo_sapiens]",
                   "NP_001054.2 serotransferrin isoform 1 precursor [Homo sapiens]",
                   "pdb|1SUV|E_Chain_E,_Serotransferrin,_C-lobe",
                   "pdb|2HAU|A_Chain_A,_Serotransferrin",
                   "pdb|3S9L_C_Chain_C,_Serotransferrin",
                   "EAW79166.1_transferrin,_isoform_CRA_c_[Homo_sapiens]",
                   "BAG64724.1_unnamed_protein_product_[Homo_sapiens]",
                   "AAA59992.1_melanotransferrin_[Homo_sapiens]",
                   "NP_005920.2 melanotransferrin isoform 1 preproprotein [Homo sapiens]",
                   "XP_047304106.1_melanotransferrin_isoform_X3_[Homo_sapiens]",
                   "XP_011511152.1_melanotransferrin_isoform_X1_[Homo_sapiens]",
                   "XP_006713706.2_melanotransferrin_isoform_X2_[Homo_sapiens]")
lt_hitname <- factor(lactotransferrin_data[,3], levels = seq_names_lt)
lt_iden <- lactotransferrin_data[,4]
lt_ev <- lactotransferrin_data[,5]

lt_data <- cbind(lt_query,lt_hit, lt_hitname, lt_iden, lt_ev)
colnames(lt_data) = colnames(lactotransferrin_data)

#diferentes símbolos para cada proteina "query", diferentes colores para cada parámetro
par(mar = c(20, 4, 1, 4) + 0.3)
plot(Percent_identity~hit_name, lt_data, col = "black", pch = 1, cex = 0.8, axes = FALSE, ylab = "% identity", xlab = "")
axis(side = 2, at = pretty(range(l_data[,4])))
axis(side = 1, at = 1:100, labels = seq_names_lt, las = 2, cex.axis = 0.5)
#legend(2.2, 34.5, title = "Cow lactotransferrin:", legend = c(""), pch = c(1, 5, 4), col = c("black", "black", "black"), cex = 0.7, box.lty = 0, bg = "transparent")
#rect(2, 27.7, 5.5, 34.8)
par(new = TRUE)
plot(E_value~hit_name, lt_data, col = "dark grey", axes = FALSE, xlab = "", ylab = "", pch = 1, cex = 0.8)
axis(side = 4, at = pretty(range(lt_data[,5])))
mtext("E-value", side = 4, line = 3)  

