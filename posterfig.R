#figure for poster

#chosen proteins: profilin, 2s albumin, HSP, argkin

par(mfrow=c(2,2), mar=c(0,3,2,0), oma =c(3,3,1,1), cex.axis=0.8) 
plot(identity ~ comparison_type, prof, 
     names = c("", ""), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL,
     main = "Profilin")

plot(identity ~ comparison_type, heatshockp, 
     names = c("", "", ""), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL,
     main = "Heat Shock Protein")

plot(identity ~ comparison_type, twoSalb, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL,
     main = "2S albumin")

plot(identity ~ comparison_type, argkin, 
     names = c("Allergenic", "Allergenic v. \n Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL,
     main = "Arginine kinase")

