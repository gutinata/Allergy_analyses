#figure for poster

#re-analyzing tropomyosin without fish data
trop <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_identity.csv", header = T, row.names = 1)
trop <- as.matrix(trop)

trop_groups[1:10,]

trop_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_groups.csv", header = T)

#delete fish data from main matrices
trop_nf <- trop[-(64:67),-(64:67)]
trop_groups_nf <- trop_groups[-(64:67),]

tropomyosin_nf <- convert_matrix(trop_nf, trop_groups_nf, group_column = "group_type", type_column = "allergenic")
head(tropomyosin_nf)
levels(tropomyosin_nf$comparison_type)
summary(tropomyosin_nf$comparison)


trop_org_order <-c ("Mollusk","Crustacean","Insect","Mite","Worm","Bird","Mammal")
trop_org <-  organism_split(tropomyosin_nf, trop_org_order)
summary(trop_org)

#reorder comparison column levels for boxplot
for (i in trop_org) {
        print(levels(droplevels(i$comparison)))
}

trop_org$Mollusk[,1] <- factor(trop_org$Mollusk[,1], levels = c("MolluskMollusk","CrustaceanMollusk","InsectMollusk","MiteMollusk","MolluskWorm","BirdMollusk","MammalMollusk"))
trop_org$Crustacean[,1] <- factor(trop_org$Crustacean[,1], levels = c("CrustaceanMollusk","CrustaceanCrustacean","CrustaceanInsect","CrustaceanMite","CrustaceanWorm","BirdCrustacean","CrustaceanMammal"))
trop_org$Insect[,1] <- factor(trop_org$Insect[,1], levels = c("InsectMollusk","CrustaceanInsect","InsectInsect","InsectMite","InsectWorm","BirdInsect","InsectMammal"))
trop_org$Mite[,1] <- factor(trop_org$Mite[,1], levels = c("MiteMollusk","CrustaceanMite","InsectMite","MiteMite","MiteWorm","BirdMite","MammalMite"))
trop_org$Worm[,1] <- factor(trop_org$Worm[,1], levels = c("MolluskWorm","CrustaceanWorm","InsectWorm","MiteWorm","WormWorm","BirdWorm","MammalWorm"))
trop_org$Bird[,1] <- factor(trop_org$Bird[,1], levels = c("BirdMollusk","BirdCrustacean","BirdInsect","BirdMite","BirdWorm","BirdBird","BirdMammal"))
trop_org$Mammal[,1] <- factor(trop_org$Mammal[,1], levels = c("MammalMollusk","CrustaceanMammal","InsectMammal","MammalMite","MammalWorm","BirdMammal","MammalMammal"))

trop_min_iden <- min(tropomyosin_nf$identity)

#plot pairwise boxplots
organism_row_plot(trop_org, trop_min_iden, trop_org_order, puntos = TRUE)

#plot comparison bloxplot
plot(identity ~ comparison_type, tropomyosin_nf, 
     names = c("Allergenic", "Allergenic v. Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL)


#general, summary boxplot
#chosen proteins: profilin, 2s albumin, HSP, tropomyosin

par(mfrow=c(2,2), mar=c(0,3,2,0), oma =c(3,3,1,1), cex.axis=0.8) 
plot(identity ~ comparison_type, prof, 
     names = c("", ""), 
     col = c("#af8dc3", "#7fbf7b"), 
     ylab = "% identity", xlab = NULL,
     ylim = c(15,100),
     main = "Profilin")
abline(h=35, col = "#d7191c", lty="solid", lwd=2)
abline(h=50, col = "#2c7bb6", lty="longdash", lwd=2)

plot(identity ~ comparison_type, heatshockp, 
     names = c("", "", ""), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL,
     ylim = c(15,100),
     main = "Heat Shock Protein")
abline(h=35, col = "#d7191c", lty="solid", lwd=2)
abline(h=50, col = "#2c7bb6", lty="longdash", lwd=2)

plot(identity ~ comparison_type, twoSalb, 
     names = c("Within groups", "Between groups"), 
     col = c("#af8dc3", "#7fbf7b"),
     ylab = "% identity", xlab = NULL,
     ylim = c(15,100),
     main = "2S albumin")
abline(h=35, col = "#d7191c", lty="solid", lwd=2)
abline(h=50, col = "#2c7bb6", lty="longdash", lwd=2)

plot(identity ~ comparison_type, tropomyosin_nf, 
     names = c("Allergenic", "Allergenic v. \n Non-allergenic", "Non-Allergenic"), 
     col = c("white", "#d8b365", "#5ab4ac"), 
     ylab = "% identity", xlab = NULL,
     ylim = c(15,100),
     main = "Tropomyosin")
abline(h=35, col = "#d7191c", lty="solid", lwd=2)
abline(h=50, col = "#2c7bb6", lty="longdash", lwd=2)



