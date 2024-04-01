## Tropomyosin epitope graphs

trop_epitopes <- data.frame(read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_epitopes.csv", header = T))

#reorganize the levels of the Type variable
trop_epitopes$Type <- factor(trop_epitopes$Type, levels = c("Crustacean","Insect","Mite","Mollusk","Worm","Fish","Bird","Mammal"))

min(trop_epitopes[,3:5])

par(mfrow = c(5, 1), mai = c(1, 0.2, 0.2, 0.2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(5,5,0.2,0.2))
plot(Epitope.1 ~ Type, trop_epitopes, ylab = "% Identity", xlab = FALSE, xaxt = "n", ylim = c(20, 100), lwd = 0.75, las = 1)
points(Epitope.1 ~ Type, trop_epitopes, pch = 20)
plot(Epitope.2 ~ Type, trop_epitopes, ylab = "% Identity", xlab = FALSE, xaxt = "n", ylim = c(20, 100), lwd = 0.75, las = 1)
points(Epitope.2 ~ Type, trop_epitopes, pch = 20)
plot(Epitope.3 ~ Type, trop_epitopes, ylab = "% Identity", xlab = FALSE, xaxt = "n", ylim = c(20, 100), lwd = 0.75, las = 1)
points(Epitope.3 ~ Type, trop_epitopes, pch = 20)
plot(Epitope.4 ~ Type, trop_epitopes, ylab = "% Identity", xlab = FALSE, xaxt = "n", ylim = c(20, 100), lwd = 0.75, las = 1)
points(Epitope.4 ~ Type, trop_epitopes, pch = 20)
plot(Epitope.5 ~ Type, trop_epitopes, ylab = "% Identity", ylim = c(20, 100), lwd = 0.75, las = 1)
points(Epitope.5 ~ Type, trop_epitopes, pch = 20)

xlab = "% Identity", ylab = FALSE, ylim = c(20, 100), lwd = 0.75, las = 1)
