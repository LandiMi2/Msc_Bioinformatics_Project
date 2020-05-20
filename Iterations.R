#create a matrix of values 
rounds <- matrix(c(00, 27, 27, 27, 27, 01, 60, 40, 91,98, 02, 111, 49, 220, 158, 03, 166, 69, 256, 191, 04,
                   191, 86, 265, 204, 05, 214, 98, 267, 222, 06, 231, 122, 271, 225, 07, 193, 144, 269,
                   234, 08, 195, 118, 274, 236, 09, 195, 118, 276, 246, 10, 183, 125, 278, 248), ncol = 5, byrow = T)

colnames(rounds) <- c("Iterations","Boran", "Freisian", "Ankole", "Ndama")
rounds_df <- as.data.frame(rounds)

# Draw a line graph 
plot(rounds_df$Boran, type = "b", lwd = 2, xaxt = "n", ylim = c(0,300), color = "black", 
     xlab = "Iterations", ylab = "Bovine breeds", main = "IgDiscover Individualized V gene database build-up in each iteration")
axis(1, at=1:length(rounds_df$Iterations), labels = rounds_df$Iterations)
lines(rounds_df$Freisian, col = "blue",  type = "b",lwd = 2)
lines(rounds_df$Ankole, col = "orange",  type = "b",lwd = 2)
lines(rounds_df$Ndama, col = "red",  type = "b",lwd = 2)
legend("topleft", legend = c("Boran", "Freisian", "Ankole", "Ndama"), lty = 1, lwd = 2, 
       pch = 21, col = c("black", "blue", "orange", "red"), ncol = 2, bty = "n", cex = 0.8,
       text.col =  c("black", "blue", "orange", "red"), inset = 0.01, horiz = T)
grid()
