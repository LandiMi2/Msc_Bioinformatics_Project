#create a matrix of values of pre-germline filter 
rounds <- matrix(c(00, 27, 27, 27, 27, 01, 60, 40, 91,98, 02, 111, 49, 220, 158, 03, 166, 69, 256, 191, 04,
                   191, 86, 265, 204, 05, 214, 98, 267, 222, 06, 231, 122, 271, 225, 07, 193, 144, 269,
                   234, 08, 195, 118, 274, 236, 09, 195, 118, 276, 246, 10, 183, 125, 278, 248, 11, 186, 139, 278, 258,
                   12, 186, 140, 279, 262, 13, 186, 145, 276, 266, 14, 186, 150, 279, 264, 15, 186, 151, 277, 273,
                   16, 186, 152, 279, 269, 17, 186, 153, 277, 276, 18, 186, 153, 279, 274, 
                   19, 186, 155, 277, 282, 20, 186, 155, 279, 281, 21,0 ,0 ,0 , 292,
                   22,0 ,0 ,0 ,298, 23,0 ,0 ,0 ,305, 24, 0 ,0 ,0 ,302, 25, 0 ,0 ,0 ,307, 26
                   ,0 ,0 ,0 , 304, 27,0 ,0 ,0 ,307, 28,0 ,0 ,0 ,304, 29,0 ,0 ,0 ,307, 30,0 ,0 ,0 ,304), ncol = 5, byrow = T)

colnames(rounds) <- c("Iterations","Boran", "Freisian", "Ankole", "Ndama")
rounds_df <- as.data.frame(rounds)
tail(rounds_df)
# Draw a line graph 
plot(rounds_df$Boran, type = "b", lwd = 2, xaxt = "n", ylim = c(0,300), color = "black", 
     xlab = "Iterations", ylab = "Number of Individualized V gene", main = "IgDiscover Individualized V gene database build-up after pre-germline filter")
axis(1, at=1:length(rounds_df$Iterations), labels = rounds_df$Iterations)
lines(rounds_df$Freisian, col = "blue",  type = "b",lwd = 2)
lines(rounds_df$Ankole, col = "orange",  type = "b",lwd = 2)
lines(rounds_df$Ndama, col = "red",  type = "b",lwd = 2)
legend("right", legend = c("Boran", "Freisian", "Ankole", "Ndama"), lty = 1, lwd = 2, 
       pch = 21, col = c("black", "blue", "orange", "red"), ncol = 1, bty = "n", cex = 0.8,
       text.col =  c("black", "blue", "orange", "red"), inset = 0.01, horiz = F)
grid()


#create a matrix of values of germline filter 

germline_rounds <- matrix(c(0,0,0,0,0,1,1,7,1,2,2,4,15,1,5,3,12,21,1,5,4,14,27,2,6,5,15,29,3,6,
                            6,16,33,3,6,7,18,40,3,6,8,18,44,2,6,9,18,44,2,6,10,19,43,2,6,11,18,49,2,6,
                            12,18,52,3,7,13,18,56,2,7,14,18,57,2,7,15,18,59,2,7,16,18,58,3,7,
                            17,18,59,2,6,18,18,58,3,6,19,18,59,2,6,20,18,58,3,6,21,0,0,0,6,22,0,0,0,6,
                            22,0,0,0,6,23,0,0,0,6,24,0,0,0,6,25,0,0,0,6,26,0,0,0,6,27,0,0,0,6,28,0,0,0,6,
                            29,0,0,0,6,30,0,0,0,6), ncol = 5, byrow = T)

colnames(germline_rounds) <- c("Iterations","Boran", "Freisian", "Ankole", "Ndama")
germline_rounds_df <- as.data.frame(germline_rounds)
tail(germline_rounds_df)
# Draw a line graph 
plot(germline_rounds_df$Boran, type = "b", lwd = 2, xaxt = "n", ylim = c(0,60), color = "black", 
     xlab = "Iterations", ylab = "Number of Individualized V gene", main = "IgDiscover Individualized V gene database build-up after germline filter")
axis(1, at=1:length(germline_rounds_df$Iterations), labels = germline_rounds_df$Iterations)
lines(germline_rounds_df$Freisian, col = "blue",  type = "b",lwd = 2)
lines(germline_rounds_df$Ankole, col = "orange",  type = "b",lwd = 2)
lines(germline_rounds_df$Ndama, col = "red",  type = "b",lwd = 2)
legend("right", legend = c("Boran", "Freisian", "Ankole", "Ndama"), lty = 1, lwd = 2, 
       pch = 21, col = c("black", "blue", "orange", "red"), ncol = 1, bty = "n", cex = 0.8,
       text.col =  c("black", "blue", "orange", "red"), inset = 0.01, horiz = F)
grid()
