n <- 5000
x <- sort(rnorm(n,0,1))
u <- pnorm(x,0,1)
u.color <- "#2184C0"
x.color <- "#A821C0"

layout_mat <- matrix(c(2, 0, 1, 3), nrow = 2, ncol = 2, byrow = TRUE)
my_lay <- layout(mat = layout_mat, heights = c(1, 2), widths = c(2, 0.25), respect =TRUE)
par(mar = c(2, 1, 0, 0))
plot(x,pnorm(x,0,1), type="l", ann=F, axes=F, xlim=range(x))
axis(1,c(-5,5),c("",""),line=-0.7)
axis(2,seq(0,1,0.1),seq(0,1,0.1),las=2, line=-0.6)
mtext("X", 1, line=0.5, cex=1.5)
mtext(expression('Y = F(x)'[X]),2, line=2.75, cex=1.5)
points(min(x),0.75,col=u.color,pch=16,cex=2)
arrows(min(x), 0.75, qnorm(0.75,0,1), 0.75, length = 0.15, angle = 30,
       code = 2, col = u.color, lty = par("lty"), lwd = 2)
points(qnorm(0.75,0,1),0,col=x.color,pch=16,cex=2)
arrows(qnorm(0.75,0,1), 0, qnorm(0.75,0,1), 0.75, length = 0.15, angle = 30,
       code = 2, col = x.color, lty = par("lty"), lwd = 2)

par(mar = c(0, 1, 0, 0))
hist(qnorm(u,0,1), ann=F, axes=F, col=x.color)

par(mar = c(2, 0.5, 0, 0))
xhist <- hist(u, plot = FALSE)
barplot(xhist$counts,space = 0, horiz=TRUE, xlab= "Counts", ylab="Bins", ann=F, axes=F, col=u.color)