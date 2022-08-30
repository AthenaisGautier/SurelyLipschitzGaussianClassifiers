u <- seq(0, 2, 0.01)
f1 <- (u<1)*(u)+(u>=1)*(2-u)
df <- data.frame(x=u, f=f1)
F <- matrix(f1, nrow=length(u), ncol=length(u))*matrix(f1, nrow=length(u), ncol=length(u), byrow=TRUE)
contour(F)
plot(diag(F))


library(ggplot2)
library(metR)
library(viridis)

plot1 <- ggplot(df, mapping=aes(x=x, y=f))+
  geom_line()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(bquote(x[i]))+
  ylab(bquote(f(x[i])))+
  ggtitle("1D function tensorized\n")

df2 <- expand.grid(seq(nrow(df)), seq(nrow(df)))
df2$x1 <- u[df2$Var1]
df2$x2 <- u[df2$Var2]
df2$f1 <- f1[df2$Var1]
df2$f2 <- f1[df2$Var2]
df2$f <- df2$f1*df2$f2
df2$Var1 <- df2$Var2 <- NULL

plot2 <- ggplot(df2, mapping=aes(x=x1, y=x2, z=f))+
  geom_contour(mapping=aes(colour = ..level..))+
  theme_bw()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_viridis(option = "magma", 
                       name = "Function value",
                       direction=-1,
                       # here we use guide_colourbar because it is still a continuous scale
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(50, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       ))+
  xlab(bquote(x[1]))+
  ylab(bquote(x[2]))+
  geom_abline(slope=1,lty=2, col="blue")+
  ggtitle("Contour of the\ntensored function")

df3 <- df2[df2$x1==df2$x2, ]
plot3 <- ggplot(df3, mapping=aes(x=x1, y=f))+
  geom_line()+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(bquote(x[1]==x[2]))+
  ylab(bquote(f(x)))+
  geom_hline(yintercept=0, col="blue", lty=2)+
  ggtitle("Tensored function\nalongside main diagonal")

library(ggpubr)
ggarrange(plot1, plot2, plot3, ncol=3, common.legend = TRUE, legend="bottom")
ggsave("./figures_misc/non-linear.png", width=8, height=4)

df4 <- df2
df4$h <- 1-pmax(abs(df4$x1-1), abs(df4$x2-1))
plot4 <- ggplot(df4, mapping=aes(x=x1, y=x2, z=h))+
  geom_contour(mapping=aes(colour = ..level..))+
  theme_bw()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_viridis(option = "magma", 
                       name = "Function value",
                       direction=-1,
                       # here we use guide_colourbar because it is still a continuous scale
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(50, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       ))+
  xlab(bquote(x[1]))+
  ylab(bquote(x[2]))+
  geom_abline(slope=1,lty=2, col="blue")+
  geom_hline(yintercept=1, col="darkgreen", lty=2)+
  ggtitle(expression(atop("Hat function considered", "")))


df5 <- df4[df4$x1==df4$x2, ]
plot5 <- ggplot(df5, mapping=aes(x=x1, y=h))+
  geom_line()+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(bquote(x[1]==x[2]))+
  ylab(bquote(f(x)))+
  geom_hline(yintercept=0, col="blue", lty=2)+
  ggtitle(expression(atop("Hat function considered", "alongside main diagonal")))

df6 <- df4[df4$x2==1, ]
plot6 <- ggplot(df6, mapping=aes(x=x1, y=h))+
  geom_line()+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(bquote(x[1]))+
  ylab(bquote(f(x)))+
  geom_hline(yintercept=0, col="darkgreen", lty=2)+
  ggtitle(expression(atop("Hat function considered", paste("alongside x")[2]==1)))
plot6

ggarrange(plot4, plot5, plot6, ncol=3, common.legend = TRUE, legend="bottom")
ggsave("./figures_misc/linear.png", width=8, height=4)
