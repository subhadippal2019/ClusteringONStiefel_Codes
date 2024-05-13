

library(easyGgplot2)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara") 


myData <- data.frame(Group1=F1[i,j,],
                     Group3=F3[i,j,])
)

dens <- apply(myData, 2, density)

#plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
#mapply(lines, dens, col=1:length(dens))

#legend("topright", legend=names(dens), fill=1:length(dens))



#library("lattice")
#dat <- data.frame(dens = cbind(Group1=F1[i,j,],F3[i,j,]), gg = rep(c("a", "b"), each = length(F1[i,j,])))

#densityplot(~dens,data=dat,groups = gg,
  #          plot.points = FALSE, ref = TRUE, 
    #        auto.key = list(space = "right"))


install.packages("gridExtra")
library(gridExtra)
grid.arrange(myplots, ncol=3, nrow = 3)
library(ggplot2)

#Sample data

myplots=list();
File=  'C:\\Users\\S0PAL001\\Dropbox\\projects\\'
pdf(paste0(File,'EstimatedDensityForF_Group_1_3_burnIN2000_1','.pdf'),width = 9, height=6)
#par(mfrow=c(2,3))
plot_count=1;
for(i in 1:3){
  for(j in 1:2){
    dat <- data.frame(dens = c(Group1=F1[i,j,],Group3=F3[i,j,]), Groups = rep(c("Group1", "Group3"), each = length(F3[i,j,])))
    temp_plot= ggplot(dat, aes(x = dens, fill = Groups,color=Groups))+ geom_density(alpha=0.1)+ scale_color_brewer(palette = "Set1")+scale_fill_manual( values = c("red","blue"))+theme(legend.position="none")
   temp_plot=temp_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #temp_plot=temp_plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
     myplots[[plot_count]]=temp_plot+labs(x = paste0(' (',i,',',j,') th component of F' ),y=" ")
    
    
    plot_count=plot_count+1;
    }
}
 
ggplot2.multiplot(plotlist = myplots,cols=3)
dev.off()




#dev.off()

