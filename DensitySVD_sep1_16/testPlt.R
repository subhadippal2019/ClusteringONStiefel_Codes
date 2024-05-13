# load libraries
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
data = read.table('density1.txt')
   
x=data[,1]
y=data[,2]
z=data[,3]
# build your data.frame
df <- data.frame(x=x, y=y, z=z)

# build color Palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

ggplot(df) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.5) 
  #+ scale_fill_distiller(palette=myPalette, na.value="white") #+ 
  #theme_bw()
# Plot

ggplot(df, aes(x,y, fill=..level..) ) + 
  stat_density_2d( bins=11, geom = "polygon") +
  scale_fill_gradientn(colours = myPalette(11)) +
  theme_minimal() +
  coord_fixed(ratio = 1)


