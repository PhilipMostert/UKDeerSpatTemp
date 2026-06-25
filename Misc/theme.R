library(ggplot2)
library(showtext)
showtext_auto()
f1 <- "Open sans"
font_add_google(f1, f1)

themePlot <- theme_classic() +
  theme(#plot.tag = element_text(angle = 90, hjust = 0.5),
    text = element_text(family=f1, size = 90),
    #plot.tag.position = c(-0.03, 0.45),
    # axis.text = element_text(size = 18),
    # axis.title = element_text(size = 24),
    # legend.text = element_text(size = 20),
    # legend.title = element_text(size = 24, face = 'bold'),
    strip.text= element_text(size = 15*4),
    axis.text.x = element_text(size = 12*4),
    axis.text.y = element_text(size = 12*4),
    axis.title = element_text(size = 15*5),
    legend.text = element_text(size = 12*4),
    legend.title = element_text(size = 15*4),
    plot.tag.position = c(0,1),
    plot.tag = element_text(face = 'bold', hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'mm')),
    plot.title = element_text(hjust = 0.5, size = 20))

