## Code to produce Figure 3 in R 
## Data manually compiled, taking the PPVs calculated for Figure 2 and the clinical classifications of the 38 variants from their recruiting clinicians 

library(ggplot2)

cbPalette4 <- c( "#D55E00","#0072B2", "#009E73", "#56B4E9", "#CC79A7", "#E69F00", "#F0E442", "#999999")



data <- read.table("38_vars_PPV_data_lab_clinclass.txt", header = T)
labelsCQ = c("EP300", "ARID1B", "EFTUD2", "STXBP1", "MBD5", "FOXP1", "KIF11", "GLI3", "SMARCB1", "MEF2C", "ARID1A", "EHMT1'", "CTNNB1", "CAMTA1", "SRCAP", "PAX3", "EP300'", "CREBBP", "DNM1", "CHD7", "TCF4", "KANSL1","THRA", "CACNA1A",  "SCN11A", "TP63", "KIF22", "SMARCE1", "STXBP1'", "SCN2A", "COL9A3",  "DDX3X", "CD96", "POLD1", "RAD21", "TSC2", "ACVR1", "SETBP1")

p <- ggplot(data, aes(gene, PPV))
p <- p+geom_point(data=data, aes(colour=clinician_class, size = 0.4, shape = clinician_class)) +

scale_x_discrete(limits=labelsCQ) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
expand_limits(x= c(0.5,41.5))+
ylim(0,1)+
scale_colour_manual(values = cbPalette4) +

expand_limits(x=-0.25) +
expand_limits(y=0.25)+
theme(axis.text=element_text(size=13), axis.title=element_text(size=20)) + theme(legend.text=element_text(size=20))
    p + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))+ theme(legend.position = "none")
 ggsave("38_genes_illustrator_without_legend_clinclass_shape.svg", width = 30, height = 12, units = "cm")
