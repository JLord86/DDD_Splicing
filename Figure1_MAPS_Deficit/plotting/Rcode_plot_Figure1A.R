## Plot for Figure1A produced in R (via RStudio) (code given for generation of Figure1A, including VEP consequences, splice region, donor site split by nucleotide and grouped positions in the PolyPy for both DDD and ExAC data):
## File used: MAPS_FS10_combined_to_CI.txt

## Calculate and output the confidence intervals
data <- read.table("MAPS_FS10_combined_to_CI.txt", header = T)

ymin <- data$ps_adjusted - 1.96*data$standard_error
ymax <- data$ps_adjusted + 1.96*data$standard_error

data_new = cbind(data, ymin, ymax)
write.table(data_new, "MAPS_FS10_combined_with_CI.txt", sep = '\t')

## Plotting 
library(ggplot2)
cbPalette4 <- c( "#D55E00","#0072B2", "#009E73", "#56B4E9", "#CC79A7", "#E69F00", "#F0E442", "#999999")

data <- read.table("MAPS_FS10_combined_with_CI.txt", header = T)

labelsCQ = c("Nonsense" , "Missense" , "Synonymous" , "", "A-25" , "A-24" , "A-23" , "A-22" , "A-21" , "A-20" , "A-19" , "A-18" , "A-17" , "A-16" , "A-15" , "A-14" , "A-13" , "A-12" , "A-11" , "A-10" , "A-9" , "A-8" , "A-7" , "A-6" , "A-5" , "A-4" , "A-3" , "A-2" , "A-1" , "A" , "A+1" , "A+2" , "A+3" , "A+4" , "A+5" , "A+6" , "A+7" , "A+8" , "A+9" , "A+10" , "", "D-10" , "D-9" , "D-8" , "D-7" , "D-6" , "D-5" , "D-4" , "D-3" , "D-2" , "D-1" , "D" , "D+1" , "D+2" , "D+3" , "D+4" , "D+5" , "D+6" , "D+7" , "D+8" , "D+9" , "D+10", "", "don_A", "don_C", "don_G", "don_T", "", "PyPu", "Other")
limits <- aes(ymin = data$ymin, ymax = data$ymax)

plot <-
   ggplot(data, aes(x=Order, y=ps_adjusted, colour = Set)) + 
   geom_point(aes(size = 0.02)) +
   ylab("MAPS") +
   xlab("Position") + ylim(-0.1,0.2) +
   geom_errorbar(limits, width=0) +
   scale_colour_manual(values = cbPalette4) +
   scale_x_discrete(limits=labelsCQ) +
   theme(axis.text=element_text(size=16), axis.title=element_text(size=16)) +
   theme(legend.text=element_text(size=16))  + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0.9))

   plot + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))+ 
    expand_limits(x=c(0,71))

 ggsave("MAPS_vars_full_CI_FS10.pdf", bg="transparent",   width = 45, height = 15, units = "cm")


