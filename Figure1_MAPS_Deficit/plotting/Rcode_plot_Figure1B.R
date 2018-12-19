## Deficit
## Calculate and output the confidence intervals
## File used: Deficit_data_to_CI.txt

cbPalette <- c("#56B4E9","#009E73", "#0072B2",  "#CC79A7",  "#D55E00", "#E69F00",  "#F0E442",  "#999999")
ratios <- read.table("Deficit_data_to_CI.txt", header = T)

## Work out confidence intervals
numCQtypes<-length(ratios[,1])
CQsingletonnum<-ratios$prop
CQsingletonlower<- rep(0,numCQtypes) #lower limit of 95% confidence interval
CQsingletonupper<- rep(0,numCQtypes) #uppper limit of 95% confidence interval

# loop to calculate 95% confidence intervals for each consequence class
for (i in seq(1,numCQtypes)) {
  
  CQsingletonlower[i]<-quantile(rbinom(10000, ratios$total[i], ratios$prop[i]), probs=0.025) / ratios$total[i]

  CQsingletonupper[i]<-quantile(rbinom(10000, ratios$total[i], ratios$prop[i]), probs=0.975) / ratios$total[i]
  
}


new_ratios <- cbind(ratios, CQsingletonlower, CQsingletonupper)
write.table(new_ratios, "deficit_combined_with_CI.txt", sep='\t')

## Plot
data <- read.table("deficit_combined_with_CI.txt", header = T)
labelsCQ = c("Nonsense" , "Missense" , "Synonymous" , "", "A-25" , "A-24" , "A-23" , "A-22" , "A-21" , "A-20" , "A-19" , "A-18" , "A-17" , "A-16" , "A-15" , "A-14" , "A-13" , "A-12" , "A-11" , "A-10" , "A-9" , "A-8" , "A-7" , "A-6" , "A-5" , "A-4" , "A-3" , "A-2" , "A-1" , "A" , "A+1" , "A+2" , "A+3" , "A+4" , "A+5" , "A+6" , "A+7" , "A+8" , "A+9" , "A+10" , "", "D-10" , "D-9" , "D-8" , "D-7" , "D-6" , "D-5" , "D-4" , "D-3" , "D-2" , "D-1" , "D" , "D+1" , "D+2" , "D+3" , "D+4" , "D+5" , "D+6" , "D+7" , "D+8" , "D+9" , "D+10", "", "don_A", "don_C", "don_G", "don_T", "", "PyPu", "Other")
limits <- aes(ymin = data$CQsingletonlower, ymax = data$CQsingletonupper)

plot <-
   ggplot(data, aes(x=Order, y=prop, colour = "#56B4E9")) + 
   geom_point(aes(size = 0.02)) +
   ylab("Proportion of parental \nvariants in high pLI genes") +
   xlab("Position") + ylim(0,0.4) +
   geom_errorbar(limits, width=0) +
   scale_colour_manual(values = cbPalette) +
   scale_x_discrete(limits=labelsCQ) +
   theme(axis.text=element_text(size=16), axis.title=element_text(size=16)) +
   theme(legend.text=element_text(size=16))  + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0.9))

   plot + theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))+ 
    expand_limits(x=c(0,71))

 ggsave("proportion_parental_vars_full_CI_noNA_FS10.pdf", bg="transparent",   width = 45, height = 15, units = "cm")
