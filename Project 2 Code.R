#################Full Simulation with Chi-squared#################

power <- vector() # declare vectors for the output
N_samplesize<- vector()
Odds_ratio_set <- vector()
cats_no_lactose_cancer_set<- vector()
cats_lactose_cancer_set <- vector()

counter<-1

dev_cancer <- 0.2  # set the baseline cancer percent
#developing cancer by drinking lactose will vary depending on the Odds Ratio
#note that the power will also vary based on the baseline cancer rate
#there is greater power when it is close to 0.5 than when it is close to 0 or 1

for (m in seq(3,100,2)){  #the for loop is to simulate multiple sample sizes
  
  samplesize <- m*10 #this is the sample size for each group
  # the total sample size is 2 * m *10
  # so the minimum total sample size is 60, the maximum is 2000
  
  
  for (k in 1:5) { #run a loop to examine 5 diferent effect sizes
    #The numbers are decimals because I expect cats that consume lactose to have a higher chance of cancer
    
    if (k==1) {OR<-0.833}
    
    if (k==2) {OR <- 0.6}
    
    if (k==3) {OR <- 0.4}
    
    if(k==4){OR<-0.25}
    
    if(k==5){OR<-0.1}
    
    #declare the vector to hold the p-values
    # we declare it here because we want it to reset each time
    p_value <- vector()
    for(i in 1:1000){    #run the simulation 1,000 times for the 
      #given Odds ratio and sample size
      
      odds_cancer <- dev_cancer/(1-dev_cancer) #recall that the odds is p/(1-p) 
      
      #we can now use some algebra to solve for male survival
      dev_cancer_lactose <- dev_cancer/((OR*(1-dev_cancer))+ dev_cancer)
      
      
      #let's get the amount of control with cancer
      cats_no_lactose<-as.data.frame(rbinom(n=samplesize, size=1, prob=dev_cancer))
      colnames(cats_no_lactose)[1]<-"cancer"
      cats_no_lactose$treatment <- "no_lactose"
      
      #let's get our amount of cats that drank lactose with cancer
      cats_lactose<-as.data.frame(rbinom(n=samplesize, size=1, prob=dev_cancer_lactose))
      colnames(cats_lactose)[1]<-"cancer"
      cats_lactose$treatment <- "lactose"
      
      #let's combine the two dataframes
      population <- rbind.data.frame(cats_no_lactose,cats_lactose)
      population
      
      #run the chi-squared and get p-value
      X2result<- chisq.test(population$treatment,population$cancer)$p.value
      
      p_value[i] <- X2result
    }
    sum(p_value<0.05)/1000
    
    #put all of the output data in their vectors
    
    power[counter] <- sum(p_value<0.05)/1000
    N_samplesize[counter] <- samplesize*2 
    Odds_ratio_set[counter] <- OR
    cats_no_lactose_cancer_set[counter] <- dev_cancer
    cats_lactose_cancer_set[counter] <- dev_cancer_lactose
    
    counter<- counter+1 #advance the counter
    
  }#end the odds ratio loop
  
}# end the sample size loop



#Now, we need to graph everything

power_simulation<- cbind.data.frame(power, N_samplesize, Odds_ratio_set, cats_no_lactose_cancer_set, cats_lactose_cancer_set)

power_simulation$Odds_ratio_set2 <- as.factor(power_simulation$Odds_ratio_set)
power_simulation$cats_no_lactose_cancer_set2 <- as.factor(power_simulation$cats_no_lactose_cancer_set)

library(ggplot2)

powerplot <- ggplot(power_simulation, aes(x=N_samplesize, y=power, fill=Odds_ratio_set2, group=Odds_ratio_set2)) + 
  geom_point(shape=24,size=3) + 
  geom_line(color="gray", size=0.5)+
  geom_hline(yintercept=0.8) + 
  scale_fill_manual(values=c("goldenrod", "firebrick", "green","blue","black"))+
  xlab("Sample Size(number of cats)") + 
  ylab("Stastical Power") + 
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) + 
  scale_x_continuous(limits=c(0,2000), breaks=seq(0,2000,200))+
  labs(fill = "Odds Ratio") + 
  ggtitle(label = "Power Analysis of Cats Developing Cancer by Drinking Lactose",
          subtitle = "With Base Cancer Rate at 0.2")+
  theme_classic() + 
  theme(plot.title = element_text(size=20,hjust = 0.5),
        plot.subtitle = element_text(size=15,hjust=0.5))

ggsave(plot=powerplot, file='PowerPlot.png', device="png", units="cm", width=30, height=20) #save it using ggsave
