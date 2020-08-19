library(reshape2)
library(tidyverse)
library(cowplot)
#read in the dataset

dataset <- read.csv("C:/Users/chard/OneDrive/Desktop/Exp_13_Rescue_Data.Csv")

#change exp number into a factor
#dataset$exp_num <- as.factor(dataset$exp_num)
dataset$Cell.num <- as.factor(dataset$Cell.num)
#dataset$identifier <- as.factor(dataset$identifier)

#melt it all
dataset_melt <- melt(dataset)

#remove the NAs
dataset_melt <- na.omit(dataset_melt)
dataset_melt$Treatment <- factor(dataset_melt$Treatment , levels = rev(levels(dataset_melt$Treatment)))

#subset out the speed, lifetime, and density
comet_speed <- subset(dataset_melt, variable == "speed")
comet_lifetime <- subset(dataset_melt, variable == "lifetime")
comet_density <- subset(dataset_melt, variable == "density")


########
comet_speed_log <- comet_speed
comet_speed_log$value <- log(comet_speed_log$value)
#means of everything
average_comet_speed <- aggregate(comet_speed_log$value, list(comet_speed_log$identifier,comet_speed$condition,comet_speed$exp_num),mean)
average_comet_lifetime <- aggregate(comet_lifetime$value, list(comet_lifetime$identifier,comet_lifetime$condition),mean)

#ggplot out the stuffs
ggplot(average_comet_speed, aes(x = Group.2,y=x)) + geom_boxplot() + geom_point(aes(color = Group.3), size = 2) + theme_cowplot()
ggplot(average_comet_lifetime, aes(x = Group.2,y=x)) + geom_boxplot() + geom_point()

ggplot(comet_speed,aes(y = value)) + geom_boxplot() + geom_point()
comet_one <- subset(comet_speed, condition == "WT" | condition == "WT_Netrin")
comet_one_log <- comet_one
comet_one_log$value <- log(comet_one_log$value)

ggplot(comet_one_log, aes(x = value, group = identifier, fill = condition))+ geom_density(alpha = 0.1) +
  geom_vline(aes(group = identifier, xintercept = median(comet_one_log$value), color = "blue")) + 
  geom_vline(aes(group = identifier, xintercept = mean(comet_one_log$value)))

write.csv(average_comet_speed,"comet_speed.csv")
write.csv(average_comet_lifetime,"comet_lifetime.csv")




########
#all ROIs put together
ggplot(comet_speed, aes(x = Genotype, y = value, fill = Treatment)) + geom_boxplot()

ggplot(comet_speed, aes(x = Genotype, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)

#split by ROI
comet_speed_gc = subset(comet_speed, ROI  == "growth cone")
comet_speed_axon = subset(comet_speed, ROI == "axon")
comet_speed_both = subset(comet_speed, ROI == "both")


ggplot(comet_speed_gc, aes(x = Genotype, y = value, fill = Treatment)) + geom_boxplot()
ggplot(comet_speed_axon, aes(x = Genotype, y = value, fill = Treatment)) + geom_boxplot()
ggplot(comet_speed_both, aes(x = Genotype, y = value, fill = Treatment)) + geom_boxplot()

ggplot(comet_speed_gc, aes(x = Genotype, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)
ggplot(comet_speed_axon, aes(x = Genotype, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)
ggplot(comet_speed_both, aes(x = Genotype, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)

#split by genotype
comet_speed_wt = subset(comet_speed, Genotype == "wt")
comet_speed_t67d = subset(comet_speed, Genotype == "t67d")
comet_speed_t67r = subset(comet_speed, Genotype == "t67r")

ggplot(comet_speed_wt, aes(x = ROI, y = value, fill = Treatment)) + geom_boxplot()
ggplot(comet_speed_t67d, aes(x = ROI, y = value, fill = Treatment)) + geom_boxplot()
ggplot(comet_speed_t67r, aes(x = ROI, y = value, fill = Treatment)) + geom_boxplot()

ggplot(comet_speed_wt, aes(x = ROI, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)
ggplot(comet_speed_t67d, aes(x = ROI, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)
ggplot(comet_speed_t67r, aes(x = ROI, y = value)) + 
  geom_violin(aes(color = Treatment), position = position_dodge(0.9)) + 
  geom_boxplot(aes(color = Treatment),outlier.shape = NA,position = position_dodge(0.9),width=0.1)


########### Wt, untreated vs netrin, paired t test
average_speeds <- comet_speed %>%
  group_by(ROI,Genotype,Cell.num,Treatment) %>%
  summarize(mean_speed = mean(value, na.rm = TRUE))

WT_both_average_speeds <- average_speeds %>%
  filter(Genotype == "wt" & ROI == "both")

t.test(WT_both_average_speeds$mean_speed ~ WT_both_average_speeds$Treatment, paired = TRUE)

############ t67, untreated vs netrin, paired t test, only 3 neurons
average_speeds <- comet_speed %>%
  group_by(ROI,Genotype,Cell.num,Treatment) %>%
  summarize(mean_speed = mean(value, na.rm = TRUE))

t67_both_average_speeds <- average_speeds %>%
  filter(Genotype == "t67d" & ROI == "both" & Cell.num != 4)

t.test(t67_both_average_speeds$mean_speed ~ t67_both_average_speeds$Treatment, paired = TRUE)

############# t test between WT and T67
average_speeds <- comet_speed %>%
  group_by(ROI,Genotype,Cell.num,Treatment) %>%
  summarize(mean_speed = mean(value, na.rm = TRUE))

wt_t67_untreated_both <- average_speeds %>%
  filter(Treatment == "none" & ROI == "both" & Cell.num != 4)

t.test(wt_t67_untreated_both$mean_speed ~ wt_t67_untreated_both$Genotype, paired = FALSE)
