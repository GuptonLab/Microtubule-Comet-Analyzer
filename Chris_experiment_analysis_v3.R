library(reshape2)
library(tidyverse)
library(cowplot)
#read in the dataset

dataset <- read.csv("C:/Users/chard/OneDrive/Desktop/New_Analysis_1.Csv")

#change exp number into a factor
dataset$exp_num <- as.factor(dataset$exp_num)
dataset$cell_num <- as.factor(dataset$cell_num)
#dataset$identifier <- as.factor(dataset$identifier)

#melt it all
dataset_melt <- melt(dataset)

#remove the NAs
dataset_melt <- na.omit(dataset_melt)
dataset_melt$condition <- factor(dataset_melt$condition , levels = rev(levels(dataset_melt$condition)))

#subset out the speed, lifetime, and density
comet_speed <- subset(dataset_melt, variable == "speed")
comet_lifetime <- subset(dataset_melt, variable == "lifetime")
comet_density <- subset(dataset_melt, variable == "density")


########
comet_speed_log <- comet_speed
comet_speed_log$value <- log(comet_speed_log$value)
#means of everything

#mean of each condition and exp num
average_comet_speed <- aggregate(comet_speed$value, list(comet_speed$exp_num, comet_speed$condition),mean)

average_comet_speed_percell <- aggregate(comet_speed$value, list(comet_speed$exp_num, comet_speed$condition,comet_speed$cell_num),mean)

average_comet_speed_percell$key_val <- paste(average_comet_speed_percell$Group.1, average_comet_speed_percell$Group.3)

wt_n_dat <- subset(average_comet_speed_percell, Group.2 == "WT" | Group.2 == "WT_Netrin")

#ggplot out the stuffs

#wt vs wt netrin graph
ggplot(wt_n_dat,aes(x = Group.2, y = x)) + geom_boxplot() + geom_point() + 
  geom_path(aes(group = key_val)) + geom_text(
    aes(label = key_val),
    nudge_x = 0,
    nudge_y = 0.1
  )

#comet speeds
ggplot(comet_speed, aes(x = condition, y = value)) + 
  geom_point(aes(alpha = 0.1), position = position_jitterdodge(dodge.width = 0.75))+
  geom_boxplot(aes(alpha = 0.4)) + theme_cowplot()

#average comet speeds
ggplot(average_comet_speed_percell, aes(x = Group.2, y = x)) + 
  geom_boxplot() + geom_point()


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
