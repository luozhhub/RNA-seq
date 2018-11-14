#!/usr/bin/R

samples <- c("823_control_1.txt", "824_control_2.txt", "825_control_3.txt", "826_2week_1.txt",
 "829_2week_2.txt", "829_2week_3.txt", "830_4week_1.txt", "831_4week_2.txt", "832_4week_3.txt",
 "835_7week_1.txt", "836_7week_2.txt", "839_7week_3.txt", "840_10week_1.txt", "841_10week_2.txt",
  "847_10week_3.txt")
first.sample <- read.delim("823_control_1.txt",header=F,row.names=1)
count.table <- data.frame(first.sample)
for(s in samples[2:length(samples)]){
column <- read.delim(s,header=F,row.names=1)
count.table <- cbind(count.table,s = column)
}
colnames(count.table) <- samples
write.table(count.table,file = "count_table.txt",sep="\t",quote=F)