library(ggplot2)
library(ggpubr)
library(dplyr)
library(plotrix)
library(scales)


df <- read.csv('../raw_data/Toms_mean_intensity_change.csv')
df$Mean <- as.numeric(gsub(',', '', as.character(df$Mean)))
df$SEM <- as.numeric(gsub(',', '', as.character(df$SEM)))
# df <- df[df$Marker!='PHF-tau', ]
# df <- df[df$Marker!='Tau', ]

df$Mean[df$Mean==0] <- NA


df_human_pre <- df[df$Group!='PSAPP' & df$PrePost=='Pre', ]
# df_human_pre$Mean[!is.na(df_human_pre$Mean)] <- -1*df_human_pre$Mean[!is.na(df_human_pre$Mean)]

## create plot
df_human_pre %>%
  ggplot(aes(x=Marker, y=Mean, fill=Group)) +
  geom_bar(position="dodge", stat="identity") +
#   scale_y_continuous(breaks = breaks_values,
                    #  labels = abs(breaks_values))+
  theme_classic() +
  scale_fill_manual(values = c("#bf812d", "#35978f")) +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) -1*10^x),
    labels = trans_format("log10", math_format(10^.x)))+
#   scale_y_log10() +
  coord_flip() 
  # geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.2,
                #  position=position_dodge(.9)) + xlab('')

dev.off()



from <- 50

gap.barplot(y=df_human_pre$Mean, gap=c(50, 16000), col=as.numeric(df_human_pre$Marker), 
            xlab="index", ylab="value", horiz=T)

axis.break(1, from, breakcol="snow", style="gap")
axis.break(1, from*(1+0.02), breakcol="black", style="slash")
axis.break(3, from*(1+0.02), breakcol="black", style="slash")
axis(1, at=from)

dev.off()

