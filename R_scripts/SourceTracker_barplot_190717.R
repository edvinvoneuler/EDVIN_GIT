library(ggplot2)
library(reshape2)

# SourceTracker analysis - filtered results, with and without human calc
# Unfiltered results, Jena samples, 17-06-19

#             unknown       soil   labcontam    gut       skin    plaque  calculus
#                grey     purple    pink       blue      green    yellow    orange
colors.st <- c("#999999","#5E4FA2","#993299","#3288BD","#ABDDA4","#FFFFBF","#FDAE61")

st.all.filt <- read.delim("/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/SourceTracker/sink_predictions.txt")
names(st.all.filt)
st.all.filt <- st.all.filt[,c(1,8,7,6,3,5,4,2)]

st.all.filt.m <- melt(st.all.filt,by="SampleID")

ggplot(st.all.filt.m, aes(SampleID,value,fill=variable))+
  theme(axis.text.x = element_text(angle = 90)) +
  geom_col(position="stack")+
  labs(x="Sample", y="Proportion of source", title="Unfiltered")+
  scale_fill_manual(values=colors.st[1:7])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))