library("survminer")
require("survival")
library(survival)
library(sjmisc)
library(ggplot2)
library(forestplot)

###################################


# forest plot

Forest1 <- read.csv("Univariate_Stroma_471_PFS_2Years_067CP_forest.csv", stringsAsFactors = F) #Yale Cohort
np <- ifelse(!is.na(Forest1$Count), paste(Forest1$Count," (",Forest1$Percent,")",sep=""), NA)
tabletext <- cbind(c("Cell Types","\n",Forest1$Variable), 
                   #c("No. of Patients(%)","\n",np),
                   c("Hazard Ratio","\n",Forest1$HR),
                   c("95% CI(from)","\n",Forest1$X95..CI.from), 
                   c("95% CI(to)","\n",Forest1$X95..CI.to), 
                   c("P","\n",Forest1$P.Value),
                   c("P-Adjust","\n",Forest1$P.Adjust..BH.))
library(forestplot)

# Basic forestplot call without attempting to insert expressions into tabletext
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,NA,Forest1$Median), 
           lower=c(NA,NA,Forest1$X95..CI.from), upper=c(NA,NA,Forest1$X95..CI.to),
           xlab="",
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 0.8),
                          title=gpar(cex = 0.8)),
           col=fpColors(box="lightpink2", lines="black", zero = "black"),
           zero=1, cex=0.9, lineheight = "auto", boxsize=1.5, colgap=unit(1,"mm"),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.4)

