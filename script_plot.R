library(ggplot2)
library(openxlsx)
library(gridExtra)
# setwd("C:/Users/valva/OneDrive/Bureau/tmp_2delete_in_august/graph/box_plot")

table=read.xlsx('metric_table_all_plot.xlsx')
table = table2[which(table2$chain=="0"),]

name1 = sprintf("CLL (N=%i)",sum(table$group==1))
name0 = sprintf("MBL (N=%i)",sum(table$group==0))

pval1=round(wilcox.test(table$convergence_score[table$group==1], table$convergence_score[table$group==0])$p.value,5)
pval2=round(wilcox.test(table$end_nodes_density[table$group==1], table$end_nodes_density[table$group==0])$p.value,5)
pval3=round(wilcox.test(table$avg_degree[table$group==1], table$avg_degree[table$group==0])$p.value,5)
pval4=round(wilcox.test(table$max_path_length[table$group==1], table$max_path_length[table$group==0])$p.value,5)
pval5=round(wilcox.test(table$avg_distance[table$group==1], table$avg_distance[table$group==0])$p.value,5)
pval6=round(wilcox.test(table$max_muts_length[table$group==1], table$max_muts_length[table$group==0])$p.value,5)

p1=ggplot(data = table, aes(x=factor(group), y=convergence_score,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("relative (reads) convergence") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval1),limits=factor(c(0,1)),labels=c(name0,name1))
p2=ggplot(data = table, aes(x=factor(group), y=end_nodes_density,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("end nodes density") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval2),limits=factor(c(0,1)),labels=c(name0,name1))
p3=ggplot(data = table, aes(x=factor(group), y=avg_degree,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("average degree") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval3),limits=factor(c(0,1)),labels=c(name0,name1))
p4=ggplot(data = table, aes(x=factor(group), y=max_path_length,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("maximal path length") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval4),limits=factor(c(0,1)),labels=c(name0,name1))
p5=ggplot(data = table, aes(x=factor(group), y=avg_distance,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("average distance") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval5),limits=factor(c(0,1)),labels=c(name0,name1))
p6=ggplot(data = table, aes(x=factor(group), y=max_muts_length,fill=as.factor(group))) + geom_boxplot(aes(group=group))+stat_summary(fun=mean,col='red',geom='point',fill="red",shape=20, size=8)+theme(legend.position="none") + scale_fill_manual(values=c("#56B4E9", "#E69F00")) + ylab("maximal mutational length") + scale_x_discrete(name =sprintf("p-val (Wilcoxon test): %g",pval6),limits=factor(c(0,1)),labels=c(name0,name1))
p_all=grid.arrange(p1,p2,p4,p6,p3,p5, nrow = 3)
ggsave('score_box_plot_heavy.pdf',plot=p_all,width=8, height=9)
p_all_2=grid.arrange(p1,p2,p4,p6, nrow = 1)
ggsave('score_box_plot2_heavy.pdf',plot=p_all_2,width=15,height=4)