meta<-read.csv('../data/metab_meta.csv',header=TRUE);
data<-read.csv('../data/metab_conc_data.csv',header=TRUE);
meta<-meta[meta[,'Variable.Name'] %in% colnames(data),];
meta<-meta[,c('Variable.Name','Final_Group')];
colnames(meta)<-c('variables','module');
write.csv(meta,'../results/conc/biological/metab_bio_group.csv', row.names=FALSE);