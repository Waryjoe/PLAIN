#!/usr/bin/env Rscript
source("scripts/plain_paths.R")
options(stringsAsFactors = FALSE)
.libPaths(c(plain_path("R_LIBRARY", ""), .libPaths()))
suppressPackageStartupMessages({library(readxl); library(openxlsx); library(jsonlite); library(digest)})

root <- plain_path("PROJECT", "汇总")
out <- plain_path("HEATMAP", "")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'Source_Data'), showWarnings=FALSE)
dir.create(file.path(out,'Reproducibility'), showWarnings=FALSE)
paths <- c(
  clinical=file.path(root,'训练验证集/训练集_70%_574例.xlsx'),
  radiomics=file.path(root,'影像组学_训练验证集_v5_GRAPE对齐/影像组学v5_训练集_70%_574例_3190特征.xlsx'),
  shap=file.path(root,'figure/figure汇总/Radiomics_SHAP_20260911/Radiomics_SHAP_SourceData_and_QC.xlsx'),
  benchmark=file.path(root,'最终模型_10特征双ROI_LightGBM_全流程评估图表代码包_20260812/01_Tables/双ROI影像组学_机器学习算法比较及临床融合_7_10_15_20_30项.xlsx'),
  final=file.path(root,'论文Figure_最终整理_10特征LightGBM_内外部验证_20260828/03_Source_Data/Final_Figure_Source_Data_10Feature_LightGBM_20260828.xlsx')
)
key <- function(x) sub('^0+','',sub('\\.0$','',trimws(as.character(x))))
checks <- data.frame(Check=character(), Value=character(), Status=character())
check <- function(name, pass, detail='') {
  checks <<- rbind(checks,data.frame(Check=name,Value=as.character(detail),Status=ifelse(pass,'PASS','FAIL')))
  if(!pass) stop(name, ': ', detail)
}
read <- function(p,sheet=1) as.data.frame(read_xlsx(p,sheet=sheet,.name_repair='unique'))
clinical <- read(paths[['clinical']], '临床数据')
cross <- read(paths[['clinical']], '病例映射')
mapping <- read(paths[['shap']], 'Feature_mapping')
mapping <- mapping[match(c(paste0('S',1:5),paste0('T',1:5)),mapping$ID),]
ids <- mapping$ID
headers <- names(read_xlsx(paths[['radiomics']],n_max=0))
types <- rep('skip',length(headers))
types[1:3] <- c('text','text','numeric')
types[match(mapping$Feature,headers)] <- 'numeric'
rad <- as.data.frame(read_xlsx(paths[['radiomics']],col_types=types))
cat('Radiomics metadata columns:',paste(names(rad)[1:3],collapse=', '),'\n')
hospital_column <- which(names(rad)%in%c('hospital_id','住院号','Hospital_ID'))
case_column <- which(names(rad)%in%c('case_id','影像组学case_id','Case_ID'))
check('Radiomics identifier schema',length(hospital_column)==1 && length(case_column)==1)
rad$key <- key(rad[[hospital_column]])
clinical$key <- key(clinical[['住院号']])
cross$key <- key(cross[['住院号']])
check('Unique training hospital identifiers',!anyDuplicated(rad$key)&&!anyDuplicated(clinical$key)&&!anyDuplicated(cross$key))
check('Identical training membership',setequal(rad$key,clinical$key)&&setequal(rad$key,cross$key))
clinical <- clinical[match(rad$key,clinical$key),]
cross <- cross[match(rad$key,cross$key),]
y <- as.integer(rad[[3]])
check('Training count and events',nrow(rad)==574 && sum(y)==44, paste(nrow(rad),sum(y),sep='/'))
check('Clinical outcome agrees by hospital ID',identical(y,as.integer(clinical[['吻合口']])) && identical(y,as.integer(cross[['吻合口']])) )
check('Case IDs agree with independent clinical crosswalk',identical(as.character(rad[[case_column]]),as.character(cross[['影像组学case_id']])))
x <- as.matrix(rad[,mapping$Feature]); colnames(x)<-ids
check('Ten finite features for every patient',identical(dim(x),c(574L,10L)) && all(is.finite(x)))
saved <- read(paths[['shap']],'Training_feature_values')
saved_y <- read(paths[['shap']],'Training_SHAP')
check('Feature matrix reproduces archived SHAP source',max(abs(x-as.matrix(saved[,ids])))<1e-8)
check('Archived training outcomes and row IDs agree',identical(y,as.integer(saved_y$Outcome)) && identical(saved$Sample,saved_y$Sample))
prob <- read(paths[['benchmark']],'Fusion_predictions')
cat('Prediction dataset labels:',paste(unique(prob$Dataset),collapse=' | '),'\n')
prob <- prob[prob$Feature_count==10 & prob$Radiomics_model=='LightGBM' & grepl('Training',prob$Dataset),]
prob$key <- key(prob$Hospital_ID)
check('574 matched frozen fusion OOF predictions',nrow(prob)==574&&!anyDuplicated(prob$key)&&setequal(rad$key,prob$key))
prob <- prob[match(rad$key,prob$key),]
check('Prediction outcomes agree by hospital ID',identical(y,as.integer(prob$Outcome)))
final <- read(paths[['final']],'Internal_predictions')
final <- final[final$Dataset=='Training OOF',]
check('Published-source training OOF scores unchanged',nrow(final)==574 &&
  max(abs(as.matrix(prob[,c('Clinical_probability','Radiomics_probability','Fusion_probability')])-
  as.matrix(final[,c('Clinical_LightGBM','Radiomics_10Feature_LightGBM','Clinical_Radiomics_Fusion')])))<1e-12)

mu<-colMeans(x); sdv<-apply(x,2,sd)
z<-sweep(sweep(x,2,mu,'-'),2,sdv,'/')
forest<-do.call(rbind,lapply(seq_along(ids),function(j){
  warnings<-character()
  fit<-withCallingHandlers(glm(y~z[,j],family=binomial(),control=glm.control(maxit=100,epsilon=1e-12)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  b<-coef(summary(fit))[2,]
  # Fixed-feature univariable associations, not coefficients of the prediction model.
  data.frame(ID=ids[j],Feature=mapping$Feature[j],ROI=mapping$ROI[j],N=length(y),Events=sum(y),
    Mean=mu[j],SD=sdv[j],Beta=b[1],SE=b[2],OR=exp(b[1]),
    Lower95=exp(b[1]-qnorm(.975)*b[2]),Upper95=exp(b[1]+qnorm(.975)*b[2]),P=b[4],
    Converged=fit$converged,Warning=paste(unique(warnings),collapse='; '))
}))
forest$Q_BH<-p.adjust(forest$P,'BH')
forest$Stable<-with(forest,Converged & is.finite(OR) & is.finite(Lower95) & is.finite(Upper95) & Lower95>0 & abs(Beta)<20 & SE<10)
cat('Univariable OR fits:\n');print(forest[,c('ID','OR','Lower95','Upper95','Stable','Warning')],row.names=FALSE)

continuous_columns<-c(Age='年龄',PT='凝血酶原时间PT',TT='凝血酶时间TT',Basophils='嗜碱性粒细胞百分比',AG_ratio='白蛋白球蛋白比值',HDL_C='高密度脂蛋白',LDH='乳酸脱氢酶LDH')
annotations<-data.frame(Sample=saved$Sample,Outcome=y,
  Fusion_probability=prob$Fusion_probability,Radiomics_probability=prob$Radiomics_probability,Clinical_probability=prob$Clinical_probability)
for(nm in names(continuous_columns)) {
  v<-suppressWarnings(as.numeric(clinical[[continuous_columns[nm]]]))
  if(nm!='Basophils')v[v<=0]<-NA_real_
  annotations[[nm]]<-v
}
annotations$Sex<-as.character(clinical[['性别']])
annotations$T_stage_code<-as.character(clinical[['T1']])
annotations$Differentiation_code<-as.character(clinical[['分化程度2']])
cat('Categorical annotation counts:\n')
print(table(annotations$Sex,useNA='always'));print(table(annotations$T_stage_code,useNA='always'));print(table(annotations$Differentiation_code,useNA='always'))
cat('Continuous annotation missing counts:\n');print(colSums(is.na(annotations[,names(continuous_columns)])))
ord<-order(annotations$Fusion_probability,annotations$Sample,method='radix')
annotations$Column_position<-match(seq_len(nrow(annotations)),ord)
private<-data.frame(Sample=saved$Sample,Hospital_ID=rad[[hospital_column]],Case_ID=rad[[case_column]],Outcome=y,Column_position=annotations$Column_position)
source_hashes<-vapply(paths,function(p)digest(file=p,algo='sha256'),character(1))
bundle<-list(x=x,z=z,y=y,mapping=mapping,annotations=annotations,order=ord,forest=forest,checks=checks,paths=paths,source_hashes=source_hashes)
saveRDS(bundle,file.path(out,'Reproducibility','verified_training_data.rds'))
write.csv(private,file.path(out,'Reproducibility','PRIVATE_patient_crosswalk.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(forest,file.path(out,'Source_Data','Univariable_OR_per_SD.csv'),row.names=FALSE)
write.csv(annotations,file.path(out,'Source_Data','Training_annotations.csv'),row.names=FALSE)
write.csv(data.frame(Sample=saved$Sample,Outcome=y,x,check.names=FALSE),file.path(out,'Source_Data','Training_10_features.csv'),row.names=FALSE)
write.csv(checks,file.path(out,'Reproducibility','Checks.csv'),row.names=FALSE)
write_json(list(N=574,events=44,nonevents=530,cohort='Training only',no_new_exclusions=TRUE,
  prediction_models_refitted=FALSE,feature_selection_repeated=FALSE,sort='Frozen fusion OOF probability ascending; Sample as tie-breaker',
  forest='Univariable logistic regression on each raw feature scaled by training SD; Wald 95% CI; descriptive post-selection associations',
  source_sha256=as.list(source_hashes)),file.path(out,'Reproducibility','manifest.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Prepared:',out,'\n')
