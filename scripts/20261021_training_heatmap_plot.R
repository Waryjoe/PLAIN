#!/usr/bin/env Rscript
source("scripts/plain_paths.R")
options(stringsAsFactors=FALSE)
.libPaths(c(plain_path("R_LIBRARY", ""),.libPaths()))
suppressPackageStartupMessages({library(grid);library(svglite);library(ragg);library(openxlsx);library(jsonlite);library(digest)})
base_out <- plain_path("HEATMAP", "")
args<-commandArgs(trailingOnly=TRUE)
t4_sensitivity <- '--t4-sensitivity' %in% args
vertical_legends <- t4_sensitivity || '--vertical-legends' %in% args
out <- if(t4_sensitivity)file.path(base_out,'T4_sensitivity_v3') else if(vertical_legends)file.path(base_out,'Vertical_legends_v2') else base_out
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(file.path(base_out,'Reproducibility','verified_training_data.rds'))
if(t4_sensitivity) {
  original_forest<-d$forest
  excluded<-d$x[,'T4']==1e6
  stopifnot(sum(excluded)==1L,d$annotations$Sample[excluded]=='PRIVATE_CASE_003')
  keep<-!excluded
  t4_mean<-mean(d$x[keep,'T4']);t4_sd<-sd(d$x[keep,'T4'])
  t4_z<-rep(NA_real_,nrow(d$x));t4_z[keep]<-(d$x[keep,'T4']-t4_mean)/t4_sd
  fit<-glm(d$y[keep]~t4_z[keep],family=binomial(),control=glm.control(maxit=100,epsilon=1e-12))
  co<-coef(summary(fit))[2,];beta<-unname(co[1]);se<-unname(co[2])
  ci<-exp(beta+c(-1,1)*qnorm(.975)*se)
  stopifnot(isTRUE(fit$converged),all(is.finite(c(exp(beta),ci))))
  j<-match('T4',d$forest$ID)
  values<-list(N=sum(keep),Events=sum(d$y[keep]),Mean=t4_mean,SD=t4_sd,
    Beta=beta,SE=se,OR=exp(beta),Lower95=ci[1],Upper95=ci[2],P=unname(co[4]),
    Converged=fit$converged,Warning='',Stable=TRUE)
  for(nm in names(values))d$forest[j,nm]<-values[[nm]]
  d$forest$Q_BH<-p.adjust(d$forest$P,method='BH')
  d$z[,'T4']<-t4_z
  d$t4_sensitivity<-list(excluded=excluded,rule='T4 equals 1000000; feature-level sensitivity analysis only',original_forest=original_forest)
}
ord<-d$order; ann<-d$annotations[ord,]; z<-d$z[ord,]; f<-d$forest
W<-if(vertical_legends)230 else 220; H<-180
ink<-'#222D35'; muted<-'#65727C'; missing<-'#D4D8DC'
excluded_color<-'#E9BD73'
pal<-list(Fusion=c('#F5F0F5','#826181'),Radiomics=c('#EFF8F6','#2A9D8F'),Clinical=c('#EEF4F8','#457B9D'))
continuous_palette<-c('#4B9DB3','#FAFAFA','#D6757E')
sex_palette<-c('女'='#9FC9D5','男'='#47768E')
t_palette<-c('1'='#F4D7B2','2'='#E9B27D','3'='#D87B53','4'='#AE433E')
grade_palette<-c('1'='#DCD6E9','2'='#A994BE','3'='#6B5385')
leak_palette<-c('0'='#EDF0F2','1'='#725082')
font<-'Arial'
colour <- function(v,cols,lim=c(0,1)) {
  idx<-1+round(255*pmin(1,pmax(0,(v-lim[1])/diff(lim))))
  result<-colorRampPalette(cols)(256)[idx]; result[is.na(v)]<-missing; result
}
text <- function(label,x,y,size=8,face='plain',col=ink,just='left',rot=0) {
  grid.text(label,x=unit(x,'mm'),y=unit(H-y,'mm'),just=just,rot=rot,
    gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col))
}
rect <- function(x,y,w,h,fill,border=NA,lwd=.3) {
  grid.rect(x=unit(x,'mm'),y=unit(H-y,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),
    just=c('left','top'),gp=gpar(fill=fill,col=border,lwd=lwd))
}
line <- function(x1,y1,x2,y2,col=ink,lwd=.5,lty=1) {
  grid.lines(x=unit(c(x1,x2),'mm'),y=unit(H-c(y1,y2),'mm'),gp=gpar(col=col,lwd=lwd,lty=lty))
}
gradient <- function(label,x,y,cols,low='0',high='1',middle=NULL,width=32) {
  text(label,x,y,size=7.6,face='bold')
  n<-128; rect(x+seq(0,width,length.out=n+1)[1:n],y+2,width/n+.002,2.2,colorRampPalette(cols)(n))
  text(low,x,y+7,size=7);text(high,x+width,y+7,size=7,just='right')
  if(!is.null(middle))text(middle,x+width/2,y+7,size=7,just='center')
}
key <- function(label,col,x,y) {rect(x,y-1.25,2.5,2.5,col);text(label,x+4,y,size=7.6)}
vertical_gradient <- function(label,x,y,cols,low='0',high='1',middle='0.5') {
  text(label,x,y,8,'bold')
  by<-y+8; bh<-24; bw<-3.2; n<-128
  rect(x,by+(0:(n-1))*bh/n,bw,bh/n+.003,rev(colorRampPalette(cols)(n)))
  rect(x,by,bw,bh,NA,border=muted,lwd=.3)
  text(high,x+bw+1.5,by,8)
  text(middle,x+bw+1.5,by+bh/2,8)
  text(low,x+bw+1.5,by+bh,8)
}
vertical_categories <- function(label,x,y,labels,cols) {
  text(label,x,y,8,'bold')
  by<-y+6; rh<-4.5; bw<-3.2
  for(i in seq_along(labels)) {
    rect(x,by+(i-1)*rh,bw,rh,unname(cols[i]))
    text(labels[i],x+bw+1.5,by+(i-.5)*rh,8)
  }
  rect(x,by,bw,length(labels)*rh,NA,border=muted,lwd=.3)
}
raw_ann <- d$annotations
cont_names<-c('Age','PT','TT','Basophils','AG_ratio','HDL_C','LDH')
ann_z<-sapply(cont_names,function(nm){v<-raw_ann[[nm]];(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)})
display_tracks<-list(
  list(label='Fusion probability',value=ann$Fusion_probability,colors=colour(ann$Fusion_probability,pal$Fusion)),
  list(label='Radiomics probability',value=ann$Radiomics_probability,colors=colour(ann$Radiomics_probability,pal$Radiomics)),
  list(label='Clinical probability',value=ann$Clinical_probability,colors=colour(ann$Clinical_probability,pal$Clinical)),
  list(label='Age',value=ann$Age,colors=colour(ann_z[ord,'Age'],continuous_palette,c(-2,2))),
  list(label='Sex',value=ann$Sex,colors=unname(sex_palette[ann$Sex])),
  list(label='T stage',value=ann$T_stage_code,colors=unname(t_palette[ann$T_stage_code])),
  list(label='Differentiation',value=ann$Differentiation_code,colors=unname(grade_palette[ann$Differentiation_code])),
  list(label='Prothrombin time',value=ann$PT,colors=colour(ann_z[ord,'PT'],continuous_palette,c(-2,2))),
  list(label='Thrombin time',value=ann$TT,colors=colour(ann_z[ord,'TT'],continuous_palette,c(-2,2))),
  list(label='Basophils (%)',value=ann$Basophils,colors=colour(ann_z[ord,'Basophils'],continuous_palette,c(-2,2))),
  list(label='A/G ratio',value=ann$AG_ratio,colors=colour(ann_z[ord,'AG_ratio'],continuous_palette,c(-2,2))),
  list(label='HDL-C',value=ann$HDL_C,colors=colour(ann_z[ord,'HDL_C'],continuous_palette,c(-2,2))),
  list(label='LDH',value=ann$LDH,colors=colour(ann_z[ord,'LDH'],continuous_palette,c(-2,2))),
  list(label='Anastomotic leakage',value=ann$Outcome,colors=unname(leak_palette[as.character(ann$Outcome)]))
)
draw <- function() {
  grid.newpage(); rect(0,0,W,H,'white')
  text(if(t4_sensitivity)'Training cohort: T4 feature-level sensitivity analysis' else 'Training cohort: clinical annotations and radiomic features',8,9,11,'bold')
  text('574 patients | 44 leakage events | Columns ordered by frozen fusion OOF probability',8,15,8.2,col=muted)
  x0<-8; hw<-96; cw<-hw/574
  xx<-x0+(0:573)*cw
  for(k in seq_along(display_tracks)) {
    yy<-23+(k-1)*3.65+ifelse(k>3,1.5,0)
    tr<-display_tracks[[k]]
    rect(xx,yy,cw+.001,3.15,tr$colors)
    text(tr$label,107,yy+1.57,7.8)
  }
  text('Selected radiomic features',8,81,8.5,'bold')
  text('OR (95% CI)',153,81,8.5,'bold',just='center')
  top<-85; rh<-5.9
  upper<-if(t4_sensitivity)8 else 4
  ticks<-if(t4_sensitivity)c(.25,.5,1,2,4,8) else c(.25,.5,1,2,4)
  fx<-function(v)136+(log(v)-log(.125))/(log(upper)-log(.125))*34
  for(tick in ticks)line(fx(tick),top,fx(tick),top+10*rh,col=ifelse(tick==1,'#8B959D','#E5E9EC'),lwd=ifelse(tick==1,.65,.4),lty=ifelse(tick==1,2,1))
  for(j in seq_len(ncol(z))) {
    yy<-top+(j-1)*rh
    cell_colors<-colour(z[,j],continuous_palette,c(-2,2))
    if(t4_sensitivity && colnames(z)[j]=='T4')cell_colors[excluded[ord]]<-excluded_color
    rect(xx,yy,cw+.001,rh-.15,cell_colors)
    label<-if(t4_sensitivity && colnames(z)[j]=='T4')'T4* (n=573)' else paste0(colnames(z)[j],ifelse(f$Stable[j],'','*'))
    text(label,107,yy+rh/2,if(t4_sensitivity && colnames(z)[j]=='T4')8 else 8.5)
    if(j==5)line(8,yy+rh,170,yy+rh,col='white',lwd=1.2)
    if(f$Stable[j]) {
      co<-if(f$OR[j]<1)'#338E99' else '#C77373'
      line(fx(f$Lower95[j]),yy+rh/2,fx(f$Upper95[j]),yy+rh/2,col=co,lwd=.9)
      for(b in c(f$Lower95[j],f$Upper95[j]))line(fx(b),yy+rh/2-.8,fx(b),yy+rh/2+.8,col=co,lwd=.65)
      rect(fx(f$OR[j])-1,yy+rh/2-1,2,2,fill='white',border=co,lwd=.8)
      rect(fx(f$OR[j])-.45,yy+rh/2-.45,.9,.9,fill=co)
    } else text('OR overflow*',153,yy+rh/2,7.2,col=muted,just='center')
  }
  bottom<-top+10*rh
  line(136,bottom,170,bottom,lwd=.65)
  for(tick in ticks) {
    line(fx(tick),bottom,fx(tick),bottom+1,lwd=.5)
    text(format(tick,trim=TRUE),fx(tick),bottom+3.4,7.5,just='center')
  }
  text('Odds ratio per 1 SD (log scale)',153,bottom+8.5,7.3,just='center')
  text('Lower fusion probability',8,bottom+4,7.4,col=muted)
  text('Higher fusion probability',104,bottom+4,7.4,col=muted,just='right')
  line(8,bottom+7,104,bottom+7,col=muted,lwd=.5)
  line(104,bottom+7,102.4,bottom+6.3,col=muted,lwd=.5)
  line(104,bottom+7,102.4,bottom+7.7,col=muted,lwd=.5)
  text(if(vertical_legends)'S1-S5: stomach features    T1-T5: tumor features (not T stage)' else 'S1-S5: stomach ROI    T1-T5: tumor ROI',8,158,8)
  text('Univariable associations of previously selected features; not adjusted or causal effects.',8,164,7.3,col=muted)
  text('Colors: row-wise z-scores, saturated at -2 / +2; probabilities remain on the 0-1 scale.',8,168.4,7.3,col=muted)
  text(if(t4_sensitivity)'*T4 only: one value of 1e6 excluded; n=573 (43 events). Other rows/scores unchanged.' else '*T4 (coarseness): extreme values cause per-SD OR overflow; all patients retained.',8,172.8,7.3,col=muted)

  if(vertical_legends) {
    vertical_gradient('Fusion\nprobability',181,25,pal$Fusion)
    vertical_gradient('Radiomics\nprobability',207,25,pal$Radiomics)
    vertical_gradient('Clinical\nprobability',181,67,pal$Clinical)
    vertical_gradient('Row z-score',207,67,continuous_palette,low='-2',high='+2',middle='0')
    vertical_categories('Sex',181,109,c('Female','Male'),sex_palette)
    vertical_categories('T stage',207,109,paste0('T',1:4),t_palette)
    vertical_categories('Differentiation\ncode',181,141,as.character(1:3),grade_palette)
    vertical_categories('Leakage',207,141,c('No','Yes'),leak_palette)
    key('Missing value',missing,181,170)
    if(t4_sensitivity)key('T4 excluded',excluded_color,207,170)
  } else {
  lx<-181
  gradient('Fusion probability',lx,24,pal$Fusion)
  gradient('Radiomics probability',lx,37,pal$Radiomics)
  gradient('Clinical probability',lx,50,pal$Clinical)
  gradient('Row z-score',lx,63,continuous_palette,low='-2',high='+2',middle='0')
  text('Sex',lx,79,8,'bold');key('Female',sex_palette['女'],lx,84);key('Male',sex_palette['男'],lx,89)
  text('T stage',lx,99,8,'bold')
  for(j in 1:4)key(paste0('T',j),t_palette[as.character(j)],lx,100+5*j)
  text('Differentiation code',lx,129,8,'bold')
  for(j in 1:3)key(as.character(j),grade_palette[as.character(j)],lx,130+5*j)
  text('Leakage',lx,154,8,'bold');key('No',leak_palette['0'],lx,159);key('Yes',leak_palette['1'],lx+17,159)
  key('Missing clinical value',missing,lx,168)
  }
}
stem<-file.path(out,'Training574_Heatmap_Forest')
svglite(paste0(stem,'.svg'),width=W/25.4,height=H/25.4,bg='white');draw();dev.off()
pdf(paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family='Helvetica',useDingbats=FALSE)
font<-'Helvetica';draw();dev.off();font<-'Arial'
agg_png(paste0(stem,'_preview.png'),width=W,height=H,units='mm',res=220,background='white');draw();dev.off()
agg_png(paste0(stem,'_600dpi.png'),width=W,height=H,units='mm',res=600,background='white');draw();dev.off()

if(t4_sensitivity) {
  dir.create(file.path(out,'Reproducibility'),showWarnings=FALSE)
  dir.create(file.path(out,'Source_Data'),showWarnings=FALSE)
  saveRDS(d,file.path(out,'Reproducibility','t4_sensitivity_data.rds'))
  t4_table<-data.frame(Sample=d$annotations$Sample,Outcome=d$y,
    T4_original=d$x[,'T4'],T4_analysis=ifelse(excluded,NA,d$x[,'T4']),
    Excluded_T4_only=excluded,T4_z_analysis=d$z[,'T4'],
    Column_position=d$annotations$Column_position)
  write.csv(t4_table,file.path(out,'Source_Data','T4_analysis_values.csv'),row.names=FALSE,na='')
  write.csv(f,file.path(out,'Source_Data','Forest_sensitivity.csv'),row.names=FALSE)
  comparison<-rbind(original_forest[original_forest$ID=='T4',],f[f$ID=='T4',])
  comparison<-data.frame(Analysis=c('Original_all_values','T4_1e6_excluded'),comparison,row.names=NULL)
  for(nm in c('OR','Lower95','Upper95'))comparison[[nm]]<-ifelse(is.finite(comparison[[nm]]),format(comparison[[nm]],digits=12),'Overflow')
  write.csv(comparison,file.path(out,'Source_Data','T4_before_after.csv'),row.names=FALSE)
  caption<-c('Figure. Training-cohort feature heatmap and univariable odds ratios: T4 feature-level sensitivity analysis.',
    'All 574 training patients (44 leakage events, 530 non-events) remain as columns in their original ascending frozen fusion meta-OOF probability order. The sole T4 value equal to 1000000 (pseudonymous sample PRIVATE_CASE_003) is excluded only from T4 calculations and its heatmap cell; the cell is marked amber. T4 therefore uses 573 patients (43 leakage events, 530 non-events). Other feature rows, clinical annotations, outcomes and archived model probabilities are unchanged.',
    'T4 is tumor exponential NGTDM Coarseness, not clinical T stage. Excluding 1e6 is an explicitly requested post-hoc sensitivity analysis, not proof of source-data error. PyRadiomics documents a 1e6 fallback when the Coarseness denominator is zero; the exact extraction condition for this case has not been reconstructed. No additional patient or value exclusions were searched.',
    'T4 mean and sample SD are recalculated from the 573 retained values, then used for its z-score heatmap and separate univariable logistic regression. Its OR per one retained-value SD and Wald 95% CI are shown. The other nine feature estimates are unchanged. The forest axis is expanded to 8 to include the complete T4 confidence interval. Heatmap colors remain saturated at z=-2/+2; model probabilities remain on their original 0-1 scales. Missing clinical values are grey. The original value is preserved in the source CSV.',
    sprintf('T4: OR %.4f (95%% CI %.4f-%.4f), N=%d, events=%d; SD %.8f. These values use a different SD than the original full-data per-SD OR.',f$OR[j],f$Lower95[j],f$Upper95[j],f$N[j],f$Events[j],f$SD[j]),
    'These are training-set descriptive associations of previously selected features, not independent validation or causal effects. Wald intervals, P-values and BH q-values do not account for prior feature selection or this post-hoc exclusion. No archived predictive model has been refit; the figure is not evidence of improved AUC, calibration or DCA. T-stage and differentiation source coding retain their original clinical-timing limitations. SHAP is omitted.',
    'PyRadiomics source: https://pyradiomics.readthedocs.io/en/latest/_modules/radiomics/ngtdm.html#RadiomicsNGTDM.getCoarsenessFeatureValue')
  writeLines(caption,file.path(out,'Figure_legend_EN.txt'))
  writeLines(c('# T4 feature-level sensitivity analysis','',
    'Native editable AI: Training574_Heatmap_Forest.ai. SVG/PDF and previews are companions.',
    'Scope: only the PRIVATE_CASE_003 T4 value of 1000000 is excluded for this sensitivity analysis. The patient remains in the full 574-column heatmap and every other feature/model score. Original raw files and prior figures are unchanged.',
    'T4 is standardized and fit on N=573, with 43 leakage events. The excluded cell is amber, not imputed. This is not an adjudicated data correction or a model retraining.',
    'Source_Data contains original and analysis T4 values, exclusion status, before/after statistics and all ten forest results. The original complete data workbook remains in ../Source_Data.',
    'Keep this version together with the original analysis; see Figure_legend_EN.txt for full methods and limitations.'),file.path(out,'README.md'))
  file.copy('scripts/20261021_training_heatmap_plot.R',file.path(out,'Reproducibility'),overwrite=TRUE)
  cat('Exported T4 sensitivity:',stem,'\n')
  print(f[f$ID=='T4',])
  quit(save='no',status=0)
}

if(vertical_legends) {
  dir.create(file.path(out,'Reproducibility'),showWarnings=FALSE)
  feature<-d$x[,'T4']; idx<-which.max(feature)
  diagnostic<-list(
    feature='tumor__exponential_ngtdm_Coarseness',feature_id='T4',
    clinical_T_stage_4_count=sum(d$annotations$T_stage_code=='4'),
    n=length(feature),missing=sum(is.na(feature)),median=median(feature),
    maximum=max(feature),maximum_sample=d$annotations$Sample[idx],
    mean=mean(feature),sd=sd(feature),log_OR=d$forest$Beta[d$forest$ID=='T4'],
    log_OR_SE=d$forest$SE[d$forest$ID=='T4'],
    implementation_reference='https://pyradiomics.readthedocs.io/en/latest/_modules/radiomics/ngtdm.html#RadiomicsNGTDM.getCoarsenessFeatureValue',
    interpretation='The 1e6 value matches the documented zero-denominator fallback in PyRadiomics NGTDM Coarseness. The exact patient-level extraction condition has not been reconstructed. The per-SD OR overflows; this is not a missing feature or clinical T-stage detection failure.',
    changes='Legend layout and clarification only; no data, fit, patients, probabilities, sorting or color limits changed.')
  write_json(diagnostic,file.path(out,'T4_diagnostic.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
  old_caption<-readLines(file.path(base_out,'Figure_legend_EN.txt'))
  writeLines(c(old_caption,'Revision: the continuous color scales and discrete category legends are vertical. The radiomic identifier T4 is distinct from clinical T stage 4. The T4 value of 1e6 is consistent with the documented PyRadiomics zero-denominator fallback, but the case-specific extraction condition has not been reconstructed.'),file.path(out,'Figure_legend_EN.txt'))
  writeLines(c('# Training574 heatmap: vertical legends v2','',
    'Open Training574_Heatmap_Forest.ai in Adobe Illustrator. SVG and PDF companions retain vector cells and text.',
    'Only legend layout and a feature-ID clarification changed. All 574 patients and archived scores, associations, sorting and color ranges are unchanged.',
    'T4 denotes tumor exponential NGTDM Coarseness, not clinical T4 stage. The clinical stage-4 category is present (396 patients).',
    'There are no missing T4 feature values. One value is 1000000 (PRIVATE_CASE_003), versus median 0.003199595. The per-SD log-OR is finite but exponentiating it overflows. The 1e6 value is consistent with a documented PyRadiomics fallback, not sufficient evidence of a bad CT or incorrect outcome.',
    'For exact source data use ../Source_Data/Training574_Figure_Source.xlsx. No duplicate patient source files were created. The original figure files were not overwritten.',
    'See T4_diagnostic.json and Figure_legend_EN.txt for provenance and limitations.'),file.path(out,'README.md'))
  file.copy('scripts/20261021_training_heatmap_plot.R',file.path(out,'Reproducibility'),overwrite=TRUE)
  cat('Exported vertical legends:',stem,'\n')
  quit(save='no',status=0)
}

feature_key<-d$mapping[,c('ID','ROI','Feature','Short_name','Training_mean','Training_SD')]
limits<-data.frame(Item=c('Cohort','Column order','Probabilities','Feature rows','Heatmap','Forest','T4 limitation','Clinical annotation missingness','Basophils zero','Inference','T stage provenance','Differentiation','Reference'),
  Description=c('All 574 training patients: 44 leakage, 530 no leakage; no exclusions.',
    'Ascending archived fusion meta-OOF probability; deterministic Sample tie-breaker; outcomes not used for sorting.',
    'Three saved OOF probability columns; no recalibration, refitting or edits.',
    'Ten previously locked dual-ROI model features, S1-S5 stomach and T1-T5 tumor; no new feature selection.',
    'Continuous clinical and radiomic values standardized separately by training observed-value mean and sample SD. Colors saturated at +/-2 only; raw values unchanged. Probabilities use 0-1.',
    'Separate univariable logistic regression for each raw radiomic feature standardized by training SD; Wald 95% CI. P and BH q across ten displayed tests included in source data, not selection-adjusted.',
    'T4 has an extreme value of 1000000; conventional per-SD OR and CI overflow and a numerical 0/1 probability warning occurs. Log-OR and SE remain in the source table; the forest row is marked OR overflow. No substitution, winsorization or patient removal.',
    'Clinical annotations use observed values only. Nonpositive age, PT, TT, A/G, HDL-C and LDH displayed as missing; missing cells grey, no imputation.',
    'Observed basophil zeros retained for annotation; this display does not alter any saved clinical-model preprocessing or predictions.',
    'Training-set post-selection descriptive associations, not independent confirmation or causal effects. Confidence intervals do not account for prior feature selection.',
    'T1 source column used as recorded (codes 1-4). The clinical table also contains pathological staging; the clinical-versus-pathological acquisition time has not been adjudicated here. Annotation is not proof of preoperative availability.',
    'Source differentiation codes 1/2/3 retained; no unverified relabelling as well/moderate/poor.',
    'Layout adapted from Huang et al., npj Digital Medicine 2026, Fig. 4A, DOI 10.1038/s41746-026-03415-6. All plotted numbers are from this gastric-cancer training cohort.'))
wb<-createWorkbook()
sheets<-list(Readme=limits,Feature_key=feature_key,OR_per_SD=f,
  Annotations_sorted=ann,Features_raw=data.frame(Sample=d$annotations$Sample,d$x),
  Feature_z_scores=data.frame(Sample=d$annotations$Sample,d$z),
  Clinical_z_scores=data.frame(Sample=d$annotations$Sample,ann_z),Checks=d$checks)
for(nm in names(sheets)) {
  addWorksheet(wb,nm);writeData(wb,nm,sheets[[nm]],withFilter=TRUE)
  freezePane(wb,nm,firstRow=TRUE);setColWidths(wb,nm,cols=seq_len(ncol(sheets[[nm]])),widths=18)
}
saveWorkbook(wb,file.path(out,'Source_Data','Training574_Figure_Source.xlsx'),overwrite=TRUE)
write.csv(feature_key,file.path(out,'Source_Data','Feature_key.csv'),row.names=FALSE)
caption<-c(
  'Figure. Clinical annotations and selected dual-ROI radiomic features in the training cohort.',
  'Each column is one of 574 training patients (44 with anastomotic leakage and 530 without), ordered by ascending archived fusion meta-OOF probability. The same order is used for every row. Upper tracks display the three frozen OOF predictions, observed clinical variables and leakage outcome. Lower tracks display the ten previously selected features (S1-S5, stomach; T1-T5, tumor). Continuous clinical and radiomic values are standardized separately using the training observed-value mean and sample SD; colors saturate at z=-2 and z=+2 for display only. Missing clinical annotations are grey. Probability tracks use their original 0-1 scale.',
  'The aligned forest plot shows univariable odds ratios per one training SD increase in each raw radiomic feature, with Wald 95% confidence intervals on a logarithmic axis. Teal indicates OR<1 and coral indicates OR>1, not causality. T4 (exponential NGTDM Coarseness) contains an extreme value of 1000000; its per-SD OR and confidence limits numerically overflow. The log-OR and SE are retained in the source table, and its forest row is marked OR overflow. All ten features and all patients remain in the heatmap. No model refitting, re-selection, outcome changes or exclusions were made.',
  'These are descriptive associations of features previously selected in the training data. Confidence intervals and P-values do not account for feature selection; the figure is not independent validation. T-stage and differentiation annotations use existing source codes and do not establish their preoperative availability. Code meanings and exact feature definitions are provided in the source workbook. SHAP is intentionally omitted.',
  'Abbreviations: OOF, out of fold; OR, odds ratio; CI, confidence interval; SD, standard deviation; ROI, region of interest; A/G, albumin/globulin; HDL-C, high-density lipoprotein cholesterol; LDH, lactate dehydrogenase.',
  'Reference layout: https://doi.org/10.1038/s41746-026-03415-6 (Fig. 4A).')
writeLines(caption,file.path(out,'Figure_legend_EN.txt'))
writeLines(c('# Training574 heatmap and OR forest',
  '', 'Final AI: Training574_Heatmap_Forest_verified.ai. Companion files: Training574_Heatmap_Forest.svg / .pdf. All contain editable vector cells, paths and live text. Preview and 600-dpi PNGs are for viewing. The earlier AI was preserved because Illustrator reported unsaved changes.',
  '', 'Source_Data includes pseudonymized patient-level plotting values and the feature key. Reproducibility/PRIVATE_patient_crosswalk.csv contains direct research identifiers: do not upload or share publicly.',
  '', 'Only the training cohort was used. No SHAP panel and no new patient exclusions. All scores match their archived source.',
  '', 'Important: T4 has an extreme value; its raw per-SD logistic OR is unstable. The row is retained and flagged, not corrected or removed.',
  '', 'T-stage and differentiation tracks retain source coding. Their clinical timing and code semantics need source-data confirmation before a preoperative interpretation.',
  '', 'See Figure_legend_EN.txt and Source_Data/Training574_Figure_Source.xlsx for methods and limits.'),file.path(out,'README.md'))
file.copy('scripts/20261021_training_heatmap_prepare.R',file.path(out,'Reproducibility'),overwrite=TRUE)
file.copy('scripts/20261021_training_heatmap_plot.R',file.path(out,'Reproducibility'),overwrite=TRUE)
cat('Exported:',stem,'\n')
