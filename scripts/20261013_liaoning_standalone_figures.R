source("scripts/plain_paths.R")
.libPaths(c(plain_path("R_LIBRARY", ""), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
src <- file.path(args[1], '03_Source_Data')
out <- args[2]
with_ci <- length(args)>=3 && args[3]=='auc-ci'
roc <- read.csv(file.path(src, 'ROC_before_after.csv'))
cal <- read.csv(file.path(src, 'calibration_before_after.csv'))
metrics <- read.csv(file.path(src, 'model_metrics.csv'))
cols <- c(Clinical='#457B9D', Radiomics='#2A9D8F', Fusion='#826181')
roc <- subset(roc, Stage=='Clinical_restored')
cal <- subset(cal, Stage=='Clinical_restored' & Model=='Fusion')
metrics <- subset(metrics, Stage=='Clinical_restored')
for (name in c('roc', 'metrics')) {
  q <- get(name)
  q$Model <- factor(q$Model, levels=names(cols))
  assign(name, q)
}
if (with_ci) {
  pred <- read.csv(file.path(src, 'probability_restoration_PRIVATE.csv'),
    colClasses=c(patient_id='character', Chart_id='character'))
  stopifnot(nrow(pred)==355, sum(pred$Outcome)==17,
    !anyDuplicated(pred$patient_id), all(pred$Outcome %in% c(0,1)))
  columns <- c(Clinical='Clinical_restored_original',
    Radiomics='Radiomics_unchanged', Fusion='Fusion_unchanged')
  ci_table <- do.call(rbind, lapply(names(cols), function(model) {
    p <- pred[[columns[[model]]]]; y <- pred$Outcome
    stopifnot(all(is.finite(p)), all(p>=0 & p<=1))
    r <- pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)
    ci <- as.numeric(pROC::ci.auc(r,method='delong',conf.level=.95))
    # Independently verify the DeLong variance from pairwise placements.
    pairs <- outer(p[y==1],p[y==0],function(a,b)(a>b)+.5*(a==b))
    variance <- var(rowMeans(pairs))/sum(y==1)+var(colMeans(pairs))/sum(y==0)
    bounds <- pmax(0,pmin(1,mean(pairs)+c(-1,1)*qnorm(.975)*sqrt(variance)))
    stopifnot(abs(mean(pairs)-ci[2])<1e-12,
      max(abs(bounds-ci[c(1,3)]))<1e-12,
      abs(metrics$AUC[metrics$Model==model]-ci[2])<1e-12)
    data.frame(Model=model,N=length(y),Events=sum(y==1),Nonevents=sum(y==0),
      AUC=ci[2],CI95_lower=ci[1],CI95_upper=ci[3],
      Method='DeLong',DeLong_SE=sqrt(variance))
  }))
  metrics$CI95_lower <- ci_table$CI95_lower[match(metrics$Model,ci_table$Model)]
  metrics$CI95_upper <- ci_table$CI95_upper[match(metrics$Model,ci_table$Model)]
  write.csv(ci_table,file.path(out,'AUC_95CI_DeLong.csv'),row.names=FALSE)
  jsonlite::write_json(list(N=355,events=17,nonevents=338,method='DeLong',
    confidence_level=.95,pROC_version=as.character(packageVersion('pROC')),
    variance_independently_verified=TRUE,
    scope='Conditional on the selected cohort and saved predictions; not selection-adjusted or independent validation.',
    reference='https://cran.r-project.org/web/packages/pROC/refman/pROC.html#ci.auc'),
    file.path(out,'method.json'),pretty=TRUE,auto_unbox=TRUE)
  print(ci_table,digits=10,row.names=FALSE)
}
theme_original <- function() theme_bw(base_family='Arial', base_size=10)+
  theme(panel.border=element_rect(colour='#E5EAEE', linewidth=.3),
    axis.line=element_line(linewidth=.3),
    panel.grid.major=element_line(colour='#E5EAEE', linewidth=.25),
    panel.grid.minor=element_blank(),
    plot.title=element_text(face='bold', size=12),
    plot.subtitle=element_text(size=9, colour='#65717C'),
    plot.caption=element_text(size=8, colour='#65717C', hjust=0, lineheight=1.15),
    plot.margin=margin(10,12,10,10))
export <- function(p, stem) {
  w <- 190/25.4; h <- 205/25.4
  ggsave(file.path(out, 'SVG', paste0(stem,'.svg')), p,
    width=w, height=h, device=svglite::svglite, bg='white')
  quartz(type='pdf', file=file.path(out, 'PDF', paste0(stem,'.pdf')),
    width=w, height=h, family='Arial', bg='white')
  print(p); dev.off()
  ggsave(file.path(out, 'Previews', paste0(stem,'.png')), p,
    width=w, height=h, dpi=200, device=ragg::agg_png, bg='white')
}
metrics$Label <- if(with_ci) sprintf('%s  %.4f (%.4f-%.4f)',
  metrics$Model,metrics$AUC,metrics$CI95_lower,metrics$CI95_upper) else
  sprintf('%s  AUC %.4f', metrics$Model, metrics$AUC)
metrics$y <- .24-(as.integer(metrics$Model)-1)*.075
metrics$x <- if(with_ci) .34 else .43
pr <- ggplot(roc, aes(FPR,TPR,colour=Model))+
  geom_abline(slope=1, intercept=0, linetype='dashed', colour='#ABB4BC', linewidth=.4)+
  geom_path(linewidth=.75)+
  geom_text(data=metrics, aes(x=x,y=y,label=Label), hjust=0,
    size=(if(with_ci)9 else 10)/ggplot2::.pt, show.legend=FALSE, family='Arial')+
  scale_colour_manual(values=cols, guide='none')+
  scale_x_continuous(limits=c(0,1), breaks=c(0,.5,1), expand=c(.005,0))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,.5,1), expand=c(.005,0))+
  coord_fixed()+theme_original()+
  labs(title=if(with_ci)'ROC curves | AUC with 95% CI' else 'ROC curves',
    subtitle='Same 355 patients: 17 events / 338 non-events.',
    x='1 - Specificity', y='Sensitivity',
    caption=paste('Clinical: archived pre-update probabilities restored.',
      'Radiomics and Fusion: previous updated scores, unchanged.',
      if(with_ci)'DeLong 95% CIs; no adjustment for prior outcome-based selection.' else
        'No new exclusions or model refitting.',
      'Post-hoc outcome-selected cohort; NOT independent validation.', sep='\n'))
export(pr, if(with_ci)'01_ROC_with_95CI' else '01_ROC')
pc <- ggplot(cal, aes(Mean_prediction,Observed_fraction))+
  geom_abline(slope=1, intercept=0, linetype='dashed', colour='#ABB4BC', linewidth=.4)+
  geom_errorbar(aes(ymin=Wilson95_lower,ymax=Wilson95_upper),
    width=.014, colour=cols['Fusion'], linewidth=.4)+
  geom_line(colour=cols['Fusion'], linewidth=.7)+
  geom_point(shape=21, fill='white', colour=cols['Fusion'], size=1.9, stroke=.5)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,.5,1), expand=c(.005,0))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,.5,1), expand=c(.005,0))+
  coord_fixed()+theme_original()+
  labs(title='Fusion calibration',
    subtitle='Same 355 patients: 17 events / 338 non-events.',
    x='Mean predicted probability', y='Observed event rate',
    caption=paste('Saved Fusion probabilities and calibration bins unchanged.',
      'Bars: Wilson 95% intervals; selection uncertainty not included.',
      'No new exclusions or model refitting.',
      'Post-hoc outcome-selected cohort; NOT independent validation.', sep='\n'))
if(!with_ci) export(pc, '02_Fusion_Calibration')
