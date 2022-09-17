##############################################################################-#
# This code can replicate the result in manuscript titled 
# "Mapping the long-term associations between air pollutants and COVID-19 
# risk and the attributable burdens in the continental United States"
##############################################################################-#
###############################################################################
############# variable explanation in adat #####################################
# doctor: active medical doctor per 1,000 people in 2019
# hospital_bed: the number of hospital beds in 2019
# medianhousevalue: median house value in 2019
# medhouseholdincome: median house income in 2019
# popdensity: population density in 2020
# mean_bmi: mean body mass index ; tem: mean temperature across 2020-2021
# education: the percent of less than high school education in 2020
# older_pecent: the percent of ≥ 65 years old in 2020
# smoke_rate: smoke rate ; pct_blk: the percent of black race in 2020
# poverty: the percent in poverty in 2020
###############################################################################

library(tmap)
library(splines)
library(RColorBrewer)
library(INLA)
library(ggplot2)
library(corrplot)
library(stringr)
library(car)
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
path <- "E:\\BaiduNetdiskWorkspace\\research\\covid\\airpollution-covid\\US-longterm\\submission\\Data and Code"
setwd(path)
source("PlotForesterFunction.R")
load("data/usmap_counties.Rdata"); mapdata_county <- mapdata
load("data/usmap_49states.Rdata"); mapdata_state <- mapdata
load("data/Analysis data in main analysis.Rdata")
apnames_all <- c("PM2.5","O3","SO2","NO2","CO")
IQR <- apply(adat[,apnames_all], 2, function(x) {a <- quantile(x);a[4]-a[2]})
a1 <- mean(adat$NO2<10); a2 <- mean(adat$PM2.5<5)
standard_level <- c(5,60,ceiling(quantile(adat$SO2,probs =a1)),10,ceiling(quantile(adat$CO,probs = a2)))
names(standard_level) <- apnames_all
rm(mapdata,a1,a2)
dir.create("results/single-pollutant",recursive = T);dir.create("results/two-pollutants",recursive = T)
###############################################################################
###########descriptive analysis ###############################################
###############################################################################
### correaltion between pollutants and incidence and morbility
adat$incidence <- adat$cases/adat$pop * 100
adat$morbility <- adat$deaths/adat$pop * 100
method <- "kendall"
corres <- matrix(NA,nrow = 5,ncol = 4,dimnames = list(apnames_all,c("Incidence","P","Morbility","P")))
a <- c("incidence","morbility")
for (i in 1:2) {
  for (j in 1:length(apnames_all)) {
    b <- cor.test(adat[,apnames_all[j]],adat[,a[i]],method = method,continuity = T)
    corres[j,c(i*2-1,i*2)] <- c(b$estimate,b$p.value)
  }
}
corres <- round(corres,3)
cordata <- adat[,apnames_all]
M <- (cor(cordata)-0.005) %>% round(2)
tiff("Corr-pollutant-matrix.tiff",width = 15,height = 15,units = "cm",res = 600)
corrplot(M, order = 'AOE', type = 'upper', tl.pos = 'd',tl.cex = 1.5,cl.cex = 1)
corrplot(M, add = TRUE, type = 'lower', method = 'number', order = 'AOE',
         diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = "black",number.cex = 1.5)
dev.off()

### sptial distribution of air pollutants
tmap_mode("plot")
for (i in apnames_all) {
  tempdata <- adat[,c("countyID",i)]
  names(tempdata)[1] <- "GEOID"
  plotdata <- merge(mapdata_county,tempdata)
  x <- round(range(plotdata[[i]],na.rm = T) + c(-0.05,0.05),1)
  v <- c(x[1],quantile(plotdata[[i]],probs = c(0.02,0.1,0.3,0.5,0.7,0.9,0.98)),x[2]) %>% round(1)
  colorlables <- NULL
  for (j in 2:length(v)) {
    colorlables <- c(colorlables,"["%+%v[j-1]%+%"-"%+%v[j]%+%")")
  }
  paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
  mapi <- tm_shape(plotdata) + 
    tm_polygons(col = i,palette=paleta, title=i,legend.show=T, border.alpha=0.1,
                legend.reverse=T, style="fixed",breaks=v, interval.closure="left",
                labels=colorlables) +
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)
  cat(i," ")
  tmap_save(mapi, "airmap-"%+%i%+%".tiff", width = 15, height = 8, units = "cm",
            dpi = 600, asp = 0, outer.margins = 0,scale = 0.9) # scale: control the magnitude of legend
}
### map for COVID-19 outcomes
adat$incidence <- adat$cases/adat$pop * 100
adat$morbility <- adat$deaths/adat$pop * 100
for (i in c("cases","incidence","deaths","morbility")) {
  tempdata <- adat[,c("countyID",i)]
  names(tempdata)[1] <- "GEOID"
  plotdata <- merge(mapdata_county,tempdata)
  if(i== "incidence"){
    x <- round(range(plotdata[[i]][adat$countyID != "48301"],na.rm = T) + c(-0.005,0.005),2) # 48301的发病数远大于人口
    v <- c(x[1],quantile(plotdata[[i]],probs = c(0.02,0.1,0.3,0.5,0.7,0.9,0.98)),Inf) %>% round(2)
    title <- "Incidence (%)"
    colorlables <- NULL
    v1 <- c(x[1],quantile(plotdata[[i]],probs = c(0.02,0.1,0.3,0.5,0.7,0.9,0.98)),x[2]) %>% round(2)
    for (j in 2:length(v1)) {
      colorlables <- c(colorlables,"["%+%v1[j-1]%+%"-"%+%v1[j]%+%")")
    }
  } else if(i== "morbility"){
    x <- round(range(plotdata[[i]],na.rm = T) + c(-0.005,0.005),2)
    v <- c(x[1],quantile(plotdata[[i]],probs = c(0.02,0.1,0.3,0.5,0.7,0.9,0.98)),x[2]) %>% round(2)
    title <- "Morbility (%)"
    colorlables <- NULL
    for (j in 2:length(v)) {
      colorlables <- c(colorlables,"["%+%v[j-1]%+%"-"%+%v[j]%+%")")
    }
  } else{
    x <- round(range(plotdata[[i]],na.rm = T) + c(-0.5,0.5))
    v <- c(x[1],quantile(plotdata[[i]],probs = c(0.02,0.1,0.3,0.5,0.7,0.9,0.98)),x[2]) %>% round()
    title <- "Cum. " %+% i
    colorlables <- NULL
    for (j in 2:length(v)) {
      colorlables <- c(colorlables,"["%+%v[j-1]%+%"-"%+%v[j]%+%")")
    }
  }
  
  paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
  mapi <- tm_shape(plotdata) + 
    tm_polygons(col = i,palette=paleta, title=title,legend.show=T, border.alpha=0.1,
                legend.reverse=T, style="fixed",breaks=v, interval.closure="left",
                labels=colorlables) +
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)
  cat(i," ")
  tmap_save(mapi, "airmap-"%+%i%+%".tiff", width = 15, height = 8, units = "cm",
            dpi = 600, asp = 0, outer.margins = 0,scale = 0.9) 
}

###############################################################################
########### model constructing: incidence######################################
###############################################################################
########### setting prior  ==================================================
ncounty <- nrow(adat)
load("data\\nb_county.Rdata")
W <- spdep::nb2mat(nb)
Wnb <- W * apply(W, 1, function(x) sum(x!=0))
R <- diag(rowSums(Wnb)) - Wnb
Cmatrix_county <- diag(ncounty) - R

nstate <- length(unique(adat$state))
load("data\\nb_state.Rdata")
W <- spdep::nb2mat(nb)
Wnb <- W * apply(W, 1, function(x) sum(x!=0))
R <- diag(rowSums(Wnb)) - Wnb
Cmatrix_state <- diag(nstate) - R

sdunif="expression:
      logdens=-log_precision/2;
      return(logdens)"
lunif = "expression:
      a = 1;
      b = 1;
      beta = exp(theta)/(1+exp(theta));
      logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
      log_jacobian = log(beta*(1-beta));
      return(logdens+log_jacobian)"
famtype <- "nbinomial" 

########### incidence single air pollutant==========================================
## LCAR for average ----
statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
fixedbeta <- NULL
for (apname in apnames_all) {
  formula <- substitute(cases ~ x + offset(log(pop)) + 
                          log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
                          education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
                          f(countylevel,model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
                            hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))),
                        list(x = as.name(apname)))
  model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.inla=list(strategy="simplified.laplace"))
  fixedbeta <- rbind(fixedbeta,model$summary.fixed[apname,c("0.025quant","0.5quant","0.975quant")] * IQR[apname]) 
  cat(apname," ")
}
fixedbeta_single <- fixedbeta
# save(fixedbeta_single,file = "Average effect-incidence-single.Rdata")
# load("Average effect-incidence-_single.Rdata")

## LCAR map + 95%CI mark ------
apnames_all <- c("PM2.5","O3", "SO2","NO2","CO")
res_attr <- NULL # 储存总归因人数
for (apname in apnames_all) {
  statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
  countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
  formula <- cases ~ offset(log(pop)) +
    log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
    + education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk + 
    f(statelevel,eval(parse(text = apname)),model="generic1", Cmatrix = Cmatrix_state, constr=FALSE,
      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
    f(countylevel, model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
  model <- inla(formula, family="nbinomial", data=adatmiss, # nbinomial
                control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.inla=list(strategy="simplified.laplace"))
  randombeta <- model$summary.random$statelevel
  beta_state <- randombeta[,c("0.025quant","0.5quant","0.975quant")]
  ### plot map RR
  plotdata <- mapdata_state 
  RR_state <- exp(beta_state * IQR[apname])
  plotdata$RR <- RR_state$`0.5quant`
  sig <- RR_state$`0.025quant` > 1 | RR_state$`0.975quant` < 1
  plotdata$mark <- NA
  plotdata$mark[sig] <- "*" ; plotdata$mark[!sig] <- "" 
  
  a <- RR_state$`0.5quant`
  x <- round(range(a,na.rm = T) + c(-0.005,0.005),2)
  if(!any(a<1)){
    v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
  } else if(!any(a>1)){
    v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
  } else {
    v <- c(x[1],quantile(a[a<1],probs = c(0.05,0.3,0.5,0.7)),1,
           quantile(a[a>1],probs = c(0.3,0.5,0.7,0.95)),x[2]) %>% round(3)
  }
  colorlables <- NULL
  for (i in 2:length(v)) {
    colorlables <- c(colorlables,"["%+%v[i-1]%+%"-"%+%v[i]%+%")")
  }
  
  paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
  mapi <- tm_shape(plotdata) + 
    tm_polygons(col = "RR",palette=paleta, title=apname, 
                legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
                breaks=v, interval.closure="left",labels=colorlables) + 
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)+
    tm_text(text = "mark",col = "black")
  tmap_save(mapi, "results/single-pollutant/RRmap-state-incidence-"%+%apname%+%".tiff", 
            width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
  cat(apname," ")
  
  ### mapping the attributable cases ---------
  standard <- standard_level[apname]
  a <- beta_state$`0.5quant`; a[a<0] <- 0
  plotdata_attr <- data.frame(state = mapdata_state$STUSPS,beta = a)
  plotdata_attr <- merge(plotdata_attr,adat[,c("state","countyID",apname,"cases")],by = "state")
  plotdata_attr$dif <- plotdata_attr[,apname] - standard
  plotdata_attr$dif[plotdata_attr$dif < 0] <- 0
  plotdata_attr$cases_attr_percent <- 1-1/exp(plotdata_attr$beta * plotdata_attr$dif)
  plotdata_attr$cases_attr <- plotdata_attr$cases_attr_percent * plotdata_attr$cases
  plotdata_attr$GEOID <- plotdata_attr$countyID
  plotdata_attr <- merge(mapdata_county,plotdata_attr)
  
  a <- sum(plotdata_attr$cases_attr); b <- a/sum(adat$cases) * 100
  res_attr <- rbind(res_attr,c(a,b))
  # the number of attributable cases
  a <- plotdata_attr$cases_attr
  x <- round(range(a,na.rm = T) + c(-0.05,0.05),1)
  v <- c(-0.1,0,1,5,quantile(a[a>5],probs = c(0.2,0.5,0.7,0.9,0.95,0.99)),x[2]) %>% ceiling()
  colorlables <- "0"
  for (i in 3:length(v)) {
    colorlables <- c(colorlables,"("%+%v[i-1]%+%"-"%+%v[i]%+%"]")
  }
  a <- length(colorlables)-1
  paleta <- c("white",brewer.pal(a,"Reds"))
  # scales::show_col(paleta)
  mapi <- tm_shape(plotdata_attr) + 
    tm_polygons(col = "cases_attr",palette=paleta, title=apname, 
                legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
                breaks=v, interval.closure="right",labels=colorlables,border.col = "grey10") + 
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)
  tmap_save(mapi, "results/single-pollutant/RRmap-state-incidence-"%+%apname%+%"_cases_attr.tiff", 
            width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
}
rownames(res_attr) <- apnames_all
write.csv(res_attr,file = "Attributed cases-single.csv")
########### incidence two air pollutants============================================
## LCAR for average ----
statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
fixedbeta <- NULL
a <- length(apnames_all)
for (i in 1:(a-1)) {
  for (j in (i+1):a) {
    apname1 <- apnames_all[i];apname2 <- apnames_all[j]
    formula <- substitute(cases ~ x1 + x2 + offset(log(pop)) + 
                            log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
                            education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
                            f(countylevel,model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
                              hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))),
                          list(x1 = as.name(apname1),x2 = as.name(apname2)))
    model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"))
    b <- c(apname1,apname2)
    d <- model$summary.fixed[b,c("0.025quant","0.5quant","0.975quant")] * IQR[b]
    rownames(d) <- c(apname1%+%"|"%+%apname2,apname2%+%"|"%+%apname1)
    fixedbeta <- rbind(fixedbeta,d)
    cat(apname1," ",apname2,"\n")
  }
}
fixedbeta_two <- fixedbeta
# save(fixedbeta_two,file = "Average effect-incidence-two.Rdata")
# load("Average effect-incidence-two.Rdata")

## LCAR map ----
apnames_all <- c("PM2.5","O3", "SO2","NO2","CO")
for (apname in apnames_all) {
  for (conap in setdiff(apnames_all,apname)) {
    statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
    countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
    controlx <- adat[,conap]
    formula <- cases ~ offset(log(pop)) + controlx +
      log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
      + education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
      f(statelevel,eval(parse(text = apname)),model="generic1", Cmatrix = Cmatrix_state, constr=FALSE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
      f(countylevel, model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
    formula <- gsub(apname%+%" +","",paste(as.character(formula)[c(2,1,3)],collapse = " ")) %>% as.formula()
    model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"))
    randombeta <- model$summary.random$statelevel
    beta_state <- randombeta[,c("0.025quant","0.5quant","0.975quant")]
    apply(randombeta, 2, mean) 
    
    ### plot map
    plotdata <- mapdata_state 
    RR_state <- exp(beta_state * IQR[apname])
    plotdata$RR <- RR_state$`0.5quant`
    sig <- RR_state$`0.025quant` > 1 | RR_state$`0.975quant` < 1
    plotdata$mark <- NA
    plotdata$mark[sig] <- "*" ; plotdata$mark[!sig] <- "" 
    
    a <- RR_state$`0.5quant`
    x <- round(range(a,na.rm = T) + c(-0.005,0.005),2)
    if(!any(a<1)){
      v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
    } else if(!any(a>1)){
      v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
    } else {
      v <- c(x[1],quantile(a[a<1],probs = c(0.05,0.3,0.5,0.7)),1,
             quantile(a[a>1],probs = c(0.3,0.5,0.7,0.95)),x[2]) %>% round(3)
    }
    colorlables <- NULL
    for (i in 2:length(v)) {
      colorlables <- c(colorlables,"["%+%v[i-1]%+%"-"%+%v[i]%+%")")
    }
    
    paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
    mapi <- tm_shape(plotdata) + 
      tm_polygons(col = "RR",palette=paleta, title=apname%+%"|"%+%conap, 
                  legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
                  breaks=v, interval.closure="left",labels=colorlables) + 
      tm_layout(legend.position = c("right","bottom"),frame = FALSE)+
      tm_text(text = "mark",col = "black")
    tmap_save(mapi, "results/two-pollutants/RRmap-state-incidence-"%+%apname%+%"-"%+%conap%+%".tiff", 
              width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
    cat(apname,"|",conap," ")
  }
}
########### incidence plot forest for average associations: single + two ===========
# load("Average effect-incidence-two.Rdata")
# load("Average effect-incidence-single.Rdata")
RR_average <- exp(rbind(fixedbeta_single,fixedbeta_two))
names(RR_average) <- c("beta0.025","beta","beta0.975")
plotRR <- as.data.frame(RR_average)
plotRR$class <- rownames(plotRR)
plotRR <- plotRR[sort(rownames(plotRR)),]

# plot forester
plotRR$Pollutant <- plotRR$class
plotRR$class <- stringr::str_split(plotRR$class,pattern = "\\|",simplify = T)[,1]
mydata <- rbind(plotRR[plotRR$class == apnames_all[1],],NA,plotRR[plotRR$class == apnames_all[2],],
                NA,plotRR[plotRR$class == apnames_all[3],],NA,plotRR[plotRR$class == apnames_all[4],],
                NA,plotRR[plotRR$class == apnames_all[5],])
xlims <- range(mydata[,1:3],na.rm = T); xlims <- c(0.90,1.10) 
mydata$col <- "black"; mydata$col[mydata$Pollutant %in% apnames_all] <- "red"
mydata$IQR <- NA; mydata$IQR[!is.na(mydata$Pollutant)] <- rep(round(IQR,1),each = length(apnames_all))
forester_update(left_side_data = mydata[,c("Pollutant","IQR")],   
                estimate = round(mydata$beta,3),    
                ci_low = mydata$beta0.025,      
                ci_high = mydata$beta0.975,   
                xlim = xlims,             
                estimate_precision = 3,
                null_line_at = 1,
                file_path = "forester_plot_incidence.tiff",
                point_sizes = 1,
                xbreaks = c(0.90,0.95,1,1.05,1.10),
                ggplot_width = 20,
                arrowss = F,               
                estimate_col_name = "RR per IQR increase",colors = mydata$col)

###############################################################################
########### model constructing: mortality ######################################
###############################################################################
########### mobility single air pollutants=========================================
## LCAR for average ----
statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
fixedbeta <- NULL
for (apname in apnames_all) {
  formula <- substitute(deaths ~ x + offset(log(pop)) + 
                   log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
                   education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
                   f(countylevel,model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
                     hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))),
                   list(x = as.name(apname)))
  model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.inla=list(strategy="simplified.laplace"))
  fixedbeta <- rbind(fixedbeta,model$summary.fixed[apname,c("0.025quant","0.5quant","0.975quant")] * IQR[apname]) 
  cat(apname," ")
}
fixedbeta_single <- fixedbeta
# save(fixedbeta_single,file = "Average effect-morbility-single.Rdata")
# load("Average effect-morbility-single.Rdata")

## LCAR map + 95%CI mark ------
apnames_all <- c("PM2.5","O3", "SO2","NO2","CO")
res_attr <- NULL 
for (apname in apnames_all) {
  statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
  statelevel1 <- statelevel
  countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
  formula <- deaths ~ offset(log(pop)) +
    log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
    + education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
    f(statelevel,eval(parse(text = apname)),model="generic1", Cmatrix = Cmatrix_state, constr=FALSE,
      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
    f(countylevel, model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
  formula <- gsub(apname%+%" +","",paste(as.character(formula)[c(2,1,3)],collapse = " ")) %>% as.formula()
  model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.inla=list(strategy="simplified.laplace"))
  randombeta <- model$summary.random$statelevel
  beta_state <- randombeta[,c("0.025quant","0.5quant","0.975quant")]
  ### plot map
  plotdata <- mapdata_state 
  RR_state <- exp(beta_state * IQR[apname])
  plotdata$RR <- RR_state$`0.5quant`
  sig <- RR_state$`0.025quant` > 1 | RR_state$`0.975quant` < 1
  plotdata$mark <- NA
  plotdata$mark[sig] <- "*" ; plotdata$mark[!sig] <- "" 
  
  a <- RR_state$`0.5quant`
  x <- round(range(a,na.rm = T) + c(-0.005,0.005),2)
  if(!any(a<1)){
    v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
  } else if(!any(a>1)){
    v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
  } else {
    v <- c(x[1],quantile(a[a<1],probs = c(0.05,0.3,0.5,0.7)),1,
           quantile(a[a>1],probs = c(0.3,0.5,0.7,0.95)),x[2]) %>% round(3)
  }
  colorlables <- NULL
  for (i in 2:length(v)) {
    colorlables <- c(colorlables,"["%+%v[i-1]%+%"-"%+%v[i]%+%")")
  }
  
  paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
  mapi <- tm_shape(plotdata) + 
    tm_polygons(col = "RR",palette=paleta, title=apname, 
                legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
                breaks=v, interval.closure="left",labels=colorlables) + 
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)+
    tm_text(text = "mark",col = "black")
  tmap_save(mapi, "results/single-pollutant/RRmap-state-morbility-"%+%apname%+%".tiff", 
            width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
  cat(apname," ")
  
  ### mapping the attributable deaths ---------
  standard <- standard_level[apname]
  a <- beta_state$`0.5quant`; a[a<0] <- 0
  plotdata_attr <- data.frame(state = mapdata_state$STUSPS,beta = a)
  plotdata_attr <- merge(plotdata_attr,adat[,c("state","countyID",apname,"deaths")],by = "state")
  plotdata_attr$dif <- plotdata_attr[,apname] - standard
  plotdata_attr$dif[plotdata_attr$dif < 0] <- 0
  plotdata_attr$cases_attr_percent <- 1-1/exp(plotdata_attr$beta * plotdata_attr$dif)
  plotdata_attr$cases_attr <- plotdata_attr$cases_attr_percent * plotdata_attr$deaths
  plotdata_attr$GEOID <- plotdata_attr$countyID
  plotdata_attr <- merge(mapdata_county,plotdata_attr)
  
  a <- sum(plotdata_attr$cases_attr); b <- a/sum(adat$deaths) * 100
  res_attr <- rbind(res_attr,c(a,b))
  # attributabble deaths
  a <- plotdata_attr$cases_attr
  x <- round(range(a,na.rm = T) + c(-0.05,0.05),1)
  v <- c(-0.1,0,1,5,quantile(a[a>5],probs = c(0.2,0.5,0.7,0.9,0.95,0.99)),x[2]) %>% ceiling()
  colorlables <- "0"
  for (i in 3:length(v)) {
    colorlables <- c(colorlables,"("%+%v[i-1]%+%"-"%+%v[i]%+%"]")
  }
  a <- length(colorlables)-1
  paleta <- c("white",brewer.pal(a,"Reds"))
  # scales::show_col(paleta)
  mapi <- tm_shape(plotdata_attr) +
    tm_polygons(col = "cases_attr",palette=paleta, title=apname,
                legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed",
                breaks=v, interval.closure="right",labels=colorlables,border.col = "grey10") +
    tm_layout(legend.position = c("right","bottom"),frame = FALSE)
  tmap_save(mapi, "results/single-pollutant/RRmap-state-mortality-"%+%apname%+%"_deaths_attr.tiff",
            width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
}
rownames(res_attr) <- apnames_all
write.csv(res_attr,file = "Attributed deaths-single.csv")
########### mobility two air pollutants=========================================
## LCAR for average ---------
statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
fixedbeta <- NULL
a <- length(apnames_all)
for (i in 1:(a-1)) {
  for (j in (i+1):a) {
    apname1 <- apnames_all[i];apname2 <- apnames_all[j]
    formula <- substitute(deaths ~ x1 + x2 + offset(log(pop)) + 
                          log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
                          education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
                          f(countylevel,model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
                            hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))),
                          list(x1 = as.name(apname1),x2 = as.name(apname2)))
    model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"))
    b <- c(apname1,apname2)
    d <- model$summary.fixed[b,c("0.025quant","0.5quant","0.975quant")] * IQR[b]
    rownames(d) <- c(apname1%+%"|"%+%apname2,apname2%+%"|"%+%apname1)
    fixedbeta <- rbind(fixedbeta,d)
    cat(apname1," ",apname2,"\n")
  }
}
fixedbeta_two <- fixedbeta
# save(fixedbeta_two,file = "Average effect-morbility-two.Rdata")
# load("Average effect-morbility-two.Rdata")

## LCAR map ----
apnames_all <- c("PM2.5","O3", "SO2","NO2","CO")
for (apname in apnames_all) {
  for (conap in setdiff(apnames_all,apname)) {
    statelevel <- factor(adat$state,levels = mapdata_state$STUSPS) %>% as.integer()
    countylevel <- factor(adat$countyID,levels = mapdata_county$GEOID) %>% as.integer()
    controlx <- adat[,conap]
    formula <- deaths ~ offset(log(pop)) + controlx +
      log(medhouseholdincome) + hospital_bed + ppp + popdensity + tem + mean_bmi +
      + education + older_pecent + smoke_rate + log(medianhousevalue) + poverty + pct_blk +
      f(statelevel,eval(parse(text = apname)),model="generic1", Cmatrix = Cmatrix_state, constr=FALSE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
      f(countylevel, model="generic1", Cmatrix = Cmatrix_county, constr=TRUE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
    formula <- gsub(apname%+%" +","",paste(as.character(formula)[c(2,1,3)],collapse = " ")) %>% as.formula()
    model <- inla(formula, family=famtype, data=adatmiss, # nbinomial
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"))
    randombeta <- model$summary.random$statelevel
    beta_state <- randombeta[,c("0.025quant","0.5quant","0.975quant")]
    apply(randombeta, 2, mean) 
    
    ### plot map
    plotdata <- mapdata_state 
    RR_state <- exp(beta_state * IQR[apname])
    plotdata$RR <- RR_state$`0.5quant`
    sig <- RR_state$`0.025quant` > 1 | RR_state$`0.975quant` < 1
    plotdata$mark <- NA
    plotdata$mark[sig] <- "*" ; plotdata$mark[!sig] <- "" 
    
    a <- RR_state$`0.5quant`
    x <- round(range(a,na.rm = T) + c(-0.005,0.005),2)
    if(!any(a<1)){
      v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
    } else if(!any(a>1)){
      v <- c(x[1],quantile(a,probs = c(0.2,0.4,0.6,0.8,0.95)),x[2]) %>% round(3)
    } else {
      v <- c(x[1],quantile(a[a<1],probs = c(0.05,0.3,0.5,0.7)),1,
             quantile(a[a>1],probs = c(0.3,0.5,0.7,0.95)),x[2]) %>% round(3)
    }
    colorlables <- NULL
    for (i in 2:length(v)) {
      colorlables <- c(colorlables,"["%+%v[i-1]%+%"-"%+%v[i]%+%")")
    }
    
    paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
    mapi <- tm_shape(plotdata) + 
      tm_polygons(col = "RR",palette=paleta, title=apname%+%"|"%+%conap, 
                  legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
                  breaks=v, interval.closure="left",labels=colorlables) + 
      tm_layout(legend.position = c("right","bottom"),frame = FALSE)+
      tm_text(text = "mark",col = "black")
    tmap_save(mapi, "results/two-pollutants/RRmap-state-morbility-"%+%apname%+%"-"%+%conap%+%".tiff", 
              width = 15, height = 8, units = "cm",dpi = 600, asp = 0, outer.margins = 0,scale = 0.9)
    cat(apname,"|",conap," ")
  }
}
########### mortality plot forest: single + two ----
# load("Average effect-morbility-two.Rdata")
# load("Average effect-morbility-single.Rdata")
RR_average <- exp(rbind(fixedbeta_single,fixedbeta_two))
names(RR_average) <- c("beta0.025","beta","beta0.975")
plotRR <- as.data.frame(RR_average)
plotRR$class <- rownames(plotRR)
plotRR <- plotRR[sort(rownames(plotRR)),]

# plot forester
plotRR$Pollutant <- plotRR$class
plotRR$class <- stringr::str_split(plotRR$class,pattern = "\\|",simplify = T)[,1]
mydata <- rbind(plotRR[plotRR$class == apnames_all[1],],NA,plotRR[plotRR$class == apnames_all[2],],
                NA,plotRR[plotRR$class == apnames_all[3],],NA,plotRR[plotRR$class == apnames_all[4],],
                NA,plotRR[plotRR$class == apnames_all[5],])
xlims <- range(mydata[,1:3],na.rm = T); xlims <- c(0.90,1.10) 
mydata$col <- "black"; mydata$col[mydata$Pollutant %in% apnames_all] <- "red"
mydata$IQR <- NA; mydata$IQR[!is.na(mydata$Pollutant)] <- rep(round(IQR,1),each = length(apnames_all))
forester_update(left_side_data = mydata[,c("Pollutant","IQR")],  
                estimate = round(mydata$beta,3),   
                ci_low = mydata$beta0.025,      
                ci_high = mydata$beta0.975,   
                xlim = xlims,            
                estimate_precision = 3,
                null_line_at = 1,
                file_path = "forester_plot_mortality.tiff",
                point_sizes = 1,
                xbreaks = c(0.90,0.95,1,1.05,1.10),
                ggplot_width = 20,
                arrowss = F,              
                estimate_col_name = "RR per IQR increase",colors = mydata$col)


