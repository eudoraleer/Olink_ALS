#########################################################################################################
# Title: Protein biomarkers in risk and prognosis of amyotrophic lateral sclerosis
# Data: Olink Proteomics
# Creator: Lu Pan, lu.pan@ki.se
# Date last modified: 2023-12-13
#########################################################################################################

#########################################################################################################
library(gee)
library(xlsx)
library(eeptools)
library(dplyr)
library(KMsurv)

data <- NULL
pheno <- NULL
data <- readRDS("Data/Olink_ALS_Data.RDS")
pheno <- readRDS("Data/Olink_ALS_Metadata.RDS")

project_name <- "OLINK_ALS"
output_dir <- paste("Output/", sep = "")
dir.create(output_dir)
setwd(output_dir)

fig1_data <- NULL
for(i in 1:length(data)){
  print(i)
  cproteins <- NULL
  cproteins <- sort(unique(data[[i]]$OlinkID))
  for(j in 1:length(cproteins)){
    current <- NULL
    current <- data[[i]][which(data[[i]]$OlinkID == cproteins[j]),]
    current$GROUP <- factor(current$GROUP, levels = c("Control","ALS"))
    cresult <- NULL
    cresult <- gee(NPX_RAW ~ GROUP + BP_HIGH + SMOKING + CLOSEST_BMI + MND_DIAGNOSIS_AGE + SEX + NUM_FREEZING_DAYS + PlateID, id = FAMILY_NUMBER, data = current)
    crp <- NULL
    crp <- 2 * pnorm(-abs(summary(cresult)$coefficients[,"Robust z"]))
    cnp <- NULL
    cnp <- 2 * pnorm(-abs(summary(cresult)$coefficients[,"Naive z"]))
    ci <- NULL
    se <- NULL
    se <- summary(cresult)$coefficients["GROUPALS","Robust S.E."]
    ci <- exp(coef(cresult)["GROUPALS"] + c(-1, 1) * se * qnorm(0.975))
    cest <- NULL
    cest <- exp(coef(cresult)["GROUPALS"])
    cfc <- NULL
    cfc <- mean(current[which(current$GROUP == "ALS"),"NPX_RAW"]) - mean(current[which(current$GROUP == "Control"),"NPX_RAW"])
    fig1_data <- rbind(fig1_data, data.frame(DATA_TYPE = unique(current$DATA_TYPE),
                                             Assay = unique(current$Assay),
                                             Panel = unique(current$Panel),
                                             OlinkID = unique(current$OlinkID),
                                             COMPARISON = paste(levels(current$GROUP), collapse = "_VS_"),
                                             ROBUST_P = crp["GROUPALS"],
                                             NAIVE_P = cnp["GROUPALS"],
                                             CONF_LOW = ci[1],
                                             CONF_HIGH = ci[2],
                                             EXP_ESTIMATE = cest,
                                             AVE_LOG2FC = cfc))
  }
}

fig1_data <- split(fig1_data, fig1_data$DATA_TYPE)
fig1_data <- lapply(fig1_data, function(x){
  x <- data.frame(x, ROBUSTPADJUST = p.adjust(x$ROBUST_P, method = "BH"))
})
fig1_data <- do.call(rbind.data.frame, fig1_data)
fig1_data <- fig1_data[order(fig1_data$ROBUST_P, decreasing = F),]
fig1_data_all <- fig1_data

fig1_data <- fig1_data[which(fig1_data$ROBUSTPADJUST < 0.05),]
colnames(fig1_data) <- gsub(" ","_", stringr::str_to_title(gsub("_"," ", colnames(fig1_data))))
colnames(fig1_data) <- gsub("Exp_estimate","Odds_Ratio",colnames(fig1_data), ignore.case = T)
colnames(fig1_data) <- gsub("Ave_log2fc","Mean_Difference_NPX",colnames(fig1_data), ignore.case = T)
colnames(fig1_data) <- gsub("Olinkid","OlinkID",colnames(fig1_data))
colnames(fig1_data) <- gsub("Robustpadjust","Robust_Padjust",colnames(fig1_data))
colnames(fig1_data)
fig1_data <- fig1_data[order(fig1_data$Mean_Difference_NPX, decreasing = T),]
fig1_data$Comparison <- "ALS_VS_Control"

xlsx::write.xlsx(fig1_data[which(fig1_data$Data_Type == "PLASMA"),],paste(output_dir,"ALS_Olink_Table_S2.xlsx", sep = ""), sheetName = "Plasma", row.names = F)
xlsx::write.xlsx(fig1_data[which(fig1_data$Data_Type == "CSF"),],paste(output_dir,"ALS_Olink_Table_S2.xlsx", sep = ""), sheetName = "CSF", row.names = F, append = T)

########################################################################################
########################################################################################
fig1_SPSI <- NULL
selected <- c("SP","SI")
for(m in 1:length(selected)){
  for(i in grep("PLASMA", names(data), ignore.case = T)){
    print(i)
    data[[i]]$SUBGROUP <- pheno[match(data[[i]]$SID, pheno$SID),"SUBGROUP"]
    data[[i]]$SUBGROUP <- gsub("Case","ALS",data[[i]]$SUBGROUP, ignore.case = T)
    
    cproteins <- NULL
    cproteins <- sort(unique(data[[i]]$OlinkID))
    for(j in 1:length(cproteins)){
      current <- NULL
      current <- data[[i]][which(data[[i]]$OlinkID == cproteins[j]),]
      current <- current[which(current$SUBGROUP %in% c("ALS",selected[m])),]
      current$SUBGROUP <- factor(current$SUBGROUP, levels = c(selected[m],"ALS"))
      cresult <- NULL
      cresult <- gee(NPX_RAW ~ SUBGROUP + BP_HIGH + SMOKING + CLOSEST_BMI + MND_DIAGNOSIS_AGE + SEX + NUM_FREEZING_DAYS + PlateID, id = FAMILY_NUMBER, data = current)
      crp <- NULL
      crp <- 2 * pnorm(-abs(summary(cresult)$coefficients[,"Robust z"]))
      cnp <- NULL
      cnp <- 2 * pnorm(-abs(summary(cresult)$coefficients[,"Naive z"]))
      ci <- NULL
      se <- NULL
      se <- summary(cresult)$coefficients["SUBGROUPALS","Robust S.E."]
      ci <- exp(coef(cresult)["SUBGROUPALS"] + c(-1, 1) * se * qnorm(0.975))
      cest <- NULL
      cest <- exp(coef(cresult)["SUBGROUPALS"])
      cfc <- NULL
      cfc <- mean(current[which(current$SUBGROUP == "ALS"),"NPX_RAW"]) - mean(current[which(current$SUBGROUP == selected[m]),"NPX_RAW"])
      fig1_SPSI <- rbind(fig1_SPSI, data.frame(DATA_TYPE = unique(current$DATA_TYPE),
                                               Assay = unique(current$Assay),
                                               OlinkID = unique(current$OlinkID),
                                               Panel = unique(current$Panel),
                                               COMPARISON = paste(levels(current$SUBGROUP), collapse = "_VS_"),
                                               ROBUST_P = crp["SUBGROUPALS"],
                                               NAIVE_P = cnp["SUBGROUPALS"],
                                               CONF_LOW = ci[1],
                                               CONF_HIGH = ci[2],
                                               EXP_ESTIMATE = cest,
                                               AVE_LOG2FC = cfc))
    }
  }
  
}

fig1_SPSI$SPLIT_TERM <- paste(fig1_SPSI$DATA_TYPE, fig1_SPSI$COMPARISON, sep = "_")
fig1_SPSI <- split(fig1_SPSI, fig1_SPSI$SPLIT_TERM)
fig1_SPSI <- lapply(fig1_SPSI, function(x){
  x <- data.frame(x, ROBUSTPADJUST = p.adjust(x$ROBUST_P, method = "BH"))
})
fig1_SPSI <- do.call(rbind.data.frame, fig1_SPSI)
fig1_SPSI <- fig1_SPSI[order(fig1_SPSI$ROBUSTPADJUST, decreasing = F),]

########################################################################################
########################################################################################
# CORRELATION TEST
plotx <- fig1_data_all
plotx <- plotx[order(plotx$ROBUSTPADJUST, decreasing = F),]
sig_markers <- plotx[which(plotx$ROBUSTPADJUST < 0.05),]
sig_overlaps <- intersect(sig_markers[which(sig_markers$DATA_TYPE == "CSF"),"OlinkID"],
                          sig_markers[which(sig_markers$DATA_TYPE == "PLASMA"),"OlinkID"])
sig_csf <- sig_markers[which(sig_markers$DATA_TYPE == "CSF"),"OlinkID"]
sig_csf <- sig_csf[which(!sig_csf %in% sig_overlaps)]
sig_plasma <- sig_markers[which(sig_markers$DATA_TYPE == "PLASMA"),"OlinkID"]
sig_plasma <- sig_plasma[which(!sig_plasma %in% sig_overlaps)]
sig_plasma_top10 <- sig_plasma[1:ifelse(length(sig_plasma) > 10, 10, length(sig_plasma))]
sig_csf_top10 <- sig_csf[1:ifelse(length(sig_csf) > 10, 10, length(sig_csf))]

subset1 <- plotx[which(plotx$DATA_TYPE == "PLASMA"),]
subset1 <- data.frame(PLASMA_ASSAY = subset1$Assay,
                      PLASMA_OlinkID = subset1$OlinkID,
                      PLASMA_OR = subset1$AVE_LOG2FC,
                      PLASMA_P = subset1$ROBUST_P,
                      PLASMA_BH = subset1$ROBUSTPADJUST,
                      PLASMA_COMPARISON = "PLASMA ALS VS CONTROL (FDR < 0.05)")
subset2 <- plotx[which(plotx$DATA_TYPE == "CSF"),]
subset2 <- data.frame(CSF_ASSAY = subset2$Assay,
                      CSF_OlinkID = subset2$OlinkID,
                      CSF_OR = subset2$AVE_LOG2FC,
                      CSF_P = subset2$ROBUST_P,
                      CSF_BH = subset2$ROBUSTPADJUST,
                      CSF_COMPARISON = "CSF ALS VS CONTROL (FDR < 0.05)")
current <- cbind(subset1, subset2[match(subset1$PLASMA_OlinkID, subset2$CSF_OlinkID),])
which(current$PLASMA_ASSAY != current$CSF_ASSAY)
current$ABS_DIFF <- as.numeric(abs(current$PLASMA_OR - current$CSF_OR))
current$ASSAY <- ifelse(current$PLASMA_OlinkID == current$PLASMA_OlinkID, current$PLASMA_ASSAY, NA)
which(is.na(current$ASSAY))

current$Label <- NULL
current$Label <- ifelse(current$PLASMA_OlinkID %in% c(sig_plasma_top10) & !(current$PLASMA_OlinkID %in% sig_overlaps), current$PLASMA_ASSAY, "")
current$GROUP <- ifelse(current$PLASMA_OlinkID %in% c(sig_plasma) & !(current$PLASMA_OlinkID %in% sig_overlaps), "PLASMA ALS VS CONTROL (FDR < 0.05)", "")
current$Label2 <- ifelse(current$CSF_OlinkID %in% c(sig_csf_top10) & !(current$CSF_OlinkID %in% sig_overlaps), current$CSF_ASSAY, "")
current$GROUP2 <- ifelse(current$CSF_OlinkID %in% c(sig_csf) & !(current$CSF_OlinkID %in% sig_overlaps), "CSF ALS VS CONTROL (FDR < 0.05)", "")
current[which(current$Label == ""),"Label"] <- current[which(current$Label == ""),"Label2"]
current[which(current$GROUP == ""),"GROUP"] <- current[which(current$GROUP == ""),"GROUP2"]
current$Label2 <- ifelse(current$CSF_OlinkID %in% c(sig_overlaps), current$CSF_ASSAY, "")
current$GROUP2 <- ifelse(current$CSF_OlinkID %in% sig_overlaps, "PLASMA & CSF  - ALS VS CONTROL (FDR < 0.05)", "")
current[which(current$Label == ""),"Label"] <- current[which(current$Label == ""),"Label2"]
current[which(current$GROUP == ""),"GROUP"] <- current[which(current$GROUP == ""),"GROUP2"]
current[which(current$GROUP == ""),"GROUP"] <- "Others"
current$GROUP <- factor(current$GROUP, levels = c("PLASMA ALS VS CONTROL (FDR < 0.05)",
                                                  "CSF ALS VS CONTROL (FDR < 0.05)",
                                                  "PLASMA & CSF  - ALS VS CONTROL (FDR < 0.05)",
                                                  "Others"))
fig3_data <- current

###################################################################
###################################################################
# SUPPLEMENTARY FIGURE 3
plotx <- fig1_SPSI
plotx$TYPE <- gsub("_VS_ALS","",plotx$COMPARISON)
plotx <- plotx[order(plotx$ROBUSTPADJUST, decreasing = F),]
sig_markers <- plotx[which(plotx$ROBUSTPADJUST < 0.05),]
sig_overlaps <- intersect(sig_markers[which(sig_markers$TYPE == "SI"),"OlinkID"],
                          sig_markers[which(sig_markers$TYPE == "SP"),"OlinkID"])
sig_g1 <- sig_markers[which(sig_markers$TYPE == "SP"),"OlinkID"]
sig_g1 <- sig_g1[which(!sig_g1 %in% sig_overlaps)]
sig_g2 <- sig_markers[which(sig_markers$TYPE == "SI"),"OlinkID"]
sig_g2 <- sig_g2[which(!sig_g2 %in% sig_overlaps)]
sig_g1_top10 <- sig_g1[1:ifelse(length(sig_g1) > 10, 10, length(sig_g1))]
sig_g2_top10 <- sig_g2[1:ifelse(length(sig_g2) > 10, 10, length(sig_g2))]

subset1 <- plotx[which(plotx$TYPE == "SP"),]
subset1 <- data.frame(SP_ASSAY = subset1$Assay,
                      SP_OlinkID = subset1$OlinkID,
                      SP_OR = subset1$AVE_LOG2FC,
                      SP_P = subset1$ROBUST_P,
                      SP_BH = subset1$ROBUSTPADJUST,
                      SP_COMPARISON = "Plasma ALS VS SP (FDR < 0.05)")
subset2 <- plotx[which(plotx$TYPE == "SI"),]
subset2 <- data.frame(SI_ASSAY = subset2$Assay,
                      SI_OlinkID = subset2$OlinkID,
                      SI_OR = subset2$AVE_LOG2FC,
                      SI_P = subset2$ROBUST_P,
                      SI_BH = subset2$ROBUSTPADJUST,
                      SI_COMPARISON = "Plasma ALS VS SI (FDR < 0.05)")
current <- cbind(subset1, subset2[match(subset1$SP_OlinkID, subset2$SI_OlinkID),])
which(current$SP_ASSAY != current$SI_ASSAY)
current$ABS_DIFF <- as.numeric(abs(current$SP_OR - current$SI_OR))
current$ASSAY <- ifelse(current$SP_OlinkID == current$SP_OlinkID, current$SP_ASSAY, NA)
which(is.na(current$ASSAY))

current$Label <- NULL
current$Label <- ifelse(current$SP_OlinkID %in% c(sig_g1_top10) & !(current$SP_OlinkID %in% sig_overlaps), current$SP_ASSAY, "")
current$GROUP <- ifelse(current$SP_OlinkID %in% c(sig_g1) & !(current$SP_OlinkID %in% sig_overlaps), "Plasma ALS VS SP (FDR < 0.05)", "")
current$Label2 <- ifelse(current$SI_OlinkID %in% c(sig_g2_top10) & !(current$SI_OlinkID %in% sig_overlaps), current$SI_ASSAY, "")
current$GROUP2 <- ifelse(current$SI_OlinkID %in% c(sig_g2) & !(current$SI_OlinkID %in% sig_overlaps), "Plasma ALS VS SI (FDR < 0.05)", "")
current[which(current$Label == ""),"Label"] <- current[which(current$Label == ""),"Label2"]
current[which(current$GROUP == ""),"GROUP"] <- current[which(current$GROUP == ""),"GROUP2"]
current$Label2 <- ifelse(current$SI_OlinkID %in% c(sig_overlaps), current$SI_ASSAY, "")
current$GROUP2 <- ifelse(current$SI_OlinkID %in% sig_overlaps, "SP & SI  - ALS VS Control (FDR < 0.05)", "")
current[which(current$Label == ""),"Label"] <- current[which(current$Label == ""),"Label2"]
current[which(current$GROUP == ""),"GROUP"] <- current[which(current$GROUP == ""),"GROUP2"]
current[which(current$GROUP == ""),"GROUP"] <- "Others"
current$GROUP <- factor(current$GROUP, levels = c("Plasma ALS VS SP (FDR < 0.05)",
                                                  "Plasma ALS VS SI (FDR < 0.05)",
                                                  "SP & SI  - ALS VS Control (FDR < 0.05)",
                                                  "Others"))
fig3_SPSI <- current
###################################################################
###################################################################
# Hazard Ratio

data.als <- do.call(rbind.data.frame, data)
data.als <- data.als[which(data.als$GROUP == "ALS"),]
data.als <- split(data.als, data.als$DATA_TYPE)

choices <- c("NPX_RAW","NPX_COV")
ratios_out <- NULL
cox_out <- NULL

for(l in 1:length(choices)){
  for(k in 1:length(data.als)){
    cname <- names(data.als)[k]
    print(cname)
    current <- data.als[[k]]
    current$PLATE <- factor(current$PlateID, levels = sort(unique(current$PlateID)))
    assays <- unique(current$OlinkID)
    for(i in 1:length(assays)){
      print(assays[i])
      temp <- current[which(current$OlinkID == assays[i]),]
      
      temp$TSTART <- as.Date(temp$TSTART, origin = "1899-12-30")
      temp$TSTOP_VENT_DEATH <- as.Date(temp$TSTOP_VENT_DEATH)
      temp$TSTOP_DEATH <- as.Date(temp$TSTOP_DEATH)
      
      temp$DEATH_VENT_DEATH <- ifelse(temp$SURVIVAL_VENT_DEATH == "DEAD", 1, 0)
      temp$DEATH_DEATH <- ifelse(temp$SURVIVAL_DEATH == "DEAD", 1, 0)
      
      temp$NEAREST_FRS_SCORE <- as.numeric(as.character(temp$NEAREST_FRS_SCORE))
      temp$BP_HIGH <- factor(temp$BP_HIGH, levels = sort(unique(as.numeric(as.character(temp$BP_HIGH)))))
      temp$SMOKING <- factor(temp$SMOKING, levels = sort(unique(as.numeric(as.character(temp$SMOKING)))))
      temp$SEX <- factor(temp$SEX, levels = c("F","M"))
      temp$ONSET_TYPE_LARGE <- factor(temp$ONSET_TYPE_LARGE, levels = c("Spinal","Bulbar","Others"))
      temp$CLOSEST_BMI <- as.numeric(as.character(temp$CLOSEST_BMI))
      temp$NPX <- temp$NPX_RAW
      temp <- temp[which(!is.na(temp$NEAREST_FRS_SCORE)),]
      
      if(choices[l] == "NPX_RAW"){
        tt <- NULL
        tt <- temp[which(temp$TSTART < temp$TSTOP_DEATH),]
        coxph_out_death <- coxph(Surv(time = as.numeric(tt$TSTART), time2 = as.numeric(as.Date(tt$TSTOP_DEATH)), DEATH_DEATH) ~ NPX, data = tt)
      }else if (choices[l] == "NPX_COV"){
        tt <- NULL
        tt <- temp[which(temp$TSTART < temp$TSTOP_DEATH),]
        coxph_out_death <- coxph(Surv(time = as.numeric(tt$TSTART), time2 = as.numeric(tt$TSTOP_DEATH), DEATH_DEATH) ~ NPX + BP_HIGH + SMOKING + CLOSEST_BMI + MND_DIAGNOSIS_AGE + SEX + NUM_FREEZING_DAYS + PLATE + ONSET_TYPE_LARGE + NEAREST_FRS_SCORE + DIAGNOSTIC_DELAY, data = tt)
      }
      
      cox_out[[length(cox_out)+1]] <- coxph_out_death
      names(cox_out)[length(cox_out)] <- paste("DEATH",cname, choices[l], assays[i], sep = "|")
      
      ratios_out <- rbind(ratios_out, data.frame(DATA_TYPE = cname,
                                                 OUTCOME = "DEATH",
                                                 CHOICE = choices[l],
                                                 OlinkID = assays[i],
                                                 Assay = unique(temp$Assay),
                                                 Panel = unique(temp$Panel),
                                                 COEF = summary(coxph_out_death)$coefficients["NPX","coef"],
                                                 EXP_COEF = summary(coxph_out_death)$coefficients["NPX","exp(coef)"],
                                                 P = summary(coxph_out_death)$coefficients["NPX","Pr(>|z|)"]))
    }
    
  }
}

ratios_out$SPLIT_TERM <- paste(ratios_out$DATA_TYPE, ratios_out$OUTCOME, ratios_out$CHOICE, sep = "_")
ratios_out <- split(ratios_out, ratios_out$SPLIT_TERM)
ratios_out <- lapply(ratios_out, function(x){
  x <- data.frame(x, PADJUST = p.adjust(x$P, method = "BH"))
})
ratios_out <- do.call(rbind.data.frame, ratios_out)
ratios_out$GROUP <- ifelse(ratios_out$P < 0.05, "p < 0.05", "p > 0.05")
ratios_out[which(ratios_out$PADJUST < 0.1),"GROUP"] <- "FDR < 0.1"
ratios_out$Label <- ""
ratios_out <- ratios_out[order(ratios_out$PADJUST, decreasing = F),]
ratios_out$Label <- ifelse(ratios_out$P < 0.05, ratios_out$Assay, "")
ratios_out <- split(ratios_out, ratios_out$SPLIT_TERM)
for(i in 1:length(ratios_out)){
  x <- ratios_out[[i]]
  if(length(which(!is.na(unique(x$Label)))) > 10){
    x$Label <- x$Assay
    x$Label[c(11:nrow(x))] <- ""
    ratios_out[[i]] <- x
  }
}

ratios_out <- do.call(rbind.data.frame, ratios_out)
ratios_out <- ratios_out[order(ratios_out$P, decreasing = F),]
fig4_HR <- ratios_out
fig4_COX <- cox_out

fig4_HR <- fig4_HR[which(fig4_HR$PADJUST < 0.1),]
colnames(fig4_HR)
colnames(fig4_HR) <- gsub(" ","_", stringr::str_to_title(gsub("_"," ", colnames(fig4_HR))))
colnames(fig4_HR) <- gsub("^Coef$","log(Hazard_Ratio)",colnames(fig4_HR), ignore.case = T)
colnames(fig4_HR) <- gsub("Exp_Coef","Hazard_Ratio",colnames(fig4_HR), ignore.case = T)
colnames(fig4_HR) <- gsub("Olinkid","OlinkID",colnames(fig4_HR))
colnames(fig4_HR)
fig4_HR <- fig4_HR[order(fig4_HR$Hazard_Ratio, decreasing = T),]
fig4_HR <- fig4_HR[,c("Data_Type","Choice","OlinkID","Assay","Panel","log(Hazard_Ratio)","Hazard_Ratio","P","Padjust")]

xlsx::write.xlsx(fig4_HR[which(fig4_HR$Data_Type == "PLASMA" & fig4_HR$Choice == "NPX_RAW"),],"ALS_Olink_Table_S3.xlsx", sheetName = "Plasma_Unadjusted", row.names = F, append = F)
xlsx::write.xlsx(fig4_HR[which(fig4_HR$Data_Type == "PLASMA" & fig4_HR$Choice == "NPX_COV"),],"ALS_Olink_Table_S3.xlsx", sheetName = "Plasma_Adjusted", row.names = F, append = T)
xlsx::write.xlsx(fig4_HR[which(fig4_HR$Data_Type == "CSF" & fig4_HR$Choice == "NPX_RAW"),],"ALS_Olink_Table_S3.xlsx", sheetName = "CSF_Unadjusted", row.names = F, append = T)
xlsx::write.xlsx(fig4_HR[which(fig4_HR$Data_Type == "CSF" & fig4_HR$Choice == "NPX_COV"),],"ALS_Olink_Table_S3.xlsx", sheetName = "CSF_Adjusted", row.names = F, append = T)

########################################################################
########################################################################
# SURVIVAL

library(survminer)
library(RTCGA)

fig5_cutoffs <- NULL
fig5_SURVOUT <- NULL

for(i in 1:length(data.als)){
  current <- data.als[[i]]
  cname <- names(data.als)[i]
  assays <- NULL
  assays <- unique(current$OlinkID)
  for(j in 1:length(assays)){
    temp <- NULL
    temp <- current[which(current$OlinkID == assays[j]),]
    cmarker <- unique(temp$Assay)
    if(nrow(temp)> 0){
      temp$NPX <- temp$NPX_RAW
      temp$DEATH <- ifelse(temp$SURVIVAL_DEATH == "DEAD", 1, 0)
      temp$TIME_SINCE_DIAGNOSIS = difftime(temp$TSTOP_DEATH,temp$TSTART, units = c("days"))
      surv_out <- NULL
      surv_out <- surv_cutpoint(
        temp,
        # minprop = 0.2,
        time = "TIME_SINCE_DIAGNOSIS",
        event = "DEATH",
        variables = c("NPX")
      )
      
      surv_cat <- surv_categorize(surv_out) 
      fig5_cutoffs <- rbind(fig5_cutoffs, data.frame(TYPE = cname, Olink_ID = assays[j], Assay = cmarker, CUT_POINT = surv_out$cutpoint$cutpoint))
      
      fig5_SURVOUT[[length(fig5_SURVOUT)+1]] <- surv_out
      names(fig5_SURVOUT)[length(fig5_SURVOUT)] <- paste(cname, assays[j], cmarker, sep = "|")
      
    }
    
  }
  
}

########################################################################
########################################################################
library(ggformula)
library(mgcv)

all_models <- NULL
frs_gam_table <- NULL

for(j in 1:length(data.als)){
  current <- data.als[[j]]
  current$NPX <- current$NPX_RAW
  current$patient_id <- pheno[match(gsub("(.*)\\|.*","\\1",current$SampleID),gsub("(.*)\\|.*","\\1",pheno$OLINK_ID)),"patient_id"]
  assays <- unique(current$OlinkID)
  for(i in 1:length(assays)){
    temp <- current[which(current$OlinkID == assays[i]),]
    ccut <- fig5_cutoffs[which(fig5_cutoffs$TYPE == names(data.als)[j] & fig5_cutoffs$Olink_ID == assays[i]),"CUT_POINT"]
    cname <- paste(unique(temp$Assay), "(", unique(temp$OlinkID), ")", sep = "")
    temp$SURVIVAL_GROUP <- ifelse(temp$NPX > ccut, "High","Low")
    frs$SURVIVAL_GROUP <- temp[match(frs$patient_id, temp$patient_id),"SURVIVAL_GROUP"]
    cfrs <- frs[which(!is.na(frs$SURVIVAL_GROUP)),]
    cfrs$SURVIVAL_GROUP <- factor(cfrs$SURVIVAL_GROUP, levels = c("Low","High"))
    cfrs$patient_id <- as.numeric(as.character(cfrs$patient_id))
    cmodel <- gam(als_alsfrs_calc_score ~ s(as.numeric(as.character(TIME_SINCE_DIAGNOSIS_MONTHS))) + SURVIVAL_GROUP + s(patient_id, bs = 're') + patient_id:TIME_SINCE_DIAGNOSIS_MONTHS, data = cfrs)
    
    all_models[[length(all_models)+1]] <- cmodel
    names(all_models)[length(all_models)] <- paste(names(data.als)[j],cname, sep = "_")
    cp <- summary(cmodel)$p.table[grep("^SURVIVAL_GROUPHigh$",row.names(summary(cmodel)$p.table)),"Pr(>|t|)"]
    frs_gam_table <- rbind(frs_gam_table, data.frame(TYPE = names(data.als)[j], Assay = unique(temp$Assay), OlinkID = unique(temp$OlinkID), P = cp, SE = summary(cmodel)$p.table[grep("^SURVIVAL_GROUPHigh$",row.names(summary(cmodel)$p.table)),"Std. Error"], ESTIMATE = summary(cmodel)$p.table[grep("^SURVIVAL_GROUPHigh$",row.names(summary(cmodel)$p.table)),"Estimate"]))
  }
}

frs_gam_table <- split(frs_gam_table, frs_gam_table$TYPE)
frs_gam_table <- lapply(frs_gam_table, function(x){
  x <- data.frame(x, PADJUST = p.adjust(x$P, method = "BH"))
})

frs_gam_table <- do.call(rbind.data.frame, frs_gam_table)
frs_gam_table <- frs_gam_table[order(frs_gam_table$P, decreasing = F),]
fig7_data <- frs_gam_table
fig7_models <- all_models

fig7_data <- fig7_data[which(fig7_data$PADJUST < 0.01),]
table(fig7_data$TYPE)
colnames(fig7_data)
colnames(fig7_data) <- gsub(" ","_", stringr::str_to_title(gsub("_"," ", colnames(fig7_data))))
colnames(fig7_data) <- gsub("^Se$","SE",colnames(fig7_data), ignore.case = T)
colnames(fig7_data) <- gsub("Estimate","Estimate(High_Expr_Group)",colnames(fig7_data), ignore.case = T)
colnames(fig7_data) <- gsub("Olinkid","OlinkID",colnames(fig7_data))
colnames(fig7_data)
fig7_data <- fig7_data[order(fig7_data$P, decreasing = F),]

xlsx::write.xlsx(fig7_data[which(fig7_data$Type == "PLASMA"),],"ALS_Olink_Table_S4.xlsx", sheetName = "Plasma", row.names = F, append = F)
xlsx::write.xlsx(fig7_data[which(fig7_data$Type == "CSF"),],"ALS_Olink_Table_S4.xlsx", sheetName = "CSF", row.names = F, append = T)

########################################################################
########################################################################
data <- readRDS(paste("Data/Olink_ALS_Data_Longitudinal.RDS", sep = ""))
pheno <- readRDS("Data/Olink_ALS_Metadata_Longitudinal.RDS")
########################################################################
# 1. Longitudinal Analysis
# Consider only those with multiple measurements
current <- unique(data[,c("DATA_TYPE","SampleID","SID","RID")]) # ,"NPX"
# Assumed that all technical replicates were from the same sampling date sample
current$SAMPLING_DATE <- pheno[match(gsub("(.*)\\|.*","\\1",current$RID), pheno$RID),c("SAMPLING_DATE")]
current$Month <- gsub(".*-(.*)-.*","\\1",current$SAMPLING_DATE)
current <- current[order(current$SAMPLING_DATE, decreasing = F),]
current <- split(current, current$SID)
current <- lapply(current, function(x){
  x <- split(x, x$DATA_TYPE)
})
current <- lapply(current, function(x){
  lapply(x, function(y){
    y <- data.frame(y, Visit = paste("Visit",1:nrow(y), sep = ""))
  })
})
current <- lapply(current, function(x){
  x <- do.call(rbind.data.frame, x)
})
current <- do.call(rbind.data.frame, current)

current$SAMPLING_DATE <- as.Date(current$SAMPLING_DATE)
current$SAMPLING_MONTHS <- as.Date(gsub(".*?-(.*)","\\1",current$SAMPLING_DATE), "%m-%d")

# With adjusted NPX
covariates <- c("BP_HIGH","SMOKING","CLOSET_BMI","MND_DIAGNOSIS_AGE","SEX","NUM_FREEZING_DAYS","PlateID","FAMILY_NUMBER")

pheno$YEARS_AFTER_DIAGNOSIS_RAW <- as.Date(pheno$SAMPLING_DATE) - as.Date(pheno$MND_DIAGNOSIS_DATE)
pheno$YEARS_AFTER_DIAGNOSIS_RAW <- round(pheno$YEARS_AFTER_DIAGNOSIS_RAW/30.417, digit=0) # 365.25
pheno$YEARS_AFTER_DIAGNOSIS <- pheno$YEARS_AFTER_DIAGNOSIS_RAW

long_sids <- pheno[which(pheno$YEARS_AFTER_DIAGNOSIS <=36),]
long_sids <- long_sids[which(long_sids$GROUP == "ALS"),"SID"]
long_data <- data[which(data$SID %in% long_sids),]
long_data$YEARS_AFTER_DIAGNOSIS <- pheno[match(gsub("(.*)\\|.*","\\1",long_data$SampleID), pheno$OLINK_ID),"YEARS_AFTER_DIAGNOSIS"]
View(long_data[which(is.na(long_data$YEARS_AFTER_DIAGNOSIS)),])
long_data$OID <- paste(long_data$Assay, long_data$OlinkID, sep = "|")

library(lme4)
library(lmerTest)
long_models <- NULL
long_stats <- NULL
proteins <- unique(long_data$OID)
datatypes <- unique(long_data$DATA_TYPE)
for(i in 1:length(proteins)){
  print(i)
  cmodel <- NULL
  x <- NULL
  plotx <- NULL
  for(j in 1:length(datatypes)){
    x <- long_data[which(long_data$OID == proteins[i] & long_data$DATA_TYPE == datatypes[j]),]
    x$CLOSEST_BMI <- as.numeric(as.character(x$CLOSEST_BMI))
    x <- x[which(x$YEARS_AFTER_DIAGNOSIS <= 36),] # quantile(x$YEARS_AFTER_DIAGNOSIS, 0.95)
    x <- x[which(x$YEARS_AFTER_DIAGNOSIS >= 0),]
    x <- x[which(x$SID %in% names(table(x$SID))[which(table(x$SID) > 1)]),]
    x$NPX <- x$NPX_RAW
    if(nrow(x) > 0){
      print("More than 0!")
      cmodel<-lmer(NPX~YEARS_AFTER_DIAGNOSIS+(1+YEARS_AFTER_DIAGNOSIS|SID)+ BP_HIGH + SMOKING + CLOSEST_BMI + MND_DIAGNOSIS_AGE + 
                     SEX + NUM_FREEZING_DAYS + PlateID, data=x)
      long_models[[length(long_models)+1]] <- cmodel
      names(long_models)[length(long_models)] <- paste(datatypes[j],proteins[i], sep = "|")
      long_stats <- rbind(long_stats, data.frame(DATA_TYPE = datatypes[j],
                                                 OID = proteins[i],
                                                 Assay = unique(x$Assay),
                                                 OlinkID = unique(x$OlinkID),
                                                 ESTIMATE = summary(cmodel)$coefficients["YEARS_AFTER_DIAGNOSIS","Estimate"],
                                                 P = summary(cmodel)$coefficients["YEARS_AFTER_DIAGNOSIS","Pr(>|t|)"]))
      plotx <- rbind(plotx,x)
    }
  }
}

fig8_long_models <- long_models

long_stats <- split(long_stats, long_stats$DATA_TYPE)
long_stats <- lapply(long_stats, function(x){
  x <- data.frame(x, PADJUST = p.adjust(x$P, method = "BH"))
})
long_stats <- do.call(rbind.data.frame, long_stats)
long_stats <- long_stats[order(long_stats$P, decreasing = F),]
fig8_long_stats <- long_stats
fig8_long_stats <- fig8_long_stats[which(fig8_long_stats$PADJUST < 0.1),]
colnames(fig8_long_stats)
colnames(fig8_long_stats) <- gsub(" ","_", stringr::str_to_title(gsub("_"," ", colnames(fig8_long_stats))))
colnames(fig8_long_stats) <- gsub("Olinkid","OlinkID",colnames(fig8_long_stats))
fig8_long_stats$Odds_Ratio <- exp(fig8_long_stats$Estimate)
colnames(fig8_long_stats)
fig8_long_stats <- fig8_long_stats[,c("Data_Type","Assay","OlinkID","Odds_Ratio","P","Padjust")]
fig8_long_stats <- fig8_long_stats[order(fig8_long_stats$Odds_Ratio, decreasing = T),]

xlsx::write.xlsx(fig8_long_stats[which(fig8_long_stats$Data_Type == "PLASMA"),],"ALS_Olink_Table_S5.xlsx", sheetName = "Plasma", row.names = F)

########################################################################################
########################################################################################
